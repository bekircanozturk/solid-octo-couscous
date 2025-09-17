# ---- ADAPTIVE LASSO → RANDOM FOREST (paper-style, fixed rolling) ------------
# Direct h-step, fixed rolling (nprev=99), 4 AR, 4 lags of X, + 4 frozen PCs
# Selection via glmnet + BIC (df = #nonzeros), then RF on selected features.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(foreach)
  library(doParallel)
  library(doRNG)
})

# ---------------- Utilities (RIDGE_REFERENCE) --------------------------------
rmse <- function(e) sqrt(mean(e^2)); mae <- function(e) mean(abs(e))
mad_ <- function(e) median(abs(e - median(e)))

compute_metrics <- function(y_real, yhat, y_all, idx_oos, h) {
  err <- y_real - yhat
  ae  <- abs(err); se <- err^2
  m_rmse <- sqrt(mean(se)); m_mae <- mean(ae)
  m_mad <- median(abs(err - median(err)))
  m_meanad <- mean(abs(err - mean(err)))
  rng <- diff(range(y_real)); m_nrmse <- m_rmse / ifelse(rng == 0, 1, rng)
  # h-step RW baseline
  y_lag_h <- y_all[idx_oos - h]
  base_err <- y_real - y_lag_h
  base_rmse <- sqrt(mean(base_err^2))
  eps <- .Machine$double.eps
  m_mrae <- mean(ae / pmax(abs(base_err), eps))
  denom_mase <- mean(abs(diff(y_all[(idx_oos[1]-1):idx_oos[length(idx_oos)]])))
  m_mase <- m_mae / pmax(denom_mase, eps)
  m_mape <- 100 * mean(ae / pmax(abs(y_real), eps))
  c(rmse=m_rmse, mae=m_mae, mad=m_mad, mean_ad=m_meanad,
    nrmse=m_nrmse, mrae=m_mrae, mase=m_mase, mape=m_mape,
    rmse_rel_rw = m_rmse / base_rmse)
}

standardize_impute <- function(X_tr, x_oos) {
  mu <- apply(X_tr, 2, function(c) mean(c, na.rm = TRUE))
  sd <- apply(X_tr, 2, function(c) stats::sd(c, na.rm = TRUE))
  sd[!is.finite(sd) | sd < 1e-8] <- 1
  Xs <- sweep(X_tr, 2, mu, "-"); Xs <- sweep(Xs, 2, sd, "/")
  xo <- sweep(x_oos, 2, mu, "-"); xo <- sweep(xo, 2, sd, "/")
  Xs[!is.finite(Xs)] <- 0; xo[!is.finite(xo)] <- 0
  list(X_tr = as.matrix(Xs), x_oos = as.matrix(xo))
}

# ---------------- Selection helpers ------------------------------------------
# BIC selection for glmnet (lasso / elastic net). df = glmnet's df (#nonzeros).
# --- REPLACE the whole function with this ---
select_lambda_bic_glmnet <- function(X, y, alpha = 1, penalty.factor = NULL) {
  X <- as.matrix(X); storage.mode(X) <- "double"
  y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  
  # Always give glmnet an explicit penalty vector
  if (is.null(penalty.factor)) {
    penalty.factor <- rep(1, p)
  } else {
    if (length(penalty.factor) != p) {
      warning(sprintf("[select_lambda_bic_glmnet] penalty.factor len=%d != ncol(X)=%d; recycling.",
                      length(penalty.factor), p))
      penalty.factor <- rep_len(penalty.factor, p)
    }
    penalty.factor[!is.finite(penalty.factor) | penalty.factor <= 0] <- 1
  }
  
  fit <- tryCatch(
    glmnet::glmnet(x = X, y = y,
                   alpha = alpha, family = "gaussian",
                   standardize = FALSE, intercept = TRUE,
                   penalty.factor = penalty.factor),
    error = function(e) {
      stop(sprintf("[glmnet] %s\n  ncol(X)=%d  len(pf)=%d",
                   conditionMessage(e), p, length(penalty.factor)))
    }
  )
  
  Yhat <- predict(fit, newx = X)
  rss  <- colSums((matrix(y, nrow=n, ncol=ncol(Yhat)) - Yhat)^2)
  eps  <- .Machine$double.eps
  bic  <- n * log(pmax(rss, eps) / n) + fit$df * log(n)
  j    <- which.min(bic)
  list(lambda = fit$lambda[j], idx=j, fit_path=fit, bic=bic, df=fit$df)
}


# Build adaptive weights w_i = 1/(|β_i| + 1/sqrt(n)) aligned to X's columns by NAME
build_adalasso_weights <- function(fit_path, s, X) {
  cn <- colnames(X)
  cf <- coef(fit_path, s = s)            # dgCMatrix (p+1 x 1)
  b  <- as.matrix(cf)                    # convert
  slope <- drop(b[-1, , drop = FALSE])   # remove intercept
  names(slope) <- rownames(b)[-1]
  
  beta_aligned <- slope[cn]              # reorder to match X
  beta_aligned[is.na(beta_aligned)] <- 0
  
  n <- nrow(X)
  w <- (abs(beta_aligned) + 1/sqrt(n))^(-1)
  w[!is.finite(w) | w <= 0] <- 1
  as.numeric(w)
}

# ---------------- Paper-style design (RIDGE_REFERENCE) ------------------------
build_xy_paper <- function(dados, yname, h, p_y=4, x_lags=4, pc_fit) {
  stopifnot(yname %in% colnames(dados))
  dates <- rownames(dados)
  y <- as.numeric(dados[, yname]); names(y) <- dates
  Xraw <- dados[, setdiff(colnames(dados), yname), drop = FALSE]
  
  # y AR lags
  ylag <- sapply(1:p_y, function(L) c(rep(NA, L), y[1:(length(y)-L)]))
  colnames(ylag) <- paste0("y_lag", 1:p_y)
  
  # X lags for all non-target variables
  make_lag <- function(v, L) c(rep(NA, L), v[1:(length(v)-L)])
  Xlags <- do.call(cbind, lapply(1:x_lags, function(L) {
    xs <- apply(Xraw, 2, make_lag, L = L)
    colnames(xs) <- paste0(colnames(Xraw), "_L", L)
    xs
  }))
  
  # 4 PCs of contemporaneous X_t with frozen loadings
  X_t_scaled <- scale(Xraw, center = pc_fit$center, scale = pc_fit$scale)
  PC <- as.matrix(X_t_scaled) %*% pc_fit$rotation[, 1:4, drop = FALSE]
  colnames(PC) <- paste0("PC", 1:4)
  
  maxdrop <- max(p_y, x_lags)
  anchors <- maxdrop:(nrow(dados) - h)  # indices of t
  
  X <- cbind(ylag, Xlags, PC)[anchors, , drop = FALSE]
  rownames(X) <- dates[anchors]
  y_target <- y[anchors + h]; names(y_target) <- dates[anchors + h]
  list(X = X, y = y_target)
}

# ---------------- Model-specific helper: one k-step forecast ------------------
# type: "adalasso" (default) or "lasso"; alpha for glmnet path (1 = lasso)
adalassorf_one_k <- function(k, Xall, yall, dates, idx_oos, h, R,
                             cap=8, type="adalasso", alpha=1,
                             min_sel=2,
                             rf_ntree=500, rf_nodesize=5, rng_seed=12345) {
  origin_k <- dates[idx_oos[k] - h]
  last_train_anchor_k <- dates[idx_oos[k] - 2*h]
  
  end_row <- which(rownames(Xall) == last_train_anchor_k)
  if (length(end_row) != 1)
    stop("Anchor not found for last_train_anchor_k at h=", h,
         " (", last_train_anchor_k, ")")
  start_row <- end_row - R + 1
  if (start_row < 1)
    stop("start_row < 1 at h=", h, " — not enough history. Reduce nprev or h.")
  
  X_tr  <- Xall[start_row:end_row, , drop = FALSE]
  y_tr  <- yall[start_row:end_row]
  
  oos_row <- which(rownames(Xall) == origin_k)
  if (length(oos_row) != 1)
    stop("Anchor not found for origin_k at h=", h, " (", origin_k, ")")
  x_oos <- Xall[oos_row, , drop = FALSE]
  
  # Standardize + impute, then winsorize |z|<=cap
  SI <- standardize_impute(X_tr, x_oos)
  Xs <- SI$X_tr; xo <- SI$x_oos
  Xs[Xs >  cap] <- cap;  Xs[Xs < -cap] <- -cap
  xo[xo >  cap] <- cap;  xo[xo < -cap] <- -cap
  
  # Guard degenerate columns
  keep_cols <- apply(Xs, 2, function(v) is.finite(sd(v)) && sd(v) > 1e-8)
  if (!all(keep_cols)) {
    Xs <- Xs[, keep_cols, drop = FALSE]
    xo <- xo[, keep_cols, drop = FALSE]
  }
  # Ensure unique, aligned names (prevents mapping issues)
  if (anyDuplicated(colnames(Xs))) {
    colnames(Xs) <- make.unique(colnames(Xs))
    colnames(xo) <- colnames(Xs)
  }
  stopifnot(is.matrix(Xs), ncol(Xs) > 0, all(is.finite(Xs)))
  p <- ncol(Xs)
  
  # ---- Stage 1: glmnet path + BIC ----
  sel1 <- select_lambda_bic_glmnet(Xs, y_tr, alpha = alpha)
  
  if (type == "adalasso") {
    # Build name-aligned adaptive weights on CURRENT Xs
    w <- build_adalasso_weights(sel1$fit_path, s = sel1$lambda, X = Xs)
    if (length(w) != p) {
      warning(sprintf("[adalassorf_one_k] weight length %d != p=%d; recycling.", length(w), p))
      w <- rep_len(w, p)
    }
    storage.mode(w) <- "double"
    
    # ---- Stage 2 (adaptive) ----
    sel2 <- select_lambda_bic_glmnet(Xs, y_tr, alpha = alpha, penalty.factor = w)
    b2   <- as.numeric(coef(sel2$fit_path, s = sel2$lambda))
    nz   <- which(abs(b2[-1]) > 0)
    lambda2 <- sel2$lambda
  } else {
    b1   <- as.numeric(coef(sel1$fit_path, s = sel1$lambda))
    nz   <- which(abs(b1[-1]) > 0)
    lambda2 <- NA
  }
  
  # Fallback: ensure at least min_sel variables
  if (length(nz) < min_sel) {
    cors <- tryCatch(cor(Xs, y_tr), error = function(e) rep(0, ncol(Xs)))
    ord  <- order(abs(cors), decreasing = TRUE)
    nz   <- unique(c(nz, ord[seq_len(min(min_sel, length(ord)))]))
  }
  
  # ---- Random Forest on selected features (trained on standardized/winsorized X) ----
  Xsub <- Xs[, nz, drop = FALSE]
  xsub <- xo[, nz, drop = FALSE]
  mtry_val <- max(1L, floor(ncol(Xsub) / 3))
  set.seed(rng_seed + h*1000 + k)
  
  rf <- randomForest(x = Xsub, y = y_tr,
                     ntree = rf_ntree, mtry = mtry_val,
                     nodesize = rf_nodesize, importance = TRUE)
  yhat <- as.numeric(predict(rf, newdata = xsub))
  
  # Optional debug payload
  dbg <- NULL
  if (h == 1 && k == length(idx_oos)) {
    imp <- tryCatch(importance(rf)[, 1], error = function(e) NULL)
    top_imp <- if (length(imp)) {
      ord <- order(imp, decreasing = TRUE)
      head(data.frame(var = colnames(Xsub)[ord], IncMSE = imp[ord]), 10)
    } else NULL
    dbg <- list(
      origin_k = origin_k,
      last_train_anchor_k = last_train_anchor_k,
      alpha = alpha, type = type,
      lambda_stage1 = sel1$lambda,
      lambda_stage2 = lambda2,
      n_selected = length(nz),
      top_importance = top_imp
    )
  }
  
  list(pred = yhat, debug = dbg)
}

# ---------------- Driver ------------------------------------------------------
model_folder <- "adalassorf"   # short folder name
type_select  <- "adalasso"     # "adalasso" or "lasso"
alpha_sel    <- 1              # 1 = LASSO; set 0.5 for EN(α=0.5) if desired

use_parallel <- TRUE
max_workers  <- 8
winsor_cap   <- 8
min_selected <- 2

load("dados.rda")  # object: dados
if (!("cpi_total" %in% colnames(dados)))
  stop("cpi_total not found in 'dados' (see dprep).")

dates     <- rownames(dados)
y_all     <- as.numeric(dados[, "cpi_total"])
nprev     <- 99
idx_oos   <- (nrow(dados) - nprev + 1):nrow(dados)
dates_oos <- dates[idx_oos]
horizons  <- 1:12

dir.create(file.path("forecasts", model_folder), recursive = TRUE, showWarnings = FALSE)
ADALASSORF_cpi_list <- vector("list", length(horizons))

# Optional parallel backend
cl <- NULL
if (use_parallel) {
  workers <- max(1L, min(max_workers, parallel::detectCores(logical = TRUE) - 1L))
  cl <- parallel::makeCluster(workers)
  doParallel::registerDoParallel(cl)
}

for (h in horizons) {
  # ----- Frozen PCs from the initial training window (ending at origin_first) -----
  origin_first <- dates[idx_oos[1] - h]
  train_rows   <- which(dates <= origin_first)
  
  other_vars <- setdiff(colnames(dados), "cpi_total")
  X_train0   <- scale(dados[train_rows, other_vars, drop = FALSE], center = TRUE, scale = TRUE)
  pc         <- prcomp(X_train0, center = FALSE, scale. = FALSE)
  pc_fit     <- list(rotation = pc$rotation,
                     center   = attr(X_train0, "scaled:center"),
                     scale    = attr(X_train0, "scaled:scale"))
  
  sup  <- build_xy_paper(dados, "cpi_total", h = h, p_y = 4, x_lags = 4, pc_fit = pc_fit)
  Xall <- sup$X; yall <- sup$y
  
  last_train_anchor_first <- dates[idx_oos[1] - 2*h]
  if (!(last_train_anchor_first %in% rownames(Xall))) {
    stop("Not enough history for h=", h, ": missing last_train_anchor_first = ", last_train_anchor_first)
  }
  keep <- rownames(Xall) <= last_train_anchor_first
  X0   <- Xall[keep, , drop = FALSE]; y0 <- yall[keep]
  R    <- nrow(X0)
  if (R <= 0) stop("Training window length R<=0 at h=", h, ". Check sample size.")
  
  # ---- Forecast loop across k (parallel-safe & reproducible) -----------------
  if (use_parallel) {
    res_list <- foreach(k = 1:nprev,
                        .packages = c("glmnet","randomForest"),
                        .export   = c("adalassorf_one_k","standardize_impute",
                                      "select_lambda_bic_glmnet","build_xy_paper",
                                      "compute_metrics","build_adalasso_weights"),
                        .options.RNG = 20240601) %dorng% {
                          adalassorf_one_k(k, Xall, yall, dates, idx_oos, h, R,
                                           cap = winsor_cap, type = type_select, alpha = alpha_sel,
                                           min_sel = min_selected,
                                           rf_ntree = 500, rf_nodesize = 5, rng_seed = 12345)
                        }
    preds <- vapply(res_list, function(x) x$pred, numeric(1))
    if (!is.null(res_list[[nprev]]$debug)) {
      dbg <- res_list[[nprev]]$debug
      cat("\n[h=1, k=last] origin_k =", dbg$origin_k,
          " last_train_anchor_k =", dbg$last_train_anchor_k, "\n")
      cat("alpha =", dbg$alpha, " type =", dbg$type,
          "  lambda1* =", dbg$lambda_stage1, "  lambda2* =", dbg$lambda_stage2,
          "  #selected =", dbg$n_selected, "\n")
      if (!is.null(dbg$top_importance)) {
        cat("Top 10 RF importance (IncMSE):\n"); print(dbg$top_importance)
      }
    }
  } else {
    preds <- numeric(nprev)
    for (k in 1:nprev) {
      outk <- adalassorf_one_k(k, Xall, yall, dates, idx_oos, h, R,
                               cap = winsor_cap, type = type_select, alpha = alpha_sel,
                               min_sel = min_selected,
                               rf_ntree = 500, rf_nodesize = 5, rng_seed = 12345)
      preds[k] <- outk$pred
      if (h == 1 && k == nprev && !is.null(outk$debug)) {
        dbg <- outk$debug
        cat("\n[h=1, k=last] origin_k =", dbg$origin_k,
            " last_train_anchor_k =", dbg$last_train_anchor_k, "\n")
        cat("alpha =", dbg$alpha, " type =", dbg$type,
            "  lambda1* =", dbg$lambda_stage1, "  lambda2* =", dbg$lambda_stage2,
            "  #selected =", dbg$n_selected, "\n")
        if (!is.null(dbg$top_importance)) {
          cat("Top 10 RF importance (IncMSE):\n"); print(dbg$top_importance)
        }
      }
    }
  }
  
  names(preds) <- dates_oos
  y_real <- y_all[idx_oos]
  mets   <- compute_metrics(y_real, preds, y_all, idx_oos, h)
  
  ADALASSORF_cpi_list[[h]] <- list(
    pred     = preds,
    errors   = mets,
    loss_abs = abs(y_real - preds),
    loss_sq  = (y_real - preds)^2
  )
}

if (!is.null(cl)) parallel::stopCluster(cl)

# ---------------- Save (same shape/naming as ridge/rw/etc.) -------------------
mat_adalassorf <- do.call(cbind, lapply(1:12, function(h) ADALASSORF_cpi_list[[h]]$pred))
colnames(mat_adalassorf) <- paste0("h", 1:12)
rownames(mat_adalassorf) <- dates_oos

outfile_csv  <- file.path("forecasts", model_folder, "ADALASSORF_cpi_total.csv")
outfile_rdat <- file.path("forecasts", model_folder, "ADALASSORF-list-cpi_total.Rdata")

write.table(mat_adalassorf, outfile_csv, sep = ";", row.names = TRUE, col.names = TRUE)
save(ADALASSORF_cpi_list, file = outfile_rdat)

cat("Saved ADALASSORF to ", file.path("forecasts", model_folder),
    " — nprev=", nprev, ", target=cpi_total\n", sep = "")
