# ==== QERT (ranger quantile + extratrees) — paper-style pipeline ====
suppressPackageStartupMessages({
  library(ranger)
  library(foreach)
  library(doParallel)
  library(doRNG)
})

# ---------------- Utilities (RIDGE_REFERENCE) ----------------
rmse <- function(e) sqrt(mean(e^2)); mae <- function(e) mean(abs(e))
mad_ <- function(e) median(abs(e - median(e)))

compute_metrics <- function(y_real, yhat, y_all, idx_oos, h) {
  err <- y_real - yhat
  ae  <- abs(err); se <- err^2
  m_rmse <- sqrt(mean(se)); m_mae <- mean(ae)
  m_mad <- median(abs(err - median(err)))
  m_meanad <- mean(abs(err - mean(err)))
  rng <- diff(range(y_real)); m_nrmse <- m_rmse / ifelse(rng == 0, 1, rng)
  y_lag_h <- y_all[idx_oos - h]   # RW(h) baseline
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

# Build paper-style design: y AR(4), X lags 1..4, + 4 PCs (frozen)
build_xy_paper <- function(dados, yname, h, p_y=4, x_lags=4, pc_fit) {
  stopifnot(yname %in% colnames(dados))
  dates <- rownames(dados)
  y <- as.numeric(dados[, yname]); names(y) <- dates
  Xraw <- dados[, setdiff(colnames(dados), yname), drop = FALSE]
  ylag <- sapply(1:p_y, function(L) c(rep(NA, L), y[1:(length(y)-L)]))
  colnames(ylag) <- paste0("y_lag", 1:p_y)
  make_lag <- function(v, L) c(rep(NA, L), v[1:(length(v)-L)])
  Xlags <- do.call(cbind, lapply(1:x_lags, function(L) {
    xs <- apply(Xraw, 2, make_lag, L = L)
    colnames(xs) <- paste0(colnames(Xraw), "_L", L); xs
  }))
  X_t_scaled <- scale(Xraw, center = pc_fit$center, scale = pc_fit$scale)
  PC <- as.matrix(X_t_scaled) %*% pc_fit$rotation[, 1:4, drop = FALSE]
  colnames(PC) <- paste0("PC", 1:4)
  maxdrop <- max(p_y, x_lags)
  anchors <- maxdrop:(nrow(dados) - h)
  X <- cbind(ylag, Xlags, PC)[anchors, , drop = FALSE]
  rownames(X) <- dates[anchors]
  y_target <- y[anchors + h]; names(y_target) <- dates[anchors + h]
  list(X = X, y = y_target)
}

# ---------- Helper: robust quantile extraction (avoids subscript errors) ----
safe_quantiles <- function(fit, xo, qs = c(0.1, 0.5, 0.9)) {
  out <- try(predict(fit, as.data.frame(xo), type = "quantiles",
                     quantiles = qs)$predictions, silent = TRUE)
  if (inherits(out, "try-error") || is.null(out)) return(rep(NA_real_, length(qs)))
  # ranger returns either a length-3 vector or a 1x3 matrix
  if (is.matrix(out)) {
    v <- suppressWarnings(as.numeric(out[1, ]))
  } else {
    v <- suppressWarnings(as.numeric(out))
  }
  # pad/truncate to length(qs)
  if (length(v) < length(qs) || any(!is.finite(v))) {
    v2 <- rep(NA_real_, length(qs))
    v2[seq_len(min(length(v), length(qs)))] <- v[seq_len(min(length(v), length(qs)))]
    return(v2)
  }
  v[seq_along(qs)]
}

# ---------------- Model helper: one (h,k) forecast with QERT -----------------
qert_one_k <- function(k, Xall, yall, dates, idx_oos, h, R,
                       cap = 8, ntree = 3000, mtry_frac = 1/3,
                       min_node = 8, num_random_splits = 5, rng_seed = 12345) {
  
  origin_k <- dates[idx_oos[k] - h]
  last_train_anchor_k <- dates[idx_oos[k] - 2*h]
  
  end_row <- which(rownames(Xall) == last_train_anchor_k)
  if (length(end_row) != 1)
    stop("Anchor not found for last_train_anchor_k at h=", h, " (", last_train_anchor_k, ")")
  start_row <- end_row - R + 1
  if (start_row < 1)
    stop("start_row < 1 at h=", h, " — not enough history. Reduce nprev or h.")
  
  X_tr  <- Xall[start_row:end_row, , drop = FALSE]
  y_tr  <- yall[start_row:end_row]
  
  oos_row <- which(rownames(Xall) == origin_k)
  if (length(oos_row) != 1)
    stop("Anchor not found for origin_k at h=", h, " (", origin_k, ")")
  x_oos <- Xall[oos_row, , drop = FALSE]
  
  # Standardize + winsorize
  SI <- standardize_impute(X_tr, x_oos)
  Xs <- SI$X_tr; xo <- SI$x_oos
  Xs[Xs >  cap] <- cap;  Xs[Xs < -cap] <- -cap
  xo[xo >  cap] <- cap;  xo[xo < -cap] <- -cap
  
  # Guard degeneracy
  keep_cols <- apply(Xs, 2, function(v) is.finite(sd(v)) && sd(v) > 1e-8)
  if (!all(keep_cols)) { Xs <- Xs[, keep_cols, drop = FALSE]; xo <- xo[, keep_cols, drop = FALSE] }
  p <- ncol(Xs)
  if (p == 0) stop("All predictors degenerate after SI at h=", h, ", k=", k)
  
  # Edge case: near-constant design → trivial quantiles (use mean)
  if (nrow(unique(Xs)) <= 1L) {
    yhat_mean <- mean(y_tr)
    return(list(pred = yhat_mean,
                varimp = setNames(numeric(p), colnames(Xs)),
                debug = if (h == 1 && k == length(idx_oos)) list(
                  origin_k = origin_k, last_train_anchor_k = last_train_anchor_k,
                  note = "Xs had <=1 unique row; mean used; quantiles degenerate") else NULL))
  }
  
  # Safe mtry with backoff if ranger complains
  mtry0 <- if (p > 1L) max(1L, min(floor(p * mtry_frac), floor(sqrt(p)), p - 1L)) else 1L
  set.seed(rng_seed + h*1000 + k)
  
  fit <- NULL; mtry_try <- mtry0; attempts <- 0; last_err <- NULL
  repeat {
    attempts <- attempts + 1
    fit <- tryCatch(
      ranger(x = as.data.frame(Xs), y = y_tr,
             num.trees = ntree,
             mtry = mtry_try,
             min.node.size = min_node,
             splitrule = "extratrees",
             num.random.splits = num_random_splits,
             replace = TRUE, sample.fraction = 1.0,
             importance = "impurity_corrected",
             quantreg = TRUE,
             num.threads = 1),
      error = function(e) e
    )
    if (!inherits(fit, "error")) break
    msg <- conditionMessage(fit); last_err <- msg
    if (grepl("sample larger than population", msg, ignore.case = TRUE) && mtry_try > 1L) {
      mtry_try <- max(1L, floor(mtry_try / 2)); next
    }
    # Bail out safely
    yhat_mean <- mean(y_tr)
    return(list(pred = yhat_mean,
                varimp = setNames(numeric(p), colnames(Xs)),
                debug = if (h == 1 && k == length(idx_oos))
                  list(origin_k = origin_k, last_train_anchor_k = last_train_anchor_k,
                       note = paste0("ranger fallback to mean. last_err: ", last_err)) else NULL))
  }
  
  # ---- Diagnostics you asked for -------------------------------------------
  has_splits <- {
    sv <- try(fit$forest$split.varIDs, silent = TRUE)
    if (inherits(sv, "try-error") || is.null(sv)) {
      NA_integer_
    } else if (is.list(sv)) {
      sum(vapply(sv, function(v) sum(v >= 0L, na.rm = TRUE), integer(1)))
    } else if (is.matrix(sv) || is.vector(sv)) {
      sum(sv >= 0L, na.rm = TRUE)
    } else NA_integer_
  }
  n_tr <- nrow(Xs)
  ok_splits <- (n_tr >= 2L * min_node)
  ranger_version <- as.character(utils::packageVersion("ranger"))
  
  # Robust quantile prediction
  qv <- safe_quantiles(fit, xo, qs = c(0.1, 0.5, 0.9))
  yhat <- if (is.finite(qv[2])) qv[2] else mean(y_tr)  # fallback to mean if median is NA
  
  # Variable importance + robust fallback
  vi <- fit$variable.importance
  if (!is.null(vi)) names(vi) <- colnames(Xs)
  if (is.null(vi) || !any(is.finite(vi)) || sum(abs(vi), na.rm = TRUE) == 0) {
    fit_vi <- ranger(
      x = as.data.frame(Xs), y = y_tr,
      num.trees = min(1000L, ntree),
      mtry = mtry_try, min.node.size = min_node,
      splitrule = "variance", replace = TRUE,
      importance = "permutation",
      quantreg = FALSE,
      num.threads = 1
    )
    vi <- fit_vi$variable.importance
    if (!is.null(vi)) names(vi) <- colnames(Xs)
  }
  
  dbg <- NULL
  if (h == 1 && k == length(idx_oos)) {
    print(utils::packageVersion("ranger"))  # (3) version
    dbg <- list(
      origin_k = origin_k, last_train_anchor_k = last_train_anchor_k,
      q10 = qv[1], q50 = qv[2], q90 = qv[3],
      mtry_used = mtry_try, mtry_init = mtry0,
      ntree = ntree, min_node = min_node,
      splitrule = "extratrees", num_random_splits = num_random_splits,
      attempts = attempts, has_splits = has_splits,
      ok_splits = ok_splits, ranger_version = ranger_version
    )
  }
  
  list(pred = yhat, varimp = vi, debug = dbg)
}

# ---------------- Driver ------------------------------------------------------
model_folder <- "qert"
use_parallel <- TRUE
max_workers  <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
winsor_cap   <- 8

load("dados.rda")
if (!("cpi_total" %in% colnames(dados))) stop("cpi_total not found in dados.")

dates     <- rownames(dados)
y_all     <- as.numeric(dados[, "cpi_total"])
nprev     <- 99
idx_oos   <- (nrow(dados) - nprev + 1):nrow(dados)
dates_oos <- dates[idx_oos]
horizons  <- 1:12

dir.create(file.path("forecasts", model_folder), recursive = TRUE, showWarnings = FALSE)
QERT_cpi_list <- vector("list", length(horizons))

cl <- NULL
if (use_parallel) {
  workers <- max(1L, min(max_workers, parallel::detectCores(logical = TRUE) - 1L))
  cl <- parallel::makeCluster(workers)
  doParallel::registerDoParallel(cl)
}

for (h in horizons) {
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
  if (!(last_train_anchor_first %in% rownames(Xall)))
    stop("Not enough history for h=", h, ": missing last_train_anchor_first = ", last_train_anchor_first)
  keep <- rownames(Xall) <= last_train_anchor_first
  X0   <- Xall[keep, , drop = FALSE]; y0 <- yall[keep]
  R    <- nrow(X0)
  if (R <= 0) stop("Training window length R<=0 at h=", h, ". Check sample size.")
  
  if (use_parallel) {
    res_list <- foreach(k = 1:nprev,
                        .packages = "ranger",
                        .export = c("qert_one_k","standardize_impute","build_xy_paper","compute_metrics","safe_quantiles"),
                        .options.RNG = 20240601) %dorng% {
                          qert_one_k(k, Xall, yall, dates, idx_oos, h, R,
                                     cap = winsor_cap, ntree = 3000, mtry_frac = 1/3,
                                     min_node = 8, num_random_splits = 5, rng_seed = 12345)
                        }
    preds <- vapply(res_list, function(x) x$pred, numeric(1))
    varimp_list <- lapply(res_list, function(x) x$varimp)
    
    if (!is.null(res_list[[nprev]]$debug)) {
      d <- res_list[[nprev]]$debug
      cat(sprintf(
        "\n[h=1, k=last] origin_k=%s last_train_anchor_k=%s (q10,q50,q90)=(%.4f,%.4f,%.4f) mtry_used=%d ntree=%d min_node=%d splitrule=%s num.random.splits=%d attempts=%d has_splits=%s ok_splits=%s ranger=%s\n",
        d$origin_k, d$last_train_anchor_k, d$q10, d$q50, d$q90,
        as.integer(d$mtry_used), as.integer(d$ntree), as.integer(d$min_node),
        "extratrees", as.integer(d$num_random_splits), as.integer(d$attempts),
        as.character(d$has_splits), as.character(d$ok_splits), d$ranger_version
      ))
    }
  } else {
    preds <- numeric(nprev); varimp_list <- vector("list", nprev)
    for (k in 1:nprev) {
      outk <- qert_one_k(k, Xall, yall, dates, idx_oos, h, R,
                         cap = winsor_cap, ntree = 3000, mtry_frac = 1/3,
                         min_node = 8, num_random_splits = 5, rng_seed = 12345)
      preds[k] <- outk$pred; varimp_list[[k]] <- outk$varimp
      if (h == 1 && k == nprev && !is.null(outk$debug)) {
        d <- outk$debug
        cat(sprintf(
          "\n[h=1, k=last] origin_k=%s last_train_anchor_k=%s (q10,q50,q90)=(%.4f,%.4f,%.4f) mtry_used=%d ntree=%d min_node=%d splitrule=%s num.random.splits=%d attempts=%d has_splits=%s ok_splits=%s ranger=%s\n",
          d$origin_k, d$last_train_anchor_k, d$q10, d$q50, d$q90,
          as.integer(d$mtry_used), as.integer(d$ntree), as.integer(d$min_node),
          "extratrees", as.integer(d$num_random_splits), as.integer(d$attempts),
          as.character(d$has_splits), as.character(d$ok_splits), d$ranger_version
        ))
      }
    }
  }
  
  names(preds) <- dates_oos
  y_real <- y_all[idx_oos]
  mets   <- compute_metrics(y_real, preds, y_all, idx_oos, h)
  
  QERT_cpi_list[[h]] <- list(
    pred     = preds,                      # point = median (q50)
    errors   = mets,
    loss_abs = abs(y_real - preds),
    loss_sq  = (y_real - preds)^2,
    varimp   = varimp_list                 # per-k VI (corrected impurity or fallback permutation)
  )
}

if (!is.null(cl)) parallel::stopCluster(cl)

# ---------------- Save (same shape/naming as ridge) --------------------------
mat_qert <- do.call(cbind, lapply(1:12, function(h) QERT_cpi_list[[h]]$pred))
colnames(mat_qert) <- paste0("h", 1:12)
rownames(mat_qert) <- dates_oos

outfile_csv  <- file.path("forecasts", model_folder, "QERT_cpi_total.csv")
outfile_rdat <- file.path("forecasts", model_folder, "QERT-list-cpi_total.Rdata")

write.table(mat_qert, outfile_csv, sep = ";", row.names = TRUE, col.names = TRUE)
save(QERT_cpi_list, file = outfile_rdat)

cat("Saved QERT (quantile extratrees; VI w/ robust fallback) to forecasts/",
    model_folder, " — nprev=", nprev, ", target=cpi_total\n", sep = "")
