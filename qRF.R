# ==== QRF (ranger quantile) — paper-style pipeline ====
# Direct h-step; fixed rolling; anchors (no leakage); frozen PCs; nprev=99
# Uses ridge's compute_metrics(), standardize_impute(), build_xy_paper()
# Saves: forecasts/qrf/QRF_cpi_total.csv and QRF-list-cpi_total.Rdata
# With importance = "impurity_corrected" saved per (h,k) in QRF_cpi_list[[h]]$varimp

suppressPackageStartupMessages({
  library(ranger)
  library(foreach)
  library(doParallel)
  library(doRNG)
})

# ---------------- Utilities (exactly as in RIDGE_REFERENCE) ----------------
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
    colnames(xs) <- paste0(colnames(Xraw), "_L", L)
    xs
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

# ---------------- Model helper: one (h,k) forecast with QRF -------------------
qrf_one_k <- function(k, Xall, yall, dates, idx_oos, h, R,
                      cap = 8, ntree = 3000, mtry_frac = 1/3,
                      min_node = 8, rng_seed = 12345) {
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
  
  # Standardize + impute + winsorize
  SI <- standardize_impute(X_tr, x_oos)
  Xs <- SI$X_tr; xo <- SI$x_oos
  Xs[Xs >  cap] <- cap;  Xs[Xs < -cap] <- -cap
  xo[xo >  cap] <- cap;  xo[xo < -cap] <- -cap
  
  # Guard degeneracy
  keep_cols <- apply(Xs, 2, function(v) is.finite(sd(v)) && sd(v) > 1e-8)
  if (!all(keep_cols)) { Xs <- Xs[, keep_cols, drop = FALSE]; xo <- xo[, keep_cols, drop = FALSE] }
  p <- ncol(Xs)
  if (p == 0) stop("All predictors degenerate after SI at h=", h, ", k=", k)
  
  # Effective splittable predictors (some columns can be constant even if sd>0 on doubles)
  p_eff <- sum(apply(Xs, 2, function(v) length(unique(v)) > 1L))
  if (p_eff == 0L) {
    yhat_mean <- mean(y_tr)
    dbg <- list(
      origin_k = origin_k, last_train_anchor_k = last_train_anchor_k,
      note = "No splittable predictors (p_eff=0); mean used",
      n_tr = nrow(Xs), p = p, p_eff = p_eff,
      mtry_used = NA, mtry_init = NA, attempts = 0L,
      has_splits = NA, ok_splits = (nrow(Xs) >= 2*min_node),
      ranger_version = as.character(utils::packageVersion("ranger")),
      q10 = NA_real_, q50 = yhat_mean, q90 = NA_real_
    )
    return(list(pred = yhat_mean, varimp = setNames(numeric(p), colnames(Xs)), debug = dbg))
  }
  
  # Safe mtry + retry backoff (cap by p_eff)
  mtry0 <- if (p_eff > 1L) max(1L, min(floor(p * mtry_frac), floor(sqrt(p_eff)), p_eff, p - 1L)) else 1L
  set.seed(rng_seed + h*1000 + k)
  
  fit <- NULL; mtry_try <- mtry0; attempts <- 0; last_err <- NULL
  repeat {
    attempts <- attempts + 1
    fit <- tryCatch(
      ranger(x = as.data.frame(Xs), y = y_tr,
             quantreg = TRUE, num.trees = ntree,
             mtry = mtry_try, min.node.size = min_node,
             splitrule = "variance",
             replace = TRUE, sample.fraction = 1.0,
             keep.inbag = TRUE,
             importance = "impurity_corrected",
             num.threads = 1),
      error = function(e) e
    )
    if (!inherits(fit, "error")) break
    msg <- conditionMessage(fit); last_err <- msg
    
    # Robust match for all common phrasings of the sampling error
    if (grepl("cannot take a sample larger than the population", msg, ignore.case = TRUE) ||
        grepl("sample larger than the population", msg, ignore.case = TRUE) ||
        grepl("sample larger than population", msg, ignore.case = TRUE)) {
      if (mtry_try > 1L) {
        mtry_try <- max(1L, min(floor(mtry_try / 2), p_eff))
        next
      }
    }
    # Irrecoverable: fallback to mean so the run continues
    yhat_mean <- mean(y_tr)
    dbg <- list(
      origin_k = origin_k, last_train_anchor_k = last_train_anchor_k,
      note = paste0("ranger fallback to mean. last_err: ", last_err),
      n_tr = nrow(Xs), p = p, p_eff = p_eff,
      mtry_used = mtry_try, mtry_init = mtry0, attempts = attempts,
      has_splits = NA, ok_splits = (nrow(Xs) >= 2*min_node),
      ranger_version = as.character(utils::packageVersion("ranger")),
      q10 = NA_real_, q50 = yhat_mean, q90 = NA_real_
    )
    return(list(pred = yhat_mean, varimp = setNames(numeric(p), colnames(Xs)), debug = dbg))
  }
  
  # ---- Diagnostics right after successful fit ----
  has_splits <- {
    sv <- try(fit$forest$split.varIDs, silent = TRUE)
    if (inherits(sv, "try-error") || is.null(sv)) NA_integer_ else sum(sv >= 0)
  }
  n_tr <- nrow(Xs)
  ok_splits <- (n_tr >= 2 * min_node)
  ranger_version <- as.character(utils::packageVersion("ranger"))
  
  # Prediction (quantiles)
  q_raw <- predict(fit, as.data.frame(xo), type = "quantiles",
                   quantiles = c(0.1, 0.5, 0.9))$predictions
  qv <- if (is.matrix(q_raw)) as.numeric(q_raw[1, ]) else as.numeric(q_raw)
  yhat <- qv[2]
  
  # Variable importance (+ permutation fallback if all zero/missing)
  vi <- fit$variable.importance
  if (!is.null(vi)) names(vi) <- colnames(Xs)
  if (is.null(vi) || !any(is.finite(vi)) || sum(abs(vi)) == 0) {
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
  
  dbg <- list(
    origin_k = origin_k,
    last_train_anchor_k = last_train_anchor_k,
    q10 = qv[1], q50 = qv[2], q90 = qv[3],
    mtry_used = mtry_try, mtry_init = mtry0,
    ntree = ntree, min_node = min_node,
    attempts = attempts,
    n_tr = n_tr, p = p, p_eff = p_eff,
    has_splits = has_splits,
    ok_splits = ok_splits,
    ranger_version = ranger_version,
    note = NULL
  )
  
  list(pred = yhat, varimp = vi, debug = dbg)
}

# ---------------- Driver ------------------------------------------------------
model_folder <- "qrf"
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
QRF_cpi_list <- vector("list", length(horizons))

# Optional parallel backend (parallelize across k; keep ranger num.threads=1)
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
                        .packages = c("ranger"),
                        .export = c("qrf_one_k","standardize_impute","build_xy_paper",
                                    "compute_metrics"),
                        .options.RNG = 20240601) %dorng% {
                          qrf_one_k(k, Xall, yall, dates, idx_oos, h, R,
                                    cap = winsor_cap, ntree = 3000, mtry_frac = 1/3,
                                    min_node = 8, rng_seed = 12345)
                        }
    preds <- vapply(res_list, function(x) x$pred, numeric(1))
    varimp_list <- lapply(res_list, function(x) x$varimp)
    
    if (!is.null(res_list[[nprev]]$debug)) {
      d <- res_list[[nprev]]$debug
      cat(sprintf(
        "\n[h=1, k=last] origin_k=%s last_train_anchor_k=%s (q10,q50,q90)=(%.4f,%.4f,%.4f) mtry_used=%s ntree=%d min_node=%d attempts=%d | n_tr=%d p=%d p_eff=%d has_splits=%s ok_splits=%s ranger=%s\n",
        d$origin_k, d$last_train_anchor_k,
        as.numeric(d$q10), as.numeric(d$q50), as.numeric(d$q90),
        as.character(d$mtry_used), as.integer(d$ntree), as.integer(d$min_node), as.integer(d$attempts),
        as.integer(d$n_tr), as.integer(d$p), as.integer(d$p_eff),
        as.character(d$has_splits), as.character(d$ok_splits),
        d$ranger_version
      ))
    }
  } else {
    preds <- numeric(nprev)
    varimp_list <- vector("list", nprev)
    for (k in 1:nprev) {
      outk <- qrf_one_k(k, Xall, yall, dates, idx_oos, h, R,
                        cap = winsor_cap, ntree = 3000, mtry_frac = 1/3,
                        min_node = 8, rng_seed = 12345)
      preds[k] <- outk$pred
      varimp_list[[k]] <- outk$varimp
      if (h == 1 && k == nprev && !is.null(outk$debug)) {
        d <- outk$debug
        cat(sprintf(
          "\n[h=1, k=last] origin_k=%s last_train_anchor_k=%s (q10,q50,q90)=(%.4f,%.4f,%.4f) mtry_used=%s ntree=%d min_node=%d attempts=%d | n_tr=%d p=%d p_eff=%d has_splits=%s ok_splits=%s ranger=%s\n",
          d$origin_k, d$last_train_anchor_k,
          as.numeric(d$q10), as.numeric(d$q50), as.numeric(d$q90),
          as.character(d$mtry_used), as.integer(d$ntree), as.integer(d$min_node), as.integer(d$attempts),
          as.integer(d$n_tr), as.integer(d$p), as.integer(d$p_eff),
          as.character(d$has_splits), as.character(d$ok_splits),
          d$ranger_version
        ))
      }
    }
  }
  
  names(preds) <- dates_oos
  y_real <- y_all[idx_oos]
  mets   <- compute_metrics(y_real, preds, y_all, idx_oos, h)
  
  QRF_cpi_list[[h]] <- list(
    pred     = preds,
    errors   = mets,
    loss_abs = abs(y_real - preds),
    loss_sq  = (y_real - preds)^2,
    varimp   = varimp_list
  )
}

if (!is.null(cl)) parallel::stopCluster(cl)

# ---------------- Save (same shape/naming as ridge) --------------------------
mat_qrf <- do.call(cbind, lapply(1:12, function(h) QRF_cpi_list[[h]]$pred))
colnames(mat_qrf) <- paste0("h", 1:12)
rownames(mat_qrf) <- dates_oos

outfile_csv  <- file.path("forecasts", model_folder, "QRF_cpi_total.csv")
outfile_rdat <- file.path("forecasts", model_folder, "QRF-list-cpi_total.Rdata")

write.table(mat_qrf, outfile_csv, sep = ";", row.names = TRUE, col.names = TRUE)
save(QRF_cpi_list, file = outfile_rdat)

cat("Saved QRF (ranger quantreg; impurity_corrected VI) to forecasts/", model_folder,
    " — nprev=", nprev, ", target=cpi_total\n", sep = "")
