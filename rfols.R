# ---- RF-OLS (paper-style) in the RIDGE_REFERENCE framework -------------------
# Model name (affects output folders/files only)
model_folder <- "rfols"
#install.packages("doRng")
suppressPackageStartupMessages({
  library(randomForest)
  library(foreach)
  library(doParallel)
  library(doRNG)
})

use_parallel <- TRUE                  # <- toggle here
max_workers  <- 8                    # cap workers to avoid RAM spikes

# ---------------- Utilities (copied from RIDGE_REFERENCE) ---------------------
rmse <- function(e) sqrt(mean(e^2)); mae <- function(e) mean(abs(e))
mad_ <- function(e) median(abs(e - median(e)))

compute_metrics <- function(y_real, yhat, y_all, idx_oos, h) {
  err <- y_real - yhat
  ae  <- abs(err); se <- err^2
  m_rmse <- sqrt(mean(se)); m_mae <- mean(ae)
  m_mad <- median(abs(err - median(err)))
  m_meanad <- mean(abs(err - mean(err)))
  rng <- diff(range(y_real)); m_nrmse <- m_rmse / ifelse(rng == 0, 1, rng)
  
  # h-step RW baseline (same convention as rwtext.txt)
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

# Build paper-style design: y AR(4), X lags 1..4 of all non-targets, + 4 PCs (fixed)
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
  
  # PCs on contemporaneous X_t using FIXED center/scale/rotation from the initial train
  X_t_scaled <- scale(Xraw, center = pc_fit$center, scale = pc_fit$scale)
  PC <- as.matrix(X_t_scaled) %*% pc_fit$rotation[, 1:4, drop = FALSE]
  colnames(PC) <- paste0("PC", 1:4)
  
  maxdrop <- max(p_y, x_lags)
  last_anchor <- nrow(dados) - h
  if (last_anchor < maxdrop) {
    stop("Not enough history for horizon h=", h,
         " (need anchor >= ", maxdrop, ", have last_anchor = ", last_anchor, ").")
  }
  anchors <- seq.int(maxdrop, last_anchor)  # indices of t
  
  X <- cbind(ylag, Xlags, PC)[anchors, , drop = FALSE]
  rownames(X) <- dates[anchors]
  y_target <- y[anchors + h]; names(y_target) <- dates[anchors + h]
  list(X = X, y = y_target)
}

# ---------------- RF-OLS core: per-tree OLS on split vars; average ------------
predict_rfols_once <- function(X_tr, y_tr, x_oos, rf_model) {
  inbag <- rf_model$inbag
  if (is.null(inbag)) stop("randomForest did not return inbag counts; set keep.inbag=TRUE.")
  
  # Normalize 'inbag' to a matrix of counts
  inbag_mat <- if (is.list(inbag)) do.call(cbind, inbag) else inbag
  ntree <- ncol(inbag_mat)
  n <- nrow(X_tr); p <- ncol(X_tr)
  
  splitcol_name <- function(tt) {
    cand <- c("split var", "splitvar")
    cand[which(cand %in% colnames(tt))[1]]
  }
  
  preds <- numeric(ntree)
  for (k in seq_len(ntree)) {
    saux <- inbag_mat[, k]; saux[is.na(saux)] <- 0L
    if (sum(saux) <= 1L) { preds[k] <- mean(y_tr); next }
    sboot <- rep.int(seq_len(n), times = saux)
    
    tr <- tryCatch(randomForest::getTree(rf_model, k, labelVar = FALSE), error = function(e) NULL)
    if (is.null(tr)) { preds[k] <- mean(y_tr[sboot]); next }
    scn <- splitcol_name(tr)
    sel <- unique(tr[, scn]); sel <- sort(sel[is.finite(sel) & sel > 0 & sel <= p])
    
    if (length(sel) == 0L) { preds[k] <- mean(y_tr[sboot]); next }
    
    Xb <- X_tr[sboot, sel, drop = FALSE]; yb <- y_tr[sboot]
    Xb_ <- cbind(`(Intercept)` = 1, Xb)
    fit <- tryCatch(stats::lm.fit(x = Xb_, y = yb), error = function(e) NULL)
    if (is.null(fit) || anyNA(fit$coefficients)) { preds[k] <- mean(yb); next }
    
    x_ <- c(1, as.numeric(x_oos[, sel, drop = TRUE]))
    if (length(x_) != length(fit$coefficients)) {
      len <- min(length(x_), length(fit$coefficients))
      preds[k] <- sum(x_[seq_len(len)] * fit$coefficients[seq_len(len)])
    } else {
      preds[k] <- sum(x_ * fit$coefficients)
    }
  }
  mean(preds)
}

# -------- helper: one k-step (used by both serial and parallel paths) ---------
rf_ols_one_k <- function(k, Xall, yall, R, dates, idx_oos, h) {
  origin_k            <- dates[idx_oos[k] - h]
  last_train_anchor_k <- dates[idx_oos[k] - 2*h]
  
  end_row <- which(rownames(Xall) == last_train_anchor_k)
  if (length(end_row) != 1)
    stop("Anchor not found for last_train_anchor_k at h=", h, " (", last_train_anchor_k, ")")
  start_row <- end_row - R + 1
  if (start_row < 1) stop("start_row < 1 at h=", h, " — not enough history. Reduce nprev or h.")
  
  X_tr  <- Xall[start_row:end_row, , drop = FALSE]
  y_tr  <- yall[start_row:end_row]
  
  oos_row <- which(rownames(Xall) == origin_k)
  if (length(oos_row) != 1) stop("Anchor not found for origin_k at h=", h, " (", origin_k, ")")
  x_oos <- Xall[oos_row, , drop = FALSE]
  
  # Standardize + impute + winsorize
  SI     <- standardize_impute(X_tr, x_oos)
  X_tr_s <- SI$X_tr; x_oos_s <- SI$x_oos
  cap <- 8
  X_tr_s[X_tr_s >  cap] <- cap;  X_tr_s[X_tr_s < -cap] <- -cap
  x_oos_s[x_oos_s >  cap] <- cap; x_oos_s[x_oos_s < -cap] <- -cap
  
  # Drop degenerate columns
  keep_cols <- which(apply(X_tr_s, 2, function(z) stats::sd(z) > 0))
  if (length(keep_cols) == 0L) stop("All columns degenerate after standardization at h=", h, ", k=", k)
  if (length(keep_cols) < ncol(X_tr_s)) {
    X_tr_s  <- X_tr_s[, keep_cols, drop = FALSE]
    x_oos_s <- x_oos_s[, keep_cols, drop = FALSE]
  }
  
  # RF with single-tree spec (k nodes), then per-tree OLS prediction
  rf_fit <- randomForest(
    x = X_tr_s, y = y_tr,
    keep.inbag = TRUE,
    ntree = 500  ,       # B
    maxnodes = 20,     # k nodes (terminal nodes) — matches paper snippet
    nodesize = 5  
  )
  predict_rfols_once(X_tr_s, y_tr, x_oos_s, rf_fit)
}

# ---------------- Driver (anchors/rolling identical to RIDGE_REFERENCE) -------
load("dados.rda")  # provides `dados`
if (!("cpi_total" %in% colnames(dados))) stop("cpi_total not found in dados.")

dates     <- rownames(dados)
y_all     <- as.numeric(dados[, "cpi_total"])
nprev     <- 99
idx_oos   <- (nrow(dados) - nprev + 1):nrow(dados)
dates_oos <- dates[idx_oos]
horizons  <- 1:12

# Output dir
out_dir <- file.path("forecasts", model_folder)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Start cluster once (reuse across horizons) for lower overhead
cl <- NULL
if (use_parallel) {
  cores <- max(1, min(parallel::detectCores() - 1, max_workers))
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  RNGkind("L'Ecuyer-CMRG")  # reproducible independent streams
  set.seed(123)             # master seed for doRNG
}

RFOLS_cpi_list <- vector("list", length(horizons))
set.seed(123)  # RF reproducibility in serial path (harmless when doRNG is used)

for (h in horizons) {
  # --------- Initial PCA window & frozen loadings (no leakage) ---------------
  origin_first <- dates[idx_oos[1] - h]                 # issue-time anchor for first OOS
  train_rows   <- which(dates <= origin_first)
  
  other_vars <- setdiff(colnames(dados), "cpi_total")
  X_train0   <- scale(dados[train_rows, other_vars, drop = FALSE], center = TRUE, scale = TRUE)
  pc         <- prcomp(X_train0, center = FALSE, scale. = FALSE)
  pc_fit     <- list(rotation = pc$rotation,
                     center   = attr(X_train0, "scaled:center"),
                     scale    = attr(X_train0, "scaled:scale"))
  
  sup  <- build_xy_paper(dados, "cpi_total", h = h, p_y = 4, x_lags = 4, pc_fit = pc_fit)
  Xall <- sup$X; yall <- sup$y
  
  # --------- Fixed rolling window length R from FIRST origin -----------------
  last_train_anchor_first <- dates[idx_oos[1] - 2*h]
  if (!(last_train_anchor_first %in% rownames(Xall))) {
    stop("Not enough history for h=", h, ": missing last_train_anchor_first = ", last_train_anchor_first)
  }
  keep <- rownames(Xall) <= last_train_anchor_first
  X0   <- Xall[keep, , drop = FALSE]; y0 <- yall[keep]
  R    <- nrow(X0)
  if (R <= 0) stop("Training window length R<=0 at h=", h, ". Check sample size.")
  
  # -------------------- Parallel/serial k-loop -------------------------------
  if (use_parallel) {
    preds <- foreach(k = 1:nprev,
                     .combine = "c",
                     .packages = c("randomForest"),
                     .export = c("rf_ols_one_k", "standardize_impute",
                                 "predict_rfols_once")) %dorng% {
                                   rf_ols_one_k(k, Xall, yall, R, dates, idx_oos, h)
                                 }
    # Optional: if you want the detailed debug for h=1,k=last, re-run once serially
    if (h == 1) {
      k <- nprev
      origin_k            <- dates[idx_oos[k] - h]
      last_train_anchor_k <- dates[idx_oos[k] - 2*h]
      cat("\n[h=1, k=last] origin_k =", origin_k,
          " last_train_anchor_k =", last_train_anchor_k, "\n")
      # You can paste your earlier pretty-print block here if needed.
    }
  } else {
    preds <- numeric(nprev)
    for (k in 1:nprev) {
      preds[k] <- rf_ols_one_k(k, Xall, yall, R, dates, idx_oos, h)
      if (h == 1 && k == nprev) {
        origin_k            <- dates[idx_oos[k] - h]
        last_train_anchor_k <- dates[idx_oos[k] - 2*h]
        cat("\n[h=1, k=last] origin_k =", origin_k,
            " last_train_anchor_k =", last_train_anchor_k, "\n")
      }
    }
  }
  
  names(preds) <- dates_oos
  y_real <- y_all[idx_oos]
  mets   <- compute_metrics(y_real, preds, y_all, idx_oos, h)
  
  RFOLS_cpi_list[[h]] <- list(
    pred     = preds,
    errors   = mets,
    loss_abs = abs(y_real - preds),
    loss_sq  = (y_real - preds)^2
  )
}

if (!is.null(cl)) parallel::stopCluster(cl)

# --------------------------- Save (like ridge) --------------------------------
mat_rfols <- do.call(cbind, lapply(1:12, function(h) RFOLS_cpi_list[[h]]$pred))
colnames(mat_rfols) <- paste0("h", 1:12)
rownames(mat_rfols) <- dates_oos

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.table(mat_rfols, file.path(out_dir, "RFOLS_cpi_total.csv"),
            sep = ";", row.names = TRUE, col.names = TRUE)
save(RFOLS_cpi_list, file = file.path(out_dir, "RFOLS-list_cpi_total.Rdata"))

cat("Saved RF-OLS (", model_folder, ") to ", out_dir,
    " — nprev=", nprev, ", target=cpi_total\n", sep = "")
