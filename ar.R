# ---- AR (paper-style, ridge-plumbing): direct h-step, fixed rolling, NO LEAKAGE ----
suppressPackageStartupMessages({
  library(stats)
})

# ---------------- Utilities copied from RIDGE_REFERENCE ----------------
rmse <- function(e) sqrt(mean(e^2)); mae <- function(e) mean(abs(e))
mad_ <- function(e) median(abs(e - median(e)))

compute_metrics <- function(y_real, yhat, y_all, idx_oos, h) {
  err <- y_real - yhat
  ae  <- abs(err); se <- err^2
  m_rmse <- sqrt(mean(se)); m_mae <- mean(ae)
  m_mad <- median(abs(err - median(err)))
  m_meanad <- mean(abs(err - mean(err)))
  rng <- diff(range(y_real)); m_nrmse <- m_rmse / ifelse(rng == 0, 1, rng)
  
  # h-step RW baseline (aligned with rwtext.R conventions)
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

# Build paper-style design: y AR(4), X lags 1..4 for all non-targets, + 4 PCs (fixed)
# NOTE: For AR we will later subset to AR columns only (ylag*), but we DO NOT remove PCs from the design.
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
  anchors <- maxdrop:(nrow(dados) - h)  # indices of t
  
  X <- cbind(ylag, Xlags, PC)[anchors, , drop = FALSE]
  rownames(X) <- dates[anchors]
  y_target <- y[anchors + h]; names(y_target) <- dates[anchors + h]
  list(X = X, y = y_target)
}

# ---------------- AR-specific helper: choose p via BIC on initial window ----------------
# We ONLY use the AR columns (ylag1..ylag4) for this selection, honoring paper-style p_y=4.
choose_p_bic_ar <- function(X_y_ar, y, p_max = ncol(X_y_ar)) {
  p_max <- min(p_max, ncol(X_y_ar))
  best_p <- 1; best_bic <- Inf
  for (p in 1:p_max) {
    df <- data.frame(y = y, X_y_ar[, 1:p, drop = FALSE])
    fit <- lm(y ~ ., data = df)
    b <- BIC(fit)
    if (is.finite(b) && b < best_bic) { best_bic <- b; best_p <- p }
  }
  best_p
}

# ---------------- Driver ----------------
# Input: dados.rda created by your dprep (matrix with date rownames, "cpi_total" present)
load("dados.rda")
if (!("cpi_total" %in% colnames(dados))) stop("cpi_total not found in dados.")

dates   <- rownames(dados)
y_all   <- as.numeric(dados[, "cpi_total"])
nprev   <- 99
idx_oos <- (nrow(dados) - nprev + 1):nrow(dados)
dates_oos <- dates[idx_oos]

horizons <- 1:12
dir.create("forecasts/ar", recursive = TRUE, showWarnings = FALSE)

AR_cpi_list <- vector("list", length(horizons))

for (h in horizons) {
  # --------- PCA fit on the initial training window (freeze loadings/center/scale) ----------
  origin_first <- dates[idx_oos[1] - h]              # issue-time anchor for the first OOS
  train_rows   <- which(dates <= origin_first)       # indices for X_t PCA fit (non-targets)
  
  other_vars <- setdiff(colnames(dados), "cpi_total")
  if (length(other_vars) == 0L) stop("Need non-target predictors for PCs; none found.")
  X_train0   <- scale(dados[train_rows, other_vars, drop = FALSE], center = TRUE, scale = TRUE)
  pc         <- prcomp(X_train0, center = FALSE, scale. = FALSE)
  pc_fit     <- list(rotation = pc$rotation,
                     center   = attr(X_train0, "scaled:center"),
                     scale    = attr(X_train0, "scaled:scale"))
  
  sup  <- build_xy_paper(dados, "cpi_total", h = h, p_y = 4, x_lags = 4, pc_fit = pc_fit)
  Xall <- sup$X; yall <- sup$y
  
  # last training anchor for the FIRST origin is (origin_first - h)
  last_train_anchor_first <- dates[idx_oos[1] - 2*h]
  if (!(last_train_anchor_first %in% rownames(Xall))) {
    stop("Not enough history for h=", h, " : missing last_train_anchor_first = ", last_train_anchor_first)
  }
  keep <- rownames(Xall) <= last_train_anchor_first
  X0   <- Xall[keep, , drop = FALSE]; y0 <- yall[keep]
  R    <- nrow(X0)
  if (R <= 0) stop("Training window length R<=0 at h=", h, ". Check sample size.")
  
  # Subset AR columns ONLY for the AR model & choose p* via BIC on initial window
  ar_cols <- grep("^y_lag[1-4]$", colnames(X0), value = TRUE)
  if (length(ar_cols) < 1) stop("AR lag columns not found in design.")
  X0_y <- X0[, ar_cols, drop = FALSE]
  p_star <- choose_p_bic_ar(X0_y, y0, p_max = ncol(X0_y))
  
  preds <- numeric(nprev)
  
  for (k in 1:nprev) {
    # --------- Fixed rolling, NO LEAKAGE (exactly ridge anchors) ----------
    origin_k <- dates[idx_oos[k] - h]               # issue-time anchor (features at t)
    last_train_anchor_k <- dates[idx_oos[k] - 2*h]  # last training anchor (features at t-h)
    
    end_row <- which(rownames(Xall) == last_train_anchor_k)
    if (length(end_row) != 1) stop("Anchor not found for last_train_anchor_k at h=", h, " (", last_train_anchor_k, ")")
    start_row <- end_row - R + 1
    if (start_row < 1) stop("start_row < 1 at h=", h, " — not enough history. Reduce nprev or h.")
    
    # Training block (AR columns only)
    X_tr  <- Xall[start_row:end_row, ar_cols, drop = FALSE]
    y_tr  <- yall[start_row:end_row]
    
    # x_oos is taken at the issue-time anchor (NOT in training block)
    oos_row <- which(rownames(Xall) == origin_k)
    if (length(oos_row) != 1) stop("Anchor not found for origin_k at h=", h, " (", origin_k, ")")
    x_oos <- Xall[oos_row, ar_cols, drop = FALSE]
    
    # standardize+impute using training window; apply same to x_oos
    SI     <- standardize_impute(X_tr, x_oos)
    X_tr_s <- SI$X_tr; x_oos_s <- SI$x_oos
    
    # winsorize z-scores at |z|<=8 for pipeline uniformity
    cap <- 8
    X_tr_s[X_tr_s >  cap] <- cap;  X_tr_s[X_tr_s < -cap] <- -cap
    x_oos_s[x_oos_s >  cap] <- cap; x_oos_s[x_oos_s < -cap] <- -cap
    
    # Fit AR(p*) by OLS on standardized features; predict at standardized x_oos
    df_fit <- data.frame(y = y_tr, X_tr_s[, 1:p_star, drop = FALSE])
    fit <- lm(y ~ ., data = df_fit)
    DEBUG_AR <- TRUE
    if (DEBUG_AR && h == 1 && k == nprev) {
      cat("\n[AR DEBUG] h=", h, "k=last\n", sep = "")
      cat("p_star =", p_star, "\n")
      cat("Regressor columns used in TRAIN:\n")
      print(colnames(X_tr)[1:p_star])
      cat("Regressor columns used in OOS x:\n")
      print(colnames(x_oos)[1:p_star])
      
      # Assert there are NO PCs or X-lag columns in use
      stopifnot(!any(grepl("^PC\\d+$", colnames(X_tr)[1:p_star])))
      stopifnot(!any(grepl("_L[1-4]$",  colnames(X_tr)[1:p_star])))
      stopifnot(all(grepl("^y_lag[1-4]$", colnames(X_tr)[1:p_star])))
    }
    # prediction
    xnew <- as.data.frame(x_oos_s[, 1:p_star, drop = FALSE])
    preds[k] <- as.numeric(predict(fit, newdata = xnew))
  }
  
  names(preds) <- dates_oos
  y_real <- y_all[idx_oos]
  mets   <- compute_metrics(y_real, preds, y_all, idx_oos, h)
  
  AR_cpi_list[[h]] <- list(
    pred     = preds,
    errors   = mets,
    loss_abs = abs(y_real - preds),
    loss_sq  = (y_real - preds)^2
  )
}

# --------------- Save exactly like ridge/rw (matrix H=12 x nprev with rownames=dates_oos) ---------------
mat_ar <- do.call(cbind, lapply(1:12, function(h) AR_cpi_list[[h]]$pred))
colnames(mat_ar) <- paste0("h", 1:12)
rownames(mat_ar) <- dates_oos

write.table(mat_ar, "forecasts/ar/AR_cpi_total.csv",
            sep = ";", row.names = TRUE, col.names = TRUE)
save(AR_cpi_list, file = "forecasts/ar/AR-list-cpi_total.Rdata")

cat("Saved AR (paper-style plumbing, BIC p*, direct h-step, fixed rolling, no leakage) to forecasts/ar — nprev=",
    nprev, ", target=cpi_total\n", sep = "")



