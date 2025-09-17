# ---- ADAPTIVE LASSO (paper-style) with BIC λ, direct h-step, fixed rolling, nprev=99 ----
suppressPackageStartupMessages(library(glmnet))

# ---------------- Utilities (canonical from RIDGE_REFERENCE) ----------------
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
  y_lag_h  <- y_all[idx_oos - h]
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
  sd[!is.finite(sd) | sd < 1e-8] <- 1  # safer floor
  Xs <- sweep(X_tr, 2, mu, "-"); Xs <- sweep(Xs, 2, sd, "/")
  xo <- sweep(x_oos, 2, mu, "-"); xo <- sweep(xo, 2, sd, "/")
  Xs[!is.finite(Xs)] <- 0; xo[!is.finite(xo)] <- 0
  list(X_tr = as.matrix(Xs), x_oos = as.matrix(xo))
}

# For glmnet models with alpha in (0,1] we use glmnet's df (nonzeros) in BIC
select_lambda_bic_glmnet <- function(X, y, alpha = 1) {
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X)
  fit <- glmnet::glmnet(x = X, y = y,
                        alpha = alpha, family = "gaussian",
                        standardize = FALSE, intercept = TRUE)
  if (length(fit$lambda) == 0) stop("glmnet path is empty.")
  Yhat <- predict(fit, newx = X)
  rss  <- colSums((matrix(y, nrow=n, ncol=ncol(Yhat)) - Yhat)^2)
  df   <- fit$df                                   # nonzeros (excl intercept)
  eps  <- .Machine$double.eps
  bic  <- n * log(pmax(rss, eps) / n) + df * log(n)
  j    <- which.min(bic)
  list(lambda = fit$lambda[j], idx = j, fit_path = fit, bic = bic, df = df)
}

# Build paper-style design: y AR(4), X lags 1..4 for all non-targets, + 4 PCs (fixed)
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
  anchors <- seq.int(maxdrop, last_anchor)  # indices of t (issue-time anchors)
  X <- cbind(ylag, Xlags, PC)[anchors, , drop = FALSE]
  rownames(X) <- dates[anchors]
  
  y_target <- y[anchors + h]; names(y_target) <- dates[anchors + h]
  list(X = X, y = y_target)
}

# ---------------- Driver (data & OOS alignment per dprep/rwtext) ----------------
load("dados.rda")  # expects matrix with date rownames and 'cpi_total'
if (!("cpi_total" %in% colnames(dados))) stop("cpi_total not found in dados.")  # dprep contract
dates   <- rownames(dados)
y_all   <- as.numeric(dados[, "cpi_total"])
nprev   <- 99
idx_oos <- (nrow(dados) - nprev + 1):nrow(dados)
dates_oos <- dates[idx_oos]

horizons <- 1:12
alpha <- 1  # LASSO base for ADAPTIVE LASSO

dir.create("forecasts/adalasso", recursive = TRUE, showWarnings = FALSE)

ADALASSO_cpi_list <- vector("list", length(horizons))

for (h in horizons) {
  # --------- Anchors & fixed-R window (canonical) ----------
  origin_first <- dates[idx_oos[1] - h]                 # issue-time anchor for first OOS
  train_rows   <- which(dates <= origin_first)
  
  # PCA on initial train window (X_t of non-targets), FREEZE loadings
  other_vars <- setdiff(colnames(dados), "cpi_total")
  X_train0   <- scale(dados[train_rows, other_vars, drop = FALSE], center = TRUE, scale = TRUE)
  pc         <- prcomp(X_train0, center = FALSE, scale. = FALSE)
  pc_fit     <- list(rotation = pc$rotation,
                     center   = attr(X_train0, "scaled:center"),
                     scale    = attr(X_train0, "scaled:scale"))
  
  sup  <- build_xy_paper(dados, "cpi_total", h = h, p_y = 4, x_lags = 4, pc_fit = pc_fit)
  Xall <- sup$X; yall <- sup$y
  
  # Training window length R fixed from the FIRST origin's last training anchor:
  last_train_anchor_first <- dates[idx_oos[1] - 2*h]
  if (!(last_train_anchor_first %in% rownames(Xall))) {
    stop("Not enough history for h=", h, ": missing last_train_anchor_first = ", last_train_anchor_first)
  }
  keep <- rownames(Xall) <= last_train_anchor_first
  X0   <- Xall[keep, , drop = FALSE]; y0 <- yall[keep]
  R    <- nrow(X0)
  if (R <= 0) stop("Training window length R<=0 at h=", h, ". Check sample size.")
  
  preds_adalasso <- numeric(nprev)
  
  for (k in 1:nprev) {
    # For each OOS k, end training at last_train_anchor_k = dates[idx_oos[k] - 2*h]
    origin_k <- dates[idx_oos[k] - h]               # issue-time anchor t
    last_train_anchor_k <- dates[idx_oos[k] - 2*h]  # last training anchor s
    
    end_row <- which(rownames(Xall) == last_train_anchor_k)
    if (length(end_row) != 1) stop("Anchor not found for last_train_anchor_k at h=", h, " (", last_train_anchor_k, ")")
    start_row <- end_row - R + 1
    if (start_row < 1) stop("start_row < 1 at h=", h, " — not enough history. Reduce nprev or h.")
    
    X_tr  <- Xall[start_row:end_row, , drop = FALSE]
    y_tr  <- yall[start_row:end_row]
    
    # x_oos must be taken at the issue-time anchor (NOT in training)
    oos_row <- which(rownames(Xall) == origin_k)
    if (length(oos_row) != 1) stop("Anchor not found for origin_k at h=", h, " (", origin_k, ")")
    x_oos <- Xall[oos_row, , drop = FALSE]
    
    # Standardize on training window + impute; apply same to x_oos
    SI     <- standardize_impute(X_tr, x_oos)
    X_tr_s <- SI$X_tr; x_oos_s <- SI$x_oos
    
    # Optional winsorization at |z| ≤ 8 to keep pipeline uniform with ridge
    cap <- 8
    X_tr_s[X_tr_s >  cap] <- cap;  X_tr_s[X_tr_s < -cap] <- -cap
    x_oos_s[x_oos_s >  cap] <- cap; x_oos_s[x_oos_s < -cap] <- -cap
    
    # ---- Stage 1: LASSO with BIC λ (glmnet df) ----
    sel1 <- select_lambda_bic_glmnet(X_tr_s, y_tr, alpha = alpha)  # standardize=FALSE inside
    fit1 <- sel1$fit_path
    beta1 <- as.numeric(coef(fit1, s = sel1$lambda))  # includes intercept
    b1 <- beta1[-1]
    
    # ---- Stage 2: Adaptive LASSO (weights = (|b_j| + 1/sqrt(n))^{-1}) ----
    n_tr <- length(y_tr)
    w <- (abs(b1) + 1/sqrt(n_tr))^(-1)
    fit2 <- glmnet::glmnet(x = X_tr_s, y = y_tr,
                           alpha = alpha, family = "gaussian",
                           penalty.factor = w,
                           standardize = FALSE, intercept = TRUE)
    # BIC selection along weighted path (glmnet df)
    Yhat2 <- predict(fit2, newx = X_tr_s)
    rss2  <- colSums((matrix(y_tr, nrow=n_tr, ncol=ncol(Yhat2)) - Yhat2)^2)
    df2   <- fit2$df
    eps   <- .Machine$double.eps
    bic2  <- n_tr * log(pmax(rss2, eps) / n_tr) + df2 * log(n_tr)
    j2    <- which.min(bic2)
    lambda2 <- fit2$lambda[j2]
    
    preds_adalasso[k] <- as.numeric(predict(fit2, newx = x_oos_s, s = lambda2))
  }
  
  names(preds_adalasso) <- dates_oos
  
  y_real <- y_all[idx_oos]
  mets_adalasso <- compute_metrics(y_real, preds_adalasso, y_all, idx_oos, h)
  
  ADALASSO_cpi_list[[h]] <- list(
    pred     = preds_adalasso,
    errors   = mets_adalasso,
    loss_abs = abs(y_real - preds_adalasso),
    loss_sq  = (y_real - preds_adalasso)^2
  )
}

# --------------- Save exactly like ridge/rw (matrix H=12 x nprev) ---------------
mat_adalasso <- do.call(cbind, lapply(1:12, function(h) ADALASSO_cpi_list[[h]]$pred))
colnames(mat_adalasso) <- paste0("h", 1:12)
rownames(mat_adalasso) <- dates_oos

write.table(mat_adalasso, "forecasts/adalasso/ADALASSO_cpi_total.csv",
            sep = ";", row.names = TRUE, col.names = TRUE)

save(ADALASSO_cpi_list, file = "forecasts/adalasso/ADALASSO-list-cpi_total.Rdata")

cat("Saved ADALASSO (two-stage, BIC λ) to forecasts/adalasso — nprev=",
    nprev, ", target=cpi_total\n", sep = "")
