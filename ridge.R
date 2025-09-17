# ---- RIDGE (paper-style) with BIC λ, direct h-step, fixed rolling, nprev=99 ----
suppressPackageStartupMessages(library(glmnet))

# ---------------- Utilities (aligned with your RW/other shrinkage scripts) ----------------
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
  sd[!is.finite(sd) | sd < 1e-8] <- 1  # was '== 0'; this is safer
  Xs <- sweep(X_tr, 2, mu, "-"); Xs <- sweep(Xs, 2, sd, "/")
  xo <- sweep(x_oos, 2, mu, "-"); xo <- sweep(x_oos, 2, sd, "/")
  Xs[!is.finite(Xs)] <- 0; xo[!is.finite(xo)] <- 0
  list(X_tr = as.matrix(Xs), x_oos = as.matrix(xo))
}



select_lambda_bic_ridge <- function(X, y) {
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X)
  
  fit <- glmnet::glmnet(x = X, y = y,
                        alpha = 0, family = "gaussian",
                        standardize = FALSE, intercept = TRUE)
  
  # RSS along the path (train fit)
  Yhat <- predict(fit, newx = X)
  rss  <- colSums((matrix(y, nrow=n, ncol=ncol(Yhat)) - Yhat)^2)
  
  # SVD of training design
  s <- svd(X, nu = 0, nv = 0)$d
  # ---- KEY FIX: n * lambda in the denominator ----
  df_eff <- sapply(fit$lambda, function(l) sum((s^2) / (s^2 + n*l)))
  
  eps <- .Machine$double.eps
  bic <- n * log(pmax(rss, eps) / n) + df_eff * log(n)
  j   <- which.min(bic)
  
  list(lambda = fit$lambda[j], idx = j, fit_path = fit,
       bic = bic, df = df_eff)
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
  anchors <- maxdrop:(nrow(dados) - h)  # indices of t
  
  # >>> RE-ENABLE PCs here <<<
  X <- cbind(ylag, Xlags, PC)[anchors, , drop = FALSE]
  
  rownames(X) <- dates[anchors]
  y_target <- y[anchors + h]; names(y_target) <- dates[anchors + h]
  list(X = X, y = y_target)
}

# ---------------- Driver ----------------
load("dados.rda")  # expects matrix with date rownames and 'cpi_total'
if (!("cpi_total" %in% colnames(dados))) stop("cpi_total not found in dados.")

dates   <- rownames(dados)
y_all   <- as.numeric(dados[, "cpi_total"])
nprev   <- 99
idx_oos <- (nrow(dados) - nprev + 1):nrow(dados)
dates_oos <- dates[idx_oos]

horizons  <- 1:12
alpha_rdg <- 0  # ridge

dir.create("forecasts/ridge", recursive = TRUE, showWarnings = FALSE)

RIDGE_cpi_list <- vector("list", length(horizons))

for (h in horizons) {
  # --------- PATCH 1: initial window must end at (origin_first - h) in anchor time ----------
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
  
  # last training anchor for the FIRST origin is (origin_first - h)
  last_train_anchor_first <- dates[idx_oos[1] - 2*h]
  if (!(last_train_anchor_first %in% rownames(Xall))) {
    stop("Not enough history for h=", h, ": missing last_train_anchor_first = ", last_train_anchor_first)
  }
  keep <- rownames(Xall) <= last_train_anchor_first
  X0   <- Xall[keep, , drop = FALSE]; y0 <- yall[keep]
  R    <- nrow(X0)
  if (R <= 0) stop("Training window length R<=0 at h=", h, ". Check sample size.")
  
  preds_ridge <- numeric(nprev)
  
  for (k in 1:nprev) {
    # --------- PATCH 2: for each OOS k, end training at (origin_k - h) ----------
    origin_k <- dates[idx_oos[k] - h]               # issue-time anchor t = T - h
    last_train_anchor_k <- dates[idx_oos[k] - 2*h]  # last training anchor s = T - 2h
    
    end_row <- which(rownames(Xall) == last_train_anchor_k)
    if (length(end_row) != 1) stop("Anchor not found for last_train_anchor_k at h=", h, " (", last_train_anchor_k, ")")
    start_row <- end_row - R + 1
    if (start_row < 1) stop("start_row < 1 at h=", h, " — not enough history. Reduce nprev or h.")
    
    X_tr  <- Xall[start_row:end_row, , drop = FALSE]
    y_tr  <- yall[start_row:end_row]
    
    # --------- PATCH 3: x_oos taken at the issue-time anchor (not in training) ----------
    oos_row <- which(rownames(Xall) == origin_k)
    if (length(oos_row) != 1) stop("Anchor not found for origin_k at h=", h, " (", origin_k, ")")
    x_oos <- Xall[oos_row, , drop = FALSE]
    
    # standardize on training window + impute; apply same to x_oos
    SI     <- standardize_impute(X_tr, x_oos)
    X_tr_s <- SI$X_tr; x_oos_s <- SI$x_oos
    
    cap <- 8
    X_tr_s[X_tr_s >  cap] <- cap;  X_tr_s[X_tr_s < -cap] <- -cap
    x_oos_s[x_oos_s >  cap] <- cap; x_oos_s[x_oos_s < -cap] <- -cap
    
    # Ridge with BIC λ
    sel <- select_lambda_bic_ridge(X_tr_s, y_tr)
    fit <- sel$fit_path
    preds_ridge[k] <- as.numeric(predict(fit, newx = x_oos_s, s = sel$lambda))
    # --- DEBUG: inspect one forecast ---
    if (h == 1 && k == nprev) {
      # 1) Report anchors
      cat("\n[h=1, k=last] origin_k =", origin_k,
          " last_train_anchor_k =", last_train_anchor_k, "\n")
      
      # 2) Selected lambda & effective df (ridge)
      sel_dbg <- sel
      cat("lambda* =", sel_dbg$lambda, "  df_eff* ≈", sel_dbg$df[sel_dbg$idx], "\n")
      
      # 3) How many NAs in raw x_oos before SI? (should be 0 ideally)
      na_raw <- colSums(is.na(x_oos))
      cat("x_oos raw NA count =", sum(na_raw), "\n")
      
      # 4) Any extreme z-scores in x_oos_s? (|z|>8 is suspicious)
      z <- as.numeric(x_oos_s)
      big <- which(abs(z) > 8)
      if (length(big)) {
        cat("Extreme |z| in x_oos_s:", length(big), "features. Top 10:\n")
        ord <- order(abs(z), decreasing=TRUE)
        print(head(data.frame(var = colnames(X_tr)[ord], z = z[ord]), 10))
      } else {
        cat("No extreme |z| > 8 in x_oos_s.\n")
      }
      
      # 5) Top coefficient contributions to the prediction
      b <- as.numeric(coef(fit, s = sel_dbg$lambda))  # [1] intercept; [-1] slopes
      contrib <- z * b[-1]
      ordc <- order(abs(contrib), decreasing=TRUE)
      cat("Top 10 contribs to yhat:\n")
      print(head(data.frame(var = colnames(X_tr)[ordc],
                            z = z[ordc], beta = b[-1][ordc],
                            contrib = contrib[ordc]), 10))
    }
    
  }
  
  names(preds_ridge) <- dates_oos
  y_real <- y_all[idx_oos]
  mets   <- compute_metrics(y_real, preds_ridge, y_all, idx_oos, h)
  
  RIDGE_cpi_list[[h]] <- list(
    pred     = preds_ridge,
    errors   = mets,
    loss_abs = abs(y_real - preds_ridge),
    loss_sq  = (y_real - preds_ridge)^2
  )
}

# --------------- Save like RW/LASSO/ELNET/ADALASSO ----------------
mat_ridge <- do.call(cbind, lapply(1:12, function(h) RIDGE_cpi_list[[h]]$pred))
colnames(mat_ridge) <- paste0("h", 1:12)
rownames(mat_ridge) <- dates_oos

write.table(mat_ridge, "forecasts/ridge/RIDGE_cpi_total.csv",
            sep = ";", row.names = TRUE, col.names = TRUE)
save(RIDGE_cpi_list, file = "forecasts/ridge/RIDGE-list-cpi_total.Rdata")

cat("Saved RIDGE (alpha=0, BIC λ) to forecasts/ridge — nprev=",
    nprev, ", target=cpi_total\n", sep = "")
