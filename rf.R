# ---- RANDOM FOREST (paper-style) direct h-step, fixed rolling, nprev=99 ----
suppressPackageStartupMessages({
  library(randomForest)
  library(parallel)
})

set.seed(123)

# ---------- Utilities ----------
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

# Train-window standardize + mean-impute; apply to x_oos
standardize_impute <- function(X_tr, x_oos) {
  mu <- apply(X_tr, 2, function(c) mean(c, na.rm = TRUE))
  sd <- apply(X_tr, 2, function(c) stats::sd(c, na.rm = TRUE))
  sd[!is.finite(sd) | sd == 0] <- 1
  Xs <- sweep(X_tr, 2, mu, "-"); Xs <- sweep(Xs, 2, sd, "/")
  xo <- sweep(x_oos, 2, mu, "-"); xo <- sweep(xo, 2, sd, "/")
  Xs[!is.finite(Xs)] <- 0; xo[!is.finite(xo)] <- 0
  list(X_tr = as.matrix(Xs), x_oos = as.matrix(xo))
}

# Paper-style design: 4 AR of y, 4 lags of all X, + 4 PCs (frozen)
build_xy_paper <- function(dados, yname, h, p_y=4, x_lags=4, pc_fit) {
  stopifnot(yname %in% colnames(dados))
  dates <- rownames(dados)
  y <- as.numeric(dados[, yname]); names(y) <- dates
  Xraw <- dados[, setdiff(colnames(dados), yname), drop = FALSE]
  
  # y AR lags
  ylag <- sapply(1:p_y, function(L) c(rep(NA, L), y[1:(length(y)-L)]))
  colnames(ylag) <- paste0("y_lag", 1:p_y)
  
  # X lags
  make_lag <- function(v, L) c(rep(NA, L), v[1:(length(v)-L)])
  Xlags <- do.call(cbind, lapply(1:x_lags, function(L) {
    xs <- apply(Xraw, 2, make_lag, L = L)
    colnames(xs) <- paste0(colnames(Xraw), "_L", L)
    xs
  }))
  
  # 4 PCs from contemporaneous X_t using fixed loadings from initial train
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

# ---------- Driver ----------
load("dados.rda")
if (!("cpi_total" %in% colnames(dados))) stop("cpi_total not found in dados.")

dates   <- rownames(dados)
y_all   <- as.numeric(dados[, "cpi_total"])
nprev   <- 99
idx_oos <- (nrow(dados) - nprev + 1):nrow(dados)
dates_oos <- dates[idx_oos]
horizons <- 1:12

# RF hyperparams
ntree_rf    <- 1000
nodesize_rf <- 5
mtry_frac   <- 1/3
ncores      <- max(1, min(parallel::detectCores()-1, 6))

dir.create("forecasts/rf", recursive = TRUE, showWarnings = FALSE)
RF_cpi_list <- vector("list", length(horizons))

for (h in horizons) {
  origin_first <- dates[idx_oos[1] - h]                 # issue-time anchor
  train_rows   <- which(dates <= origin_first)
  
  # PCs on initial train window (X_t of non-targets), freeze loadings
  other_vars <- setdiff(colnames(dados), "cpi_total")
  X_train0   <- scale(dados[train_rows, other_vars, drop = FALSE], center = TRUE, scale = TRUE)
  pc         <- prcomp(X_train0, center = FALSE, scale. = FALSE)
  pc_fit     <- list(rotation = pc$rotation,
                     center   = attr(X_train0, "scaled:center"),
                     scale    = attr(X_train0, "scaled:scale"))
  
  sup  <- build_xy_paper(dados, "cpi_total", h = h, p_y = 4, x_lags = 4, pc_fit = pc_fit)
  Xall <- sup$X; yall <- sup$y
  
  # >>> FIX: initial training window ends at (origin_first - h)
  last_train_anchor_first <- dates[idx_oos[1] - 2*h]
  if (!(last_train_anchor_first %in% rownames(Xall))) {
    stop("Not enough history for h=", h, " (missing last_train_anchor_first).")
  }
  keep <- rownames(Xall) <= last_train_anchor_first
  X0   <- Xall[keep, , drop = FALSE]; y0 <- yall[keep]
  R    <- nrow(X0)
  if (R <= 0) stop("R <= 0 at h=", h)
  
  rf_one <- function(k) {
    origin_k <- dates[idx_oos[k] - h]              # issue-time anchor t
    last_train_anchor_k <- dates[idx_oos[k] - 2*h] # last training anchor t-h
    
    end_row <- which(rownames(Xall) == last_train_anchor_k)
    if (length(end_row) != 1) stop("No row for last_train_anchor_k at h=", h)
    start_row <- end_row - R + 1
    if (start_row < 1) stop("start_row < 1 at h=", h)
    
    X_tr_raw  <- Xall[start_row:end_row, , drop = FALSE]
    y_tr      <- yall[start_row:end_row]
    
    # OOS features at issue time (not in training)
    oos_row <- which(rownames(Xall) == origin_k)
    if (length(oos_row) != 1) stop("No row for origin_k at h=", h)
    x_oos_raw <- Xall[oos_row, , drop = FALSE]
    
    # Standardize + impute (optional for trees, but keeps pipeline uniform)
    SI      <- standardize_impute(X_tr_raw, x_oos_raw)
    X_tr_s  <- SI$X_tr
    x_oos_s <- SI$x_oos
    
    p <- ncol(X_tr_s)
    mtry_val <- max(1L, floor(p * mtry_frac))
    
    fit <- randomForest(
      x = as.data.frame(X_tr_s),
      y = as.numeric(y_tr),
      xtest = as.data.frame(x_oos_s),
      ntree = ntree_rf,
      mtry  = mtry_val,
      nodesize = nodesize_rf,
      keep.forest = TRUE,
      importance = FALSE
    )
    as.numeric(fit$test$predicted)
  }
  
  # Parallel or serial loop
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("Xall","yall","dates","idx_oos","R","h",
                  "standardize_impute","mtry_frac",
                  "ntree_rf","nodesize_rf","randomForest"),
      envir = environment()
    )
    preds <- parallel::parSapply(cl, 1:nprev, rf_one)
  } else {
    preds <- sapply(1:nprev, rf_one)
  }
  
  names(preds) <- dates_oos
  y_real <- y_all[idx_oos]
  mets   <- compute_metrics(y_real, preds, y_all, idx_oos, h)
  
  RF_cpi_list[[h]] <- list(
    pred     = preds,
    errors   = mets,
    loss_abs = abs(y_real - preds),
    loss_sq  = (y_real - preds)^2
  )
}

# --------------- Save ----------------
mat_rf <- do.call(cbind, lapply(1:12, function(h) RF_cpi_list[[h]]$pred))
colnames(mat_rf) <- paste0("h", 1:12)
rownames(mat_rf) <- dates_oos

write.table(mat_rf, "forecasts/rf/RF_cpi_total.csv",
            sep = ";", row.names = TRUE, col.names = TRUE)
save(RF_cpi_list, file = "forecasts/rf/RF-list-cpi_total.Rdata")

cat("Saved RF (paper-style, ntree=", ntree_rf, ", mtry≈p*", mtry_frac,
    ") to forecasts/rf — nprev=", nprev, ", target=cpi_total\n", sep = "")
