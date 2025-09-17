# ---- RANDOM WALK (paper-style) direct h-step, fixed rolling, nprev=99 ----
suppressPackageStartupMessages(library(forecast))  # keep RW estimator parity with legacy

# ---------------- Utilities (identical to RIDGE reference) ----------------
rmse <- function(e) sqrt(mean(e^2)); mae <- function(e) mean(abs(e))
mad_ <- function(e) median(abs(e - median(e)))

compute_metrics <- function(y_real, yhat, y_all, idx_oos, h) {
  err <- y_real - yhat
  ae  <- abs(err); se <- err^2
  m_rmse <- sqrt(mean(se)); m_mae <- mean(ae)
  m_mad <- median(abs(err - median(err)))
  m_meanad <- mean(abs(err - mean(err)))
  rng <- diff(range(y_real)); m_nrmse <- m_rmse / ifelse(rng == 0, 1, rng)
  
  # h-step RW baseline — same convention across models
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

# (Extracted for parity with the ridge pipeline; unused for RW)
standardize_impute <- function(X_tr, x_oos) {
  mu <- apply(X_tr, 2, function(c) mean(c, na.rm = TRUE))
  sd <- apply(X_tr, 2, function(c) stats::sd(c, na.rm = TRUE))
  sd[!is.finite(sd) | sd < 1e-8] <- 1
  Xs <- sweep(X_tr, 2, mu, "-"); Xs <- sweep(Xs, 2, sd, "/")
  xo <- sweep(x_oos, 2, mu, "-"); xo <- sweep(x_oos, 2, sd, "/")
  Xs[!is.finite(Xs)] <- 0; xo[!is.finite(xo)] <- 0
  list(X_tr = as.matrix(Xs), x_oos = as.matrix(xo))
}

# ---------------- Driver ----------------
# expects matrix `dados` with date rownames; produced by your dprep
load("dados.rda")
if (!("cpi_total" %in% colnames(dados))) stop("cpi_total not found in dados.")

dates    <- rownames(dados)
y_all    <- as.numeric(dados[, "cpi_total"])
nprev    <- 99
idx_oos  <- (nrow(dados) - nprev + 1):nrow(dados)   # target indices (t)
dates_oos <- dates[idx_oos]
horizons <- 1:12

dir.create("forecasts/rw", recursive = TRUE, showWarnings = FALSE)

RW_cpi_list <- vector("list", length(horizons))

# Small helper: RW h-step using data up to an origin (issue-time) anchor
runRW <- function(y_up_to_origin, h) {
  # forecast::rwf returns the last observed value as the mean path (no drift)
  out <- forecast::rwf(y_up_to_origin, h = h, drift = FALSE, fan = FALSE)
  as.numeric(out$mean[h])
}

for (h in horizons) {
  # ---- Anchors & fixed rolling length R (ridge methodology) ----
  # First OOS origin (issue-time) t = dates[idx_oos[1] - h]
  origin_first <- dates[idx_oos[1] - h]
  # Last training anchor for THIS horizon (first origin): t - h
  last_train_anchor_first <- dates[idx_oos[1] - 2*h]
  
  # Sanity checks
  if ((idx_oos[1] - h) < 1) {
    stop("Not enough history for h=", h, ": origin_first would be before start of sample.")
  }
  if ((idx_oos[1] - 2*h) < 1) {
    stop("Not enough history for h=", h, ": last_train_anchor_first would be before start of sample.")
  }
  
  # R is fixed from the first origin: number of observations up to last_train_anchor_first
  last_train_anchor_first_idx <- which(dates == last_train_anchor_first)
  if (length(last_train_anchor_first_idx) != 1) {
    stop("Missing last_train_anchor_first (", last_train_anchor_first, ") for h=", h)
  }
  R <- last_train_anchor_first_idx
  if (R <= 0) stop("Training window length R<=0 at h=", h)
  
  preds <- numeric(nprev)
  
  for (k in 1:nprev) {
    # Issue-time anchor for k-th OOS target
    origin_k <- dates[idx_oos[k] - h]
    # Last training anchor for k-th origin
    last_train_anchor_k <- dates[idx_oos[k] - 2*h]
    
    # Locate indices
    origin_k_idx <- which(dates == origin_k)
    last_train_anchor_k_idx <- which(dates == last_train_anchor_k)
    
    if (length(origin_k_idx) != 1) {
      stop("Anchor not found for origin_k at h=", h, " (", origin_k, ")")
    }
    if (length(last_train_anchor_k_idx) != 1) {
      stop("Anchor not found for last_train_anchor_k at h=", h, " (", last_train_anchor_k, ")")
    }
    
    # Fixed rolling: start_row ensures a constant window length R
    start_row <- last_train_anchor_k_idx - R + 1
    if (start_row < 1) {
      stop("start_row < 1 at h=", h, " — not enough history. Reduce nprev or h.")
    }
    
    # ---- Fit RW using ONLY data up to the origin (no leakage) ----
    y_train_up_to_origin <- y_all[1:origin_k_idx]
    preds[k] <- runRW(y_train_up_to_origin, h = h)
    
    # (Note: For RW, this equals y_all[origin_k_idx] but we keep runRW for estimator parity)
  }
  
  names(preds) <- dates_oos
  y_real <- y_all[idx_oos]
  
  # Metrics (identical to ridge)
  mets <- compute_metrics(y_real, preds, y_all, idx_oos, h)
  
  RW_cpi_list[[h]] <- list(
    pred     = preds,
    errors   = mets,
    loss_abs = abs(y_real - preds),
    loss_sq  = (y_real - preds)^2
  )
  
  # Optional debug (mirror ridge’s one-shot print)
  if (h == 1) {
    cat("\n[h=1] origin_first =", origin_first,
        " last_train_anchor_first =", last_train_anchor_first,
        "  R =", R, "\n")
  }
}

# --------------- Save (identical format to ridge) ----------------
mat_rw <- do.call(cbind, lapply(1:12, function(h) RW_cpi_list[[h]]$pred))
colnames(mat_rw) <- paste0("h", 1:12)
rownames(mat_rw) <- dates_oos

write.table(mat_rw, "forecasts/rw/RW_cpi_total.csv",
            sep = ";", row.names = TRUE, col.names = TRUE)
save(RW_cpi_list, file = "forecasts/rw/RW-list-cpi_total.Rdata")

cat("Saved RW (direct h-step, fixed rolling anchors) to forecasts/rw — nprev=",
    nprev, ", target=cpi_total\n", sep = "")
