# ==== RLASSO (hdm) — paper-style pipeline =====================================
# Direct h-step; fixed rolling; anchors (NO LEAKAGE); frozen PCs; nprev = 99
# Uses canonical: compute_metrics(), standardize_impute(), build_xy_paper()
# Saves: forecasts/rlasso/RLASSO_cpi_total.csv and RLASSO-list-cpi_total.Rdata
# ==============================================================================

suppressPackageStartupMessages({
  library(hdm)        # rlasso
  library(foreach)
  library(doParallel)
  library(doRNG)
})

# ---------------- Utilities (copied from RIDGE_REFERENCE) ----------------------
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
  
  X <- cbind(ylag, Xlags, PC)[anchors, , drop = FALSE]
  rownames(X) <- dates[anchors]
  y_target <- y[anchors + h]; names(y_target) <- dates[anchors + h]
  list(X = X, y = y_target)
}

# ---------------- Driver -------------------------------------------------------
model_folder <- "rlasso"  # short name for output folder & filenames

# Load your cleaned matrix with date rownames (see dprep.txt) -------------------
load("dados.rda")  # provides `dados` (matrix with date rownames)
if (!is.matrix(dados) || is.null(rownames(dados))) {
  stop("Expected 'dados' as a matrix with date rownames from dprep.")
}
if (!("cpi_total" %in% colnames(dados))) {
  stop("Column 'cpi_total' not found in `dados`.")
}

dates    <- rownames(dados)
y_all    <- as.numeric(dados[, "cpi_total"])
nprev    <- 99
idx_oos  <- (nrow(dados) - nprev + 1):nrow(dados)
dates_oos <- dates[idx_oos]
horizons <- 1:12

# Output dir (same style as RW/RIDGE) ------------------------------------------
dir.create(file.path("forecasts", model_folder), recursive = TRUE, showWarnings = FALSE)

# Optional parallelism
use_parallel <- TRUE
max_workers  <- max(1, parallel::detectCores(logical = TRUE) - 1)
if (use_parallel) {
  cl <- parallel::makeCluster(max_workers)
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(123)
}

RLASSO_cpi_list <- vector("list", length(horizons))

for (h in horizons) {
  
  # ---- Anchors & initial PC fit (NO LEAKAGE) ---------------------------------
  origin_first <- dates[idx_oos[1] - h]           # issue-time anchor for first OOS
  last_train_anchor_first <- dates[idx_oos[1] - 2*h]
  
  # Freeze PCA from the *initial training window* ending at last_train_anchor_first
  other_vars <- setdiff(colnames(dados), "cpi_total")
  if (!(last_train_anchor_first %in% dates)) {
    stop("Not enough history for h=", h, ": missing last_train_anchor_first = ", last_train_anchor_first)
  }
  pc_rows <- which(dates <= last_train_anchor_first)
  if (!length(pc_rows)) stop("Empty PC window at h=", h)
  X_train0   <- scale(dados[pc_rows, other_vars, drop = FALSE], center = TRUE, scale = TRUE)
  pc         <- prcomp(X_train0, center = FALSE, scale. = FALSE)
  pc_fit     <- list(rotation = pc$rotation,
                     center   = attr(X_train0, "scaled:center"),
                     scale    = attr(X_train0, "scaled:scale"))
  
  # Build paper-style design (AR4, 4 lags for all X, +4 PCs with frozen loadings)
  sup  <- build_xy_paper(dados, "cpi_total", h = h, p_y = 4, x_lags = 4, pc_fit = pc_fit)
  Xall <- sup$X
  yall <- sup$y
  
  # Set fixed window length R from FIRST origin (training ends at last_train_anchor_first)
  end_row_first <- which(rownames(Xall) == last_train_anchor_first)
  if (length(end_row_first) != 1) {
    stop("Anchor not found in design for last_train_anchor_first at h=", h, " (", last_train_anchor_first, ")")
  }
  X0 <- Xall[1:end_row_first, , drop = FALSE]
  y0 <- yall[1:end_row_first]
  R  <- nrow(X0)
  if (R <= 0) stop("Training window length R<=0 at h=", h, ". Check sample size.")
  
  # Storage for predictions
  preds <- numeric(nprev)
  
  # ---- Rolling evaluation for nprev OOS origins ------------------------------
  Kidx <- 1:nprev
  
  # Optionally run in parallel over k
  pred_list <- foreach(k = Kidx, .export = c("standardize_impute"),
                       .combine = "c") %dorng% {
                         
                         origin_k <- dates[idx_oos[k] - h]               # predict using x_t at t = T - h
                         last_train_anchor_k <- dates[idx_oos[k] - 2*h]  # training ends at t - h
                         
                         end_row <- which(rownames(Xall) == last_train_anchor_k)
                         if (length(end_row) != 1) {
                           stop("Anchor not found for last_train_anchor_k at h=", h, " (", last_train_anchor_k, ")")
                         }
                         start_row <- end_row - R + 1
                         if (start_row < 1) stop("start_row < 1 at h=", h, " — not enough history. Reduce nprev or h.")
                         
                         X_tr <- Xall[start_row:end_row, , drop = FALSE]
                         y_tr <- yall[start_row:end_row]
                         
                         # x_oos taken at issue-time anchor (not in training)
                         oos_row <- which(rownames(Xall) == origin_k)
                         if (length(oos_row) != 1) stop("Anchor not found for origin_k at h=", h, " (", origin_k, ")")
                         x_oos <- Xall[oos_row, , drop = FALSE]
                         
                         # Standardize+impute on training, apply to x_oos (canonical)
                         SI     <- standardize_impute(X_tr, x_oos)
                         X_tr_s <- SI$X_tr; x_oos_s <- SI$x_oos
                         
                         # Winsorize |z| at 8 (canonical)
                         cap <- 8
                         X_tr_s[X_tr_s >  cap] <- cap;  X_tr_s[X_tr_s < -cap] <- -cap
                         x_oos_s[x_oos_s >  cap] <- cap; x_oos_s[x_oos_s < -cap] <- -cap
                         
                         # Guard degenerate columns (near-constant). Keep those with sample sd > 0.
                         # (After standardization, constant cols become zeros.)
                         keep <- which(apply(X_tr_s, 2, function(v) sd(as.numeric(v)) > 0))
                         if (length(keep) == 0) {
                           # fallback: intercept-only prediction = mean(y_tr)
                           yhat_k <- mean(y_tr)
                         } else {
                           X_tr_s2 <- X_tr_s[, keep, drop = FALSE]
                           x_oos_s2 <- x_oos_s[, keep, drop = FALSE]
                           
                           # ---- NEW_MODEL core: square-root Lasso (hdm::rlasso) -------------------
                           # post = FALSE (pure penalized estimator); intercept included by default.
                           fit <- hdm::rlasso(x = X_tr_s2, y = y_tr, post = FALSE, intercept = TRUE)
                           yhat_k <- as.numeric(predict(fit, newdata = x_oos_s2))
                         }
                         
                         # DEBUG (only once): h=1,k=last
                         if (h == 1 && k == nprev) {
                           cat("\n[h=1, k=last] origin_k =", origin_k,
                               " last_train_anchor_k =", last_train_anchor_k, "\n")
                           na_raw <- colSums(is.na(x_oos))
                           cat("x_oos raw NA count =", sum(na_raw), "\n")
                           z <- as.numeric(x_oos_s)
                           big <- which(abs(z) > 8)
                           if (length(big)) {
                             cat("Extreme |z| in x_oos_s:", length(big), "features. Top 10 by |z|:\n")
                             ord <- order(abs(z), decreasing=TRUE)
                             print(head(data.frame(var = colnames(X_tr)[ord], z = z[ord]), 10))
                           } else {
                             cat("No extreme |z| > 8 in x_oos_s.\n")
                           }
                           # If we have a model, show top abs(beta*z) contributions on kept cols
                           if (exists("fit") && length(keep) > 0) {
                             b <- as.numeric(coef(fit))  # includes intercept as first
                             # Map coefficients to kept columns (skip intercept)
                             beta <- b[-1]
                             contrib <- as.numeric(x_oos_s2) * beta
                             ordc <- order(abs(contrib), decreasing = TRUE)
                             cat("Top 10 contribs to yhat (|z*beta|):\n")
                             print(head(data.frame(var = colnames(X_tr_s2)[ordc],
                                                   z = as.numeric(x_oos_s2)[ordc],
                                                   beta = beta[ordc],
                                                   contrib = contrib[ordc]), 10))
                           }
                         }
                         
                         yhat_k
                       } # foreach k
  
  preds[] <- pred_list
  
  names(preds) <- dates_oos
  y_real <- y_all[idx_oos]
  mets   <- compute_metrics(y_real, preds, y_all, idx_oos, h)
  
  RLASSO_cpi_list[[h]] <- list(
    pred     = preds,
    errors   = mets,
    loss_abs = abs(y_real - preds),
    loss_sq  = (y_real - preds)^2
  )
}

if (use_parallel) parallel::stopCluster(cl)

# --------------- Save (same format as RW/RIDGE) -------------------------------
mat_pred <- do.call(cbind, lapply(1:12, function(h) RLASSO_cpi_list[[h]]$pred))
colnames(mat_pred) <- paste0("h", 1:12)
rownames(mat_pred) <- dates_oos

write.table(mat_pred,
            file = file.path("forecasts", model_folder, "RLASSO_cpi_total.csv"),
            sep = ";", row.names = TRUE, col.names = TRUE)

save(RLASSO_cpi_list,
     file = file.path("forecasts", model_folder, "RLASSO-list-cpi_total.Rdata"))

cat("Saved RLASSO (hdm) to forecasts/", model_folder,
    " — nprev=", nprev, ", target=cpi_total\n", sep = "")
