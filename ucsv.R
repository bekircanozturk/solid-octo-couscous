# ---- UCSV (paper-style) with train-to-origin, NO LEAKAGE, nprev=99 ----
# Parallelized across forecast origins (k) with parLapplyLB.
# - Monthly horizons: h = 1..12, forecast at t+h equals last tauhat at origin t
# - Accumulated horizons: h ∈ {3,6,12}, dependent = trailing h-month sum of CPI
# - Outputs (same files, both variants combined):
#     forecasts/ucsv/UCSV_cpi_total.csv
#     forecasts/ucsv/UCSV-list-cpi_total.Rdata

suppressPackageStartupMessages({
  library(parallel)
  library(stats)
})

# ============================== UCSV CORE (unchanged) ==============================
ucsv <- function(y, display = FALSE){
  T  <- length(y)
  nloop <- 4000; burnin <- 1000
  Vtau <- 0.12; Vh <- 0.12
  atau <- 10;   ltau <- 0.04*9
  ah   <- 10;   lh   <- 0.03*9
  
  omega2tau <- ltau/(atau-1)
  omega2h   <- lh/(ah-1)
  h <- as.numeric(log(var(y)*.8)) * matrix(rep(1,T))
  Sp <- diag(T); d <- matrix(rep(0,T),1); sparse <- rbind(d,Sp); sparse <- sparse[-(T+1),]
  H <- diag(T) - sparse
  
  store_omega2tau <- matrix(0, nloop-burnin, 1)
  store_omega2h   <- matrix(0, nloop-burnin, 1)
  store_tau       <- matrix(0, nloop-burnin, T)
  store_h         <- matrix(0, nloop-burnin, T)
  newatau <- (T-1)/2 + atau; newah <- (T-1)/2 + ah
  
  for (loop in 1:nloop){
    invOmegatau <- diag(T)*c(1/Vtau, 1/omega2tau*rep(1,T-1))
    invSigy     <- diag(T)*c(exp(-h))
    Ktau <- t(H) %*% invOmegatau %*% H + invSigy
    Ctau <- t(chol(Ktau))
    tauhat <- solve(Ktau,(invSigy %*% y))
    tau <- tauhat + solve(t(Ctau), matrix(rnorm(T),,1))
    
    ystar <- log((y - tau)^2 + .0001)
    result <- SVRW(ystar,h,omega2h,Vh)
    h <- result[[1]]
    
    newltau <- ltau + sum((tau[2:nrow(tau)]-tau[1:(nrow(tau)-1)])^2)/2
    omega2tau <- 1/gamrand(newatau, newltau)
    newlh <- lh + sum((h[2:nrow(h)]-h[1:(nrow(h)-1)])^2)/2
    omega2h <- 1/gamrand(newah, newlh)
    
    if (loop>burnin){
      i <- loop-burnin
      store_tau[i,]       <- t(tau)
      store_h[i,]         <- t(h)
      store_omega2tau[i,] <- omega2tau
      store_omega2h[i,]   <- omega2h
    }
  }
  tauhat <- matrix(rowMeans(t(store_tau)))
  hhat   <- matrix(colMeans(store_h))
  
  if(display){
    plot(y,type="l"); lines(tauhat,col="blue")
  }
  list(tauhat = tauhat, hhat = hhat, store_tau = store_tau, store_h = store_h)
}

SVRW <- function(ystar,h,omega2h,Vh){
  T <- length(h)
  pi <- c(0.0073, .10556, .00002, .04395, .34001, .24566, .2575)
  mui <- c(-10.12999, -3.97281, -8.56686, 2.77786, .61942, 1.79518,-1.08819) - 1.2704
  sig2i <- c(5.79596, 2.61369, 5.17950, .16735, .64009, .34023, 1.26261)
  sigi <- sqrt(sig2i)
  
  temprand <- matrix(runif(T))
  q <- matrix(rep(pi,each=T),T) * dnorm(matrix(rep(ystar,each=7),byrow=TRUE,,7),
                                        matrix(rep(h,each=7),byrow=TRUE,,7) + matrix(rep(mui,each=T),T),
                                        matrix(rep(sigi,each=T),T))
  q <- q/matrix(rep(rowSums(q),each=7),byrow=TRUE,,7)
  S <- matrix(7 - rowSums(matrix(rep(temprand,each=7),byrow=TRUE,,7) < t(apply(q,1,cumsum))) + 1)
  
  Sp <- diag(T); d <- matrix(rep(0,T),1); sparse <- rbind(d,Sp); sparse <- sparse[-(T+1),]
  H <- diag(T) - sparse
  invOmegah <- diag(T)*c(1/Vh, 1/omega2h*rep(1,T-1))
  d <- matrix(mui[S]); invSigystar <- diag(T)*c(1/sig2i[S])
  Kh <- t(H) %*% invOmegah %*% H + invSigystar
  Ch <- t(chol(Kh))
  hhat <- solve(Kh,(invSigystar %*% (ystar - d)))
  h <- hhat + solve(t(Ch),matrix(rnorm(T)))
  list(h,S)
}

gamrand <- function(alpha, lambda){
  if (alpha>1){
    d <- alpha-1/3; c <- 1/sqrt(9*d)
    repeat{
      Z <- rnorm(1); if (Z > (-1/c)){
        V <- (1+c*Z)^3; U <- runif(1)
        if (log(U) <= (0.5*Z^2 + d - d*V + d*log(V))) break
      }
    }
    x <- d*V/lambda
  } else {
    x <- gamrand(alpha+1,lambda)
    x <- x*runif(1)^(1/alpha)
  }
  x
}

# ======================= Canonical utilities (ridge-style) =======================
rmse <- function(e) sqrt(mean(e^2)); mae <- function(e) mean(abs(e))
mad_ <- function(e) median(abs(e - median(e)))

compute_metrics <- function(y_real, yhat, y_all, idx_oos, h) {
  err <- y_real - yhat
  ae  <- abs(err); se <- err^2
  m_rmse <- sqrt(mean(se)); m_mae <- mean(ae)
  m_mad <- median(abs(err - median(err)))
  m_meanad <- mean(abs(err - mean(err)))
  rng <- diff(range(y_real)); m_nrmse <- m_rmse / ifelse(rng == 0, 1, rng)
  
  # h-step RW baseline on the SAME target series y_all (monthly or accumulated)
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

# trailing h-month accumulator (produces NAs for first h-1 obs)
accum_series <- function(y, h) as.numeric(stats::filter(y, rep(1, h), sides = 1))

# helper: last R finite indices ≤ end_idx (for accumulated targets)
last_R_finite <- function(finite_mask, end_idx, R) {
  cand <- which(finite_mask & seq_along(finite_mask) <= end_idx)
  if (length(cand) < R) stop("Not enough finite observations to build accumulated window.")
  cand[(length(cand)-R+1):length(cand)]
}

# =============================== Driver (train-to-origin) ===============================
load("dados.rda")  # expects 'dados' (matrix with date rownames) and 'cpi_total' column
if (!("cpi_total" %in% colnames(dados))) stop("cpi_total not found in dados.")

dates   <- rownames(dados)
y_all   <- as.numeric(dados[, "cpi_total"]); names(y_all) <- dates

nprev     <- 99
horizons  <- 1:12
acc_set   <- c(3,6,12)
idx_oos   <- (length(y_all) - nprev + 1):length(y_all)
dates_oos <- dates[idx_oos]

dir.create("forecasts/ucsv", recursive = TRUE, showWarnings = FALSE)

# -------- cluster setup (parallel across k) --------
cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
cl <- parallel::makeCluster(cores)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
parallel::clusterSetRNGStream(cl, iseed = 123)

UCSV_cpi_list <- list()
pred_cols     <- list()

# --------- A) MONTHLY horizons: train up to origin_k, predict tauhat_t ----------
for (h in horizons) {
  origin_first <- dates[idx_oos[1] - h]                     # issue-time anchor for first OOS
  end_idx_y_first <- which(dates == origin_first)
  if (length(end_idx_y_first) != 1) stop("origin_first not found for h=", h)
  
  # Fixed window length R_y: use full history up to the first origin (contiguous)
  R_y <- end_idx_y_first
  if (R_y <= 10) stop("Too short history at h=", h, " (R_y=", R_y, ")")
  
  # export per-h symbols once
  parallel::clusterExport(cl,
                          varlist = c("ucsv","SVRW","gamrand","idx_oos","dates","y_all","h","R_y"),
                          envir = environment())
  
  preds <- unlist(parallel::parLapplyLB(cl, 1:nprev, function(k) {
    origin_k  <- dates[idx_oos[k] - h]
    end_idx_y <- which(dates == origin_k)
    start_idx_y <- end_idx_y - R_y + 1
    if (start_idx_y < 1) stop("start_row < 1 at h=", h, " — not enough history")
    y_tr <- y_all[start_idx_y:end_idx_y]
    if (any(!is.finite(y_tr))) stop("Non-finite values in monthly training segment at h=", h, ", k=", k)
    fit <- ucsv(y_tr, display = FALSE)
    as.numeric(tail(fit$tauhat, 1))  # π̂_{t+h|t} = τ̂_{t|t}
  }), use.names = FALSE)
  
  names(preds) <- dates_oos
  y_real <- y_all[idx_oos]
  mets   <- compute_metrics(y_real, preds, y_all, idx_oos, h)
  
  UCSV_cpi_list[[paste0("UCSV_h", h)]] <- list(
    pred     = preds,
    errors   = mets,
    loss_abs = abs(y_real - preds),
    loss_sq  = (y_real - preds)^2
  )
  pred_cols[[paste0("h", h)]] <- preds
  
  if (h == 1) {
    cat("[MONTHLY] h=1 done with", cores, "cores. Window R_y =", R_y,
        "  last origin:", dates[idx_oos[nprev] - 1], "\n")
  }
}

# --------- B) ACCUMULATED horizons: h ∈ {3,6,12} ----------
for (h in acc_set) {
  y_acc <- accum_series(y_all, h)
  finite_mask <- is.finite(y_acc)
  
  origin_first <- dates[idx_oos[1] - h]
  end_idx_first <- which(dates == origin_first)
  if (length(end_idx_first) != 1) stop("origin_first not found for accum h=", h)
  
  # Fixed R for accumulated: number of finite points up to origin_first
  R_acc <- sum(finite_mask[1:end_idx_first])
  if (R_acc <= 10) stop("Too short accumulated history at h=", h, " (R=", R_acc, ")")
  
  parallel::clusterExport(cl,
                          varlist = c("ucsv","SVRW","gamrand","idx_oos","dates","y_acc","finite_mask","R_acc","h"),
                          envir = environment())
  
  preds <- unlist(parallel::parLapplyLB(cl, 1:nprev, function(k) {
    origin_k  <- dates[idx_oos[k] - h]
    end_idx_k <- which(dates == origin_k)
    cand <- which(finite_mask & seq_along(finite_mask) <= end_idx_k)
    if (length(cand) < R_acc) stop("Not enough finite observations at accum h=", h, ", k=", k)
    idx_train <- cand[(length(cand)-R_acc+1):length(cand)]
    y_tr_acc  <- y_acc[idx_train]
    if (any(!is.finite(y_tr_acc))) stop("Non-finite values in accumulated training at h=", h, ", k=", k)
    fit <- ucsv(y_tr_acc, display = FALSE)
    as.numeric(tail(fit$tauhat, 1))  # accumulated forecast: τ̂_{t|t} on y_acc
  }), use.names = FALSE)
  
  names(preds) <- dates_oos
  y_real_acc   <- y_acc[idx_oos]
  mets_acc     <- compute_metrics(y_real_acc, preds, y_acc, idx_oos, h)
  
  UCSV_cpi_list[[paste0("UCSV_h", h, "_acc")]] <- list(
    pred     = preds,
    errors   = mets_acc,
    loss_abs = abs(y_real_acc - preds),
    loss_sq  = (y_real_acc - preds)^2
  )
  pred_cols[[paste0("h", h, "_acc")]] <- preds
  
  if (h == 3) {
    cat("[ACCUM] h=3 done with", cores, "cores. R_acc =", R_acc,
        "  last origin:", dates[idx_oos[nprev] - 3], "\n")
  }
}

# --------------- Save (monthly + accumulated in the same files) ----------------
all_pred_mat <- do.call(cbind, pred_cols)
all_pred_mat <- as.matrix(all_pred_mat)
rownames(all_pred_mat) <- dates_oos
ord_names <- c(paste0("h", 1:12), paste0("h", acc_set, "_acc"))
ord_names <- intersect(ord_names, colnames(all_pred_mat))
all_pred_mat <- all_pred_mat[, ord_names, drop = FALSE]

dir.create("forecasts/ucsv", recursive = TRUE, showWarnings = FALSE)

write.table(all_pred_mat, "forecasts/ucsv/UCSV_cpi_total.csv",
            sep = ";", row.names = TRUE, col.names = TRUE)

save(UCSV_cpi_list, file = "forecasts/ucsv/UCSV-list-cpi_total.Rdata")

cat("\nSaved UCSV monthly (h=1..12) + accumulated (h=3,6,12) to forecasts/ucsv\n",
    "Cores used: ", cores, "\n",
    "Columns in CSV: ", paste(colnames(all_pred_mat), collapse = ", "), "\n", sep = "")
