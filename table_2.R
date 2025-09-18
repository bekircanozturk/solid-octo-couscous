# ======================= make_tables.R  —  Tables 1,2,4,5 =======================
# Project layout expected:
#   forecasts/<model>/*-list-<series>.RData|.Rdata|.rds   (primary; list of 12 elems w/ loss_abs, loss_sq)
#   forecasts/<model>/<model>_<series>.csv                (optional; preds only, not required here)
# Outputs (LaTeX & TXT): tables_out2/<series>/*

suppressPackageStartupMessages({
  library(tidyverse)   # purrr/dplyr/readr/tibble
  library(data.table)
  library(stringr)
  library(xtable)
  library(zoo)
  library(boot)        # RNG; we implement our own stationary bootstrap indices
  library(sandwich)    # HAC (QS) + Andrews bandwidth
  library(MCS)         # Model Confidence Set (CRAN)
})

# ---------------------------- USER CONFIG ---------------------------------------
root_dir        <- "forecasts"            # where model folders live
series          <- "cpi_total"            # target series name key used in filenames
out_dir         <- file.path("tables_out2", series)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Which model names in folder match these roles:
baseline_rw_key <- "rw"                   # RW folder name contains this key
rf_key          <- "rf"                   # RF folder name contains this key
bench_keys      <- c("rw","ar","ucsv")    # used in Table 4 (RF vs each of these)

# -------------------- Multi-Horizon SPA (Quaedvlieg, 2017) ----------------------
ensure_MHSPA <- function(){
  if (!requireNamespace("MultiHorizonSPA", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools")
    devtools::install_github("lucabarbaglia/MultiHorizonSPA")
  }
}
try(ensure_MHSPA(), silent = TRUE)

# -------------------------- DISCOVERY / LOAD ------------------------------------
discover_models <- function(root_dir, series, strict_series = FALSE){
  if (!dir.exists(root_dir)) stop("Root dir not found: ", root_dir)
  mdl_dirs <- list.dirs(root_dir, full.names = FALSE, recursive = FALSE)
  if (length(mdl_dirs) == 0) stop("No model subfolders inside: ", root_dir)
  
  find_candidate <- function(model_dir, series){
    dirpath <- file.path(root_dir, model_dir)
    rds_pat <- paste0("(-list-)", stringr::fixed(series), "\\.(rds|RDS|RData|Rdata)$")
    rds_cand <- list.files(dirpath, full.names = TRUE, ignore.case = TRUE, pattern = rds_pat)
    if (!length(rds_cand) && !strict_series) {
      rds_cand <- list.files(dirpath, full.names = TRUE, ignore.case = TRUE,
                             pattern = paste0(stringr::fixed(series), ".*\\.(rds|RDS|RData|Rdata)$"))
    }
    csv_cand <- list.files(
      dirpath, full.names = TRUE, ignore.case = TRUE,
      pattern = paste0("(^|/)", stringr::fixed(model_dir), "(_|-)", stringr::fixed(series), "\\.csv$")
    )
    if (!length(csv_cand) && !strict_series) {
      csv_cand <- list.files(dirpath, full.names = TRUE, ignore.case = TRUE,
                             pattern = paste0(stringr::fixed(series), ".*\\.csv$"))
    }
    tibble::tibble(
      model   = model_dir,
      rds_path= if (length(rds_cand)) rds_cand[1] else NA_character_,
      csv_path= if (length(csv_cand)) csv_cand[1] else NA_character_
    )
  }
  
  res <- purrr::map_dfr(mdl_dirs, ~find_candidate(.x, series)) %>%
    dplyr::filter(!is.na(rds_path) | !is.na(csv_path))
  
  if (nrow(res) == 0) {
    cat("\n--- Debug listing ---\n")
    for (d in mdl_dirs) {
      cat("==", file.path(root_dir,d), "==\n")
      print(list.files(file.path(root_dir,d), full.names = TRUE))
    }
    stop("No RDS/RData candidates found for series='", series, "'.")
  }
  
  cat("\nFound model files:\n")
  print(res %>% dplyr::select(model, rds_path, csv_path))
  res
}

# Reads a model object from .rds OR .RData/.Rdata and returns loss matrices
# Expect a list of length 12 (monthly) or 15 (monthly + acc3/6/12); each element has $loss_abs/$loss_sq
read_model_file <- function(path){
  obj <- NULL
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
  } else if (grepl("\\.rdata$", path, ignore.case = TRUE) ||
             grepl("\\.rdata$", tolower(path))) {
    env <- new.env(parent = emptyenv())
    nms <- load(path, envir = env)
    picks <- purrr::keep(nms, \(nm){
      x <- get(nm, envir = env)
      is.list(x) &&
        length(x) >= 12 &&
        all(vapply(x[seq_len(12)], function(z) is.list(z) && (!is.null(z$loss_abs) || !is.null(z$loss_sq) || !is.null(z$pred)), TRUE))
    })
    if (length(picks) == 0) {
      stop("No suitable 12/15-horizon list found in ", basename(path),
           ". Objects present: ", paste(nms, collapse = ", "))
    }
    obj <- get(picks[[1]], envir = env)
  } else {
    stop("Unsupported file type: ", path)
  }
  
  stopifnot(is.list(obj))
  H <- length(obj)
  
  len_loss_abs <- sapply(obj, function(z) if (!is.null(z$loss_abs)) length(z$loss_abs) else NA_integer_)
  len_loss_sq  <- sapply(obj, function(z) if (!is.null(z$loss_sq )) length(z$loss_sq ) else NA_integer_)
  len_pred     <- sapply(obj, function(z) if (!is.null(z$pred    )) length(z$pred    ) else NA_integer_)
  nprev <- max(c(len_loss_abs, len_loss_sq, len_pred), na.rm = TRUE)
  
  LA_mat <- matrix(NA_real_, nrow = nprev, ncol = H)
  LS_mat <- matrix(NA_real_, nrow = nprev, ncol = H)
  
  if (H == 15) {
    coln <- c(paste0("h", 1:12), "acc3", "acc6", "acc12")
  } else {
    coln <- paste0("h", seq_len(H))
  }
  colnames(LA_mat) <- colnames(LS_mat) <- coln
  
  for (h in seq_len(H)) {
    if (!is.null(obj[[h]]$loss_abs)) LA_mat[seq_along(obj[[h]]$loss_abs), h] <- obj[[h]]$loss_abs
    if (!is.null(obj[[h]]$loss_sq )) LS_mat[seq_along(obj[[h]]$loss_sq ), h] <- obj[[h]]$loss_sq
  }
  list(y_real = NULL, P = NULL, loss_abs = LA_mat, loss_sq = LS_mat)
}

load_model <- function(rds_path, csv_path){
  if (!is.na(rds_path) && file.exists(rds_path)) return(read_model_file(rds_path))
  stop("No RDS/RData with losses found for model. CSV without y_real cannot be used to compute losses.")
}

mdl_tbl <- discover_models(root_dir, series, strict_series = FALSE)
stopifnot(nrow(mdl_tbl) > 0)
models <- mdl_tbl$model

loaded <- mdl_tbl |>
  mutate(data = purrr::pmap(list(rds_path, csv_path), load_model))

# Identify RW & RF models by folder key
pick_first_by_key <- function(key){
  m <- models[str_detect(models, regex(key, ignore_case = TRUE))]
  if (length(m) == 0) NA_character_ else m[1]
}
rw_model <- pick_first_by_key(baseline_rw_key)
rf_model <- pick_first_by_key(rf_key)
stopifnot(!is.na(rw_model))
message("Baseline RW model = ", rw_model)
if (is.na(rf_model)) warning("RF model not found by key '", rf_key, "'. GW/SPA RF vs benchmarks will be limited.")

# Build loss arrays (per horizon) for all models
get_by_model <- function(field){
  setNames(lapply(loaded$data, `[[`, field), loaded$model)
}
LOSS_A <- get_by_model("loss_abs")  # list of matrices nprev x (12 or 15)
LOSS_S <- get_by_model("loss_sq")

# Use RW to set base dims for monthly horizons
stopifnot(!is.null(LOSS_S[[rw_model]]), !is.null(LOSS_A[[rw_model]]))
nprev <- nrow(LOSS_S[[rw_model]])
horizons_monthly <- paste0("h", 1:12)
acc_horizons     <- c(3,6,12)  # accumulated targets

# -------------------------- METRICS & RATIOS (for Tables 1,2,5) ------------------
take_monthly <- function(M) {
  keep <- intersect(horizons_monthly, colnames(M))
  M[, keep, drop = FALSE]
}
rmse_from_sq  <- function(LS) apply(LS, 2, function(col) sqrt(mean(col, na.rm = TRUE)))
mae_from_abs  <- function(LA) apply(LA, 2, function(col) mean(col, na.rm = TRUE))
mad_from_abs  <- function(LA){ apply(LA, 2, function(col) median(col, na.rm = TRUE)) }

metrics_by_model <- function(m){
  LSm <- take_monthly(LOSS_S[[m]])
  LAm <- take_monthly(LOSS_A[[m]])
  list(
    rmse = rmse_from_sq(LSm),
    mae  = mae_from_abs(LAm),
    mad  = mad_from_abs(LAm)
  )
}
METRICS <- setNames(lapply(models, metrics_by_model), models)

ratios_vs_rw <- function(mets_m, mets_rw){
  common <- intersect(names(mets_m$rmse), names(mets_rw$rmse))
  list(
    rmse = mets_m$rmse[common] / mets_rw$rmse[common],
    mae  = mets_m$mae [common] / mets_rw$mae [common],
    mad  = mets_m$mad [common] / mets_rw$mad [common]
  )
}
METS_RW <- METRICS[[rw_model]]
RATIOS  <- setNames(lapply(models, \(m) ratios_vs_rw(METRICS[[m]], METS_RW)), models)

# Accumulated horizons: prefer explicit acc3/6/12 columns; else fall back to h=3/6/12
acc_loss_vec <- function(LS_or_LA_m, k){
  cn <- colnames(LS_or_LA_m)
  accname <- paste0("acc", k)
  if (accname %in% cn) return(LS_or_LA_m[, accname])
  hname <- paste0("h", k)
  if (hname %in% cn) return(LS_or_LA_m[, hname])
  rep(NA_real_, nrow(LS_or_LA_m))
}

ACC_METS <- setNames(lapply(models, function(m){
  list(
    rmse = vapply(acc_horizons, function(k) sqrt(mean(acc_loss_vec(LOSS_S[[m]], k), na.rm = TRUE)), 0.0),
    mae  = vapply(acc_horizons, function(k) mean(acc_loss_vec(LOSS_A[[m]], k), na.rm = TRUE), 0.0),
    mad  = vapply(acc_horizons, function(k) median(acc_loss_vec(LOSS_A[[m]], k), na.rm = TRUE), 0.0)
  )
}), models)

ACC_RATIOS <- setNames(lapply(models, function(m){
  list(
    rmse = ACC_METS[[m]]$rmse / ACC_METS[[rw_model]]$rmse,
    mae  = ACC_METS[[m]]$mae  / ACC_METS[[rw_model]]$mae,
    mad  = ACC_METS[[m]]$mad  / ACC_METS[[rw_model]]$mad
  )
}), models)

# -------------------------- SPA (Hansen 2005) per horizon ------------------------
spa_pvalue_one_h <- function(loss_mat, b, B = 3000, geom_p = 1/10){
  L0 <- loss_mat[, b]
  La <- loss_mat[, -b, drop = FALSE]
  if (ncol(La) == 0) return(NA_real_)
  D  <- sweep(La, 1, L0, FUN = "-")
  Tn <- nrow(D)
  mu_hat <- colMeans(D, na.rm = TRUE)
  mu_plus <- pmax(mu_hat, 0)
  sd_hat <- apply(D, 2, sd)
  sd_hat[!is.finite(sd_hat) | sd_hat == 0] <- Inf
  t_obs  <- sqrt(Tn) * (mu_hat - mu_plus) / sd_hat
  t_obs[!is.finite(t_obs)] <- -Inf
  T_obs  <- max(t_obs, na.rm = TRUE)
  
  gen_idx <- function(Tn, p) {
    lens <- rgeom(10*Tn, prob = p) + 1L
    idx  <- integer(0)
    while (length(idx) < Tn) {
      start <- sample.int(Tn, 1)
      len   <- lens[sample.int(length(lens),1)]
      run   <- seq.int(start, length.out = len)
      run[run > Tn] <- run[run > Tn] - Tn
      idx   <- c(idx, run)
    }
    idx[seq_len(Tn)]
  }
  T_boot <- numeric(B)
  for (r in seq_len(B)) {
    idx <- gen_idx(Tn, geom_p)
    Db  <- D[idx, , drop = FALSE]
    mu_b <- colMeans(Db)
    mu_b_plus <- pmax(mu_hat, 0)
    sd_b <- apply(Db, 2, sd)
    sd_b[!is.finite(sd_b) | sd_b == 0] <- Inf
    t_b  <- sqrt(Tn) * (mu_b - mu_b_plus) / sd_b
    t_b[!is.finite(t_b)] <- -Inf
    T_boot[r] <- max(t_b, na.rm = TRUE)
  }
  mean(T_boot >= T_obs, na.rm = TRUE)
}

spa_avg_pvals <- function(loss_array, include_acc = c(3,6,12), B = 3000, geom_p = 1/10){
  Ms <- models
  out <- tibble(model = Ms, p_spa_avg = NA_real_)
  for (j in seq_along(Ms)) {
    pvals <- c()
    for (h in 1:12) {
      Lh <- sapply(Ms, function(m) loss_array[[m]][, paste0("h", h)])
      p_h <- try(spa_pvalue_one_h(Lh, b = j, B = B, geom_p = geom_p), silent = TRUE)
      if (!inherits(p_h, "try-error")) pvals <- c(pvals, p_h)
    }
    for (k in include_acc) {
      Lk <- sapply(Ms, function(m) acc_loss_vec(loss_array[[m]], k))
      keep <- which(colSums(is.finite(as.matrix(Lk))) > 0)
      if (length(keep) > 1) {
        bj <- which(keep == j); if (length(bj) == 0) next
        p_h <- try(spa_pvalue_one_h(as.matrix(Lk)[, keep, drop = FALSE], b = bj, B = B, geom_p = geom_p), silent = TRUE)
        if (!inherits(p_h, "try-error")) pvals <- c(pvals, p_h)
      }
    }
    out$p_spa_avg[j] <- mean(pvals, na.rm = TRUE)
  }
  out
}
SPA_sq <- spa_avg_pvals(LOSS_S, include_acc = acc_horizons, B = 3000)
SPA_ab <- spa_avg_pvals(LOSS_A, include_acc = acc_horizons, B = 3000)

# ----------------------- Giacomini–White (GW) & Multi-Horizon SPA (Table 4) -----
gw_pval <- function(loss_A, loss_B, alternative = c("two.sided","less","greater")){
  alternative <- match.arg(alternative)
  
  a <- as.numeric(loss_A)
  b <- as.numeric(loss_B)
  ok <- is.finite(a) & is.finite(b)          # use only periods observed for BOTH models
  if (!any(ok)) return(NA_real_)
  
  d <- a[ok] - b[ok]                          # >0 => A worse than B
  n <- length(d)
  if (n < 8) return(NA_real_)                 # sample too small -> NA
  
  fit <- try(stats::lm(d ~ 1), silent = TRUE) # mean(d) test
  if (inherits(fit, "try-error")) return(NA_real_)
  
  # HAC SE for the intercept:
  V <- try(
    sandwich::kernHAC(fit, kernel = "Quadratic Spectral", bw = bwAndrews, prewhite = FALSE),
    silent = TRUE
  )
  if (inherits(V, "try-error")) {
    V <- try(sandwich::NeweyWest(fit, prewhite = FALSE, adjust = TRUE), silent = TRUE)
    if (inherits(V, "try-error")) return(NA_real_)
  }
  se <- sqrt(V[1, 1])
  if (!is.finite(se) || se <= 0) return(NA_real_)
  
  t  <- coef(fit)[1] / se
  if (alternative == "two.sided") 2 * (1 - stats::pnorm(abs(t)))
  else if (alternative == "less")  stats::pnorm(t)
  else                             (1 - stats::pnorm(t))
}


build_table4 <- function(){
  if (is.na(rf_model)) return(NULL)
  present_bench <- models[ Reduce(`|`, lapply(bench_keys, \(k) str_detect(models, regex(k, TRUE)))) ]
  present_bench <- unique(setdiff(present_bench, rf_model))
  if (length(present_bench) == 0) return(NULL)
  
  tbl <- lapply(present_bench, function(bm){
    p_sq <- vapply(1:12, function(h) gw_pval(LOSS_S[[rf_model]][, paste0("h",h)], LOSS_S[[bm]][, paste0("h",h)], "two.sided"), 0.0)
    p_ab <- vapply(1:12, function(h) gw_pval(LOSS_A[[rf_model]][, paste0("h",h)], LOSS_A[[bm]][, paste0("h",h)], "two.sided"), 0.0)
    p_sq_acc <- vapply(acc_horizons, function(k)
      gw_pval(acc_loss_vec(LOSS_S[[rf_model]], k), acc_loss_vec(LOSS_S[[bm]], k), "two.sided"), 0.0)
    p_ab_acc <- vapply(acc_horizons, function(k)
      gw_pval(acc_loss_vec(LOSS_A[[rf_model]], k), acc_loss_vec(LOSS_A[[bm]], k), "two.sided"), 0.0)
    
    mh_u <- mh_a <- NA_real_
    if (requireNamespace("MultiHorizonSPA", quietly = TRUE)) {
      good_h <- which(vapply(1:12, function(h)
        all(is.finite(LOSS_S[[rf_model]][, paste0("h",h)])) &&
          all(is.finite(LOSS_S[[bm      ]][, paste0("h",h)])), TRUE))
      if (length(good_h) >= 2) {
        LD <- sapply(good_h, function(h)
          LOSS_S[[rf_model]][, paste0("h",h)] - LOSS_S[[bm]][, paste0("h",h)])
        mh_u <- try(MultiHorizonSPA::Test_uSPA(LossDiff = LD, L = 10, B = 999)$p_value, silent = TRUE)
        mh_a <- try(MultiHorizonSPA::Test_aSPA(LossDiff = LD, weights = rep(1/ncol(LD), ncol(LD)), L = 10, B = 999)$p_value, silent = TRUE)
        if (inherits(mh_u,"try-error")) mh_u <- NA_real_
        if (inherits(mh_a,"try-error")) mh_a <- NA_real_
      }
    }
    tibble(
      model = bm,
      horizon = c(paste0("h",1:12), paste0("acc",acc_horizons), "MH-uSPA", "MH-aSPA"),
      squared = c(p_sq, p_sq_acc, mh_u, mh_a),
      absolute= c(p_ab, p_ab_acc, NA_real_, NA_real_)
    )
  }) |> dplyr::bind_rows()
  tbl
}
TABLE4 <- build_table4()

# -------------------------- Helpers: 2-decimal text dumps ------------------------
round_df_to_text <- function(x, digits = 2, skip = character()){
  df <- as.data.frame(x)
  for (nm in names(df)) {
    if (is.numeric(df[[nm]]) && !(nm %in% skip)) {
      df[[nm]] <- formatC(round(df[[nm]], digits), digits = digits, format = "f")
    }
  }
  df
}
write_txt <- function(x, path, include_rownames = TRUE){
  if (isTRUE(include_rownames)) {
    write.table(x, file = path, sep = "\t", quote = FALSE,
                col.names = NA, row.names = TRUE)
  } else {
    write.table(x, file = path, sep = "\t", quote = FALSE,
                col.names = TRUE, row.names = FALSE)
  }
}

# ---------------------------- Write Table 4 (LaTeX + TXT) -----------------------
write_table4_tex_and_txt <- function(TABLE4){
  if (is.null(TABLE4) || nrow(TABLE4) == 0) return(invisible())
  rows <- c(paste0("h",1:12), paste0("acc",acc_horizons), "MH-uSPA", "MH-aSPA")
  cols <- unique(TABLE4$model)
  
  M_sq <- matrix(NA_real_, nrow = length(rows), ncol = length(cols),
                 dimnames = list(rows, cols))
  for (bm in cols) {
    tb <- dplyr::filter(TABLE4, model == bm)
    vals <- tb$squared; names(vals) <- tb$horizon
    M_sq[intersect(rows, names(vals)), bm] <- vals[intersect(rows, names(vals))]
  }
  X_sq <- as.data.frame(apply(M_sq, 2, function(c) sprintf("%.3f", c)))
  print(xtable(X_sq, align = c("l", rep("r", ncol(X_sq)))),
        file = file.path(out_dir, "Table4_GW_squared.tex"),
        include.rownames = TRUE, sanitize.text.function = identity)
  # TXT (2 decimals)
  X_sq_txt <- round_df_to_text(M_sq, digits = 2)
  write_txt(cbind(Horizon = rownames(X_sq_txt), X_sq_txt),
            file.path(out_dir, "Table4_GW_squared.txt"), include_rownames = FALSE)
  
  M_ab <- matrix(NA_real_, nrow = length(rows), ncol = length(cols),
                 dimnames = list(rows, cols))
  for (bm in cols) {
    tb <- dplyr::filter(TABLE4, model == bm)
    vals <- tb$absolute; names(vals) <- tb$horizon
    M_ab[intersect(rows, names(vals)), bm] <- vals[intersect(rows, names(vals))]
  }
  X_ab <- as.data.frame(apply(M_ab, 2, function(c) ifelse(is.na(c), "", sprintf("%.3f", c))))
  print(xtable(X_ab, align = c("l", rep("r", ncol(X_ab)))),
        file = file.path(out_dir, "Table4_GW_absolute.tex"),
        include.rownames = TRUE, sanitize.text.function = identity)
  # TXT (2 decimals; keep blanks for NA)
  X_ab_txt <- as.data.frame(M_ab)
  X_ab_txt[] <- lapply(X_ab_txt, function(col) {
    ifelse(is.na(col), "", formatC(round(col, 2), digits = 2, format = "f"))
  })
  write_txt(cbind(Horizon = rownames(M_ab), X_ab_txt),
            file.path(out_dir, "Table4_GW_absolute.txt"), include_rownames = FALSE)
}
write_table4_tex_and_txt(TABLE4)

# ----------------------------- MCS per horizon (for Table 5) --------------------
loss_matrix_at_h <- function(LS_or_LA, h, min_T = 10, eps_var = 1e-12){
  M <- sapply(models, function(m) {
    v <- try(LS_or_LA[[m]][, paste0("h",h)], silent = TRUE)
    if (inherits(v, "try-error") || is.null(v)) rep(NA_real_, nprev) else v
  })
  M <- as.matrix(M)
  have_T  <- colSums(is.finite(M)) >= min_T
  var_ok  <- apply(M, 2, function(col){
    x <- col[is.finite(col)]
    length(x) >= min_T && length(unique(x)) >= 2 && is.finite(stats::var(x)) && stats::var(x) > eps_var
  })
  good <- have_T & var_ok
  list(M = M[, good, drop = FALSE], keep_cols = good)
}

extract_mcs_membership <- function(obj, which_col = c("auto","MCS","MCS_M","MCS_R")){
  which_col <- match.arg(which_col)
  res <- NULL
  if (methods::is(obj, "MCS")) {
    sn <- methods::slotNames(obj)
    if ("MCS" %in% sn) {
      MCSslot <- methods::slot(obj, "MCS")
      if (is.data.frame(MCSslot) || is.matrix(MCSslot)) {
        rn <- rownames(MCSslot)
        cols <- colnames(MCSslot)
        pick_col <- if (which_col == "auto") intersect(c("MCS","MCS_M","MCS_R"), cols)[1] else intersect(which_col, cols)[1]
        if (!is.na(pick_col) && length(pick_col)) {
          keep_models <- rn[ MCSslot[, pick_col, drop = TRUE] == 1 ]
          res <- keep_models
        } else {
          res <- rn[ apply(MCSslot, 1, function(row) any(row == 1, na.rm = TRUE)) ]
        }
      }
    }
    if (is.null(res)) {
      for (snm in c("SSM","Superior","models","Model")) {
        if (snm %in% sn) {
          obj2 <- methods::slot(obj, snm)
          if (is.character(obj2)) { res <- obj2; break }
        }
      }
    }
  }
  res
}

mcs_membership <- function(loss_mat_full, type = c("sq","ab")){
  type <- match.arg(type)
  M <- loss_mat_full$M
  keep_cols <- loss_mat_full$keep_cols
  
  # Drop any rows with NA across kept models before running MCS
  if (ncol(M) >= 2) {
    M <- M[stats::complete.cases(M), , drop = FALSE]
  }
  # Ensure unique, stable colnames (avoid check.names renaming)
  if (anyDuplicated(colnames(M))) {
    colnames(M) <- make.unique(colnames(M), sep = "_dup")
  }
  memb <- rep(NA, ncol(M))
  
  if (ncol(M) >= 2 && nrow(M) >= 10) {
    out <- try({
      suppressWarnings(
        capture.output({
          MCS::MCSprocedure(Loss = as.data.frame(M, check.names = FALSE),
                            alpha = 0.5, B = 7500, statistic = "Tmax", cl = NULL)
        }, file = NULL)
      )
    }, silent = TRUE)
    
    if (!inherits(out, "try-error") && !is.null(out)) {
      pcol <- if (type == "sq") "MCS_M" else "MCS_R"
      keep_models <- extract_mcs_membership(out, which_col = pcol)
      if (is.null(keep_models)) keep_models <- extract_mcs_membership(out, which_col = "auto")
      if (!is.null(keep_models)) memb <- colnames(M) %in% keep_models
    }
  }
  
  # Expand back to full set aligned with 'models'
  full_memb <- rep(NA, length(keep_cols))
  if (sum(keep_cols) > 0) {
    full_memb[keep_cols] <- memb
    names(full_memb) <- models
  } else {
    names(full_memb) <- models
  }
  full_memb
}

MCS_SQ <- lapply(1:12, function(h) {
  loss_matrix <- loss_matrix_at_h(LOSS_S, h)
  mcs_membership(loss_matrix, type = "sq")
})
MCS_AB <- lapply(1:12, function(h) {
  loss_matrix <- loss_matrix_at_h(LOSS_A, h)
  mcs_membership(loss_matrix, type = "ab")
})
names(MCS_SQ) <- names(MCS_AB) <- paste0("h",1:12)

# ----------------------------- TABLE 2 (ratios) ---------------------------------
make_table2 <- function(){
  make_block <- function(which){
    # monthly rows h1..h12 (force alignment) + acc3/6/12
    mat <- sapply(models, function(m) {
      v <- RATIOS[[m]][[which]]
      v[paste0("h", 1:12)]
    })
    rownames(mat) <- paste0("h",1:12)
    acc <- sapply(models, function(m) {
      v <- ACC_RATIOS[[m]][[which]]       # length 3 numeric
      names(v) <- paste0("acc", acc_horizons)
      v[paste0("acc", acc_horizons)]
    })
    rownames(acc) <- paste0("acc", acc_horizons)
    rbind(mat, acc)
  }
  blk_rmse <- make_block("rmse"); blk_mae <- make_block("mae"); blk_mad <- make_block("mad")
  
  T2 <- array(NA_real_, dim = c(nrow(blk_rmse), 3*length(models)))
  coln <- c(); j <- 1
  for (m in models) {
    T2[, j]   <- blk_rmse[, m]
    T2[, j+1] <- blk_mae[,  m]
    T2[, j+2] <- blk_mad[,  m]
    coln <- c(coln, paste0(m," RMSE"), paste0(m," MAE"), paste0(m," MAD"))
    j <- j + 3
  }
  T2 <- as.data.frame(T2); rownames(T2) <- c(paste0("h",1:12), paste0("acc",acc_horizons)); colnames(T2) <- coln
  
  boldfmt <- function(x, cols_idx){
    mins <- apply(x[, cols_idx, drop=FALSE], 1, min, na.rm = TRUE)
    for (r in seq_len(nrow(x))) {
      for (c in cols_idx) {
        if (!is.finite(x[r,c])) next
        if (abs(x[r,c] - mins[r]) < 1e-12)
          x[r,c] <- sprintf("\\textbf{%.3f}", x[r,c])
        else
          x[r,c] <- sprintf("%.3f", x[r,c])
      }
    }
    x
  }
  X <- T2
  for (block in seq(1, ncol(X), by = 3)) X <- boldfmt(X, block:(block+2))
  
  print(xtable(X, align = c("l", rep("r", ncol(X)))),
        file = file.path(out_dir, "Table2_ratios_vs_RW.tex"),
        include.rownames = TRUE, sanitize.text.function = identity)
  
  # TXT (2 decimals, numeric only)
  T2_txt <- as.data.frame(T2)
  for (nm in names(T2_txt)) {
    T2_txt[[nm]] <- as.numeric(T2_txt[[nm]])
  }
  T2_txt <- round_df_to_text(T2_txt, digits = 2)
  T2_txt <- cbind(Horizon = rownames(T2_txt), T2_txt)
  write_txt(T2_txt, file.path(out_dir, "Table2_ratios_vs_RW.txt"), include_rownames = FALSE)
}
make_table2()

# ----------------------------- TABLE 5 (grid + MCS shading) ---------------------
make_table5 <- function(){
  rm <- sapply(models, function(m) c(
    RATIOS[[m]]$rmse[paste0("h",1:12)],
    ACC_RATIOS[[m]]$rmse)) |> t()
  ma <- sapply(models, function(m) c(
    RATIOS[[m]]$mae [paste0("h",1:12)],
    ACC_RATIOS[[m]]$mae )) |> t()
  colnames(rm) <- colnames(ma) <- c(paste0("h",1:12), paste0("acc",acc_horizons))
  rownames(rm) <- rownames(ma) <- models
  
  bold_rm <- apply(rm, 2, function(col){
    mn <- min(col, na.rm = TRUE)
    vapply(col, function(x) if (is.finite(x) && abs(x - mn) < 1e-12) sprintf("\\textbf{%.3f}", x) else sprintf("%.3f", x), "")
  })
  par_ma <- apply(ma, 2, function(col) vapply(col, function(x) sprintf("(%.3f)", x), ""))
  cells <- matrix(paste0(bold_rm, "", par_ma), nrow = nrow(rm), dimnames = dimnames(rm))
  
  mcs_sq_mat <- do.call(cbind, MCS_SQ); rownames(mcs_sq_mat) <- models
  mcs_ab_mat <- do.call(cbind, MCS_AB); rownames(mcs_ab_mat) <- models
  
  shade_cell <- function(txt, in_sq, in_ab){
    if (isTRUE(in_sq) && isTRUE(in_ab)) sprintf("\\cellcolor{green!20}%s", txt)
    else if (isTRUE(in_sq))             sprintf("\\cellcolor{gray!20}%s",  txt)
    else if (isTRUE(in_ab))             sprintf("\\cellcolor{blue!15}%s",  txt)
    else                                txt
  }
  for (i in seq_len(nrow(cells))) {
    for (h in 1:12) {
      cells[i, h] <- shade_cell(cells[i,h], mcs_sq_mat[i,h], mcs_ab_mat[i,h])
    }
  }
  
  mincount <- paste0(rowSums(mcs_sq_mat, na.rm = TRUE)," | ", rowSums(mcs_ab_mat, na.rm = TRUE))
  T5 <- cbind(cells, `MCS count (sq|ab)` = mincount)
  
  print(xtable(T5, align = c("l", rep("r", ncol(T5)))),
        file = file.path(out_dir, "Table5_grid_MCS_shaded.tex"),
        include.rownames = TRUE, sanitize.text.function = identity)
  
  # TXT (2 decimals): separate clean numeric dumps for RMSE and MAE ratios
  rm_txt <- round_df_to_text(as.data.frame(rm), digits = 2)
  rm_txt <- cbind(Model = rownames(rm_txt), rm_txt)
  write_txt(rm_txt, file.path(out_dir, "Table5_rmse_ratios.txt"), include_rownames = FALSE)
  ma_txt <- round_df_to_text(as.data.frame(ma), digits = 2)
  ma_txt <- cbind(Model = rownames(ma_txt), ma_txt)
  write_txt(ma_txt, file.path(out_dir, "Table5_mae_ratios.txt"), include_rownames = FALSE)
}
make_table5()

# ----------------------------- TABLE 1 (overview + SPA) -------------------------
make_table1 <- function(){
  stack_stats <- function(which){
    mat_h <- sapply(models, function(m) {
      v <- RATIOS[[m]][[which]]
      v[paste0("h",1:12)]
    }) |> t()
    colnames(mat_h) <- paste0("h",1:12)
    mat_acc <- sapply(models, function(m) ACC_RATIOS[[m]][[which]]) |> t()
    colnames(mat_acc) <- paste0("acc",acc_horizons)
    M <- cbind(mat_h, mat_acc)  # n_models x 15
    
    # Count how many horizons each model achieves the row-min
    col_mins <- apply(M, 2, min, na.rm = TRUE)
    best_count <- rowSums(sweep(M, 2, col_mins, FUN = "=="), na.rm = TRUE)
    
    tibble(
      model = rownames(M),
      avg = rowMeans(M, na.rm = TRUE),
      max = apply(M, 1, max, na.rm = TRUE),
      min = apply(M, 1, min, na.rm = TRUE),
      best_count = best_count
    )
  }
  S_rmse <- stack_stats("rmse")
  S_mae  <- stack_stats("mae")
  S_mad  <- stack_stats("mad")
  
  spa_sq <- SPA_sq |> dplyr::rename(spa_avg_sq = p_spa_avg)
  spa_ab <- SPA_ab |> dplyr::rename(spa_avg_ab = p_spa_avg)
  
  mcs_avg_sq <- tibble(model = models, mcs_avg_sq = rowMeans(do.call(cbind, MCS_SQ), na.rm = TRUE))
  mcs_avg_ab <- tibble(model = models, mcs_avg_ab = rowMeans(do.call(cbind, MCS_AB), na.rm = TRUE))
  
  mh_mcs <- try({
    if (requireNamespace("MultiHorizonSPA", quietly = TRUE)) {
      tibble(model = models, mh_mcs_sq = vapply(models, function(m){
        good_h <- which(vapply(1:12, function(h)
          all(is.finite(LOSS_S[[m]][, paste0("h",h)])) &&
            all(is.finite(LOSS_S[[rw_model]][, paste0("h",h)])), TRUE))
        if (length(good_h) < 2) return(NA_real_)
        LD <- sapply(good_h, function(h) LOSS_S[[m]][, paste0("h",h)] - LOSS_S[[rw_model]][, paste0("h",h)])
        out <- try(MultiHorizonSPA::Test_uSPA(LossDiff = LD, L = 10, B = 999)$p_value, silent = TRUE)
        if (inherits(out,"try-error")) NA_real_ else out
      }, 0.0))
    } else NULL
  }, silent = TRUE)
  if (inherits(mh_mcs, "try-error") || is.null(mh_mcs)) {
    mh_mcs <- tibble(model = models, mh_mcs_sq = NA_real_)
    message("MultiHorizonSPA not available; 'Multi-h MCS p (sq)' will be NA.")
  }
  
  T1 <- S_rmse |>
    dplyr::select(model, avg_rmse = avg, max_rmse = max, min_rmse = min, best_rmse = best_count) |>
    dplyr::left_join(S_mae |> dplyr::select(model, avg_mae=avg, max_mae=max, min_mae=min, best_mae=best_count), by="model") |>
    dplyr::left_join(S_mad |> dplyr::select(model, avg_mad=avg, max_mad=max, min_mad=min, best_mad=best_count), by="model") |>
    dplyr::left_join(spa_sq, by="model") |>
    dplyr::left_join(spa_ab, by="model") |>
    dplyr::left_join(mcs_avg_sq, by="model") |>
    dplyr::left_join(mcs_avg_ab, by="model") |>
    dplyr::left_join(mh_mcs, by="model")
  
  T1f <- T1 |>
    dplyr::transmute(
      Model = model,
      `Avg RMSE`=avg_rmse, `Avg MAE`=avg_mae, `Avg MAD`=avg_mad,
      `Max RMSE`=max_rmse, `Max MAE`=max_mae, `Max MAD`=max_mad,
      `Min RMSE`=min_rmse, `Min MAE`=min_mae, `Min MAD`=min_mad,
      `Best RMSE`=best_rmse, `Best MAE`=best_mae, `Best MAD`=best_mad,
      `SPA avg (sq)`=spa_avg_sq,
      `SPA avg (ab)`=spa_avg_ab,
      `MCS avg p (sq)`=mcs_avg_sq,
      `MCS avg p (ab)`=mcs_avg_ab,
      `Multi-h MCS p (sq)`=mh_mcs_sq
    )
  
  # LaTeX (3 decimals, keep counts as integers)
  fmt <- function(x, is_count = FALSE){
    if (is_count) as.character(x) else sprintf("%.3f", x)
  }
  X <- T1f
  num_cols <- setdiff(names(X), c("Model","Best RMSE","Best MAE","Best MAD"))
  for (cn in num_cols) X[[cn]] <- fmt(X[[cn]])
  X[["Best RMSE"]] <- fmt(X[["Best RMSE"]], TRUE)
  X[["Best MAE" ]] <- fmt(X[["Best MAE" ]], TRUE)
  X[["Best MAD" ]] <- fmt(X[["Best MAD" ]], TRUE)
  
  print(xtable(X, align = c("l", rep("r", ncol(X)))),
        file = file.path(out_dir, "Table1_overview.tex"),
        include.rownames = FALSE, sanitize.text.function = identity)
  
  # TXT (2 decimals)
  T1_txt <- round_df_to_text(T1f, digits = 2, skip = c("Model","Best RMSE","Best MAE","Best MAD"))
  write_txt(T1_txt, file.path(out_dir, "Table1_overview.txt"), include_rownames = FALSE)
}
make_table1()

message("All tables written to: ", normalizePath(out_dir, winslash = "/"))
# -------------------------- Single-horizon MCS helper ---------------------------
# modelconf_single_mcs() replicates the iterative elimination logic of the Model
# Confidence Set procedure for a single loss matrix.  The implementation keeps
# the bootstrap logic deliberately simple (resampling rows with replacement) to
# avoid any dependency on unexported MCS internals while preserving the
# elimination order used by the reference `modelconf` utilities.  The main
# difference compared with a naïve implementation is that the model removed at
# each step is determined from the statistics returned by the testing stage
# (d-bar for "Tmax" or the range metric for "TR") rather than the raw column
# means of the filtered loss matrix.
modelconf_single_mcs <- function(loss_matrix,
                                 alpha = 0.25,
                                 B = 2000,
                                 statistic = c("Tmax", "TR"),
                                 seed = NULL,
                                 quiet = TRUE) {
  statistic <- match.arg(statistic)

  loss_matrix <- as.matrix(loss_matrix)
  if (!nrow(loss_matrix) || !ncol(loss_matrix)) {
    return(list(models = colnames(loss_matrix),
                removal = tibble::tibble(step = integer(), model = character(), p_value = numeric()),
                statistic = statistic))
  }

  # Drop rows where every model is NA – they carry no information for the MCS
  # comparisons.
  row_keep <- rowSums(is.finite(loss_matrix)) > 0
  loss_matrix <- loss_matrix[row_keep, , drop = FALSE]

  if (!nrow(loss_matrix) || !ncol(loss_matrix)) {
    return(list(models = colnames(loss_matrix),
                removal = tibble::tibble(step = integer(), model = character(), p_value = numeric()),
                statistic = statistic))
  }

  model_names <- colnames(loss_matrix)
  if (is.null(model_names)) {
    model_names <- paste0("model", seq_len(ncol(loss_matrix)))
    colnames(loss_matrix) <- model_names
  }

  # Ensure we restore the RNG state if a seed is provided.
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit({
      if (isTRUE(had_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  calc_statistics <- function(mat) {
    nm <- colnames(mat)
    m  <- ncol(mat)
    if (m <= 1L) {
      return(list(dbar = setNames(rep(0, m), nm), tr = setNames(rep(0, m), nm),
                  dbar_matrix = matrix(0, nrow = m, ncol = m, dimnames = list(nm, nm))))
    }

    dbar_matrix <- matrix(NA_real_, nrow = m, ncol = m, dimnames = list(nm, nm))
    for (i in seq_len(m)) {
      for (j in seq_len(m)) {
        if (i == j) next
        dif <- mat[, i] - mat[, j]
        if (!any(is.finite(dif))) next
        dbar_matrix[i, j] <- mean(dif, na.rm = TRUE)
      }
    }

    dbar_vals <- rowMeans(dbar_matrix, na.rm = TRUE)
    dbar_vals[!is.finite(dbar_vals)] <- NA_real_
    tr_vals <- apply(abs(dbar_matrix), 1, function(row) {
      if (all(!is.finite(row))) return(NA_real_)
      max(row, na.rm = TRUE)
    })
    tr_vals[!is.finite(tr_vals)] <- NA_real_

    list(dbar = dbar_vals, tr = tr_vals, dbar_matrix = dbar_matrix)
  }

  compute_test_stat <- function(stat_obj, stat_type) {
    if (stat_type == "TR") {
      dm <- stat_obj$dbar_matrix
      if (length(dm) == 0L) return(NA_real_)
      upper <- dm[upper.tri(dm, diag = FALSE)]
      upper <- abs(upper[is.finite(upper)])
      if (!length(upper)) return(NA_real_)
      max(upper)
    } else {
      vals <- stat_obj$dbar
      vals <- vals[is.finite(vals)]
      if (!length(vals)) return(NA_real_)
      max(vals)
    }
  }

  bootstrap_stat <- function(mat, stat_type, B) {
    n <- nrow(mat)
    if (B <= 0 || n <= 1) return(numeric(0))
    res <- numeric(B)
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      boot <- mat[idx, , drop = FALSE]
      stat_obj <- calc_statistics(boot)
      res[b] <- compute_test_stat(stat_obj, stat_type)
    }
    res[is.na(res)] <- NA_real_
    res
  }

  survivors <- model_names
  removal_log <- tibble::tibble(step = integer(), model = character(), p_value = numeric())
  current <- loss_matrix
  step <- 0L

  repeat {
    # Remove models that are entirely NA in the current subset.
    valid_cols <- colSums(is.finite(current)) > 0
    current <- current[, valid_cols, drop = FALSE]
    survivors <- survivors[survivors %in% colnames(current)]

    if (ncol(current) <= 1L) break

    stat_obj <- calc_statistics(current)
    obs_stat <- compute_test_stat(stat_obj, statistic)
    if (!is.finite(obs_stat)) break

    boot_stats <- bootstrap_stat(current, statistic, B)
    if (length(boot_stats)) {
      p_val <- mean(boot_stats >= obs_stat, na.rm = TRUE)
      if (is.nan(p_val)) p_val <- NA_real_
    } else {
      p_val <- NA_real_
    }

    if (is.na(p_val) || (!is.na(alpha) && p_val > alpha)) break

    elim_metric <- if (statistic == "TR") stat_obj$tr else stat_obj$dbar
    if (is.null(elim_metric) || !length(elim_metric)) break
    elim_metric <- elim_metric[colnames(current)]
    if (all(is.na(elim_metric))) break

    finite_idx <- which(is.finite(elim_metric))
    if (!length(finite_idx)) break

    drop_idx <- finite_idx[which.max(elim_metric[finite_idx])]
    drop_model <- names(elim_metric)[drop_idx]

    step <- step + 1L
    if (!quiet) {
      msg <- sprintf("Eliminating %s (step %d) with p-value %.4f", drop_model, step,
                     ifelse(is.na(p_val), NaN, p_val))
      message(msg)
    }
    removal_log <- dplyr::bind_rows(removal_log,
      tibble::tibble(step = step, model = drop_model, p_value = p_val))

    survivors <- setdiff(survivors, drop_model)
    current <- current[, setdiff(colnames(current), drop_model), drop = FALSE]
  }

  list(models = survivors, removal = removal_log, statistic = statistic)
}

