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
  P_mat  <- matrix(NA_real_, nrow = nprev, ncol = H)
  Y_mat  <- matrix(NA_real_, nrow = nprev, ncol = H)

  if (H == 15) {
    coln <- c(paste0("h", 1:12), "acc3", "acc6", "acc12")
  } else {
    coln <- paste0("h", seq_len(H))
  }
  colnames(LA_mat) <- colnames(LS_mat) <- coln
  colnames(P_mat)  <- colnames(Y_mat)  <- coln

  for (h in seq_len(H)) {
    if (!is.null(obj[[h]]$loss_abs)) LA_mat[seq_along(obj[[h]]$loss_abs), h] <- obj[[h]]$loss_abs
    if (!is.null(obj[[h]]$loss_sq )) LS_mat[seq_along(obj[[h]]$loss_sq ), h] <- obj[[h]]$loss_sq

    if (!is.null(obj[[h]]$pred)) {
      P_mat[seq_along(obj[[h]]$pred), h] <- obj[[h]]$pred
    }

    if (is.null(obj[[h]]$pred) && !is.null(obj[[h]]$forecast)) {
      P_mat[seq_along(obj[[h]]$forecast), h] <- obj[[h]]$forecast
    }

    get_actual <- NULL
    for (cand in c("y", "y_real", "real", "actual", "obs", "target", "Y", "Real")) {
      if (!is.null(obj[[h]][[cand]])) {
        get_actual <- obj[[h]][[cand]]
        break
      }
    }
    if (!is.null(get_actual)) {
      Y_mat[seq_along(get_actual), h] <- get_actual
    }
  }

  if (all(is.na(P_mat))) P_mat <- NULL
  if (all(is.na(Y_mat))) Y_mat <- NULL

  list(y_real = Y_mat, pred = P_mat, loss_abs = LA_mat, loss_sq = LS_mat)
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
PRED   <- get_by_model("pred")
Y_REAL <- get_by_model("y_real")

# Use RW to set base dims for monthly horizons
stopifnot(!is.null(LOSS_S[[rw_model]]), !is.null(LOSS_A[[rw_model]]))
nprev <- nrow(LOSS_S[[rw_model]])
horizons_monthly <- paste0("h", 1:12)
acc_horizons     <- c(3,6,12)  # accumulated targets

# -------------------------- Ensemble (mean/median/trim) -------------------------
normalize_pred_matrix <- function(mat, reference_cols){
  if (is.null(mat)) return(NULL)
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (is.null(colnames(mat)) && length(reference_cols) == ncol(mat)) {
    colnames(mat) <- reference_cols
  }
  mat
}

compute_ensembles <- function(actual_ref, ref_cols, base_models, pred_list, nprev){
  if (is.null(actual_ref) || is.null(ref_cols)) return(list())
  row_stat <- function(mat, fn) apply(mat, 1, fn)
  mean_fun <- function(preds){
    res <- rowMeans(preds, na.rm = TRUE)
    res[is.nan(res)] <- NA_real_
    res
  }
  median_fun <- function(preds){
    row_stat(preds, function(x){
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      stats::median(x)
    })
  }
  trim_fun <- function(preds){
    row_stat(preds, function(x){
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      mean(x, trim = 0.1)
    })
  }
  build_ensemble <- function(fun){
    pred_mat <- matrix(NA_real_, nrow = nprev, ncol = length(ref_cols),
                       dimnames = list(NULL, ref_cols))
    loss_sq  <- pred_mat
    loss_ab  <- pred_mat
    for (col in ref_cols) {
      preds <- sapply(base_models, function(m) {
        mat <- pred_list[[m]]
        if (is.null(mat)) return(rep(NA_real_, nprev))
        if (!(col %in% colnames(mat))) return(rep(NA_real_, nprev))
        vec <- mat[, col]
        if (length(vec) < nprev) vec <- c(vec, rep(NA_real_, nprev - length(vec)))
        else if (length(vec) > nprev) vec <- vec[seq_len(nprev)]
        vec
      })
      if (!is.matrix(preds)) preds <- matrix(preds, ncol = 1)
      actual <- actual_ref[, col]
      if (all(!is.finite(preds)) || all(!is.finite(actual))) {
        next
      }
      pred_col <- fun(preds)
      if (length(pred_col) != nprev) pred_col <- rep(NA_real_, nprev)
      err <- actual - pred_col
      pred_mat[, col] <- pred_col
      loss_sq[, col]  <- err^2
      loss_ab[, col]  <- abs(err)
    }
    list(pred = pred_mat, loss_sq = loss_sq, loss_abs = loss_ab)
  }
  list(
    `mean` = build_ensemble(mean_fun),
    `median` = build_ensemble(median_fun),
    `trimmed mean` = build_ensemble(trim_fun)
  )
}

extract_reference_actual <- function(){
  for (nm in models) {
    mat <- Y_REAL[[nm]]
    if (is.matrix(mat) && any(is.finite(mat))) {
      return(list(actual = mat, source = nm))
    }
  }
  NULL
}

load_actual_from_dados <- function(series, nprev, ref_cols,
                                    start_date = "1999-03-01",
                                    end_date = "2023-12-01"){
  if (is.null(ref_cols)) return(NULL)
  candidate_paths <- c("dados.rda","dados.RData","dados.Rdata","dados.rds","dados.RDS")
  for (path in candidate_paths) {
    if (!file.exists(path)) next
    env <- new.env(parent = emptyenv())
    loaded_objs <- try(load(path, envir = env), silent = TRUE)
    if (inherits(loaded_objs, "try-error")) next
    obj_names <- as.character(loaded_objs)
    for (nm in obj_names) {
      obj <- get(nm, envir = env)
      if (!(is.matrix(obj) || is.data.frame(obj))) next
      if (!(series %in% colnames(obj))) next
      vec <- as.numeric(obj[, series])
      if (!length(vec)) next

      row_ids <- rownames(obj)
      used_idx <- seq_along(vec)
      used_dates <- NULL
      if (!is.null(row_ids)) {
        date_vals <- suppressWarnings(as.Date(row_ids))
        if (all(is.na(date_vals))) {
          sel <- which(row_ids >= start_date & row_ids <= end_date)
        } else {
          start_dt <- try(as.Date(start_date), silent = TRUE)
          end_dt   <- try(as.Date(end_date),   silent = TRUE)
          if (inherits(start_dt, "try-error")) start_dt <- as.Date(start_date)
          if (inherits(end_dt, "try-error"))   end_dt   <- as.Date(end_date)
          sel <- which(!is.na(date_vals) & date_vals >= start_dt & date_vals <= end_dt)
        }
        if (length(sel)) {
          used_idx <- sel
          used_dates <- row_ids[sel]
        }
      }

      if (length(used_idx) < nprev && length(vec) >= nprev) {
        used_idx <- tail(seq_along(vec), nprev)
        if (!is.null(row_ids)) used_dates <- row_ids[used_idx]
      }
      if (length(used_idx) < nprev) next

      used_idx <- tail(used_idx, nprev)
      actual_vec <- vec[used_idx]

      mat <- matrix(NA_real_, nrow = nprev, ncol = length(ref_cols),
                    dimnames = list(NULL, ref_cols))
      if (!is.null(used_dates) && length(used_dates) >= nprev) {
        rownames(mat) <- tail(used_dates, nprev)
      }

      month_cols <- intersect(paste0("h", seq_len(12)), ref_cols)
      for (col in month_cols) mat[, col] <- actual_vec

      acc_cols <- intersect(paste0("acc", acc_horizons), ref_cols)
      if (length(acc_cols)) {
        for (col in acc_cols) {
          k <- as.integer(sub("acc", "", col))
          if (!is.na(k) && k > 0 && length(vec) >= k) {
            acc_vals <- zoo::rollapply(vec, width = k, align = "right",
                                      FUN = function(x) mean(x, na.rm = TRUE), fill = NA)
            mat[, col] <- acc_vals[used_idx]
          }
        }
      }
      date_range <- NULL
      if (!is.null(rownames(mat))) {
        rn <- rownames(mat)
        date_range <- c(head(rn, 1), tail(rn, 1))
      }
      return(list(actual = mat, source = nm, path = path, range = date_range))
    }
  }
  NULL
}

ref_info <- extract_reference_actual()

actual_reference <- NULL
actual_reference_source <- NULL
actual_reference_path <- NULL
actual_reference_range <- NULL

if (!is.null(ref_info)) {
  actual_reference <- ref_info$actual
  actual_reference_source <- ref_info$source
}

if (is.null(actual_reference)) {
  dados_info <- load_actual_from_dados(series, nprev, colnames(LOSS_S[[rw_model]]))
  if (!is.null(dados_info)) {
    actual_reference <- dados_info$actual
    actual_reference_source <- dados_info$source
    actual_reference_path <- dados_info$path
    actual_reference_range <- dados_info$range
  }
}

if (!is.null(actual_reference)) {
  if (!is.matrix(actual_reference)) actual_reference <- as.matrix(actual_reference)
  ref_cols <- colnames(actual_reference)
  if (is.null(ref_cols) || !length(ref_cols)) {
    fallback_cols <- colnames(LOSS_S[[rw_model]])
    if (!is.null(fallback_cols)) {
      colnames(actual_reference) <- fallback_cols
      ref_cols <- fallback_cols
    }
  }
  if (is.null(ref_cols)) ref_cols <- paste0("h", seq_len(ncol(actual_reference)))

  if (!is.null(actual_reference) && nrow(actual_reference) < nprev) {
    actual_reference <- rbind(actual_reference,
                              matrix(NA_real_, nprev - nrow(actual_reference), ncol(actual_reference)))
  }
  if (!is.null(actual_reference) && nrow(actual_reference) > nprev) {
    actual_reference <- actual_reference[seq_len(nprev), , drop = FALSE]
  }

  base_models <- models
  for (nm in base_models) {
    PRED[[nm]] <- normalize_pred_matrix(PRED[[nm]], ref_cols)
  }

  ensembles <- compute_ensembles(actual_reference, ref_cols, base_models, PRED, nprev)

  added_ensembles <- character(0)
  for (nm in names(ensembles)) {
    ens <- ensembles[[nm]]
    if (is.null(ens)) next
    if (all(is.na(ens$loss_sq)) && all(is.na(ens$loss_abs))) next
    if (nm %in% models) next
    LOSS_S[[nm]] <- ens$loss_sq
    LOSS_A[[nm]] <- ens$loss_abs
    PRED[[nm]]   <- ens$pred
    added_ensembles <- c(added_ensembles, nm)
  }

  if (length(added_ensembles)) {
    models <- c(models, added_ensembles)
    if (!is.null(actual_reference_path)) {
      range_txt <- if (!is.null(actual_reference_range) && length(actual_reference_range) == 2)
        paste0(" [", actual_reference_range[1], " to ", actual_reference_range[2], "]") else ""
      message("Actuals sourced from ", actual_reference_path, " (object '", actual_reference_source,
              "')", range_txt, "; added ensembles: ", paste(added_ensembles, collapse = ", "))
    } else {
      message("Added ensemble forecasts: ", paste(added_ensembles, collapse = ", "))
    }
  } else if (!is.null(actual_reference_path)) {
    range_txt <- if (!is.null(actual_reference_range) && length(actual_reference_range) == 2)
      paste0(" [", actual_reference_range[1], " to ", actual_reference_range[2], "]") else ""
    warning("Actual series loaded from ", actual_reference_path, " (object '", actual_reference_source,
            "')", range_txt, " but ensembles could not be formed (missing predictions).")
  }
} else {
  warning("No actual series found among models or dados.* files; ensemble rows will be skipped.")
}

ACTUAL_REFERENCE <- actual_reference
ACTUAL_REFERENCE_SOURCE <- actual_reference_source
ACTUAL_REFERENCE_PATH <- actual_reference_path
ACTUAL_REFERENCE_RANGE <- actual_reference_range

models <- unique(models)

ensemble_models <- intersect(c("mean","median","trimmed mean"), models)
model_labels <- setNames(ifelse(models %in% ensemble_models,
                                paste0(models, " (ens)"),
                                models),
                         models)

# -------------------------- METRICS & RATIOS (for Tables 1,2,5) ------------------
take_monthly <- function(M) {
  if (is.null(M)) return(NULL)
  keep <- intersect(horizons_monthly, colnames(M))
  if (!length(keep)) return(NULL)
  M[, keep, drop = FALSE]
}

metrics_by_model <- function(m){
  monthly_cols <- horizons_monthly
  rmse_vals <- setNames(rep(NA_real_, length(monthly_cols)), monthly_cols)
  mae_vals  <- rmse_vals
  mad_vals  <- rmse_vals

  LSm <- take_monthly(LOSS_S[[m]])
  LAm <- take_monthly(LOSS_A[[m]])
  preds <- take_monthly(PRED[[m]])
  actual_m <- take_monthly(Y_REAL[[m]])
  if ((is.null(actual_m) || all(!is.finite(actual_m))) && !is.null(ACTUAL_REFERENCE)) {
    actual_m <- take_monthly(ACTUAL_REFERENCE)
  }

  for (col in monthly_cols) {
    pred_col <- if (!is.null(preds) && col %in% colnames(preds)) preds[, col] else NULL
    actual_col <- if (!is.null(actual_m) && col %in% colnames(actual_m)) actual_m[, col] else NULL
    if (!is.null(pred_col) && !is.null(actual_col)) {
      err <- actual_col - pred_col
      err <- err[is.finite(err)]
      if (length(err)) {
        rmse_vals[col] <- sqrt(mean(err^2))
        mae_vals[col]  <- mean(abs(err))
        mad_vals[col]  <- stats::median(abs(err))
      }
    }

    if (is.na(rmse_vals[col]) && !is.null(LSm) && col %in% colnames(LSm)) {
      se <- LSm[, col]
      if (any(is.finite(se))) rmse_vals[col] <- sqrt(mean(se, na.rm = TRUE))
    }

    if (is.na(mae_vals[col]) && !is.null(LAm) && col %in% colnames(LAm)) {
      ae <- LAm[, col]
      if (any(is.finite(ae))) {
        mae_vals[col] <- mean(ae, na.rm = TRUE)
        mad_vals[col] <- stats::median(ae, na.rm = TRUE)
      }
    } else if (!is.na(mae_vals[col]) && is.na(mad_vals[col]) && !is.null(LAm) && col %in% colnames(LAm)) {
      ae <- LAm[, col]
      if (any(is.finite(ae))) mad_vals[col] <- stats::median(ae, na.rm = TRUE)
    }
  }

  list(rmse = rmse_vals, mae = mae_vals, mad = mad_vals)
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
      diff_info <- lapply(1:12, function(h){
        hname <- paste0("h", h)
        rf_col <- LOSS_S[[rf_model]][, hname]
        bm_col <- LOSS_S[[bm]][, hname]
        ok <- is.finite(rf_col) & is.finite(bm_col)
        diff <- rf_col - bm_col
        diff[!ok] <- NA_real_
        list(name = hname, diff = diff, count = sum(ok))
      })
      keep <- Filter(function(x) is.list(x) && !is.null(x$count) && x$count >= 10, diff_info)
      if (length(keep) >= 2) {
        LD <- do.call(cbind, lapply(keep, `[[`, "diff"))
        colnames(LD) <- vapply(keep, `[[`, character(1), "name")
        LD <- LD[stats::complete.cases(LD), , drop = FALSE]
        if (nrow(LD) >= 10 && ncol(LD) >= 2) {
          mh_u <- try(MultiHorizonSPA::Test_uSPA(LossDiff = LD, L = 10, B = 999)$p_value, silent = TRUE)
          mh_a <- try(MultiHorizonSPA::Test_aSPA(LossDiff = LD, weights = rep(1/ncol(LD), ncol(LD)), L = 10, B = 999)$p_value, silent = TRUE)
          if (inherits(mh_u,"try-error")) mh_u <- NA_real_
          if (inherits(mh_a,"try-error")) mh_a <- NA_real_
        }
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
  colnames(X_sq) <- model_labels[colnames(X_sq)]
  print(xtable(X_sq, align = c("l", rep("r", ncol(X_sq)))),
        file = file.path(out_dir, "Table4_GW_squared.tex"),
        include.rownames = TRUE, sanitize.text.function = identity)
  # TXT (2 decimals)
  X_sq_txt <- round_df_to_text(M_sq, digits = 2)
  colnames(X_sq_txt) <- model_labels[colnames(X_sq_txt)]
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
  colnames(X_ab) <- model_labels[colnames(X_ab)]
  print(xtable(X_ab, align = c("l", rep("r", ncol(X_ab)))),
        file = file.path(out_dir, "Table4_GW_absolute.tex"),
        include.rownames = TRUE, sanitize.text.function = identity)
  # TXT (2 decimals; keep blanks for NA)
  X_ab_txt <- as.data.frame(M_ab)
  X_ab_txt[] <- lapply(X_ab_txt, function(col) {
    ifelse(is.na(col), "", formatC(round(col, 2), digits = 2, format = "f"))
  })
  colnames(X_ab_txt) <- model_labels[colnames(X_ab_txt)]
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

stationary_bootstrap_indices <- function(Tn, B, p = 0.1){
  if (Tn <= 0) stop("Sample size must be positive for stationary bootstrap.")
  idx <- matrix(1L, nrow = Tn, ncol = B)
  if (Tn == 1L) return(idx)
  for (b in seq_len(B)) {
    t <- 1L
    while (t <= Tn) {
      start <- sample.int(Tn, 1L)
      block_len <- rgeom(1L, p) + 1L
      for (k in seq_len(block_len)) {
        if (t > Tn) break
        idx[t, b] <- ((start - 1L + (k - 1L)) %% Tn) + 1L
        t <- t + 1L
      }
    }
  }
  idx
}

modelconf_d_series <- function(M){
  m <- ncol(M)
  out <- matrix(NA_real_, nrow = nrow(M), ncol = m, dimnames = list(NULL, colnames(M)))
  for (i in seq_len(m)) {
    others <- setdiff(seq_len(m), i)
    if (!length(others)) next
    diffs <- sweep(M[, others, drop = FALSE], 1, M[, i], FUN = "-")
    out[, i] <- rowMeans(diffs)
  }
  out
}

modelconf_pairwise_series <- function(M){
  if (ncol(M) < 2) {
    return(list(series = matrix(NA_real_, nrow = nrow(M), ncol = 0), combos = matrix(integer(0), nrow = 2)))
  }
  cmb <- utils::combn(seq_len(ncol(M)), 2)
  coln <- apply(cmb, 2, function(idx) paste0(colnames(M)[idx], collapse = "::"))
  arr <- matrix(NA_real_, nrow = nrow(M), ncol = ncol(cmb), dimnames = list(NULL, coln))
  for (k in seq_len(ncol(cmb))) {
    i <- cmb[1, k]; j <- cmb[2, k]
    arr[, k] <- M[, i] - M[, j]
  }
  list(series = arr, combos = cmb)
}

modelconf_statistic <- function(M, statistic = c("Tmax","TR")){
  statistic <- match.arg(statistic)
  Tn <- nrow(M)
  if (statistic == "Tmax") {
    d_series <- modelconf_d_series(M)
    dbar <- colMeans(d_series)
    sdv  <- apply(d_series, 2, stats::sd)
    valid <- is.finite(dbar) & is.finite(sdv) & sdv > 0
    t_obs <- rep(NA_real_, length(dbar))
    if (any(valid)) t_obs[valid] <- sqrt(Tn) * dbar[valid] / sdv[valid]
    T_obs <- if (any(valid)) max(t_obs[valid], na.rm = TRUE) else NA_real_
    list(statistic = "Tmax", series = d_series, dbar = dbar, sdv = sdv, t_obs = t_obs, T_obs = T_obs)
  } else {
    pw <- modelconf_pairwise_series(M)
    dbar <- colMeans(pw$series)
    sdv  <- apply(pw$series, 2, stats::sd)
    valid <- is.finite(dbar) & is.finite(sdv) & sdv > 0
    t_obs <- rep(NA_real_, length(dbar))
    if (any(valid)) t_obs[valid] <- sqrt(Tn) * dbar[valid] / sdv[valid]
    T_obs <- if (any(valid)) max(abs(t_obs[valid]), na.rm = TRUE) else NA_real_
    list(statistic = "TR", series = pw$series, dbar = dbar, sdv = sdv, t_obs = t_obs, T_obs = T_obs, combos = pw$combos)
  }
}

modelconf_bootstrap_stat <- function(stat_obj, idx){
  series_b <- stat_obj$series[idx, , drop = FALSE]
  dbar_b <- if (ncol(series_b)) colMeans(series_b) else numeric(0)
  if (!length(stat_obj$sdv)) return(NA_real_)
  valid <- is.finite(stat_obj$sdv) & stat_obj$sdv > 0
  if (!length(dbar_b) || !any(valid)) return(NA_real_)
  Tn <- nrow(stat_obj$series)
  t_b <- rep(NA_real_, length(dbar_b))
  t_b[valid] <- sqrt(Tn) * (dbar_b[valid] - stat_obj$dbar[valid]) / stat_obj$sdv[valid]
  if (stat_obj$statistic == "Tmax") max(t_b[valid], na.rm = TRUE) else max(abs(t_b[valid]), na.rm = TRUE)
}

modelconf_single_mcs <- function(M, alpha = 0.1, B = 5000, statistic = "Tmax", boot_p = 0.1){
  M <- as.matrix(M)
  if (is.null(colnames(M))) colnames(M) <- paste0("model", seq_len(ncol(M)))
  all_models <- colnames(M)
  membership <- setNames(rep(FALSE, length(all_models)), all_models)
  pvals      <- setNames(rep(NA_real_, length(all_models)), all_models)
  last_p <- NA_real_
  current <- M

  while (ncol(current) >= 1) {
    if (ncol(current) == 1) {
      survivors <- colnames(current)
      membership[survivors] <- TRUE
      if (is.na(pvals[survivors])) pvals[survivors] <- ifelse(is.na(last_p), 1, last_p)
      break
    }

    stat_obj <- modelconf_statistic(current, statistic)
    if (!is.finite(stat_obj$T_obs)) {
      membership[colnames(current)] <- TRUE
      break
    }

    boot_idx <- stationary_bootstrap_indices(nrow(current), B, boot_p)
    boot_vals <- apply(boot_idx, 2, function(id) modelconf_bootstrap_stat(stat_obj, id))
    boot_vals <- boot_vals[is.finite(boot_vals)]
    if (!length(boot_vals)) {
      membership[colnames(current)] <- TRUE
      break
    }

    pval <- mean(boot_vals >= stat_obj$T_obs)
    last_p <- pval

    if (!is.finite(pval) || is.na(pval) || pval > alpha) {
      survivors <- colnames(current)
      membership[survivors] <- TRUE
      pvals[survivors][is.na(pvals[survivors])] <- pval
      break
    }

    worst_idx <- which.max(colMeans(current))
    worst_model <- colnames(current)[worst_idx]
    pvals[worst_model] <- pval
    membership[worst_model] <- FALSE
    current <- current[, -worst_idx, drop = FALSE]
    if (!ncol(current)) break
  }

  survivors <- names(membership)[membership]
  if (length(survivors)) {
    fill <- if (is.na(last_p)) 1 else last_p
    pvals[survivors][is.na(pvals[survivors])] <- fill
  }

  list(membership = membership, pvalues = pvals, last_pvalue = last_p)
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
  pvals <- rep(NA_real_, ncol(M))

  if (ncol(M) >= 2 && nrow(M) >= 10) {
    res <- modelconf_single_mcs(M, alpha = 0.1, B = 5000, statistic = "Tmax", boot_p = 0.1)
    memb <- res$membership[colnames(M)]
    pvals <- res$pvalues[colnames(M)]
  } else if (ncol(M) == 1) {
    memb <- setNames(rep(TRUE, ncol(M)), colnames(M))
    pvals <- setNames(rep(1, ncol(M)), colnames(M))
  }

  # Expand back to full set aligned with 'models'
  full_memb <- rep(NA, length(keep_cols))
  full_pvals <- rep(NA_real_, length(keep_cols))
  names(full_memb) <- names(full_pvals) <- names(keep_cols)
  if (sum(keep_cols) > 0) {
    full_memb[keep_cols] <- memb
    full_pvals[keep_cols] <- pvals
  }

  list(membership = full_memb, pvalues = full_pvals)
}

MCS_SQ <- list()
MCS_AB <- list()
MCS_P_sq <- list()
MCS_P_ab <- list()
for (h in 1:12) {
  loss_matrix_sq <- loss_matrix_at_h(LOSS_S, h)
  res_sq <- mcs_membership(loss_matrix_sq, type = "sq")
  MCS_SQ[[h]] <- res_sq$membership
  MCS_P_sq[[h]] <- res_sq$pvalues

  loss_matrix_ab <- loss_matrix_at_h(LOSS_A, h)
  res_ab <- mcs_membership(loss_matrix_ab, type = "ab")
  MCS_AB[[h]] <- res_ab$membership
  MCS_P_ab[[h]] <- res_ab$pvalues
}
names(MCS_SQ) <- names(MCS_AB) <- names(MCS_P_sq) <- names(MCS_P_ab) <- paste0("h",1:12)

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
    label <- model_labels[[m]]
    T2[, j]   <- blk_rmse[, m]
    T2[, j+1] <- blk_mae[,  m]
    T2[, j+2] <- blk_mad[,  m]
    coln <- c(coln, paste0(label," RMSE"), paste0(label," MAE"), paste0(label," MAD"))
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
  orig_rows <- rownames(cells)
  names(mincount) <- model_labels[orig_rows]
  rownames(cells) <- model_labels[orig_rows]
  T5 <- cbind(cells, `MCS count (sq|ab)` = mincount)

  print(xtable(T5, align = c("l", rep("r", ncol(T5)))),
        file = file.path(out_dir, "Table5_grid_MCS_shaded.tex"),
        include.rownames = TRUE, sanitize.text.function = identity)

  # TXT (2 decimals): separate clean numeric dumps for RMSE and MAE ratios
  rm_txt <- round_df_to_text(as.data.frame(rm), digits = 2)
  rownames(rm_txt) <- model_labels[rownames(rm_txt)]
  rm_txt <- cbind(Model = rownames(rm_txt), rm_txt)
  write_txt(rm_txt, file.path(out_dir, "Table5_rmse_ratios.txt"), include_rownames = FALSE)
  ma_txt <- round_df_to_text(as.data.frame(ma), digits = 2)
  rownames(ma_txt) <- model_labels[rownames(ma_txt)]
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
  
  mcs_p_sq_mat <- NULL
  if (length(MCS_P_sq)) {
    tmp <- try(do.call(cbind, MCS_P_sq), silent = TRUE)
    if (!inherits(tmp, "try-error")) mcs_p_sq_mat <- tmp
  }
  if (!is.null(dim(mcs_p_sq_mat))) {
    mcs_p_sq_mat <- mcs_p_sq_mat[models, , drop = FALSE]
  }
  mcs_p_ab_mat <- NULL
  if (length(MCS_P_ab)) {
    tmp <- try(do.call(cbind, MCS_P_ab), silent = TRUE)
    if (!inherits(tmp, "try-error")) mcs_p_ab_mat <- tmp
  }
  if (!is.null(dim(mcs_p_ab_mat))) {
    mcs_p_ab_mat <- mcs_p_ab_mat[models, , drop = FALSE]
  }
  mcs_avg_sq <- tibble(model = models,
                       mcs_avg_sq = if (is.null(dim(mcs_p_sq_mat))) NA_real_ else rowMeans(mcs_p_sq_mat, na.rm = TRUE))
  mcs_avg_ab <- tibble(model = models,
                       mcs_avg_ab = if (is.null(dim(mcs_p_ab_mat))) NA_real_ else rowMeans(mcs_p_ab_mat, na.rm = TRUE))
  
  mh_mcs <- try({
    if (requireNamespace("MultiHorizonSPA", quietly = TRUE)) {
      tibble(model = models, mh_mcs_sq = vapply(models, function(m){
        diff_info <- lapply(1:12, function(h){
          hname <- paste0("h", h)
          mod_col <- LOSS_S[[m]][, hname]
          rw_col  <- LOSS_S[[rw_model]][, hname]
          ok <- is.finite(mod_col) & is.finite(rw_col)
          diff <- mod_col - rw_col
          diff[!ok] <- NA_real_
          list(name = hname, diff = diff, count = sum(ok))
        })
        keep <- Filter(function(x) is.list(x) && !is.null(x$count) && x$count >= 10, diff_info)
        if (length(keep) < 2) return(NA_real_)
        LD <- do.call(cbind, lapply(keep, `[[`, "diff"))
        LD <- LD[stats::complete.cases(LD), , drop = FALSE]
        if (nrow(LD) < 10 || ncol(LD) < 2) return(NA_real_)
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
      Model = model_labels[model],
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
