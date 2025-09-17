library(readr)
library(dplyr)

in_csv        <- "dataset_tcoded.csv"
out_data_rda  <- "data.rda"      # contains object: data (tibble with a 'date' column)
out_dados_rda <- "dados.rda"     # contains object: dados (matrix with date rownames)

# read
data <- read_csv(in_csv, show_col_types = FALSE)

# find/standardize the date column name
dcol <- names(data)[grepl("^date$", names(data), ignore.case = TRUE)]
if (length(dcol) != 1) stop("Could not find a unique 'date' column.")
if (dcol != "date") data <- rename(data, date = all_of(dcol))

# ensure date is character and append "-01" when it's YYYY-MM
data <- data %>%
  mutate(date = as.character(date),
         date = if_else(grepl("^\\d{4}-\\d{2}$", date),
                        paste0(date, "-01"),
                        date))

# drop variables (case-insensitive)
drop_custom <- c("cpi_health","cpi_durables",
                 "labor_3_ans_ou_plus","labor_anciennete_moyenne",
                 "ie_euribor_1m","cpi_core")

keep_mask <- !(tolower(names(data)) %in% tolower(drop_custom))
data <- data[, keep_mask]

# ensure CPI_ALL exists and put it right after date
if (!("cpi_total" %in% names(data))) stop("cpi_total not found.")
other_cols <- setdiff(names(data), c("date","cpi_total"))
data <- data %>% select(date, cpi_total, all_of(other_cols))

# coerce non-date columns to numeric
data <- data %>% mutate(across(-date, ~ suppressWarnings(as.numeric(.))))

# --- DROP FIRST 2 MONTHS ---
if (nrow(data) < 3) stop("Not enough rows to drop the first 2.")
data <- data %>% slice(-(1:2))

dates <- data$date
X <- data %>% select(-date) %>% as.matrix()
rownames(X) <- as.character(dates)

dados <- X
print(colnames(dados))
cat("Rows:", nrow(dados), "  Cols:", ncol(dados), "\n")

# ---- save both ----
save(data,  file = out_data_rda)
save(dados, file = out_dados_rda)



# ---------------- NA CONTROL BLOCK ----------------
na_total <- sum(is.na(dados))

if (na_total > 0) {
  cat("\n*** NA CHECK ***\n")
  cat("Found", na_total, "NA cell(s) in 'dados'.\n")
  
  # per-variable NA counts (only those with NAs)
  na_by_var <- colSums(is.na(dados))
  na_by_var <- na_by_var[na_by_var > 0]
  if (length(na_by_var)) {
    cat("\nNA counts by variable:\n")
    print(sort(na_by_var, decreasing = TRUE))
  }
  
  # per-date NA counts (list dates that have any NA)
  na_by_date <- rowSums(is.na(dados))
  if (any(na_by_date > 0)) {
    cat("\nDates with NAs (first 20 shown):\n")
    na_dates <- data.frame(date = rownames(dados), n_na = na_by_date)
    print(head(na_dates[na_dates$n_na > 0, , drop = FALSE], 20))
  }
  
  # detailed list of (date, variable) NA positions → save to CSV
  na_idx <- which(is.na(dados), arr.ind = TRUE)
  na_cells <- data.frame(
    date     = rownames(dados)[na_idx[, "row"]],
    variable = colnames(dados)[na_idx[, "col"]],
    stringsAsFactors = FALSE
  )
  na_cells <- na_cells[order(na_cells$date, na_cells$variable), ]
  write.csv(na_cells, "na_cells.csv", row.names = FALSE)
  cat("\nDetailed NA cell list written to 'na_cells.csv'.\n\n")
} else {
  cat("\n*** NA CHECK ***\nNo NAs found in 'dados'.\n\n")
}
# -------------- END NA CONTROL BLOCK -----------