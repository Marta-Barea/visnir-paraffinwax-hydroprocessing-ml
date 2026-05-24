suppressPackageStartupMessages({
  library(readxl)
})

default_workbook <- path.expand(
  "~/Documents/Doctorado/Tesis Doctoral/Investigación Cepsa/Vis-NIR/XDS-NIR_FOSS/Estudio según Tipo de Parafina e Hidrotratamiento/NIRS_HT_PW.xlsx"
)

args <- commandArgs(trailingOnly = TRUE)
workbook_path <- if (length(args) >= 1) path.expand(args[[1]]) else default_workbook

source_sheet <- "Spectra 19.10.21"
target_sheet <- "HT_Class"
tol_sq_error <- 1e-12

normalize_wavelengths <- function(x) {
  as.numeric(gsub(",", ".", as.character(x), fixed = TRUE))
}

read_source_spectra <- function(path, sheet) {
  raw_sheet <- read_excel(path, sheet = sheet, col_names = FALSE)

  sample_names <- trimws(as.character(unlist(raw_sheet[3, -1], use.names = FALSE)))
  wavelengths <- normalize_wavelengths(unlist(raw_sheet[-c(1:4), 1], use.names = FALSE))

  numeric_block <- as.data.frame(raw_sheet[-c(1:4), -1], stringsAsFactors = FALSE)
  numeric_block[] <- lapply(numeric_block, as.numeric)

  spectra_matrix <- t(as.matrix(numeric_block))
  mode(spectra_matrix) <- "numeric"

  colnames(spectra_matrix) <- format(wavelengths, nsmall = 2, trim = TRUE)
  rownames(spectra_matrix) <- sample_names

  list(
    sample_names = sample_names,
    wavelengths = wavelengths,
    spectra = spectra_matrix
  )
}

read_target_spectra <- function(path, sheet) {
  target <- read_excel(path, sheet = sheet)
  target_matrix <- as.matrix(target[, -1])
  mode(target_matrix) <- "numeric"

  list(
    classes = as.character(target[[1]]),
    wavelengths = normalize_wavelengths(names(target)[-1]),
    spectra = target_matrix
  )
}

nearest_source_match <- function(target_spectrum, source_matrix) {
  squared_errors <- rowSums((source_matrix - matrix(
    target_spectrum,
    nrow = nrow(source_matrix),
    ncol = ncol(source_matrix),
    byrow = TRUE
  ))^2)

  best_index <- which.min(squared_errors)

  list(
    index = best_index,
    sq_error = squared_errors[[best_index]],
    rmse = sqrt(squared_errors[[best_index]] / ncol(source_matrix))
  )
}

summarize_within_sample_variance <- function(base_sample, source_names, source_matrix, selected_names) {
  base_selector <- sub("_R[0-9]+$", "", source_names) == base_sample
  replicate_selector <- grepl("_R[0-9]+$", source_names)
  replicate_ids <- which(base_selector & replicate_selector)
  averaged_ids <- which(base_selector & !replicate_selector)

  if (length(replicate_ids) >= 2) {
    replicate_matrix <- source_matrix[replicate_ids, , drop = FALSE]
    wavelength_variance <- apply(replicate_matrix, 2, var)
    pairwise_rmse <- utils::combn(
      seq_len(nrow(replicate_matrix)),
      2,
      FUN = function(idx) {
        sqrt(mean((replicate_matrix[idx[1], ] - replicate_matrix[idx[2], ])^2))
      }
    )
  } else {
    wavelength_variance <- NA_real_
    pairwise_rmse <- NA_real_
  }

  selected_base_names <- selected_names[sub("_R[0-9]+$", "", selected_names) == base_sample]

  data.frame(
    base_sample = base_sample,
    n_source_replicates = length(replicate_ids),
    n_source_average_like = length(averaged_ids),
    selected_in_ht_class = sum(sub("_R[0-9]+$", "", selected_names) == base_sample),
    selected_names = paste(selected_base_names, collapse = "; "),
    mean_wavelength_variance = if (all(is.na(wavelength_variance))) NA_real_ else mean(wavelength_variance),
    median_wavelength_variance = if (all(is.na(wavelength_variance))) NA_real_ else median(wavelength_variance),
    max_wavelength_variance = if (all(is.na(wavelength_variance))) NA_real_ else max(wavelength_variance),
    mean_pairwise_rmse = if (all(is.na(pairwise_rmse))) NA_real_ else mean(pairwise_rmse),
    max_pairwise_rmse = if (all(is.na(pairwise_rmse))) NA_real_ else max(pairwise_rmse),
    stringsAsFactors = FALSE
  )
}

source_data <- read_source_spectra(workbook_path, source_sheet)
target_data <- read_target_spectra(workbook_path, target_sheet)

if (!identical(round(source_data$wavelengths, 2), round(target_data$wavelengths, 2))) {
  stop("Wavelength grids do not match between source and target sheets.")
}

source_matrix <- source_data$spectra
target_matrix <- target_data$spectra

matches <- lapply(seq_len(nrow(target_matrix)), function(i) {
  best <- nearest_source_match(target_matrix[i, ], source_matrix)
  source_name <- source_data$sample_names[[best$index]]

  data.frame(
    ht_row = i,
    ht_class = target_data$classes[[i]],
    source_name = source_name,
    base_sample = sub("_R[0-9]+$", "", source_name),
    source_type = if (grepl("_R[0-9]+$", source_name)) "replicate" else "mean_or_single",
    squared_error = best$sq_error,
    rmse = best$rmse,
    match_quality = if (best$sq_error <= tol_sq_error) "exact" else "approximate",
    stringsAsFactors = FALSE
  )
})

mapping <- do.call(rbind, matches)

mean_candidates <- unique(mapping$base_sample)
mean_checks <- do.call(
  rbind,
  lapply(mean_candidates, function(base_sample) {
    same_base <- sub("_R[0-9]+$", "", source_data$sample_names) == base_sample
    replicate_ids <- which(same_base & grepl("_R[0-9]+$", source_data$sample_names))
    average_ids <- which(same_base & !grepl("_R[0-9]+$", source_data$sample_names))

    if (length(replicate_ids) < 2 || length(average_ids) == 0) {
      return(NULL)
    }

    replicate_mean <- colMeans(source_matrix[replicate_ids, , drop = FALSE])

    do.call(rbind, lapply(average_ids, function(avg_id) {
      average_name <- source_data$sample_names[[avg_id]]
      rmse_to_mean <- sqrt(mean((source_matrix[avg_id, ] - replicate_mean)^2))
      data.frame(
        base_sample = base_sample,
        average_name = average_name,
        n_replicates = length(replicate_ids),
        replicate_names = paste(source_data$sample_names[replicate_ids], collapse = "; "),
        rmse_to_replicate_mean = rmse_to_mean,
        mean_match_quality = if (rmse_to_mean <= 1e-10) "exact_mean" else "not_exact_mean",
        stringsAsFactors = FALSE
      )
    }))
  })
)

variance_summary <- do.call(
  rbind,
  lapply(unique(mapping$base_sample), function(base_sample) {
    summarize_within_sample_variance(
      base_sample = base_sample,
      source_names = source_data$sample_names,
      source_matrix = source_matrix,
      selected_names = mapping$source_name
    )
  })
)

global_mean_variance <- mean(variance_summary$mean_wavelength_variance, na.rm = TRUE)
global_sd <- sqrt(global_mean_variance)
global_mean_rmse <- mean(variance_summary$mean_pairwise_rmse, na.rm = TRUE)
valid_samples <- sum(!is.na(variance_summary$mean_wavelength_variance))

cat("Workbook:", workbook_path, "\n")
cat("\n")

cat("Match quality summary:\n")
print(table(mapping$match_quality, useNA = "ifany"))
cat("Max squared error:", max(mapping$squared_error), "\n\n")

cat("Selected source types in HT_Class:\n")
print(table(mapping$source_type))
cat("\n")

cat("Average-column validation summary:\n")
print(table(mean_checks$mean_match_quality, useNA = "ifany"))
cat("\n")

cat("Global within-sample variance summary:\n")
cat("Valid base samples with replicates:", valid_samples, "\n")
cat("Mean within-sample variance:", format(global_mean_variance, scientific = TRUE, digits = 4), "\n")
cat("Global standard deviation:", format(global_sd, scientific = TRUE, digits = 4), "\n")
cat("Mean pairwise RMSE:", format(global_mean_rmse, scientific = TRUE, digits = 4), "\n\n")

cat(
  "La repetibilidad global de los espectros fue alta, ya que la varianza intra-muestra media calculada entre réplicas a lo largo de todas las longitudes de onda fue ",
  format(global_mean_variance, scientific = TRUE, digits = 3),
  ", correspondiente a una desviación estándar global de ",
  format(global_sd, scientific = TRUE, digits = 3),
  ".\n",
  sep = ""
)