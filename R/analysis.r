#' @title Calculate Comprehensive Thermal Statistics
#' @description Computes a detailed set of summary statistics from the thermal matrix.
#'              Metrics include ROI size (Pixels), central tendency (Mean, Median, Peak Density), dispersion (SD, IQR, CV),
#'              and range (Min, Max, Quantiles). NA values (background) are automatically excluded.
#'
#' @details The function calculates the following metrics:
#'          \itemize{
#'            \item \strong{Pixels:} Number of valid non-NA pixels included in the analysis.
#'            \item \strong{Min/Max:} Extremities of the temperature distribution.
#'            \item \strong{Mean/Median:} Measures of central tendency.
#'            \item \strong{SD (Standard Deviation):} Absolute measure of spread.
#'            \item \strong{Q25/Q75:} The 25th and 75th percentiles of the temperature distribution.
#'            \item \strong{IQR (Interquartile Range):} Robust measure of spread (Q75 - Q25).
#'            \item \strong{CV (Coefficient of Variation):} Relative measure of spread (SD / Mean), useful for assessing thermal heterogeneity.
#'            \item \strong{Peak_Density:} The temperature value corresponding to the peak of the kernel density estimate.
#'          }
#'
#' @param img_obj A 'BioThermR' object.
#' @param use_processed Logical. If \code{TRUE} (default), calculates statistics on the 'processed' matrix
#'                      (where background is likely masked as NA). If \code{FALSE}, uses the 'raw' matrix.
#'
#' @return A 'BioThermR' object with the \code{stats} slot updated containing a data frame of results.
#' @importFrom stats sd median quantile density IQR
#' @export
#' @examples
#' img_obj <- system.file("extdata", "C05.raw", package = "BioThermR")
#' img <- read_thermal_raw(img_obj)
#' img <- analyze_thermal_stats(img)
analyze_thermal_stats <- function(img_obj, use_processed = TRUE) {

  if (!inherits(img_obj, "BioThermR")) {
    stop("Error: Input must be a 'BioThermR' object.")
  }

  mat <- if (use_processed) img_obj$processed else img_obj$raw

  # Extract valid pixels (remove NA background)
  vals <- as.vector(mat)
  vals <- vals[!is.na(vals)]

  if (length(vals) == 0) {
    warning("Warning: No valid pixels found for analysis (Matrix is all NA). Stats will be NA.")
    stats_df <- data.frame(
      Metric = c("Pixels", "Min", "Max", "Mean", "Median", "SD", "Q25", "Q75", "IQR", "CV", "Peak_Density"),
      Value = c(0, rep(NA_real_, 10))
    )
  } else {
    # Calculate Density Peak
    peak_val <- if (length(vals) < 2) {
      NA_real_
    } else {
      d <- density(vals)
      d$x[which.max(d$y)]
    }

    # Create Statistics Data Frame
    stats_df <- data.frame(
      Metric = c("Pixels","Min", "Max", "Mean", "Median", "SD", "Q25", "Q75", "IQR", "CV","Peak_Density"),
      Value = c(
        length(vals),
        min(vals),
        max(vals),
        mean(vals),
        median(vals),
        sd(vals),
        unname(quantile(vals, 0.25)),
        unname(quantile(vals, 0.75)),
        IQR(vals),
        sd(vals) / mean(vals),
        peak_val
      )
    )
  }

  # Update the object
  img_obj$stats <- stats_df

  return(img_obj)
}


#' @title Compile Batch Statistics into a Tidy Data Frame
#' @description Iterates through a list of 'BioThermR' objects, extracts the pre-calculated statistics
#'              (from \code{analyze_thermal_stats}), and aggregates them into a single summary data frame.
#'              This function transforms the nested list structure into a flat, tabular format (Tidy Data)
#'              suitable for downstream statistical analysis (ANOVA, t-test) or visualization.
#'
#' @param img_list A list of 'BioThermR' objects (typically the output of a batch processing workflow).
#'                 Note: \code{\link{analyze_thermal_stats}} must be run on these objects first to populate the 'stats' slot.
#'
#' @return A \code{data.frame} where:
#'         \itemize{
#'           \item \strong{Rows} represent individual images.
#'           \item \strong{Columns} include 'Filename' and all metrics computed by \code{analyze_thermal_stats}
#'                 (e.g., Min, Max, Mean, Median, SD, IQR, CV, Peak_Density).
#'         }
#' @export
#' @examples
#' \donttest{
#' # 1. Import and Process
#' img_obj_list <- system.file("extdata",package = "BioThermR")
#' img_list <- read_thermal_batch(img_obj_list)
#' img_list <- lapply(img_list, analyze_thermal_stats)
#'
#' # 2. Compile Results
#' df_results <- compile_batch_stats(img_list)
#' }
compile_batch_stats <- function(img_list) {

  # 1. Validation
  if (!is.list(img_list)) {
    stop("Error: Input must be a list of BioThermR objects.")
  }

  message(paste("Compiling statistics for", length(img_list), "images..."))

  # 2. Extract and Transpose
  # We use lapply to iterate through the list and format each object's stats into a single row
  rows_list <- lapply(names(img_list), function(img_name) {

    obj <- img_list[[img_name]]

    # Check if valid object
    if (!inherits(obj, "BioThermR")) return(NULL)

    # Check if stats exist
    if (is.null(obj$stats)) {
      warning(paste("Warning: No stats found for", img_name, "- Skipping."))
      return(NULL)
    }

    # We convert it to a named vector, then transpose to a 1-row data frame
    stats_vec <- setNames(obj$stats$Value, obj$stats$Metric)
    row_df <- as.data.frame(t(stats_vec))

    # Add Filename as the first column
    # We prefer using the list name or metadata filename
    fname <- if (!is.null(obj$meta$filename)) obj$meta$filename else img_name
    row_df <- cbind(Filename = fname, row_df)

    return(row_df)
  })

  # 3. Combine all rows
  # Remove NULLs (failed or skipped items)
  rows_list <- rows_list[!sapply(rows_list, is.null)]

  if (length(rows_list) == 0) {
    warning("No valid statistics found to compile.")
    return(data.frame())
  }

  # Bind rows together
  final_df <- do.call(rbind, rows_list)

  # Reset row names to simple numbers
  rownames(final_df) <- NULL

  return(final_df)
}


#' @title Merge Thermal Stats with Clinical Data
#' @description Merges the compiled thermal statistics with an external clinical dataset based on filenames/IDs.
#'              Automatically handles column name conflicts (removes .x/.y suffixes).
#'
#' @param thermal_df Data frame. The output from compile_batch_stats().
#' @param clinical_df Data frame. Your external clinical data.
#' @param thermal_id String. Column name in thermal_df to merge on. Default is "Filename".
#' @param clinical_id String. Column name in clinical_df to merge on. Default is "Filename".
#' @param clean_ids Logical. If TRUE, automatically removes file paths and extensions for matching. Default is TRUE.
#'
#' @return A merged data frame.
#' @export
merge_clinical_data <- function(thermal_df, clinical_df,
                                thermal_id = "Filename",
                                clinical_id = "Filename",
                                clean_ids = TRUE) {

  # 1. Validation
  if (!is.data.frame(thermal_df) || !is.data.frame(clinical_df)) {
    stop("Error: Both inputs must be data frames.")
  }

  if (!(thermal_id %in% names(thermal_df))) {
    stop(paste("Error: Column", thermal_id, "not found in thermal_df."))
  }

  if (!(clinical_id %in% names(clinical_df))) {
    stop(paste("Error: Column", clinical_id, "not found in clinical_df."))
  }

  message("Merging thermal data with clinical records...")

  # 2. Prepare IDs for matching
  t_id_vec <- as.character(thermal_df[[thermal_id]])
  c_id_vec <- as.character(clinical_df[[clinical_id]])

  if (clean_ids) {
    t_id_vec <- tools::file_path_sans_ext(basename(t_id_vec))
    c_id_vec <- tools::file_path_sans_ext(basename(c_id_vec))
  }

  # Add temporary keys
  thermal_df$join_key_temp <- t_id_vec
  clinical_df$join_key_temp <- c_id_vec

  # 3. Perform Merge (Left Join)
  merged_df <- merge(thermal_df, clinical_df,
                     by = "join_key_temp",
                     all.x = TRUE,
                     sort = FALSE)

  # 4. Cleanup & Fix Duplicate Columns (.x, .y)
  merged_df$join_key_temp <- NULL # Remove the temporary key

  # Check if we have "Filename.x" and "Filename.y"
  conflict_col_x <- paste0(thermal_id, ".x")
  conflict_col_y <- paste0(thermal_id, ".y") # Assuming thermal_id name is same in clinical_df

  # If the thermal ID column got renamed to .x, rename it back
  if (conflict_col_x %in% names(merged_df)) {
    # 1. Rename .x back to original name (e.g., "Filename.x" -> "Filename")
    names(merged_df)[names(merged_df) == conflict_col_x] <- thermal_id

    message(paste("Duplicate column detected. Kept", thermal_id, "from thermal data."))
  }

  # If there is a matching .y column (from clinical), remove it to avoid confusion
  # (But only if the clinical_id name was the same as thermal_id)
  if (thermal_id == clinical_id) {
    if (conflict_col_y %in% names(merged_df)) {
      merged_df[[conflict_col_y]] <- NULL
    }
  } else {
    # If names were different (e.g. "Filename" vs "ID"), merge usually keeps both.
    # But if clinical_id also got a .y suffix (rare in this logic unless names clashed elsewhere), handle it.
    # Here we typically don't need to do anything if names are different.
  }

  # 5. Final Report
  missing_clinical <- sum(is.na(merged_df[[clinical_id]]))
  # Note: if clinical_id was same as thermal_id, we removed it, so we check using the kept column or match logic
  # Simpler check: check for a known clinical column, or just report success

  message("Merge completed successfully.")
  return(merged_df)
}


#' @title Aggregate Technical Replicates for Statistical Rigor
#' @description Collapses technical replicates (e.g., multiple thermal images taken of the same animal)
#'              into a single biological data point per subject. This step is crucial for avoiding
#'              pseudoreplication in downstream statistical analyses (e.g., t-tests, ANOVA).
#'
#' @param data A data frame. Typically the output from \code{\link{compile_batch_stats}} or
#'             \code{\link{merge_clinical_data}}.
#' @param id_col String. The column name representing the unique Biological Subject ID (e.g., "MouseID", "Subject_No").
#'               Rows sharing this ID will be condensed into one.
#' @param method String. The mathematical function used for aggregation: either \code{"mean"} (default)
#'               or \code{"median"}. Median is often more robust to outliers (e.g., one blurry image).
#' @param keep_cols Vector of strings. Names of non-numeric metadata columns to preserve in the final output
#'                  (e.g., "Group", "Genotype", "Sex", "Treatment").
#'
#' @return A data frame with exactly one row per unique ID. The column order is reorganized to place
#'         ID and metadata first, followed by the aggregated thermal statistics and the \code{n_replicates} count.
#' @export
#' @examples
#' # Create a toy dataset with repeated measurements
#' df_raw <- data.frame(
#'   SampleID = rep(paste0("M", 1:3), each = 3),
#'   Group = rep(c("ND", "HFD", "ND"), each = 3),
#'   Sex = rep("M", 9),
#'   Median = runif(9, 33, 36),
#'   IQR = runif(9, 0.5, 1.5)
#' )
#'
#' df <- aggregate_replicates(
#'   data = df_raw,
#'   id_col = "SampleID",
#'   method = "median",
#'   keep_cols = c("Group", "Sex")
#' )
aggregate_replicates <- function(data, id_col, method = "mean", keep_cols = NULL) {

  # 1. Validation
  if (!is.data.frame(data)) stop("Error: 'data' must be a data frame.")
  if (!id_col %in% names(data)) stop(paste("Error: Column", id_col, "not found in data."))

  message(paste0("Aggregating replicates by '", id_col, "' using ", method, "..."))

  # 2. Identify Numeric Columns to Aggregate
  # We exclude the ID col and keep_cols from calculation
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c(id_col, keep_cols))

  if (length(numeric_cols) == 0) warning("No numeric columns found to aggregate.")

  # 3. Define Aggregation Function
  # Use na.rm = TRUE to handle occasional bad pixels/images
  agg_fun <- function(x) {
    if (method == "median") {
      return(median(x, na.rm = TRUE))
    } else {
      return(mean(x, na.rm = TRUE))
    }
  }

  # 4. Perform Aggregation
  # Split data by ID
  data_split <- split(data, data[[id_col]])

  # Loop through each mouse/subject
  result_list <- lapply(data_split, function(sub_df) {

    # A. Calculate Stats for numeric columns
    if (length(numeric_cols) > 0) {
      # apply returns a vector, transpose it to a 1-row data frame
      stats_vec <- apply(sub_df[, numeric_cols, drop = FALSE], 2, agg_fun)
      res_row <- as.data.frame(t(stats_vec))
    } else {
      res_row <- data.frame()
    }

    # B. Add ID column
    res_row[[id_col]] <- sub_df[[id_col]][1]

    # C. Add Kept Metadata Columns (Preserve the first value found)
    if (!is.null(keep_cols)) {
      for (col in keep_cols) {
        if (col %in% names(sub_df)) {
          res_row[[col]] <- sub_df[[col]][1]
        }
      }
    }

    # D. Optional: Add replicate count (how many images were merged)
    res_row$n_replicates <- nrow(sub_df)

    return(res_row)
  })

  # 5. Combine back to Data Frame
  final_df <- do.call(rbind, result_list)

  # Reorder columns: ID, Metadata, Stats...
  first_cols <- c(id_col, keep_cols, "n_replicates")
  other_cols <- setdiff(names(final_df), first_cols)
  final_df <- final_df[, c(first_cols, other_cols)]

  rownames(final_df) <- NULL
  message(paste("Aggregation complete. Reduced from", nrow(data), "rows to", nrow(final_df), "subjects."))

  return(final_df)
}

#' @title Pairwise Correlation Analysis for Thermal and External Traits
#'
#' @description Computes pairwise correlations between a set of thermal variables
#' and a set of external variables. It handles missing data gracefully, performs
#' multiple testing corrections, and returns detailed statistical summaries along
#' with correlation and p-value matrices.
#'
#' @param data A `data.frame` containing the variables to be analyzed.
#' @param thermal_vars A character vector specifying the column names of the
#'   thermal variables in `data`.
#' @param external_vars A character vector specifying the column names of the
#'   external variables in `data`.
#' @param method A character string indicating which correlation coefficient to
#'   compute. Can be either `"spearman"` (default) or `"pearson"`.
#' @param adjust_method A character string specifying the method for multiple
#'   testing correction. Options are `"holm"`, `"hochberg"`, `"hommel"`,
#'   `"bonferroni"`, `"BH"` (default), `"BY"`, `"fdr"`, or `"none"`.
#' @param use A character string indicating how to handle missing values.
#'   Can be `"complete.obs"` (default, uses only rows with complete data across
#'   all selected variables) or `"pairwise.complete.obs"` (uses complete cases
#'   on a pair-by-pair basis).
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{results}{A data frame with detailed pairwise statistics including variable names, sample size (`n`), correlation coefficient (`cor`), raw p-value (`p`), adjusted p-value (`p_adj`), direction, significance label, and any warnings.}
#'   \item{cor_matrix}{A numeric matrix of the calculated correlation coefficients.}
#'   \item{p_matrix}{A numeric matrix of the raw p-values.}
#'   \item{padj_matrix}{A numeric matrix of the adjusted p-values.}
#'   \item{method}{The correlation method used.}
#'   \item{adjust_method}{The p-value adjustment method used.}
#'   \item{use}{The missing data handling method used.}
#'   \item{use_note}{A descriptive note regarding how missing data was actually handled.}
#' }
#'
#' @importFrom stats complete.cases cor.test p.adjust sd
#' @export
#'
#' @examples
#' set.seed(1234)
#' df <- data.frame(
#'   Max = rnorm(50, 35, 2),
#'   Min = rnorm(50, 15, 3),
#'   Weight = rnorm(50, 20, 5),
#'   Glu = runif(50, 5, 8)
#' )
#'
#' # Introduce a few NA values to test missing data handling
#' df$Weight[c(2, 10)] <- NA
#'
#' # Run the correlation analysis
#' res <- correlate_thermal_traits(
#'        data = df,
#'        thermal_vars = c("Max", "Min"),
#'        external_vars = c("Weight", "Glu"),
#'        method = "spearman",
#'        adjust_method = "BH",
#'        use = "pairwise.complete.obs"
#'        )
#' # View the detailed results data frame
#' head(res$results)
#'
#' # View the correlation matrix
#' print(res$cor_matrix)
correlate_thermal_traits <- function(
    data,
    thermal_vars,
    external_vars,
    method = "spearman",
    adjust_method = "BH",
    use = "complete.obs"
) {
  # -----------------------------
  # Input validation
  # -----------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  if (missing(thermal_vars) || length(thermal_vars) == 0) {
    stop("`thermal_vars` must be a non-empty character vector.", call. = FALSE)
  }

  if (missing(external_vars) || length(external_vars) == 0) {
    stop("`external_vars` must be a non-empty character vector.", call. = FALSE)
  }

  if (!is.character(thermal_vars)) {
    stop("`thermal_vars` must be a character vector of column names.", call. = FALSE)
  }

  if (!is.character(external_vars)) {
    stop("`external_vars` must be a character vector of column names.", call. = FALSE)
  }

  method <- match.arg(method, choices = c("spearman", "pearson"))
  use <- match.arg(use, choices = c("complete.obs", "pairwise.complete.obs"))

  missing_thermal <- setdiff(thermal_vars, colnames(data))
  if (length(missing_thermal) > 0) {
    stop(
      sprintf(
        "The following `thermal_vars` are not found in `data`: %s",
        paste(missing_thermal, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing_external <- setdiff(external_vars, colnames(data))
  if (length(missing_external) > 0) {
    stop(
      sprintf(
        "The following `external_vars` are not found in `data`: %s",
        paste(missing_external, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  duplicated_vars <- intersect(thermal_vars, external_vars)
  if (length(duplicated_vars) > 0) {
    warning(
      sprintf(
        "Some variables appear in both `thermal_vars` and `external_vars`: %s",
        paste(duplicated_vars, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  selected_vars <- unique(c(thermal_vars, external_vars))
  non_numeric <- selected_vars[!vapply(data[selected_vars], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(
      sprintf(
        "All selected variables must be numeric. Non-numeric columns detected: %s",
        paste(non_numeric, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  valid_adjust_methods <- c(
    "holm", "hochberg", "hommel", "bonferroni",
    "BH", "BY", "fdr", "none"
  )
  if (!adjust_method %in% valid_adjust_methods) {
    stop(
      sprintf(
        "`adjust_method` must be one of: %s",
        paste(valid_adjust_methods, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # -----------------------------
  # Handle missing data based on 'use' parameter
  # -----------------------------
  if (use == "complete.obs") {
    data <- data[stats::complete.cases(data[selected_vars]), ]
    use_note <- "Only observations complete across all selected variables were used."
  } else {
    use_note <- "Pairwise complete cases were used for each variable pair."
  }

  # -----------------------------
  # Initialize containers
  # -----------------------------
  n_thermal <- length(thermal_vars)
  n_external <- length(external_vars)

  cor_matrix <- matrix(
    NA_real_,
    nrow = n_thermal,
    ncol = n_external,
    dimnames = list(thermal_vars, external_vars)
  )

  p_matrix <- matrix(
    NA_real_,
    nrow = n_thermal,
    ncol = n_external,
    dimnames = list(thermal_vars, external_vars)
  )

  results_list <- vector("list", length = n_thermal * n_external)
  idx <- 1L

  # -----------------------------
  # Pairwise correlation analysis
  # -----------------------------
  for (tv in thermal_vars) {
    for (ev in external_vars) {
      x <- data[[tv]]
      y <- data[[ev]]

      # Pairwise complete-case handling (safe for both methods after pre-filtering)
      keep <- stats::complete.cases(x, y)

      x_sub <- x[keep]
      y_sub <- y[keep]
      n_used <- length(x_sub)

      cor_est <- NA_real_
      p_val <- NA_real_
      warn_msg <- NA_character_

      if (n_used < 3) {
        warn_msg <- "Too few complete observations (< 3) for correlation testing."
      } else if (stats::sd(x_sub) == 0 || stats::sd(y_sub) == 0) {
        warn_msg <- "At least one variable has zero variance; correlation is undefined."
      } else {
        test_res <- tryCatch(
          stats::cor.test(
            x = x_sub,
            y = y_sub,
            method = method,
            exact = if (method == "spearman") FALSE else NULL
          ),
          error = function(e) e
        )

        if (inherits(test_res, "error")) {
          warn_msg <- conditionMessage(test_res)
        } else {
          cor_est <- unname(test_res$estimate)
          p_val <- unname(test_res$p.value)
        }
      }

      cor_matrix[tv, ev] <- cor_est
      p_matrix[tv, ev] <- p_val

      direction <- if (is.na(cor_est)) {
        NA_character_
      } else if (cor_est > 0) {
        "positive"
      } else if (cor_est < 0) {
        "negative"
      } else {
        "zero"
      }

      results_list[[idx]] <- data.frame(
        thermal_var = tv,
        external_var = ev,
        method = method,
        n = n_used,
        cor = cor_est,
        p = p_val,
        direction = direction,
        warning = warn_msg,
        stringsAsFactors = FALSE
      )

      idx <- idx + 1L
    }
  }

  results_df <- do.call(rbind, results_list)

  # -----------------------------
  # Multiple testing correction
  # -----------------------------
  valid_p <- !is.na(results_df$p)
  results_df$p_adj <- NA_real_

  if (any(valid_p)) {
    results_df$p_adj[valid_p] <- stats::p.adjust(
      p = results_df$p[valid_p],
      method = adjust_method
    )
  }

  # Significance labels
  results_df$significance <- vapply(
    results_df$p_adj,
    FUN.VALUE = character(1),
    FUN = function(pv) {
      if (is.na(pv)) {
        ""
      } else if (pv < 0.001) {
        "***"
      } else if (pv < 0.01) {
        "**"
      } else if (pv < 0.05) {
        "*"
      } else {
        "ns"
      }
    }
  )

  # Build adjusted p-value matrix efficiently
  padj_matrix <- matrix(
    results_df$p_adj,
    nrow = n_thermal,
    ncol = n_external,
    byrow = TRUE,
    dimnames = list(thermal_vars, external_vars)
  )

  results_df <- results_df[, c(
    "thermal_var", "external_var", "method", "n",
    "cor", "p", "p_adj", "direction", "significance", "warning"
  )]

  rownames(results_df) <- NULL

  out <- list(
    results = results_df,
    cor_matrix = cor_matrix,
    p_matrix = p_matrix,
    padj_matrix = padj_matrix,
    method = method,
    adjust_method = adjust_method,
    use = use,
    use_note = use_note
  )
  class(out) <- c("thermal_correlation_result", class(out))
  return(out)
}

#' @title Assess replicate consistency and repeatability
#' @description Evaluates the statistical consistency and reliability of repeated thermal measurements.
#' The function calculates the Intraclass Correlation Coefficient (ICC) and performs
#' variance decomposition using Linear Mixed Models (LMM) to assess the proportion of
#' variance attributable to between-subject differences versus measurement error.
#'
#' @param data A data.frame containing the thermal metrics and identifier columns.
#' @param id_col A single character string specifying the column name for subject IDs.
#' @param metrics A character vector specifying the column names of the thermal metrics to be assessed.
#' @param replicate_col A single character string specifying the column name for the replicate index. Default is \code{NULL}. If \code{NULL}, replicate indices are inferred from row order.
#' @param sort_cols A character vector specifying the columns to sort the data by before inferring replicate indices. Default is \code{NULL}.
#' @param methods A character vector specifying the statistical methods to apply. Valid options include \code{"icc"}, \code{"variance"}, and \code{"lmm"}. Default is \code{"icc"}.
#' @param return_models Logical. If \code{TRUE} and \code{"lmm"} is in \code{methods}, the fitted linear mixed models are returned. Default is \code{FALSE}.
#' @param quiet Logical. If \code{TRUE}, console messages are suppressed. Default is \code{FALSE}.
#'
#' @return An object of class \code{BioThermR_replicate_assessment}, which is a list containing:
#' \item{settings}{A list of the input parameters used for the assessment.}
#' \item{icc}{A data.frame containing ICC results (if \code{"icc"} is selected).}
#' \item{variance}{A data.frame containing variance decomposition results (if \code{"variance"} is selected).}
#' \item{lmm}{A data.frame containing fixed effects and variance components from LMM (if \code{"lmm"} is selected).}
#' \item{models}{A list of fitted \code{lmer} model objects (if \code{return_models = TRUE}).}
#'
#' @details
#' This function requires the \code{psych} package for ICC calculations and the \code{lme4} package for variance decomposition and LMM fitting. Ensure these packages are installed before selecting the respective methods.
#' @importFrom stats ave reshape
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming df is a data.frame with columns: Sample, Rep, Mean
#' res <- assess_replicates(
#'   data = df,
#'   id_col = "Sample",
#'   metrics = "Mean",
#'   replicate_col = "Rep",
#'   methods = c("icc", "variance")
#' )
#' print(res$icc)
#' }
assess_replicates <- function(data,
                              id_col,
                              metrics,
                              replicate_col = NULL,
                              sort_cols = NULL,
                              methods = "icc",
                              return_models = FALSE,
                              quiet = FALSE) {

  # ----------------------------
  # Basic checks
  # ----------------------------
  if (!is.data.frame(data)) {
    stop("Error: 'data' must be a data.frame.")
  }

  if (!is.character(id_col) || length(id_col) != 1) {
    stop("Error: 'id_col' must be a single character string.")
  }

  if (!id_col %in% colnames(data)) {
    stop(sprintf("Error: id_col '%s' not found in data.", id_col))
  }

  if (!is.character(metrics) || length(metrics) < 1) {
    stop("Error: 'metrics' must be a non-empty character vector.")
  }

  missing_metrics <- setdiff(metrics, colnames(data))
  if (length(missing_metrics) > 0) {
    stop(sprintf(
      "Error: The following metrics were not found in data: %s",
      paste(missing_metrics, collapse = ", ")
    ))
  }

  if (!is.null(replicate_col) && !replicate_col %in% colnames(data)) {
    stop(sprintf("Error: replicate_col '%s' not found in data.", replicate_col))
  }

  if (!is.null(sort_cols)) {
    missing_sort <- setdiff(sort_cols, colnames(data))
    if (length(missing_sort) > 0) {
      stop(sprintf(
        "Error: The following sort_cols were not found in data: %s",
        paste(missing_sort, collapse = ", ")
      ))
    }
  }

  valid_methods <- c("icc", "variance", "lmm")
  methods <- unique(methods)

  if (!all(methods %in% valid_methods)) {
    stop(sprintf(
      "Error: 'methods' must be chosen from: %s",
      paste(valid_methods, collapse = ", ")
    ))
  }

  if ("icc" %in% methods && !requireNamespace("psych", quietly = TRUE)) {
    stop("Error: Package 'psych' is required for method = 'icc'. Please install it first.")
  }

  if (any(c("variance", "lmm") %in% methods) &&
      !requireNamespace("lme4", quietly = TRUE)) {
    stop("Error: Package 'lme4' is required for methods 'variance' or 'lmm'. Please install it first.")
  }

  # ----------------------------
  # Helper: prepare one metric
  # ----------------------------
  prepare_metric_data <- function(metric_name) {
    cols_needed <- unique(c(id_col, replicate_col, sort_cols, metric_name))
    df <- data[, cols_needed, drop = FALSE]

    names(df)[names(df) == id_col] <- ".id"
    if (!is.null(replicate_col)) {
      names(df)[names(df) == replicate_col] <- ".rep"
    }
    names(df)[names(df) == metric_name] <- ".value"

    df <- df[!is.na(df$.id) & !is.na(df$.value), , drop = FALSE]

    if (nrow(df) == 0) {
      stop(sprintf("Error: No valid rows available for metric '%s'.", metric_name))
    }

    if (!is.null(sort_cols)) {
      sort_cols_internal <- sort_cols
      sort_cols_internal[sort_cols_internal == id_col] <- ".id"
      if (!is.null(replicate_col)) {
        sort_cols_internal[sort_cols_internal == replicate_col] <- ".rep"
      }

      ord <- do.call(order, df[, sort_cols_internal, drop = FALSE])
      df <- df[ord, , drop = FALSE]
    }

    if (is.null(replicate_col)) {
      df$.rep <- ave(seq_len(nrow(df)), df$.id, FUN = seq_along)
    }

    dup_idx <- duplicated(df[, c(".id", ".rep")])
    if (any(dup_idx)) {
      stop(sprintf(
        "Error: Duplicate subject-replicate combinations found for metric '%s'.",
        metric_name
      ))
    }

    wide <- reshape(
      df[, c(".id", ".rep", ".value"), drop = FALSE],
      idvar = ".id",
      timevar = ".rep",
      direction = "wide"
    )

    rep_cols <- grep("^\\.value\\.", colnames(wide), value = TRUE)

    if (length(rep_cols) < 2) {
      stop(sprintf(
        "Error: Metric '%s' must have at least two replicate columns for assessment.",
        metric_name
      ))
    }

    rep_ids <- sub("^\\.value\\.", "", rep_cols)
    rep_order <- suppressWarnings(order(as.numeric(rep_ids), rep_ids))
    rep_cols <- rep_cols[rep_order]

    mat <- as.matrix(wide[, rep_cols, drop = FALSE])
    rownames(mat) <- wide$.id

    list(
      long = df,
      wide = wide,
      mat = mat,
      rep_cols = rep_cols
    )
  }

  # ----------------------------
  # Containers
  # ----------------------------
  icc_out <- list()
  variance_out <- list()
  lmm_out <- list()
  model_out <- list()

  # ----------------------------
  # Main loop over metrics
  # ----------------------------
  for (metric_name in metrics) {
    if (!quiet) {
      message(sprintf("Assessing replicates for metric: %s", metric_name))
    }

    prep <- prepare_metric_data(metric_name)
    df_long <- prep$long
    mat <- prep$mat

    # ------------------------
    # ICC
    # ------------------------
    if ("icc" %in% methods) {
      icc_res <- psych::ICC(mat)
      icc_tbl <- icc_res$results

      # Preserve row labels if present
      icc_rowname <- rownames(icc_tbl)
      if (is.null(icc_rowname)) {
        icc_rowname <- seq_len(nrow(icc_tbl))
      }

      icc_tbl$ICC_label <- icc_rowname
      icc_tbl$Metric <- metric_name

      # Standardize column names
      colnames(icc_tbl) <- gsub(" ", "_", colnames(icc_tbl), fixed = TRUE)
      colnames(icc_tbl)[colnames(icc_tbl) == "lower_bound"] <- "CI_lower"
      colnames(icc_tbl)[colnames(icc_tbl) == "upper_bound"] <- "CI_upper"

      # Reorder columns
      keep_cols <- c(
        "Metric", "ICC_label", "type", "ICC", "F", "df1", "df2", "p",
        "CI_lower", "CI_upper"
      )
      keep_cols <- keep_cols[keep_cols %in% colnames(icc_tbl)]

      icc_tbl <- icc_tbl[, keep_cols, drop = FALSE]
      rownames(icc_tbl) <- NULL

      icc_out[[metric_name]] <- icc_tbl
    }

    # ------------------------
    # Variance decomposition / LMM
    # ------------------------
    if (any(c("variance", "lmm") %in% methods)) {
      df_model <- df_long[, c(".id", ".value"), drop = FALSE]
      df_model$.id <- as.factor(df_model$.id)

      fit <- lme4::lmer(.value ~ 1 + (1 | .id), data = df_model, REML = TRUE)

      vc <- as.data.frame(lme4::VarCorr(fit))
      between_var <- vc$vcov[vc$grp == ".id"]
      residual_var <- vc$vcov[vc$grp == "Residual"]

      if (length(between_var) == 0) between_var <- NA_real_
      if (length(residual_var) == 0) residual_var <- NA_real_

      total_var <- between_var + residual_var
      repeatability <- between_var / total_var

      if ("variance" %in% methods) {
        variance_tbl <- data.frame(
          Metric = metric_name,
          Between_subject_variance = between_var,
          Residual_variance = residual_var,
          Total_variance = total_var,
          Repeatability = repeatability,
          N_subjects = length(unique(df_model$.id)),
          N_observations = nrow(df_model),
          stringsAsFactors = FALSE
        )

        variance_out[[metric_name]] <- variance_tbl
      }

      if ("lmm" %in% methods) {
        fixef_est <- unname(lme4::fixef(fit)[1])
        fixef_se <- sqrt(diag(as.matrix(stats::vcov(fit))))[1]

        lmm_tbl <- data.frame(
          Metric = metric_name,
          Intercept = fixef_est,
          Intercept_SE = fixef_se,
          Between_subject_variance = between_var,
          Residual_variance = residual_var,
          Repeatability = repeatability,
          N_subjects = length(unique(df_model$.id)),
          N_observations = nrow(df_model),
          stringsAsFactors = FALSE
        )

        lmm_out[[metric_name]] <- lmm_tbl

        if (isTRUE(return_models)) {
          model_out[[metric_name]] <- fit
        }
      }
    }
  }

  # ----------------------------
  # Assemble output
  # ----------------------------
  out <- list(
    settings = list(
      id_col = id_col,
      replicate_col = replicate_col,
      sort_cols = sort_cols,
      metrics = metrics,
      methods = methods
    )
  )

  if ("icc" %in% methods) {
    out$icc <- do.call(rbind, icc_out)
    rownames(out$icc) <- NULL
  }

  if ("variance" %in% methods) {
    out$variance <- do.call(rbind, variance_out)
    rownames(out$variance) <- NULL
  }

  if ("lmm" %in% methods) {
    out$lmm <- do.call(rbind, lmm_out)
    rownames(out$lmm) <- NULL
  }

  if (isTRUE(return_models) && length(model_out) > 0) {
    out$models <- model_out
  }

  class(out) <- c("BioThermR_replicate_assessment", class(out))
  return(out)
}

