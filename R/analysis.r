#' @title Calculate Comprehensive Thermal Statistics
#' @description Computes a detailed set of summary statistics from the thermal matrix.
#'              Metrics include central tendency (Mean, Median, Peak Density), dispersion (SD, IQR, CV),
#'              and range (Min, Max, Quantiles). NA values (background) are automatically excluded.
#'
#' @details The function calculates the following metrics:
#'          \itemize{
#'            \item \strong{Min/Max:} Extremities of the temperature distribution.
#'            \item \strong{Mean/Median:} Measures of central tendency.
#'            \item \strong{SD (Standard Deviation):} Absolute measure of spread.
#'            \item \strong{IQR (Interquartile Range):} Robust measure of spread (Q75 - Q25).
#'            \item \strong{CV (Coefficient of Variation):} Relative measure of spread (SD / Mean), useful for assessing thermal heterogeneity.
#'            \item \strong{Peak_Density:} The temperature value corresponding to the peak of the kernel density estimate (Mode).
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
      Metric = c("Min", "Max", "Mean", "Median", "SD", "Q25", "Q75", "Peak_Density"),
      Value = NA
    )
  } else {
    # Calculate Density Peak (Mode)
    d <- density(vals)
    peak_val <- d$x[which.max(d$y)]

    # Create Statistics Data Frame
    stats_df <- data.frame(
      Metric = c("Min", "Max", "Mean", "Median", "SD", "Q25", "Q75", "IQR", "CV","Peak_Density"),
      Value = c(
        min(vals),
        max(vals),
        mean(vals),
        median(vals),
        sd(vals),
        quantile(vals, 0.25),
        quantile(vals, 0.75),
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

