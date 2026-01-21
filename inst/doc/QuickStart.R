## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(BioThermR)
library(ggplot2)

## ----import_real--------------------------------------------------------------
# 1. Locate the example file provided with BioThermR
raw_path <- system.file("extdata", "C05.raw", package = "BioThermR")

# 2. Read the thermal data using read_thermal_raw
obj <- read_thermal_raw(raw_path) 

# 3. Quick sanity check
# This confirms the data is loaded and has the correct dimensions
print(paste("Filename:", obj$meta$filename))
print(paste("Dimensions:", obj$meta$dims[1], "x", obj$meta$dims[2]))

## ----visualize_raw_data-------------------------------------------------------
plot_thermal_heatmap(obj, use_processed = FALSE) + 
  ggtitle("Raw Thermal Image")

## ----segment_auto-------------------------------------------------------------
# Apply automated segmentation
# method = "otsu": Automatically calculates the best threshold
# keep_largest = TRUE: Removes small noise artifacts, keeping only the animal
obj <- roi_segment_ebimage(obj, method = "otsu", keep_largest = TRUE)

# Visualize the result (Processed Data)
# Note how the background is now clean (NA)
plot_thermal_heatmap(obj, use_processed = TRUE) + 
  ggtitle("Auto-Segmented Image")

## ----calculate_statistics-----------------------------------------------------
# Compute Min, Max, Mean, Median, etc.
obj <- analyze_thermal_stats(obj)

# Show the results
print(obj$stats)

## ----batch_process, fig.width=8, fig.height=6---------------------------------
# 1. Get the path to the folder containing the 30 raw files
data_folder <- system.file("extdata", package = "BioThermR")

print(paste("Reading batch from:", data_folder))

# 2. Read the entire batch
# Note: We use pattern = ".raw" to ensure we only load the raw thermal files
batch_list <- read_thermal_batch(data_folder, pattern = "\\.raw$")

# 3. Batch Segmentation (Automated)
# We use lapply to apply the 'roi_segment_ebimage' function to every image in the list
# This automatically removes the background for all 30 mice
batch_list_clean <- lapply(batch_list, roi_segment_ebimage)

# 4. Visualization A: Gap-Free Montage
p1 <- plot_thermal_montage(batch_list_clean, ncol = 5, padding = 2, text_size = 3) + 
  ggtitle("Montage View: 30 Mice")
print(p1)

# 5. Visualization B: Thermal Cloud
p2 <- plot_thermal_cloud(batch_list_clean, spread_factor = 1.5, jitter_factor = 0.5, show_labels = TRUE) + 
  ggtitle("Thermal Cloud: Population Overview")
print(p2)

## ----batch_stats--------------------------------------------------------------
# 1. Calculate stats for each image in the list
# This adds a 'stats' slot to each BioThermR object
batch_list_stats <- lapply(batch_list_clean, analyze_thermal_stats)

# 2. Compile into a tidy data frame
# Rows = Images, Columns = Metrics
df_results <- compile_batch_stats(batch_list_stats)

# 3. View the first few rows
head(df_results)

