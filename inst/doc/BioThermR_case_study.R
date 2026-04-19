## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## ----install_github, eval=FALSE-----------------------------------------------
# # Install the full version with data from GitHub
# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("RightSZ/BioThermR")

## ----install_dependencies,eval=FALSE------------------------------------------
# # Install dependencies
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("EBImage")

## ----lib_packages-------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(BioThermR)

## ----load_data----------------------------------------------------------------
# 1. Get the path to the 30 images included with the package
data_folder <- system.file("extdata", package = "BioThermR")

# 2. Read using the batch reading function
obj_list <- read_thermal_batch(data_folder)
length(obj_list)

# 3. Use the lapply function to batch call roi_segment_ebimage for segmentation
# Otsu's method is used by default
obj_list <- lapply(obj_list, roi_segment_ebimage)

## ----viz_raw_heatmap----------------------------------------------------------
# Display the thermal heatmap of one BioThermR object
p1 <- plot_thermal_heatmap(obj_list[[1]], use_processed = FALSE)
p1

## ----viz_roi_heatmap----------------------------------------------------------
# Display the segmented image of one BioThermR object
p2 <- plot_thermal_heatmap(obj_list[[1]], use_processed = TRUE)
p2

## ----viz_3d-------------------------------------------------------------------
# Display the 3d image of one segmented BioThermR object
plot_thermal_3d(obj_list[[1]])

## ----viz_density--------------------------------------------------------------
# Display the ROI density plot of one BioThermR object
p3 <- plot_thermal_density(obj_list[[1]])
p3

## ----viz_montage--------------------------------------------------------------
# Display the overall display montage plot
p4 <- plot_thermal_montage(obj_list, ncol = 5, padding = 2, text_size = 3)
p4

## ----viz_cloud----------------------------------------------------------------
# Display the overall display cloud plot
p5 <- plot_thermal_cloud(obj_list) 
p5

## ----viz_ggpubr---------------------------------------------------------------
ggpubr::ggarrange(p1,p2,p3,p4,p5,labels = c("A","B","C","D","E"))

## ----statistical_Extraction---------------------------------------------------
# Batch calculate statistics for each BioThermR object
obj_list <- lapply(obj_list, analyze_thermal_stats)

# Compile statistics into a data frame
df <- compile_batch_stats(obj_list)

# Read grouping data and merge
pd_path <- system.file("extdata", "group.csv", package = "BioThermR")
pd <- read.csv(pd_path)
df <- merge_clinical_data(df,pd)

# Aggregate technical replicates (e.g., 3 photos per mouse)
df_new <- aggregate_replicates(df,
                               method = "mean", # use mean method
                               id_col = "Sample", # merge by sample
                               keep_cols = c("Group")
                               )
head(df_new)

## ----comparison---------------------------------------------------------------
# Display grouping data plots, grouped by Group
my_comparisons <- list( c("ND4", "ND"))
# Mean
s1 <- viz_thermal_boxplot(data = df_new,
                          y_var = "Mean",
                          x_var = "Group") +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "t.test",
                             label = "p.signif",
                             hide.ns = FALSE
                             )

# Median
s2 <- viz_thermal_boxplot(data = df_new,
                          y_var = "Median",
                          x_var = "Group")+
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "t.test",
                             label = "p.signif",
                             hide.ns = FALSE
  )

# IQR
s3 <- viz_thermal_boxplot(data = df_new,
                          y_var = "IQR",
                          x_var = "Group")+
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "t.test",
                             label = "p.signif",
                             hide.ns = FALSE
  )
# Peak Density
s4 <- viz_thermal_boxplot(data = df_new,
                          y_var = "Peak_Density",
                          x_var = "Group")+
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             method = "t.test",
                             label = "p.signif",
                             hide.ns = FALSE
  )

## ----correlation analysis-----------------------------------------------------
# In this example, selected thermal metrics are tested against body weight and blood glucose using Spearman correlation. The results are visualized as a correlation heatmap and a representative scatter plot.
res <- correlate_thermal_traits(
  df_new,
  thermal_vars =  c("Max","Min","Mean","Median","IQR","Peak_Density"),
  external_vars = c("Weight","Glucose"),
  method = "spearman")

# Correlation heatmap summarizing the associations between thermal metrics and external traits.
s5 <- viz_cor_heatmap(res)

# Representative scatter plot illustrating the relationship between one selected thermal feature (IQR) and an external trait (Weight).
s6 <- viz_cor_scatter(df_new,
                x_col = "IQR",
                x_label = "Thermal Metric: IQR",
                y_col = "Weight"
                )

## ----Figure Assembly----------------------------------------------------------
ggarrange(
  s1, s2, s3,
  s4, s5, s6,
  ncol = 3, nrow = 2,
  labels = c("A", "B", "C", "D", "E", "F"),
  widths = c(1, 1, 1),
  heights = c(1, 1.1)
)

