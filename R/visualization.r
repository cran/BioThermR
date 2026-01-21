#' @title Visualize Thermal Matrix as 2D Heatmap
#' @description Generates a high-quality raster plot of the thermal data using the 'ggplot2' framework.
#'              This function allows for quick visualization of raw or processed matrices with
#'              customizable perceptually uniform color scales (viridis).
#'
#' @details The function performs the following steps:
#'          \itemize{
#'            \item Converts the thermal matrix into a long-format data frame suitable for ggplot.
#'            \item Renders the image using \code{geom_raster}.
#'            \item Maps temperature to color using the specified 'viridis' palette.
#'            \item Ensures the aspect ratio is preserved (\code{coord_fixed}) so the image does not appear distorted.
#'            \item Sets \code{NA} values (masked background) to transparent.
#'          }
#'          Since the output is a standard ggplot object, layers can be added subsequently (e.g., new titles or annotations).
#'
#' @param img_obj A 'BioThermR' object.
#' @param use_processed Logical. If \code{TRUE} (default), plots the 'processed' matrix (showing masking effects).
#'                      If \code{FALSE}, plots the original 'raw' data.
#' @param palette String. The color map option from the 'viridis' package.
#'                Options include: \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"viridis"}, \code{"cividis"}.
#'                Default is \code{"inferno"}.
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @export
#' @examples
#' \donttest{
#' # Load raw data
#' img_obj <- system.file("extdata", "C05.raw", package = "BioThermR")
#' img <- read_thermal_raw(img_obj)
#'
#' # Apply automated segmentation
#' img <- roi_segment_ebimage(img, keep_largest = TRUE)
#'
#' # Standard plot
#' plot_thermal_heatmap(img)
#' }
plot_thermal_heatmap <- function(img_obj, use_processed = TRUE, palette = "inferno") {

  # Check object class
  if (!inherits(img_obj, "BioThermR")) {
    stop("Error: Input must be a 'BioThermR' object.")
  }

  mat <- if (use_processed) img_obj$processed else img_obj$raw

  if (is.null(mat)) {
    stop("Error: Selected matrix is empty.")
  }

  df <- expand.grid(
    row = seq_len(nrow(mat)),
    col = seq_len(ncol(mat))
  )
  df$val <- c(mat)

  p <- ggplot(df, aes(x = col, y = row, fill = val)) +
    geom_raster() +
    scale_fill_viridis_c(option = palette, name = "Temp (\u00B0C)", na.value = "transparent") +
    coord_fixed() +
    theme_void()+
    labs(
      title = paste("Thermal Heatmap:", img_obj$meta$filename),
      subtitle = paste("Source:", if(use_processed) "Processed Data" else "Raw Data")
    )

  return(p)
}

#' @title Visualize Temperature Distribution (Density Plot)
#' @description Generates a probability density plot of the temperature values within the image.
#'              This visualization is critical for assessing the homogeneity of the subject's temperature
#'              and identifying potential artifacts (e.g., bimodal distributions often indicate poor background removal).
#'
#' @details The function computes the kernel density estimate of the valid pixels (ignoring NAs).
#'          It can optionally annotate key statistical landmarks:
#'          \itemize{
#'            \item \strong{Peak:} The mode of the distribution (most frequent temperature).
#'            \item \strong{Max/Min:} The hottest and coldest points in the ROI.
#'          }
#'          Text labels are automatically repelled using 'ggrepel' to ensure they do not overlap.
#'
#' @param img_obj A 'BioThermR' object.
#' @param use_processed Logical. If \code{TRUE} (default), uses the 'processed' matrix (masked data).
#'                      If \code{FALSE}, uses the 'raw' matrix.
#' @param show_peak Logical. If \code{TRUE}, highlights and labels the peak density value (Mode). Default is \code{TRUE}.
#' @param show_max Logical. If \code{TRUE}, highlights and labels the maximum temperature value. Default is \code{TRUE}.
#' @param show_min Logical. If \code{TRUE}, highlights and labels the minimum temperature value. Default is \code{TRUE}.
#' @param digits Integer. Number of decimal places to round the labels to. Default is 2.
#' @param color String. Fill color for the density area curve. Default is "skyblue".
#' @param point_size Numeric. Size of the points marking Peak/Min/Max. Default is 2.
#' @param point_color String. Color of the points marking Peak/Min/Max. Default is "red".
#' @param point_label_color String. Color of the text labels. Default is "black".
#' @param point_label_size Numeric. Size of the text labels. Default is 2.
#'
#' @return A \code{ggplot} object. Layers can be added subsequently.
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom stats density
#' @export
#' @examples
#' \donttest{
#' # Load raw data
#' img_obj <- system.file("extdata", "C05.raw", package = "BioThermR")
#' img <- read_thermal_raw(img_obj)
#'
#' # Apply automated segmentation
#' img <- roi_segment_ebimage(img, keep_largest = TRUE)
#'
#' # Density plot
#' plot_thermal_density(img)
#' }
plot_thermal_density <- function(img_obj, use_processed = TRUE, show_peak = TRUE,
                                 show_max = TRUE, show_min = TRUE, digits = 2,color = "skyblue",
                                 point_size = 2, point_color = "red", point_label_color = "black",
                                 point_label_size = 2) {

  if (!inherits(img_obj, "BioThermR")) {
    stop("Error: Input must be a 'BioThermR' object.")
  }

  mat <- if (use_processed) img_obj$processed else img_obj$raw
  val <- as.vector(mat)
  val <- val[!is.na(val)]

  if (length(val) == 0) {
    stop("Error: No valid temperature data found (matrix might be all NA).")
  }

  df <- data.frame(temp = val)

  p <- ggplot(df, aes(x = temp)) +
    geom_density(fill = color, alpha = 0.5) +
    labs(
      title = "Temperature Density Distribution",
      x = "Temperature (\u00B0C)",
      y = "Density"
    ) +
    theme_minimal()

  label_data <- data.frame(x = numeric(), y = numeric(), label = character())
  dens <- density(val)
  if (show_peak) {
    peak_idx <- which.max(dens$y)
    label_data <- rbind(label_data, data.frame(
      x = dens$x[peak_idx],
      y = dens$y[peak_idx],
      label = paste0("Peak: ", round(dens$x[peak_idx], digits))
    ))
  }
  if (show_max) {
    max_val <- max(val)
    max_idx <- which.max(dens$x)
    label_data <- rbind(label_data, data.frame(
      x = max_val,
      y = dens$y[max_idx],
      label = paste0("Max: ", round(max_val, digits))
    ))
  }
  if (show_min) {
    min_val <- min(val)
    min_idx <- which.min(dens$x)
    label_data <- rbind(label_data, data.frame(
      x = min_val,
      y = dens$y[min_idx],
      label = paste0("Min: ", round(min_val, digits))
    ))
  }
  if (nrow(label_data) > 0) {
    p <- p +
      geom_point(data = label_data, aes(x = x, y = y), size = point_size, color = point_color) +
      ggrepel::geom_label_repel(
        data = label_data,
        aes(x = x, y = y, label = label),
        size = point_label_size,
        color = point_label_color
      )
  }

  return(p)
}

#' @title Interactive 3D Thermal Surface Plot
#' @description Generates a rotatable, interactive 3D surface plot using the 'plotly' engine.
#'              This visualization maps temperature to the Z-axis, allowing users to intuitively
#'              explore the thermal topology, gradients, and intensity of hotspots.
#'
#' @details 3D visualization is particularly powerful for:
#'          \itemize{
#'            \item \strong{Quality Control:} Quickly identifying noise spikes or "cold" artifacts that flat heatmaps might hide.
#'            \item \strong{Gradient Analysis:} Visualizing how heat dissipates from a central source (e.g., a tumor or inflammation site).
#'            \item \strong{Presentation:} Creating engaging, interactive figures for HTML reports or Shiny dashboards.
#'          }
#'          The output is an HTML widget that allows zooming, panning, and hovering to see specific pixel values.
#'
#' @param img_obj A 'BioThermR' object.
#' @param use_processed Logical. If \code{TRUE} (default), uses the 'processed' matrix (where background is likely \code{NA}).
#'                      If \code{FALSE}, uses the 'raw' sensor data.
#'
#' @return A \code{plotly} object (HTML widget).
#' @importFrom plotly plot_ly add_surface layout
#' @export
#' @examples
#' \donttest{
#' # Load raw data
#' img_obj <- system.file("extdata", "C05.raw", package = "BioThermR")
#' img <- read_thermal_raw(img_obj)
#'
#' # Apply automated segmentation
#' img <- roi_segment_ebimage(img, keep_largest = TRUE)
#'
#' # 3d plot
#' plot_thermal_3d(img)
#' }
plot_thermal_3d <- function(img_obj, use_processed = TRUE) {

  if (!inherits(img_obj, "BioThermR")) {
    stop("Error: Input must be a 'BioThermR' object.")
  }

  mat <- if (use_processed) img_obj$processed else img_obj$raw

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Error: Package 'plotly' is required for this function. Please install it.")
  }

  p <- plotly::plot_ly(z = ~mat) |>
    plotly::add_surface() |>
    plotly::layout(
      title = list(text = paste("3D Thermal View:", img_obj$meta$filename)),
      scene = list(
        xaxis = list(title = "Width (px)"),
        yaxis = list(title = "Height (px)"),
        zaxis = list(title = "Temp (\u00B0C)")
      )
    )

  return(p)
}


#' @title Generate Publication-Ready Comparative Barplots
#' @description Creates a high-quality bar plot to compare thermal metrics across experimental groups.
#'              The visualization includes bars representing the mean, customizable error bars (SD or SE),
#'              and overlaid individual data points to show biological variation.
#'
#' @details This function is designed to produce figures that are immediately suitable for scientific manuscripts.
#'          Key features include:
#'          \itemize{
#'            \item \strong{Automatic Statistics:} Calculates Mean and SD/SE internally using \code{stat_summary}.
#'            \item \strong{Smart Coloring:} Supports scientific palettes ("npg", "jco"). If the number of groups exceeds
#'                  the palette size, it automatically interpolates to generate distinct colors.
#'            \item \strong{Layout:} Uses \code{theme_classic()} and automatically expands the Y-axis limit
#'                  by 15\% to ensure error bars and significance annotations fit comfortably.
#'          }
#'
#' @param data Data frame. The merged dataset (e.g., output from \code{\link{aggregate_replicates}}).
#' @param y_var String. The name of the numeric column to plot (e.g., "Max_Temp", "Mean_Temp").
#' @param x_var String. The name of the categorical column for the X-axis groupings (e.g., "Treatment", "Genotype").
#' @param fill_var String. The name of the variable used for fill colors. Default is \code{NULL}, which uses \code{x_var}.
#' @param error_bar String. The type of error bar to display. Options:
#'                  \itemize{
#'                    \item \code{"mean_sd"}: Mean +/- Standard Deviation (shows spread of data).
#'                    \item \code{"mean_se"}: Mean +/- Standard Error of the Mean (shows precision of the mean).
#'                  }
#' @param add_points Logical. If \code{TRUE} (default), overlays individual data points using \code{geom_jitter}.
#'                   This is highly recommended for small sample sizes (n < 20) to maintain transparency.
#' @param point_size Numeric. The size of the individual jitter points. Default is 1.5.
#' @param point_alpha Numeric. The transparency of the jitter points (0 = transparent, 1 = opaque).
#'                    Default is 0.6 to handle overlapping points.
#' @param palette String or Vector.
#'                \itemize{
#'                  \item If a string: Pre-defined scientific palettes (\code{"npg"} for Nature Publishing Group, \code{"jco"} for Journal of Clinical Oncology).
#'                  \item If a character vector: A custom list of hex codes (e.g., \code{c("#FF0000", "#0000FF")}).
#'                }
#'
#' @return A \code{ggplot} object. Can be further customized with standard ggplot2 functions.
#' @import ggplot2
#' @importFrom ggsci pal_npg pal_jco
#' @importFrom grDevices colorRampPalette
#' @export
#' @examples
#' df_bio <- data.frame(
#'   Treatment = rep(c("ND", "HFD"), each = 5),
#'   Mean = c(runif(5, 33, 35), runif(5, 34, 36))
#' )
#'
#' # Boxplot with individual points
#' p <- viz_thermal_barplot(df_bio,y_var="Mean",x_var="Treatment",error_bar = "mean_se")
#' p
viz_thermal_barplot <- function(data, y_var, x_var,
                                fill_var = NULL,
                                error_bar = "mean_sd",
                                add_points = TRUE,
                                point_size = 1.5,
                                point_alpha = 0.6,
                                palette = "npg") {

  # 1. Validation
  if (!is.data.frame(data)) stop("Error: 'data' must be a data frame.")
  if (!y_var %in% names(data)) stop(paste("Error: Column", y_var, "not found in data."))
  if (!x_var %in% names(data)) stop(paste("Error: Column", x_var, "not found in data."))

  if (is.null(fill_var)) fill_var <- x_var
  data[[x_var]] <- as.factor(data[[x_var]])

  # 2. Color Logic
  n_colors_needed <- length(unique(data[[fill_var]]))
  if (length(palette) > 1) {
    my_colors <- palette
  } else if (palette == "npg" && n_colors_needed <= 9) {
    my_colors <- ggsci::pal_npg()(n_colors_needed)
  } else if (palette == "jco" && n_colors_needed <= 10) {
    my_colors <- ggsci::pal_jco()(n_colors_needed)
  } else {
    base_colors <- ggsci::pal_npg()(9)
    my_colors <- colorRampPalette(base_colors)(n_colors_needed)
  }

  # 3. Base Plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[fill_var]]))

  # 4. Add Bar (Mean)
  p <- p + stat_summary(fun = mean, geom = "bar",
                        position = position_dodge(width = 0.8),
                        width = 0.7, color = "black", alpha = 0.8)

  # 5. Add Error Bars
  p <- p + stat_summary(fun.data = error_bar, geom = "errorbar",
                        position = position_dodge(width = 0.8),
                        width = 0.2, linewidth = 0.8)

  # 6. Add Individual Points (Customizable)
  if (add_points) {
    p <- p + geom_jitter(
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      size = point_size,
      alpha = point_alpha,
      shape = 21
    )
  }

  # 7. Aesthetics
  p <- p +
    scale_fill_manual(values = my_colors) +
    labs(y = paste("Thermal Metric:", y_var), x = x_var) +
    theme_classic()

  p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

  return(p)
}


#' @title Generate Publication-Ready Comparative Boxplots
#' @description Creates a high-quality box-and-whisker plot to visualize the distribution of thermal metrics across groups.
#'              This function is ideal for displaying median values, quartiles, and range, while optionally overlaying
#'              individual data points to reveal the underlying sample distribution.
#'
#' @details This function includes several automated optimizations for scientific reporting:
#'          \itemize{
#'            \item \strong{Smart Outlier Handling:} If \code{add_points} is \code{TRUE}, the function automatically hides
#'                  the standard boxplot outliers (\code{outlier.shape = NA}) to avoid plotting the same data point twice
#'                  (once as an outlier, once as a jittered point).
#'            \item \strong{Palette Expansion:} Like the barplot function, it automatically interpolates colors if the number
#'                  of experimental groups exceeds the palette's limit.
#'            \item \strong{Layout:} Uses \code{theme_classic()} for a clean, academic look.
#'          }
#'
#' @param data Data frame. The merged dataset (e.g., output from \code{\link{aggregate_replicates}}).
#' @param y_var String. The name of the numeric column to plot (e.g., "Max", "Mean").
#' @param x_var String. The name of the categorical column for the X-axis groupings.
#' @param fill_var String. The name of the variable used for fill colors. Default is \code{NULL}, which uses \code{x_var}.
#' @param add_points Logical. If \code{TRUE} (default), overlays individual data points using \code{geom_jitter}.
#'                   Highly recommended to show sample size and distribution density.
#' @param point_size Numeric. The size of the individual jitter points. Default is 1.5.
#' @param point_alpha Numeric. The transparency of the jitter points (0 to 1). Default is 0.6.
#' @param palette String or Vector.
#'                \itemize{
#'                  \item If a string: Pre-defined scientific palettes (\code{"npg"}, \code{"jco"}).
#'                  \item If a character vector: A custom list of hex codes.
#'                }
#'
#' @return A \code{ggplot} object. Can be further customized with standard ggplot2 functions (e.g., \code{+ ylim(20, 40)}).
#' @import ggplot2
#' @importFrom ggsci pal_npg pal_jco
#' @importFrom grDevices colorRampPalette
#' @export
#' @examples
#' df_bio <- data.frame(
#'   Treatment = rep(c("ND", "HFD"), each = 5),
#'   Mean = c(runif(5, 33, 35), runif(5, 34, 36))
#' )
#'
#' # Boxplot with individual points
#' p <- viz_thermal_boxplot(df_bio, y_var = "Mean", x_var = "Treatment")
#' p
viz_thermal_boxplot <- function(data, y_var, x_var,
                                 fill_var = NULL,
                                 add_points = TRUE,
                                 point_size = 1.5,
                                 point_alpha = 0.6,
                                 palette = "npg") {

  # 1. Validation
  if (!is.data.frame(data)) stop("Error: 'data' must be a data frame.")
  if (!y_var %in% names(data)) stop(paste("Error: Column", y_var, "not found in data."))
  if (!x_var %in% names(data)) stop(paste("Error: Column", x_var, "not found in data."))

  if (is.null(fill_var)) fill_var <- x_var
  data[[x_var]] <- as.factor(data[[x_var]])

  n_colors_needed <- length(unique(data[[fill_var]]))
  if (length(palette) > 1) {
    my_colors <- palette
  } else if (palette == "npg" && n_colors_needed <= 9) {
    my_colors <- ggsci::pal_npg()(n_colors_needed)
  } else if (palette == "jco" && n_colors_needed <= 10) {
    my_colors <- ggsci::pal_jco()(n_colors_needed)
  } else {
    base_colors <- ggsci::pal_npg()(9)
    my_colors <- colorRampPalette(base_colors)(n_colors_needed)
  }

  # 3. Base Plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[fill_var]]))

  # 4. Add Boxplot
  # Logic: If we show points, hide the outliers in the boxplot (to avoid duplicates).
  # If we don't show points, show outliers as standard dots.
  outlier_shape <- if (add_points) NA else 19

  p <- p + geom_boxplot(width = 0.7, alpha = 0.8, outlier.shape = outlier_shape)

  # 5. Add Individual Points (Customizable)
  if (add_points) {
    p <- p + geom_jitter(
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
      size = point_size,
      alpha = point_alpha,
      shape = 21
    )
  }

  # 6. Aesthetics
  p <- p +
    scale_fill_manual(values = my_colors) +
    labs(y = paste("Thermal Metric:", y_var), x = x_var) +
    theme_classic()

  # 7. Adjust Y axis
  p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

  return(p)
}



#' @title Create a Thermal Image Montage
#' @description Combines multiple BioThermR objects into a single plot (a montage/collage).
#'              Images are laid out in a grid, and filenames are labeled above each image.
#'              Uses the 'processed' matrix (NA values are transparent).
#'
#' @param img_list A list of 'BioThermR' objects.
#' @param ncol Integer. Number of columns in the grid. If NULL, tries to create a roughly square grid.
#' @param padding Integer. Pixel gap between images. Default is 10.
#' @param palette String. Color palette for temperature. Default is "inferno".
#' @param text_color String. Color of the filename labels. Default is "white" (good contrast on dark backgrounds).
#' @param text_size Integer. Font size of the labels. Default is 4.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @export
#' @examples
#' \donttest{
#' # Load a batch of images
#' img_obj_list <- system.file("extdata",package = "BioThermR")
#' batch <- read_thermal_batch(img_obj_list)
#'
#' # Create a montage2 with 4 columns
#' p <- viz_thermal_montage2(batch, ncol = 4, padding = 20)
#' }
viz_thermal_montage2 <- function(img_list, ncol = NULL, padding = 10,
                                palette = "inferno", text_color = "white", text_size = 4) {

  # 1. Validation & Setup
  if (!is.list(img_list) || length(img_list) == 0) {
    stop("Error: Input must be a non-empty list of BioThermR objects.")
  }

  n_imgs <- length(img_list)

  # Determine grid dimensions
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_imgs))
  }
  nrow <- ceiling(n_imgs / ncol)

  message(paste("Creating montage layout:", nrow, "rows x", ncol, "columns"))

  # Find max dimensions to define grid cell size
  max_h <- 0
  max_w <- 0
  valid_imgs <- list()
  img_names <- c()

  # First pass: collect valid images and find max dims
  for (i in seq_along(img_list)) {
    obj <- img_list[[i]]
    if (inherits(obj, "BioThermR") && !is.null(obj$processed)) {
      valid_imgs[[length(valid_imgs)+1]] <- obj
      img_names <- c(img_names, if(!is.null(obj$meta$filename)) obj$meta$filename else paste("Img", i))
      max_h <- max(max_h, nrow(obj$processed))
      max_w <- max(max_w, ncol(obj$processed))
    }
  }

  if (length(valid_imgs) == 0) stop("No valid processed images found in list.")

  cell_w <- max_w + padding
  cell_h <- max_h + padding

  # 2. Data Transformation Loop
  combined_df_list <- list()
  labels_df_list <- list()

  for (i in seq_along(valid_imgs)) {
    obj <- valid_imgs[[i]]
    mat <- obj$processed
    fname <- img_names[i]

    # Calculate grid position (0-indexed)
    grid_r <- (i - 1) %/% ncol
    grid_c <- (i - 1) %% ncol

    # Calculate pixel offsets (Top-Left origin strategy)
    offset_x <- grid_c * cell_w
    offset_y <- grid_r * cell_h

    df <- expand.grid(
      orig_row = seq_len(nrow(mat)),
      orig_col = seq_len(ncol(mat))
    )
    df$val <- as.vector(mat)

    # Filter out NAs
    df <- df[!is.na(df$val), ]

    if (nrow(df) > 0) {
      df$x_shifted <- df$orig_col + offset_x
      df$y_shifted <- df$orig_row + offset_y

      combined_df_list[[i]] <- df

      # Calculate label position (Centered above the image)
      labels_df_list[[i]] <- data.frame(
        x = offset_x + ncol(mat) / 2,
        y = offset_y, # Label at the "top" (lowest y value before reversal, or handling via vjust)
        label = tools::file_path_sans_ext(fname)
      )
    }
  }

  # 3. Combine Data
  big_df <- do.call(rbind, combined_df_list)
  labels_df <- do.call(rbind, labels_df_list)

  if (is.null(big_df)) {
    stop("Error: No data to plot (all images might be fully masked/NA).")
  }

  # 4. Plot
  # scale_y_reverse is used because matrix row 1 is usually "top", but in Cartesian plot y=1 is "bottom".
  p <- ggplot(big_df, aes(x = x_shifted, y = y_shifted, fill = val)) +
    geom_raster() +
    scale_fill_viridis_c(option = palette, name = "Temp (\u00B0C)", na.value = "transparent") +
    # Filename Labels
    geom_text(data = labels_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, color = text_color, size = text_size, vjust = 1.2, fontface = "bold") +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "black", color = NA), # Dark background
      legend.position = "right",
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white")
    ) +
    labs(title = NULL)

  return(p)
}


#' @title Create a Gap-Free Thermal Image Montage
#' @description Combines a list of 'BioThermR' objects into a single, high-resolution grid plot (montage).
#'              Unlike standard faceting, this function uses a custom integer-coordinate system to
#'              center each subject within a uniform grid cell, ensuring pixel-perfect alignment
#'              without the visual artifacts (white gaps) often seen in R raster plots.
#'
#' @details The function executes a two-pass algorithm:
#'          \enumerate{
#'            \item \strong{Scan Pass:} Iterates through all objects to calculate the bounding box dimensions
#'                  of the largest subject. This defines the uniform cell size for the grid.
#'            \item \strong{Layout Pass:} Calculates the integer offsets required to center each smaller subject
#'                  within the master cell. It merges all pixel data into a single master data frame.
#'          }
#'          The result is rendered using \code{geom_tile(width=1, height=1)} to guarantee continuous, gap-free visualization.
#'
#' @param img_list A list of 'BioThermR' objects (e.g., output from \code{read_thermal_batch} or \code{roi_filter_interactive}).
#' @param ncol Integer. The number of columns in the grid. If \code{NULL} (default), it is automatically calculated
#'             based on the square root of the number of images to create a roughly square layout.
#' @param padding Integer. The size of the whitespace gap (in pixels) between grid cells. Default is 10.
#' @param palette String. The color palette to use (from 'viridis' package). Default is "inferno".
#' @param text_color String. Color of the filename labels. Default is "black".
#' @param text_size Integer. Font size for the filename labels. Default is 4.
#'
#' @return A \code{ggplot} object representing the combined montage.
#' @import ggplot2
#' @importFrom tools file_path_sans_ext
#' @export
#' @examples
#' \donttest{
#' # Load a batch of images
#' img_obj_list <- system.file("extdata",package = "BioThermR")
#' batch <- read_thermal_batch(img_obj_list)
#'
#' # Create a montage with 4 columns
#' p <- plot_thermal_montage(batch, ncol = 4, padding = 20)
#' }
plot_thermal_montage <- function(img_list, ncol = NULL, padding = 10,
                                palette = "inferno", text_color = "black", text_size = 4) {

  # 1. Validation & Setup
  if (!is.list(img_list) || length(img_list) == 0) {
    stop("Error: Input must be a non-empty list of BioThermR objects.")
  }

  # 2. First Pass: Extract Data & Find Max Bounding Box
  valid_items <- list()
  global_max_h <- 0
  global_max_w <- 0

  message("Scanning objects for bounding box dimensions...")

  for (i in seq_along(img_list)) {
    obj <- img_list[[i]]
    if (inherits(obj, "BioThermR") && !is.null(obj$processed)) {

      mat <- obj$processed
      fname <- if(!is.null(obj$meta$filename)) obj$meta$filename else paste("Img", i)

      # Convert matrix to DF (Row first)
      df <- expand.grid(
        orig_row = seq_len(nrow(mat)),
        orig_col = seq_len(ncol(mat))
      )
      df$val <- as.vector(mat)

      # Filter NAs
      df <- df[!is.na(df$val), ]

      if (nrow(df) > 0) {
        # Calculate Bounding Box
        min_r <- min(df$orig_row)
        max_r <- max(df$orig_row)
        min_c <- min(df$orig_col)
        max_c <- max(df$orig_col)

        # Dimensions of the ROI
        obj_h <- max_r - min_r + 1
        obj_w <- max_c - min_c + 1

        # Update Global Max Dimensions
        global_max_h <- max(global_max_h, obj_h)
        global_max_w <- max(global_max_w, obj_w)

        valid_items[[length(valid_items)+1]] <- list(
          df = df,
          fname = fname,
          bounds = c(min_r, min_c, obj_w, obj_h)
        )
      }
    }
  }

  n_valid <- length(valid_items)
  if (n_valid == 0) stop("Error: No valid data found.")

  # 3. Determine Grid Layout
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_valid))
  }
  nrow <- ceiling(n_valid / ncol)

  cell_w <- global_max_w + padding
  cell_h <- global_max_h + padding

  message(paste("Montage Layout:", nrow, "x", ncol, "| Cell Size:", cell_w, "x", cell_h))

  # 4. Second Pass: Calculate Offsets
  combined_df_list <- list()
  labels_df_list <- list()

  for (i in seq_along(valid_items)) {
    item <- valid_items[[i]]
    df <- item$df
    fname <- item$fname

    # Unpack bounds
    min_r <- item$bounds[1]
    min_c <- item$bounds[2]
    w <- item$bounds[3]
    h <- item$bounds[4]

    # Grid Position
    grid_r <- (i - 1) %/% ncol
    grid_c <- (i - 1) %% ncol

    # Base Offset
    base_offset_x <- grid_c * cell_w
    base_offset_y <- grid_r * cell_h

    # Centering Offset (Must be INTEGER to avoid gaps!)
    center_pad_x <- floor((global_max_w - w) / 2)
    center_pad_y <- floor((global_max_h - h) / 2)

    # Transform Coordinates
    df$x_shifted <- (df$orig_col - min_c) + center_pad_x + base_offset_x
    df$y_shifted <- (df$orig_row - min_r) + center_pad_y + base_offset_y

    combined_df_list[[i]] <- df

    # Label Position
    labels_df_list[[i]] <- data.frame(
      x = base_offset_x + (global_max_w / 2),
      y = base_offset_y,
      label = tools::file_path_sans_ext(fname)
    )
  }

  big_df <- do.call(rbind, combined_df_list)
  labels_df <- do.call(rbind, labels_df_list)

  # 5. Plot
  p <- ggplot(big_df, aes(x = x_shifted, y = y_shifted, fill = val)) +
    # Use geom_tile with explicit width/height = 1 to ensure solid pixels
    geom_tile(width = 1, height = 1) +
    scale_fill_viridis_c(option = palette, name = "Temp (\u00B0C)", na.value = "transparent") +
    geom_text(data = labels_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, color = text_color, size = text_size, vjust = 1.2, fontface = "bold") +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black")
    ) +
    labs(title = NULL)

  return(p)
}





#' @title Generate a "Thermal Cloud" Visualization (Phyllotaxis Layout)
#' @description Arranges a collection of segmented thermal objects into an organic, spiral "cloud" formation.
#'              Unlike a rigid grid, this layout uses a golden-angle spiral (phyllotaxis) algorithm to cluster
#'              subjects efficiently. This is particularly effective for visualizing the diversity of thermal
#'              phenotypes in large datasets or creating artistic figures for presentations and covers.
#'
#' @details The function performs object-centric rendering:
#'          \enumerate{
#'            \item \strong{Extraction:} For each image, it extracts only the valid foreground pixels (non-NA),
#'                  ignoring the original frame dimensions.
#'            \item \strong{Re-centering:} Each object is mathematically centered at (0,0) relative to its own coordinate system.
#'            \item \strong{Placement:} Objects are placed along a spiral path defined by the Golden Angle (~137.5 degrees).
#'          }
#'          The spacing and randomness of the spiral can be tuned using \code{spread_factor} and \code{jitter_factor}.
#'
#' @param img_list A list of 'BioThermR' objects. For best results, these should be pre-processed
#'                 (e.g., background removed via \code{\link{roi_filter_threshold}}).
#' @param spread_factor Numeric. Multiplier for the distance between objects.
#'                      Values > 1.0 increase spacing (airier cloud), values < 1.0 pack objects tighter. Default is 1.1.
#' @param jitter_factor Numeric. Introduces random noise to the placement coordinates to break perfect symmetry
#'                      and create a more natural look. Default is 0.5.
#' @param palette String. The color palette from the 'viridis' package. Default is "inferno".
#' @param text_color String. Color of the filename labels. Default is "black".
#' @param text_size Integer. Font size for the labels. Default is 3.
#' @param show_labels Logical. If \code{TRUE}, displays the filename below each object. Default is \code{TRUE}.
#'
#' @return A \code{ggplot} object with a white background and void theme.
#' @import ggplot2
#' @importFrom tools file_path_sans_ext
#' @export
#' @examples
#' \donttest{
#' # Load a batch of images
#' img_obj_list <- system.file("extdata",package = "BioThermR")
#' batch <- read_thermal_batch(img_obj_list)
#' batch <- lapply(batch, roi_segment_ebimage)
#'
#' # Create an artistic thermal cloud
#' p_cloud <- plot_thermal_cloud(batch, spread_factor = 1.5, jitter_factor = 2.0)
#' }
plot_thermal_cloud <- function(img_list, spread_factor = 1.1, jitter_factor = 0.5,
                              palette = "inferno", text_color = "black", text_size = 3,
                              show_labels = TRUE) {

  # 1. Check Input
  if (!is.list(img_list) || length(img_list) == 0) {
    stop("Error: Input must be a non-empty list of BioThermR objects.")
  }

  # 2. Extract Data & Calculate Relative Coordinates
  valid_objects <- list()
  max_obj_radius <- 0

  message("Processing objects...")

  for (i in seq_along(img_list)) {
    obj <- img_list[[i]]
    if (inherits(obj, "BioThermR") && !is.null(obj$processed)) {

      mat <- obj$processed

      # Step A: Matrix to DataFrame (Row first!)
      df <- expand.grid(
        orig_row = seq_len(nrow(mat)),
        orig_col = seq_len(ncol(mat))
      )
      df$val <- as.vector(mat)

      # Step B: Filter NAs (Background)
      df <- df[!is.na(df$val), ]

      if (nrow(df) > 0) {
        # Step C: Re-center the object
        center_row <- (min(df$orig_row) + max(df$orig_row)) / 2
        center_col <- (min(df$orig_col) + max(df$orig_col)) / 2

        # Create relative coordinates (centered at 0,0)
        df$rel_x <- df$orig_col - center_col
        df$rel_y <- df$orig_row - center_row

        # Calculate approximate radius
        obj_radius <- max(max(abs(df$rel_x)), max(abs(df$rel_y)))
        max_obj_radius <- max(max_obj_radius, obj_radius)

        valid_objects[[length(valid_objects) + 1]] <- list(
          df = df,
          fname = if(!is.null(obj$meta$filename)) obj$meta$filename else paste("Img", i)
        )
      }
    }
  }

  n <- length(valid_objects)
  if (n == 0) stop("No valid data found.")

  message(paste0("Generating layout for ", n, " objects..."))

  # 3. Calculate Spiral Layout (Phyllotaxis)
  theta <- 2.39996 # Golden Angle
  step_size <- (max_obj_radius * 2) * spread_factor

  final_df_list <- list()
  labels_df_list <- list()

  for (i in seq_along(valid_objects)) {
    item <- valid_objects[[i]]
    df <- item$df

    # --- Spiral Math ---
    r <- step_size * sqrt(i)
    a <- i * theta

    cx <- r * cos(a)
    cy <- r * sin(a)

    # Add Jitter
    jitter_amt <- step_size * jitter_factor
    cx <- cx + runif(1, -jitter_amt, jitter_amt)
    cy <- cy + runif(1, -jitter_amt, jitter_amt)

    # --- Merge Coordinates ---
    df$final_x <- cx + df$rel_x
    df$final_y <- cy + df$rel_y

    final_df_list[[i]] <- df

    # Store Label Position
    if (show_labels) {
      labels_df_list[[i]] <- data.frame(
        x = cx,
        y = cy - max_obj_radius * 1.2,
        label = tools::file_path_sans_ext(item$fname)
      )
    }
  }

  # 4. Combine All Data
  big_df <- do.call(rbind, final_df_list)

  # 5. Plot (White Background)
  p <- ggplot(big_df, aes(x = final_x, y = final_y, fill = val)) +
    # Force 1x1 pixel size to prevent gaps or invisibility
    geom_tile(width = 1, height = 1) +
    scale_fill_viridis_c(option = palette, name = "Temp (\u00B0C)", na.value = "transparent") +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.text = element_text(color = "black"), # Black text
      legend.title = element_text(color = "black")
    )

  # Add Labels
  if (show_labels && length(labels_df_list) > 0) {
    labels_df <- do.call(rbind, labels_df_list)
    p <- p + geom_text(data = labels_df, aes(x = x, y = y, label = label),
                       inherit.aes = FALSE,
                       color = text_color, size = text_size, fontface = "bold")
  }

  return(p)
}

#' @title Visualize ROI Overlap and Dice Coefficient
#' @description Visually compares two segmentation masks (ROIs) by overlaying them on the original thermal image.
#'              This function is primarily used for validation purposes: to compare an automated segmentation
#'              (filled layer) against a manual ground truth (contour line).
#'              It automatically calculates and displays the \strong{Dice Similarity Coefficient (DSC)} in the title.
#'
#' @details The visualization consists of three layers:
#'          \enumerate{
#'            \item \strong{Background:} The raw thermal image from \code{img_obj1}.
#'            \item \strong{Prediction (img_obj1):} The processed mask from the first object, rendered as a semi-transparent filled raster (default green).
#'            \item \strong{Ground Truth (img_obj2):} The processed mask from the second object, rendered as a contour outline (default white).
#'          }
#'          The function calculates the Dice Similarity Coefficient (DSC) using the formula:
#'          \deqn{DSC = \frac{2 \times |X \cap Y|}{|X| + |Y|}}
#'          where X and Y are the set of pixels in the two masks. A DSC of 1 indicates perfect overlap.
#'
#' @param img_obj1 A 'BioThermR' object. Typically the \strong{Automated/Predicted} segmentation.
#'                 This object must contain the raw thermal matrix. Its mask will be plotted as a filled area.
#' @param img_obj2 A 'BioThermR' object. Typically the \strong{Manual/Ground Truth} segmentation.
#'                 Its mask will be plotted as a contour outline. Dimensions must match \code{img_obj1}.
#' @param title String. Custom title for the plot. If \code{NULL} (default), the title shows "ROI overlap (DICE: X.XXX)".
#' @param color String. Fill color for the \code{img_obj1} mask. Default is "green".
#' @param alpha Numeric. Transparency level for the \code{img_obj1} mask (0 to 1). Default is 0.5.
#' @param line_color String. Line color for the \code{img_obj2} contour. Default is "white".
#' @param palette String. Color palette for the background thermal image (passed to \code{scale_fill_viridis_c}).
#'                Default is "inferno".
#'
#' @return A \code{ggplot} object showing the overlay.
#' @import ggplot2
#' @export
#' @examples
#' \donttest{
#' #' # Load raw data
#' img_obj <- system.file("extdata", "C05.raw", package = "BioThermR")
#' img <- read_thermal_raw(img_obj)
#' # Apply automated segmentation
#' img1 <- roi_segment_ebimage(img, keep_largest = TRUE)
#'
#' # Simple background removal: Keep everything above 24 degrees
#' img2 <- roi_filter_threshold(img, threshold = c(33, Inf))
#'
#' # Compare them
#' plot_roi_overlap(img_obj1 = img1,
#'                  img_obj2 = img2)
#' }
plot_roi_overlap <- function(img_obj1, img_obj2,
                             title = NULL,
                             color = "green",
                             alpha = 0.5,
                             line_color = "white",
                             palette = "inferno") {
  if (!inherits(img_obj1, "BioThermR")) {
    stop("Error: Input must be a 'BioThermR' object.")
  }
  if (!inherits(img_obj2, "BioThermR")) {
    stop("Error: Input must be a 'BioThermR' object.")
  }
  if (!is.null(img_obj1$raw)) base_temp_mat <- img_obj1$raw
  auto_mat <- img_obj1$processed
  manual_mat <- img_obj2$processed

  stopifnot(all(dim(base_temp_mat) == dim(manual_mat)),
            all(dim(base_temp_mat) == dim(auto_mat)))
  dice_coef <- function(a, b) {
    a <- as.logical(a)
    b <- as.logical(b)
    if (sum(a) + sum(b) == 0) return(0)
    2 * sum(a & b) / (sum(a) + sum(b))
  }
  df_base <- expand.grid(
    x = seq_len(nrow(base_temp_mat)),
    y = seq_len(ncol(base_temp_mat))
  )
  df_base$temp <- c(base_temp_mat)

  auto_mask <- !is.na(auto_mat)
  df_auto <- df_base
  df_auto$auto <- c(auto_mask)

  manual_mask <- !is.na(manual_mat)
  df_manual <- df_base
  df_manual$manual <- c(manual_mask)
  dice<-dice_coef(auto_mask,manual_mask)
  if(is.null(title))title <- paste("ROI overlap (DICE:",round(dice,3),")")

  p<-ggplot() +
    geom_raster(data = df_base, aes(x = y, y = x, fill = temp)) +
    scale_fill_viridis_c(option = palette, name = "Temp (\u00B0C)", na.value = "transparent") +

    geom_raster(data = subset(df_auto, auto),
                aes(x = y, y = x),
                fill = color, alpha = alpha) +
    geom_contour(data = df_manual,
                 aes(x = y, y = x, z = as.numeric(manual)),
                 breaks = 0.5, color = line_color, linewidth = 0.5) +
    coord_fixed(expand = FALSE) +
    labs(title = title) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  return(p)
}


