#' @title Automated ROI Segmentation via EBImage
#' @description Performs automated background removal using a hybrid pipeline of global thresholding
#'              (Otsu's method), morphological operations, and connected component analysis.
#'
#' @details The segmentation pipeline consists of four steps:
#'          \enumerate{
#'            \item \strong{Normalization:} The temperature matrix is scaled to [0, 1] to be compatible with 'EBImage'.
#'            \item \strong{Thresholding:} Otsu's method is used to calculate an optimal global threshold
#'                  that separates the foreground (subject) from the background based on histogram bimodality.
#'            \item \strong{Morphology:} A disc-shaped brush (size=5) is used for 'Opening' (to remove salt noise)
#'                  and 'Closing' (to fill small holes inside the subject).
#'            \item \strong{Component Filter:} If \code{keep_largest} is \code{TRUE}, the function labels all
#'                  connected regions and retains only the largest one (assuming the animal is the largest heat source),
#'                  effectively removing smaller artifacts like bedding noise or reflections.
#'          }
#'
#' @param img_obj A 'BioThermR' object.
#' @param method String. The thresholding algorithm. Currently, only \code{"otsu"} is supported.
#' @param keep_largest Logical. If \code{TRUE} (default), keeps only the largest connected object
#'                     and converts all other clusters to background (\code{NA}). This is a powerful
#'                     denoising step for animal experiments.
#' @param morphology Logical. If \code{TRUE} (default), applies morphological opening and closing
#'                   operations to smooth edges and reduce noise.
#'
#' @return A 'BioThermR' object where the \code{processed} matrix has been masked
#'         (background pixels are set to \code{NA}).
#' @importFrom EBImage Image otsu makeBrush opening closing bwlabel
#' @export
#' @examples
#' # Load raw data
#' img_obj <- system.file("extdata", "C05.raw", package = "BioThermR")
#' img <- read_thermal_raw(img_obj)
#'
#' # Apply automated segmentation
#' img <- roi_segment_ebimage(img, keep_largest = TRUE)
roi_segment_ebimage <- function(img_obj, method = "otsu", keep_largest = TRUE, morphology = TRUE) {

  if (!inherits(img_obj, "BioThermR")) {
    stop("Error: Input must be a 'BioThermR' object.")
  }

  # Check for EBImage dependency
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop("Error: Package 'EBImage' is required. Please install it via BiocManager::install('EBImage').")
  }

  mat <- img_obj$processed

  # Create a temporary normalized image
  range_val <- range(mat, na.rm = TRUE)
  mat_norm <- (mat - range_val[1]) / (range_val[2] - range_val[1])

  # Convert to EBImage object
  eb_img <- EBImage::Image(mat_norm)

  # 1. Thresholding
  if (method == "otsu") {
    # Calculate Otsu's threshold
    thr <- EBImage::otsu(eb_img)
    mask <- mat_norm > thr
  } else {
    stop("Error: Only 'otsu' method is currently supported.")
  }

  # 2. Morphological Operations (Optional)
  if (morphology) {
    kern <- EBImage::makeBrush(5, shape = "disc")
    mask <- EBImage::opening(mask, kern)
    mask <- EBImage::closing(mask, kern)
  }

  # 3. Connected Components Labeling (Find objects)
  # Converts boolean mask to integer labels (0=bg, 1=obj1, 2=obj2...)
  labels <- EBImage::bwlabel(mask)

  # 4. Keep Largest Object
  if (keep_largest) {
    # Count pixels per label
    tbl <- table(labels)
    # Remove background (label 0)
    tbl <- tbl[names(tbl) != "0"]

    if (length(tbl) > 0) {
      # Find ID of the largest object
      largest_id <- names(which.max(tbl))
      # Update mask to only include this object
      mask <- (labels == as.integer(largest_id))
      message(paste("Auto-Segmentation: Kept largest object (", max(tbl), "pixels )"))
    } else {
      warning("Warning: No object detected after thresholding.")
      mask <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))
    }
  }

  # 5. Apply Mask to Original Matrix
  # Convert EBImage mask back to logical matrix
  final_mask <- as.matrix(mask)

  # Set background to NA
  img_obj$processed[!final_mask] <- NA
  img_obj$stats <- NULL # Reset stats

  return(img_obj)
}

#' @title Filter Thermal Image by Temperature Thresholds
#' @description Masks the thermal image by retaining only pixels that fall within a specified
#'              temperature window. Pixels outside this range are set to \code{NA} (background).
#'
#' @details This is the most fundamental segmentation method for thermal images.
#'          Since animals are typically warmer than their environment, a simple low-pass filter
#'          (e.g., keep > 22 degrees Celsius) is often sufficient to separate the subject from the cage.
#'          Open boundaries can be defined using \code{Inf} or \code{-Inf}.
#'
#' @param img_obj A 'BioThermR' object.
#' @param threshold Numeric vector of length 2, e.g., \code{c(22, 38)}.
#'                  Defines the inclusive temperature range \code{[min, max]} to keep.
#'                  Use \code{Inf} for open upper bounds (e.g., \code{c(25, Inf)} keeps everything above 25 degrees Celsius).
#' @param use_processed Logical.
#'                  \itemize{
#'                    \item If \code{FALSE} (default): The filter is applied to the \strong{raw} matrix, discarding any previous masks.
#'                    \item If \code{TRUE}: The filter is applied to the \strong{processed} matrix, allowing for cumulative masking (e.g., refining an Otsu segmentation).
#'                  }
#'
#' @return A 'BioThermR' object with the \code{processed} matrix updated.
#' @export
#' @examples
#' # Load raw data
#' img_obj <- system.file("extdata", "C05.raw", package = "BioThermR")
#' img <- read_thermal_raw(img_obj)
#'
#' # Simple background removal: Keep everything above 24 degrees
#' img <- roi_filter_threshold(img, threshold = c(24, Inf))
roi_filter_threshold <- function(img_obj, threshold, use_processed = FALSE) {

  # 1. Validate Input
  if (!inherits(img_obj, "BioThermR")) {
    stop("Error: Input must be a 'BioThermR' object.")
  }

  if (!is.numeric(threshold) || length(threshold) != 2) {
    stop("Error: 'threshold' must be a numeric vector of length 2, e.g., c(20, 35).")
  }

  # Ensure the range is sorted (min, max)
  threshold <- sort(threshold)
  min_val <- threshold[1]
  max_val <- threshold[2]

  mat <- if (use_processed) img_obj$processed else img_obj$raw

  # 2. Create Mask
  # Keep values that are >= min AND <= max
  # We wrap in () to ensure order of operations
  mask <- (mat >= min_val) & (mat <= max_val)

  # Handle existing NAs (NA comparison results in NA, we want FALSE)
  mask[is.na(mask)] <- FALSE

  # 3. Apply Mask
  # Set everything NOT in the mask to NA
  mat[!mask] <- NA

  # 4. Update Object
  img_obj$processed <- mat
  img_obj$stats <- NULL

  message(paste0("ROI Filter applied: Keeping range [", min_val, ", ", max_val, "]"))

  return(img_obj)
}



#' @title Interactive ROI Selector with Denoising (Shiny App)
#' @description Launches a Shiny GUI application that allows users to interactively refine regions of interest (ROI)
#'              for a batch of thermal images. Key features include dynamic temperature thresholding and
#'              an automated "Clean Noise" tool that removes artifacts by retaining only the largest connected object (e.g., the mouse).
#'
#' @details This function opens a local web server (Shiny App). The workflow is as follows:
#'          \enumerate{
#'            \item \strong{Navigate:} Use "Prev/Next" buttons to browse the image batch.
#'            \item \strong{Threshold:} Adjust the slider to set the min/max temperature range. Real-time feedback is shown on the plot.
#'            \item \strong{Denoise:} Toggle the "Remove Noise" button. This applies a topological filter (Connected Component Analysis)
#'                  that identifies the largest contiguous heat source and removes all smaller isolated islands (e.g., urine spots, reflections).
#'            \item \strong{Confirm:} Click "Confirm" to lock in the settings for the current image and auto-advance.
#'            \item \strong{Export:} Click "Finish & Export Data" to close the app and return the processed object list to the R console.
#'          }
#'
#' @param img_input A single 'BioThermR' object OR a list of 'BioThermR' objects (e.g., from \code{read_thermal_batch}).
#' @param start_index Integer. The index of the image to start viewing. Default is 1.
#'                    Useful for resuming work on a large batch.
#' @param use_processed Logical. If \code{TRUE} (default), the interactive filter is applied to the
#'                      already 'processed' matrix (e.g., refining an auto-segmented result).
#'                      If \code{FALSE}, it starts from the raw data.
#'
#' @return The modified 'BioThermR' object (or list of objects). \strong{Note:} You must assign the result
#'         of this function to a variable (e.g., \code{data <- roi_filter_interactive(data)}) to save the changes.
#' @import shiny
#' @import ggplot2
#' @importFrom EBImage bwlabel
#' @export
roi_filter_interactive <- function(img_input, start_index = 1, use_processed = TRUE) {

  # --- 1. Data Normalization ---
  is_list <- FALSE
  if (inherits(img_input, "BioThermR")) {
    img_list <- list(img_input)
    is_list <- FALSE
  } else if (is.list(img_input) && inherits(img_input[[1]], "BioThermR")) {
    img_list <- img_input
    is_list <- TRUE
  } else {
    stop("Error: Input must be a 'BioThermR' object or a list of 'BioThermR' objects.")
  }

  if (!requireNamespace("shiny", quietly = TRUE)) stop("Package 'shiny' is required.")
  if (!requireNamespace("EBImage", quietly = TRUE)) stop("Package 'EBImage' is required.")

  n_total <- length(img_list)
  if (start_index < 1 || start_index > n_total) start_index <- 1

  # --- 2. Pre-calculate Ranges ---
  img_ranges <- lapply(img_list, function(obj) {
    raw_vals <- as.vector(if (use_processed) obj$processed else obj$raw)
    c(floor(min(raw_vals, na.rm = TRUE)), ceiling(max(raw_vals, na.rm = TRUE)))
  })

  # --- 3. UI Definition ---
  ui <- shiny::fluidPage(
    shiny::titlePanel("BioThermR: Interactive ROI Selector"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,

        # Navigation
        shiny::h4("Navigation"),
        shiny::div(style="text-align: center;",
                   shiny::actionButton("prev_btn", "Prev", icon = shiny::icon("arrow-left")),
                   shiny::span(" "),
                   shiny::textOutput("page_info", inline = TRUE),
                   shiny::span(" "),
                   shiny::actionButton("next_btn", "Next", icon = shiny::icon("arrow-right"))
        ),
        shiny::hr(),

        # Status
        shiny::div(style="text-align: center; padding: 5px; background-color: #f0f0f0;",
                   shiny::strong("Status: "), shiny::uiOutput("status_badge", inline = TRUE)
        ),
        shiny::hr(),

        # Threshold Slider
        shiny::h4("Step 1: Threshold"),
        shiny::uiOutput("slider_ui"),

        # Denoise Button (New Feature)
        shiny::h4("Step 2: Clean Up"),
        shiny::actionButton("clean_btn", "Remove Noise (Keep Largest)", class = "btn-warning", width = "100%", icon = shiny::icon("eraser")),
        shiny::helpText("Removes isolated pixels, keeping only the largest object."),
        shiny::uiOutput("clean_status"),

        shiny::hr(),

        # Action Buttons
        shiny::actionButton("confirm_btn", "Confirm / Save Current", class = "btn-primary", width = "100%", icon = shiny::icon("check")),
        shiny::br(), shiny::br(),
        shiny::actionButton("finish_btn", "Finish & Export Data", class = "btn-success", width = "100%", icon = shiny::icon("file-export"))
      ),

      shiny::mainPanel(
        width = 9,
        shiny::h4(shiny::textOutput("filename_display")),
        shiny::plotOutput("thermal_plot", height = "600px")
      )
    )
  )

  # --- 4. Server Logic ---
  server <- function(input, output, session) {

    # State storage
    storage <- shiny::reactiveValues(
      index = start_index,
      thresholds = vector("list", n_total),
      clean_flags = rep(FALSE, n_total),
      modified_flags = rep(FALSE, n_total)
    )

    # Navigation
    shiny::observeEvent(input$prev_btn, { if (storage$index > 1) storage$index <- storage$index - 1 })
    shiny::observeEvent(input$next_btn, { if (storage$index < n_total) storage$index <- storage$index + 1 })

    # Outputs
    output$page_info <- shiny::renderText({ paste(storage$index, "/", n_total) })
    output$filename_display <- shiny::renderText({ img_list[[storage$index]]$meta$filename })

    output$status_badge <- shiny::renderUI({
      if (storage$modified_flags[storage$index]) shiny::span("Saved", style="color:green;font-weight:bold")
      else shiny::span("Unsaved", style="color:grey")
    })

    output$clean_status <- shiny::renderUI({
      if (storage$clean_flags[storage$index]) shiny::div("Noise Removed", style="color:orange; font-size:0.8em;")
      else shiny::div("")
    })

    # Slider UI
    output$slider_ui <- shiny::renderUI({
      idx <- storage$index
      limits <- img_ranges[[idx]]
      val <- if (storage$modified_flags[idx]) storage$thresholds[[idx]] else limits
      shiny::sliderInput("temp_range", "Range (\u00B0C):", min = limits[1], max = limits[2], value = val, step = 0.1)
    })

    # Helper: Apply Cleaning Logic (Largest Component)
    apply_cleaning <- function(mat) {
      # 1. Binarize (NA is background, numbers are foreground)
      binary_mask <- !is.na(mat)

      # 2. Label components
      labeled <- EBImage::bwlabel(binary_mask)

      # 3. Find largest component
      dims <- dim(labeled)
      # Tabulate frequencies (0 is background, ignore it)
      counts <- table(labeled)
      counts <- counts[names(counts) != "0"]

      if (length(counts) == 0) return(mat) # No object found

      largest_label <- as.numeric(names(counts)[which.max(counts)])

      # 4. Filter matrix
      # Set anything that is NOT the largest label to NA
      mat[labeled != largest_label] <- NA
      return(mat)
    }

    # Reactive: Current processed matrix (for plotting)
    current_mat <- shiny::reactive({
      shiny::req(input$temp_range)
      idx <- storage$index
      obj <- img_list[[idx]]

      # 1. Apply Threshold
      temp_obj <- roi_filter_threshold(obj, input$temp_range, use_processed = use_processed)
      mat <- temp_obj$processed

      # 2. Apply Cleaning (if flag is set)
      if (storage$clean_flags[idx]) {
        mat <- apply_cleaning(mat)
      }
      return(mat)
    })

    # Handle "Clean" Button Click
    shiny::observeEvent(input$clean_btn, {
      storage$clean_flags[storage$index] <- !storage$clean_flags[storage$index]
    })

    # Plotting
    output$thermal_plot <- shiny::renderPlot({
      mat <- current_mat()

      df <- expand.grid(row = seq_len(nrow(mat)), col = seq_len(ncol(mat)))
      df$val <- c(mat)

      ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = val)) +
        ggplot2::geom_raster() +
        ggplot2::scale_fill_viridis_c(option = "inferno", na.value = "grey90", name = "Temp") +
        ggplot2::coord_fixed() +
        ggplot2::theme_void()
    })

    # Save Individual
    shiny::observeEvent(input$confirm_btn, {
      idx <- storage$index
      storage$thresholds[[idx]] <- input$temp_range
      # Clean flag is already in storage$clean_flags
      storage$modified_flags[idx] <- TRUE

      if (idx < n_total) {
        storage$index <- idx + 1
        shiny::showNotification(paste("Image", idx, "saved. Next."), duration = 1)
      } else {
        shiny::showNotification("Saved. Last image.", type = "message")
      }
    })

    # Finish
    shiny::observeEvent(input$finish_btn, {
      shiny::showNotification("Processing...", type = "message")

      final_list <- lapply(seq_len(n_total), function(i) {
        if (storage$modified_flags[i]) {
          # 1. Threshold
          obj <- roi_filter_threshold(img_list[[i]], storage$thresholds[[i]], use_processed = use_processed)
          # 2. Clean (if selected)
          if (storage$clean_flags[i]) {
            obj$processed <- apply_cleaning(obj$processed)
          }
          return(obj)
        } else {
          return(img_list[[i]])
        }
      })
      names(final_list) <- names(img_list)
      if (is_list) shiny::stopApp(final_list) else shiny::stopApp(final_list[[1]])
    })
  }

  message(paste("Launching Interactive Selector. Start:", start_index))
  shiny::runApp(shiny::shinyApp(ui, server))
}



