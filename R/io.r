#' @title Read a Single .raw Thermal Image File
#' @description Reads a binary .raw file (typically 32-bit floating point), converts it into a matrix,
#'              and constructs a 'BioThermR' object containing raw data and metadata.
#'
#' @param file_path String. The full path to the .raw file.
#' @param width Integer. The width of the thermal sensor (number of columns). Default is 160.
#' @param height Integer. The height of the thermal sensor (number of rows). Default is 120.
#' @param rotate Logical. Whether to rotate the image 90 degrees counter-clockwise.
#'               Default is TRUE (corrects orientation for many standard sensor exports).
#'
#' @return A list object of class "BioThermR" containing:
#' \item{raw}{The original temperature matrix (numeric).}
#' \item{processed}{A copy of the raw matrix, intended for subsequent ROI filtering or masking.}
#' \item{meta}{A list containing metadata: \code{filename}, \code{fullpath}, and \code{dims}.}
#'
#' @export
#' @examples
#' img_obj <- system.file("extdata", "C05.raw", package = "BioThermR")
#' img <- read_thermal_raw(img_obj)
read_thermal_raw <- function(file_path, width = 160, height = 120, rotate = TRUE) {

  # 1. Check if file exists
  if (!file.exists(file_path)) {
    stop(paste0("Error: File does not exist -> ", file_path))
  }

  # 2. Read binary data
  # Use on.exit to ensure the file connection closes even if errors occur
  file_con <- file(file_path, "rb")
  on.exit(close(file_con), add = TRUE)

  raw_vector <- readBin(
    con    = file_con,
    what   = "numeric",
    n      = width * height,
    size   = 4,
    endian = "little"
  )

  # Check data length
  if (length(raw_vector) != width * height) {
    warning(paste("Warning: Data length read (", length(raw_vector),
                  ") does not match expected dimensions (", width * height, ")."))
  }

  # 3. Matrix conversion and rotation
  mat <- matrix(raw_vector, nrow = height, ncol = width, byrow = TRUE)
  mat <- t(mat)

  if (rotate) {
    rotate_90_ccw <- function(x) apply(t(x), 2, rev)
    mat <- rotate_90_ccw(mat)
  }

  # 4. Construct Object
  obj <- structure(
    list(
      raw = mat,
      processed = mat,
      meta = list(
        filename = basename(file_path),
        fullpath = normalizePath(file_path, mustWork = FALSE),
        dims = dim(mat)
      )
    ),
    class = "BioThermR"
  )

  return(obj)
}


#' @title Read FLIR Radiometric JPG File
#' @description Reads a FLIR radiometric JPG file, parses the embedded metadata (emissivity,
#'              distance, Planck constants, atmospheric parameters), converts the raw sensor data
#'              to temperature (Celsius), and constructs a 'BioThermR' object.
#'
#' @details This function relies on the 'Thermimage' package and the external tool 'ExifTool'.
#'          It automatically extracts calibration constants (Planck R1, B, F, O, R2) and
#'          environmental parameters to ensure accurate physical temperature conversion.
#'          Please ensure 'ExifTool' is installed and added to your system PATH.
#'
#' @param file_path String. The full path to the FLIR radiometric .jpg file.
#' @param exiftoolpath String. Path to the ExifTool executable. Default is "installed"
#'                     (assumes it is available in the system PATH).
#'
#' @return A list object of class "BioThermR" containing:
#' \item{raw}{The calculated temperature matrix in degrees Celsius.}
#' \item{processed}{A copy of the raw matrix, intended for subsequent ROI filtering or masking.}
#' \item{meta}{A list containing metadata: \code{filename}, \code{fullpath}, and \code{dims}.}
#'
#' @importFrom Thermimage flirsettings readflirJPG raw2temp
#' @export
#' @examples
#' \donttest{
#' # Example using a flir thermal file
#' img_obj <- system.file("extdata", "IR_2412.jpg", package = "Thermimage")
#' img <- read_thermal_flir(img_obj)
#' }
read_thermal_flir <- function(file_path, exiftoolpath = "installed") {

  # 1. Check if file exists
  if (!file.exists(file_path)) {
    stop(paste0("Error: File does not exist -> ", file_path))
  }

  # 2. Check dependencies
  if (!requireNamespace("Thermimage", quietly = TRUE)) {
    stop("Package 'Thermimage' is required.")
  }

  # 3. Extract Metadata using flirsettings
  cams <- Thermimage::flirsettings(file_path, exiftoolpath = exiftoolpath, camvals = "")

  # 4. Extract Parameters (Direct mapping from your snippet)
  # Basic Vars
  ObjectEmissivity <- cams$Info$Emissivity
  OD               <- cams$Info$ObjectDistance
  ReflT            <- cams$Info$ReflectedApparentTemperature
  AtmosT           <- cams$Info$AtmosphericTemperature
  RH               <- cams$Info$RelativeHumidity

  # Window Vars
  IRWinT           <- cams$Info$IRWindowTemperature
  IRWinTran        <- cams$Info$IRWindowTransmission

  # Planck Constants
  PlanckR1 <- cams$Info$PlanckR1
  PlanckB  <- cams$Info$PlanckB
  PlanckF  <- cams$Info$PlanckF
  PlanckO  <- cams$Info$PlanckO
  PlanckR2 <- cams$Info$PlanckR2

  # Atmospheric Constants
  ATA1 <- cams$Info$AtmosphericTransAlpha1
  ATA2 <- cams$Info$AtmosphericTransAlpha2
  ATB1 <- cams$Info$AtmosphericTransBeta1
  ATB2 <- cams$Info$AtmosphericTransBeta2
  ATX  <- cams$Info$AtmosphericTransX

  # 5. Read Raw Sensor Data
  img <- Thermimage::readflirJPG(file_path, exiftoolpath = exiftoolpath)

  # 6. Convert to Temperature
  mat <- Thermimage::raw2temp(
    raw = img,
    E = ObjectEmissivity,
    OD = OD,
    RTemp = ReflT,
    ATemp = AtmosT,
    IRWTemp = IRWinT,
    IRT = IRWinTran,
    RH = RH,
    PR1 = PlanckR1,
    PB = PlanckB,
    PF = PlanckF,
    PO = PlanckO,
    PR2 = PlanckR2,
    ATA1 = ATA1,
    ATA2 = ATA2,
    ATB1 = ATB1,
    ATB2 = ATB2,
    ATX = ATX
  )

  # 7. Construct Object
  obj <- structure(
    list(
      raw = mat,
      processed = mat,
      meta = list(
        filename = basename(file_path),
        fullpath = normalizePath(file_path, mustWork = FALSE),
        dims = dim(mat)
      )
    ),
    class = "BioThermR"
  )

  return(obj)
}


#' @title Batch Read .raw Files
#' @description Scans a folder and imports all matching .raw files into a list.
#'
#' @param folder_path String. Path to the folder.
#' @param pattern String. Regex pattern. Default is "\\.raw$".
#' @param recursive Logical. Default is FALSE.
#' @param ... Additional arguments passed to \code{read_thermal_raw}.
#'
#' @return A named list of "BioThermR" objects.
#' @export
#' @examples
#' \donttest{
#' # Example using raw thermal files
#' img_obj_list <- system.file("extdata",package = "BioThermR")
#' img_list <- read_thermal_batch(img_obj_list)
#' }
read_thermal_batch <- function(folder_path, pattern = "\\.raw$", recursive = FALSE, ...) {

  if (!dir.exists(folder_path)) {
    stop(paste("Error: Directory does not exist ->", folder_path))
  }

  # Get file list
  files <- list.files(folder_path, pattern = pattern, full.names = TRUE,
                      recursive = recursive, ignore.case = TRUE)

  if (length(files) == 0) {
    stop(paste("Error: No matching files found in", folder_path))
  }

  message(paste("Reading", length(files), "files..."))

  img_list <- lapply(files, function(f) {
    tryCatch({
      return(read_thermal_raw(f, ...))
    }, error = function(e) {
      warning(paste("Failed to read:", basename(f)))
      return(NULL)
    })
  })

  # Remove NULLs
  valid_indices <- !sapply(img_list, is.null)
  img_list <- img_list[valid_indices]

  # Assign names to the list elements based on filenames
  if (length(img_list) > 0) {
    names(img_list) <- basename(files[valid_indices])
  }

  message(paste("Batch read completed. Imported", length(img_list), "files."))
  return(img_list)
}

#' @title Create BioThermR Object from Data
#' @description Manually creates a BioThermR object from a numeric matrix or data frame.
#'              Useful for converting data loaded from other formats (Excel, CSV, txt) or simulation data.
#'
#' @param data A numeric matrix or data frame representing the thermal image grid (rows x cols).
#' @param name String. An identifier for this sample (e.g., filename or sample ID). Default is "Sample".
#'
#' @return A 'BioThermR' object.
#' @export
#' @examples
#' mat <- matrix(runif(160*120, 20, 40), nrow = 120, ncol = 160)
#' obj <- create_BioThermR(mat, name = "Simulation_01")
#' plot_thermal_heatmap(obj)
create_BioThermR <- function(data, name = "Sample") {

  # 1. Input Validation & Conversion
  mat <- data

  # If data frame, convert to matrix
  if (is.data.frame(data)) {
    # Ensure all columns are numeric
    if (!all(sapply(data, is.numeric))) {
      stop("Error: Data frame must contain only numeric values.")
    }
    mat <- as.matrix(data)
  }

  # Check if it is a matrix now
  if (!is.matrix(mat)) {
    stop("Error: Input 'data' must be a matrix or a data frame.")
  }

  # Check if numeric
  if (!is.numeric(mat)) {
    stop("Error: Matrix must contain numeric temperature values.")
  }

  # 2. Construct BioThermR Object
  obj <- structure(
    list(
      raw = mat,
      processed = mat,
      meta = list(
        filename = name,
        fullpath = NA,
        dims = dim(mat)
      )
    ),
    class = "BioThermR"
  )

  message(paste0("BioThermR object '", name, "' created. Dimensions: ", nrow(mat), "x", ncol(mat)))

  return(obj)
}


#' @title Convert BioThermR Object to EBImage Object
#' @description Extracts the thermal matrix (raw or processed) from a 'BioThermR' object and
#'              converts it into an 'EBImage' class object. This conversion is essential for
#'              applying advanced morphological operations (e.g., thresholding, watershed,
#'              labeling) provided by the 'EBImage' package.
#'
#' @param img_obj A 'BioThermR' object.
#' @param use_processed Logical. If \code{TRUE} (default), uses the 'processed' matrix (which might already have masks applied).
#'                      If \code{FALSE}, uses the 'raw' temperature matrix.
#' @param replace_na Numeric. The value to replace \code{NA}s with, as 'EBImage' does not support missing values.
#'                   Default is 0 (typically treated as background).
#' @param normalize Logical. If \code{TRUE}, scales the matrix values linearly to the range [0, 1].
#'                  This is highly recommended for visualization or standard thresholding algorithms
#'                  (like Otsu) in 'EBImage'. Default is \code{FALSE} (preserves actual temperature values).
#'
#' @return An \code{Image} object (defined in the 'EBImage' package).
#' @importFrom EBImage Image
#' @export
#' @examples
#' \donttest{
#' # Load data
#' mat <- matrix(runif(160*120, 20, 40), nrow = 120, ncol = 160)
#' obj <- create_BioThermR(mat, name = "Simulation_01")
#'
#' # Convert to EBImage format with normalization (for thresholding)
#' eb_norm <- as_EBImage(obj, normalize = TRUE)
#'
#' # Convert preserving temperature values (for calculation)
#' eb_temp <- as_EBImage(obj, normalize = FALSE)
#' }
as_EBImage <- function(img_obj, use_processed = TRUE, replace_na = 0, normalize = FALSE) {

  if (!inherits(img_obj, "BioThermR")) {
    stop("Error: Input must be a 'BioThermR' object.")
  }

  # Check dependency
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop("Error: Package 'EBImage' is required. Please install it via BiocManager::install('EBImage').")
  }

  # 1. Select Matrix
  mat <- if (use_processed) img_obj$processed else img_obj$raw

  if (is.null(mat)) stop("Error: Matrix is empty.")

  # 2. Handle NAs
  if (any(is.na(mat))) {
    mat[is.na(mat)] <- replace_na
  }

  # 3. Normalize (Optional)
  if (normalize) {
    rng <- range(mat, na.rm = TRUE)
    if (diff(rng) != 0) {
      mat <- (mat - rng[1]) / (rng[2] - rng[1])
    }
  }

  # 4. Convert to EBImage
  eb_img <- EBImage::Image(mat)

  return(eb_img)
}


#' @title Import EBImage Object into BioThermR
#' @description Converts an 'EBImage' class object back into the 'BioThermR' framework.
#'              This function is designed to integrate results from external morphological
#'              operations (e.g., watershed segmentation, filtering) back into the standard
#'              BioThermR analytical workflow.
#'
#' @details This function supports two modes:
#'          \itemize{
#'            \item \strong{Update Mode:} If \code{template_obj} is provided, the input image replaces
#'                  the \code{processed} matrix of the template. Metadata (filename, path) is preserved.
#'                  Statistics are reset to \code{NULL} as the data has changed.
#'            \item \strong{Create Mode:} If \code{template_obj} is \code{NULL}, a new 'BioThermR' object is created
#'                  from scratch using the input matrix.
#'          }
#'          If the input \code{eb_img} has more than 2 dimensions (e.g., color/RGB), only the first channel
#'          is utilized, and a warning is issued.
#'
#' @param eb_img An \code{Image} object (from the 'EBImage' package).
#' @param template_obj Optional. An existing 'BioThermR' object to update.
#'                     Providing this ensures that original metadata (like filenames and original raw data)
#'                     is retained. Default is \code{NULL}.
#' @param name String. The name assigned to the new object (used only if \code{template_obj} is \code{NULL}).
#'             Default is "Imported_EBImage".
#' @param mask_zero Logical. If \code{TRUE}, converts values of 0 in the imported matrix to \code{NA}.
#'                  This is particularly useful when importing binary masks where 0 represents the background,
#'                  as BioThermR uses \code{NA} to exclude background pixels from statistical calculations.
#'                  Default is \code{FALSE}.
#'
#' @return A 'BioThermR' object with the imported matrix stored in the \code{processed} slot.
#' @importFrom EBImage Image
#' @export
#' @examples
#' \donttest{
#' # Load data
#' mat <- matrix(runif(160*120, 20, 40), nrow = 120, ncol = 160)
#' obj <- create_BioThermR(mat, name = "Simulation_01")
#'
#' # Convert to EBImage format with normalization
#' eb_obj <- as_EBImage(obj, normalize = TRUE)
#'
#' # Convert to BioThermR from EBImage
#' new_obj <- from_EBImage(eb_obj)
#' }
from_EBImage <- function(eb_img, template_obj = NULL, name = "Imported_EBImage", mask_zero = FALSE) {

  # Check dependency
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop("Error: Package 'EBImage' is required.")
  }

  if (!inherits(eb_img, "Image")) {
    stop("Error: Input must be an EBImage 'Image' object.")
  }

  # 1. Convert EBImage to Matrix
  mat <- as.matrix(eb_img)

  # If it's a color image (3D), take the first channel or warn
  if (length(dim(mat)) > 2) {
    warning("Warning: Input image has >2 dimensions. Using the first channel.")
    mat <- mat[,,1]
  }

  # 2. Handle Masking (0 -> NA)
  if (mask_zero) {
    mat[mat == 0] <- NA
  }

  # 3. Create or Update Object
  if (!is.null(template_obj)) {
    if (!inherits(template_obj, "BioThermR")) {
      stop("Error: template_obj must be a 'BioThermR' object.")
    }

    # Check dimensions mismatch warning
    if (any(dim(mat) != dim(template_obj$processed))) {
      warning("Note: The dimensions of the imported image differ from the template object.")
      template_obj$meta$dims <- dim(mat)
    }

    # Update processed matrix
    template_obj$processed <- mat
    template_obj$stats <- NULL # Invalidate stats

    return(template_obj)

  } else {
    # -- New Object Mode --
    obj <- structure(
      list(
        raw = mat,
        processed = mat,
        stats = NULL,
        meta = list(
          filename = name,
          fullpath = NA,
          dims = dim(mat)
        )
      ),
      class = "BioThermR"
    )
    return(obj)
  }
}


#' @title Save BioThermR Data Objects to Disk
#' @description Serializes and saves a single 'BioThermR' object or a list of objects to a compressed
#'              .rds file. This ensures that all components—raw temperature matrices, processed
#'              masks, metadata, and calculation stats—are preserved accurately.
#' @param img_input A single 'BioThermR' class object or a list of 'BioThermR' objects
#'                  (e.g., the output from \code{read_thermal_batch}).
#' @param file_path String. The destination path (e.g., "results/obj.rds").
#'                  If the directory structure does not exist, it will be created automatically.
#'
#' @return None (invisible \code{NULL}). Prints a success message to the console upon completion.
#' @export
#' @examples
#' \donttest{
#' # Load data
#' mat <- matrix(runif(160*120, 20, 40), nrow = 120, ncol = 160)
#' obj <- create_BioThermR(mat, name = "Simulation_01")
#'
#' # Save a single object to a temporary directory
#' out_file <- file.path(tempdir(), "mouse_01.rds")
#' save_biothermr(obj, out_file)
#' }
save_biothermr <- function(img_input, file_path) {

  # 1. Validation
  is_single <- inherits(img_input, "BioThermR")
  is_list <- is.list(img_input) && inherits(img_input[[1]], "BioThermR")

  if (!is_single && !is_list) {
    stop("Error: Input must be a 'BioThermR' object or a list of them.")
  }

  # 2. Handle File Extension
  if (!grepl("\\.rds$", file_path, ignore.case = TRUE)) {
    file_path <- paste0(file_path, ".rds")
  }

  # 3. Create Directory if needed
  dir_name <- dirname(file_path)
  if (!dir.exists(dir_name) && dir_name != ".") {
    dir.create(dir_name, recursive = TRUE)
  }

  # 4. Save
  saveRDS(img_input, file = file_path)

  message(paste("Successfully saved BioThermR data to:", file_path))
}

#' @title Load BioThermR Data from Disk
#' @description Restores a previously saved 'BioThermR' object or a list of objects from a .rds file.
#'              This function is the counterpart to \code{\link{save_biothermr}} and allows you to
#'              resume analysis from a saved checkpoint.
#'
#' @details Upon loading, the function performs an automatic validation check to ensure the
#'          file contains a valid 'BioThermR' class instance (or a list of them).
#'          It provides feedback to the console regarding the type and quantity of objects loaded.
#'
#' @param file_path String. The full path to the .rds file (e.g., "results/experiment_data.rds").
#'
#' @return A single 'BioThermR' object or a list of 'BioThermR' objects, depending on the structure
#'         of the saved data.
#' @seealso \code{\link{save_biothermr}}
#' @export
load_biothermr <- function(file_path) {

  if (!file.exists(file_path)) {
    stop(paste("Error: File not found ->", file_path))
  }

  # Load data
  data <- readRDS(file_path)

  # Validation check
  is_single <- inherits(data, "BioThermR")
  is_list <- is.list(data) && length(data) > 0 && inherits(data[[1]], "BioThermR")

  if (is_single) {
    message(paste("Loaded 1 BioThermR object from", basename(file_path)))
  } else if (is_list) {
    message(paste("Loaded a list of", length(data), "BioThermR objects from", basename(file_path)))
  } else {
    warning("Loaded object does not appear to be a valid BioThermR object.")
  }

  return(data)
}
