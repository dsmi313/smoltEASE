#' @title Read and filter LGR smolt passage FPC data from Excel
#'
#' @description Reads the seasonal smolt passage collections spreadsheet
#'   (\code{SampleEndDate}, \code{SampleRate}, \code{SampleCount}, etc.) and
#'   applies the standard species-specific filters used to prepare the
#'   \code{passageData} input for \code{\link{SCRAPI}} / \code{\link{SCRAPI2}}.
#'
#'   Filters applied for both species:
#'   \itemize{
#'     \item \code{GearCode == gear_code} (default \code{"GC"})
#'     \item \code{Clip == "N"} (unclipped fish only)
#'     \item \code{cwt != "Y"} or \code{cwt} is \code{NA} (exclude CWT fish)
#'   }
#'
#'   Additional filter for \strong{steelhead} (\code{species = "sthd"}):
#'   \itemize{
#'     \item \code{SpecialSpeciesCode != "EF"} or is \code{NA}
#'       (exclude electrofished fish; other special codes retained)
#'   }
#'
#'   Additional filter for \strong{Chinook} (\code{species = "chnk"}):
#'   \itemize{
#'     \item \code{is.na(SpecialSpeciesCode)} (exclude all special-code fish)
#'   }
#'
#' @param path path to the Excel workbook
#'   (e.g. \code{"2025SMPCollectionsAtLGR_CH1&ST.xlsx"}).
#' @param species one of \code{"chnk"} or \code{"sthd"}.
#' @param sheet sheet name or index passed to \code{readxl::read_excel}.
#'   Default \code{1} (first sheet).
#' @param gear_code gear code to retain. Default \code{"GC"}.
#'
#' @return A filtered data frame with columns: \code{SampleEndDate},
#'   \code{batch}, \code{SubBatch}, \code{GearCode}, \code{SampleRate},
#'   \code{Species}, \code{SpecialSpeciesCode}, \code{Clip}, \code{cwt},
#'   \code{SampleCount}, \code{CollectionCount}.
#'
#' @details Requires the \pkg{readxl} package
#'   (\code{install.packages("readxl")}).
#'
#' @examples
#' \dontrun{
#' chnk_fpc <- read_fpc("2025SMPCollectionsAtLGR_CH1&ST.xlsx", species = "chnk")
#' sthd_fpc <- read_fpc("2025SMPCollectionsAtLGR_CH1&ST.xlsx", species = "sthd")
#' }
#'
#' @export
read_fpc <- function(path, species, sheet = 1, gear_code = "GC") {

  species <- match.arg(species, c("chnk", "sthd"))

  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Package 'readxl' is required. Install with install.packages('readxl').")

  raw <- readxl::read_excel(path, sheet = sheet)

  if (species == "sthd") {
    out <- raw[
      raw$Species == "ST" &
      raw$GearCode == gear_code &
      (is.na(raw$SpecialSpeciesCode) | raw$SpecialSpeciesCode != "EF") &
      raw$Clip == "N" &
      (is.na(raw$cwt) | raw$cwt != "Y"),
    ]
  } else {
    out <- raw[
      raw$Species == "CH1" &
      raw$GearCode == gear_code &
      is.na(raw$SpecialSpeciesCode) &
      raw$Clip == "N" &
      (is.na(raw$cwt) | raw$cwt != "Y"),
    ]
  }

  out
}
