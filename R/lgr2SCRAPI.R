#' @title Format LGD trap data for SCRAPI / SCRAPI2
#'
#' @description Reads raw biological data from the LGD trapping database (as
#'   returned by \code{\link{query_lgd_trap}} or a saved CSV) and prepares it
#'   for use as the \code{smoltData} argument of \code{\link{SCRAPI}} or
#'   \code{\link{SCRAPI2}}.
#'
#'   Key transformations applied:
#'   \itemize{
#'     \item \code{Rear} is created from \code{LGDRear} (W / H; used by SCRAPI
#'       to filter to wild fish via \code{RTYPE = "W"}).
#'     \item \code{MPG} is populated from \code{CHNMPG} (chinook) or
#'       \code{STHDMPG} (steelhead).
#'     \item For steelhead, \code{fwAge} is derived as the integer freshwater
#'       years from \code{BioScaleFinalAge} (the portion before the first
#'       \code{.} in the age code, e.g. \code{"2.1"} → \code{2}).
#'     \item \code{GenStock} codes are corrected to the 6-character format
#'       expected by SCRAPI.
#'   }
#'
#' @param input a data frame of raw trap records (from
#'   \code{\link{query_lgd_trap}}), or a path to a CSV file of the same.
#' @param species one of \code{"chnk"} or \code{"sthd"}.
#' @param exportFile optional output file path (with or without \code{.csv}
#'   extension). When supplied the formatted data frame is written to disk.
#' @param check_codes logical. When \code{TRUE} (default) prints a sorted
#'   frequency table of \code{GenStock} values after correction so codes can
#'   be verified before running SCRAPI.
#'
#' @return A data frame with the columns needed by \code{\link{SCRAPI2}}:
#'   \code{MasterID}, \code{CollectionDate}, \code{WeekNumber}, \code{SRR},
#'   \code{Rear}, \code{GenRear}, \code{GenStock}, \code{GenSex},
#'   \code{BioScaleFinalAge}, \code{BioSamplesID}, \code{SpawnYear},
#'   \code{LGDFLmm}, \code{LGDLifeStage}, \code{GenRun}, \code{MPG},
#'   \code{GenStockProb}, \code{GenParentHatchery}, \code{GenBY},
#'   \code{LGDMarkAD}, and (steelhead only) \code{fwAge}.
#'
#' @details
#' \strong{GenStock code corrections applied automatically:}
#'
#' The LGD database occasionally returns codes longer than 6 characters.
#' SCOBI / SCRAPI expects exactly 6. The following substitutions are applied:
#'
#' \tabular{ll}{
#'   Raw code \tab Corrected \cr
#'   LOWSALM  \tab LOSALM \cr
#'   LOWCLWR  \tab LOCLWR \cr
#'   LOWGRAN  \tab LOGRAN \cr
#'   LOWSNAK  \tab LOSNAK \cr
#' }
#'
#' Any code that is not exactly 6 characters after correction triggers a
#' warning listing the offending values — fix these before running SCRAPI.
#'
#' @examples
#' \dontrun{
#' # From a saved CSV
#' chnk <- lgr2SCRAPI("MY2025.CHNKsmolts_trapBioData_allVariables.csv",
#'                     species    = "chnk",
#'                     exportFile = "MY2025CHNK.trapData_formatted")
#'
#' # Directly from query result
#' sthd_raw <- query_lgd_trap("MY2025", srr_prefix = "3", life_stage = "JV")
#' sthd     <- lgr2SCRAPI(sthd_raw, species = "sthd",
#'                         exportFile = "MY2025STHD.trapData_formatted")
#' }
#'
#' @export
lgr2SCRAPI <- function(input,
                        species,
                        exportFile  = NULL,
                        check_codes = TRUE) {

  species <- match.arg(species, c("chnk", "sthd"))

  # ---- read input ----------------------------------------------------------
  if (is.character(input)) {
    dat <- read.csv(input, header = TRUE, stringsAsFactors = FALSE)
  } else {
    dat <- as.data.frame(input, stringsAsFactors = FALSE)
  }

  # ---- derive Rear from LGDRear --------------------------------------------
  # LGDRear has clean W / H values; SCRAPI filters on this via RTYPE= ("W").
  dat$Rear <- dat$LGDRear

  # ---- derive MPG from species-specific column -----------------------------
  if (species == "chnk" && "CHNMPG" %in% names(dat)) {
    dat$MPG <- dat$CHNMPG
  } else if (species == "sthd" && "STHDMPG" %in% names(dat)) {
    dat$MPG <- dat$STHDMPG
  }

  # ---- derive fwAge for steelhead ------------------------------------------
  # BioScaleFinalAge uses the format "fw:sw" (e.g. "2:0" = 2 freshwater years,
  # 0 saltwater years). fwAge is the integer portion before the ":".
  # "N:A" is the database missing-value code -> NA.
  if (species == "sthd" && "BioScaleFinalAge" %in% names(dat)) {
    age_raw <- as.character(dat$BioScaleFinalAge)
    age_raw[age_raw %in% c("N:A", "NA", "")] <- NA
    dat$fwAge <- suppressWarnings(
      as.integer(sub(":.*$", "", age_raw))
    )
  }

  # ---- standardise GenStock codes ------------------------------------------
  corrections <- c(
    LOWSALM = "LOSALM",
    LOWCLWR = "LOCLWR",
    LOWGRAN = "LOGRAN",
    LOWSNAK = "LOSNAK",
    NG      = "NA"       # no-genetics fish excluded by SCRAPI (Primary != "NA")
  )

  # Known valid codes that are legitimately != 6 characters
  length_exempt <- c("FALL")

  if ("GenStock" %in% names(dat)) {
    bad <- dat$GenStock %in% names(corrections)
    if (any(bad, na.rm = TRUE))
      dat$GenStock[bad] <- corrections[dat$GenStock[bad]]

    non_na  <- dat$GenStock[!is.na(dat$GenStock) &
                             dat$GenStock != "NA" &
                             !dat$GenStock %in% length_exempt]
    bad_len <- unique(non_na[nchar(non_na) != 6])
    if (length(bad_len) > 0)
      warning("GenStock codes with length != 6 after correction: ",
              paste(bad_len, collapse = ", "),
              "\nFix these before running SCRAPI.")
  }

  # ---- select output columns -----------------------------------------------
  core_cols <- c(
    "MasterID", "CollectionDate", "WeekNumber", "SRR",
    "Rear", "GenRear",
    "GenStock", "GenSex",
    "BioScaleFinalAge", "BioSamplesID",
    "SpawnYear", "LGDFLmm", "LGDLifeStage",
    "GenRun", "MPG",
    "GenStockProb", "GenParentHatchery", "GenBY",
    "LGDMarkAD"
  )
  if (species == "sthd") core_cols <- c(core_cols, "fwAge")

  keep <- intersect(core_cols, names(dat))
  out  <- dat[, keep, drop = FALSE]

  # ---- print GenStock check ------------------------------------------------
  if (check_codes) {
    cat("\nGenStock frequency after correction (", species, "):\n", sep = "")
    print(sort(table(out$GenStock), decreasing = TRUE))
  }

  # ---- export --------------------------------------------------------------
  if (!is.null(exportFile)) {
    path <- if (grepl("\\.csv$", exportFile)) exportFile else paste0(exportFile, ".csv")
    write.csv(out, path, row.names = FALSE)
    cat("\nFormatted data written to:", path, "\n")
  }

  invisible(out)
}
