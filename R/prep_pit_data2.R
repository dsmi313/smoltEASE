#' Format PIT detections for the six-cell guidance-efficiency model
#'
#' @description Standardizes the PTAGIS CSV columns used in the MY2025
#'   steelhead analysis. The main detection history and optional supplemental
#'   transport detections are returned in one event table. The source column
#'   keeps supplemental transport records from being treated as ordinary
#'   detection histories by \code{prep_ge_data2()}.
#'
#' @param pit Main PTAGIS detection-history data frame.
#' @param transport_data Optional supplemental transport detection data frame
#'   with the same PTAGIS columns as \code{pit}.
#' @param species Species name retained from \code{Species Name}. Default
#'   \code{"Steelhead"}. Set to \code{NULL} to skip this filter.
#' @param rear_type Character vector of rear-type codes retained from
#'   \code{Rear Type Code}. Default \code{"W"}. Set to \code{NULL} to skip
#'   this filter entirely. Known hatchery fish (\code{"H"}) are included only
#'   when explicitly supplied here.
#' @param include_unknown Logical. When \code{TRUE}, add rear-type code
#'   \code{"U"} to the codes selected by \code{rear_type}. Thus the defaults
#'   give a W-only analysis, while \code{include_unknown = TRUE} gives W+U and
#'   still excludes known hatchery fish. Default \code{FALSE}.
#' @param tz Time zone used to parse \code{Obs Time Value}. Default
#'   \code{"UTC"}.
#'
#' @return A data frame with columns \code{tag}, \code{site},
#'   \code{det_time}, \code{det_date}, \code{antenna},
#'   \code{release_site}, \code{rear_type}, and \code{source}.
#'
#' @export
prep_pit_data2 <- function(pit,
                           transport_data = NULL,
                           species = "Steelhead",
                           rear_type = "W",
                           include_unknown = FALSE,
                           tz = "UTC") {
  if (!is.logical(include_unknown) || length(include_unknown) != 1L ||
      is.na(include_unknown)) {
    stop("include_unknown must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(rear_type)) {
    if (!is.character(rear_type) || !length(rear_type) ||
        anyNA(rear_type) || any(!nzchar(trimws(rear_type)))) {
      stop("rear_type must be NULL or a nonempty character vector.",
           call. = FALSE)
    }
    rear_type <- unique(trimws(rear_type))
    if (include_unknown) rear_type <- unique(c(rear_type, "U"))
  }

  required <- c(
    "Tag Code", "Species Name", "Rear Type Code", "Site Name",
    "Obs Time Value", "Antenna ID", "Release Site Name"
  )
  missing_main <- setdiff(required, names(pit))
  if (length(missing_main)) {
    stop("pit is missing: ", paste(missing_main, collapse = ", "),
         call. = FALSE)
  }

  pit$.source <- "main"
  dat <- pit[, c(required, ".source"), drop = FALSE]

  if (!is.null(transport_data)) {
    missing_transport <- setdiff(required, names(transport_data))
    if (length(missing_transport)) {
      stop("transport_data is missing: ",
           paste(missing_transport, collapse = ", "), call. = FALSE)
    }
    transport_data$.source <- "transport"
    dat <- rbind(dat, transport_data[, c(required, ".source"), drop = FALSE])
  }

  if (!is.null(species)) {
    dat <- dat[!is.na(dat[["Species Name"]]) &
                 dat[["Species Name"]] == species, , drop = FALSE]
  }
  if (!is.null(rear_type)) {
    dat <- dat[!is.na(dat[["Rear Type Code"]]) &
                 trimws(as.character(dat[["Rear Type Code"]])) %in%
                   rear_type, , drop = FALSE]
  }
  if (!nrow(dat)) stop("No PIT records remain after filtering.", call. = FALSE)

  raw_time <- trimws(as.character(dat[["Obs Time Value"]]))
  det_time <- as.POSIXct(rep(NA_character_, length(raw_time)), tz = tz)
  for (fmt in c("%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M",
                "%m/%d/%Y %I:%M:%S %p", "%m/%d/%Y %I:%M %p",
                "%Y-%m-%d %H:%M:%S")) {
    unresolved <- which(is.na(det_time) & !is.na(raw_time) & nzchar(raw_time))
    if (!length(unresolved)) break
    det_time[unresolved] <- as.POSIXct(raw_time[unresolved], format = fmt,
                                      tz = tz)
  }
  if (anyNA(det_time)) {
    stop("Obs Time Value contains missing or unparseable timestamps.",
         call. = FALSE)
  }

  out <- data.frame(
    tag = trimws(as.character(dat[["Tag Code"]])),
    site = substr(trimws(as.character(dat[["Site Name"]])), 1L, 3L),
    det_time = det_time,
    det_date = as.Date(det_time),
    antenna = trimws(as.character(dat[["Antenna ID"]])),
    release_site = substr(trimws(as.character(dat[["Release Site Name"]])),
                          1L, 6L),
    rear_type = trimws(as.character(dat[["Rear Type Code"]])),
    source = dat$.source,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$tag) & nzchar(out$tag), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "rear_types_included") <- if (is.null(rear_type)) {
    sort(unique(out$rear_type))
  } else {
    rear_type
  }
  out
}
