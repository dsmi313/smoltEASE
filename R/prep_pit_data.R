#' @title Parse PTAGIS Complete Tag History data for use in prep_ge_data()
#'
#' @description Converts a PTAGIS Complete Tag History (CTH) query result into
#'   the \code{dat_up} format required by \code{\link{prep_ge_data}}: one row
#'   per detection event with columns \code{tag}, \code{site}, and
#'   \code{det_date}.
#'
#'   Site codes are extracted from the \code{Event.Site.Name} field, which
#'   PTAGIS formats as \code{"CODE - Full Site Name"}.
#'
#' @param pit data frame of PTAGIS CTH output with columns \code{Tag.Code},
#'   \code{Event.Type.Name}, \code{Event.Site.Name},
#'   \code{Event.Date.MMDDYYYY}, \code{Mark.Species.Name},
#'   \code{Migration.Year.YYYY}, \code{Mark.Rear.Type.Name},
#'   \code{Event.Life.Stage.Value}.
#' @param species character. Keep only records matching this species name in
#'   \code{Mark.Species.Name}. Default \code{"Steelhead"}. Set to \code{NULL}
#'   to skip species filter.
#' @param mig_year integer or character. Keep only records for this migration
#'   year. Default \code{NULL} (no filter).
#' @param rear_type character. Keep only records matching this value in
#'   \code{Mark.Rear.Type.Name}. Default \code{"Wild Fish or Natural Production"}.
#'   Set to \code{NULL} to skip.
#' @param life_stage character. Keep only Mark events matching this life stage
#'   in \code{Event.Life.Stage.Value}, used to identify fish tagged as juveniles
#'   upstream. E.g. \code{"Juvenile"} or \code{"Parr"}. Set to \code{NULL} to
#'   keep all tags regardless of life stage at marking. Default \code{NULL}.
#' @param event_types character vector of \code{Event.Type.Name} values to
#'   retain. Default \code{c("Mark", "Observation")} — both are needed so
#'   downstream code can reconstruct each fish's LGR route.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{tag}{PIT tag code (\code{Tag.Code})}
#'     \item{site}{3-letter PTAGIS site code extracted from
#'       \code{Event.Site.Name}}
#'     \item{det_date}{detection date as \code{Date}}
#'     \item{mark_rkm}{numeric RKM of the mark site (from
#'       \code{Mark.Site.RKM.Total} or \code{Mark.Site.RKM.Value}); \code{NA}
#'       if no RKM column is present. Used by \code{\link{prep_ge_data}} to
#'       restrict Pool A to upstream-tagged fish.}
#'   }
#'
#' @details The function only keeps fish that appear in \code{event_types}.
#'   For \code{\link{prep_ge_data}}, \code{dat_up} should contain all
#'   detection events (Mark + Observation) so that each fish's first LGR
#'   route and first downstream detection can be identified.
#'
#' @examples
#' \dontrun{
#' dat_up <- prep_pit_data(MY25.pit, species = "Steelhead", mig_year = 2025)
#' ge_data <- prep_ge_data(dat_up, strat_assign = sthd_strata,
#'                          spill_data = spill, species = "sthd")
#' }
#'
#' @export
prep_pit_data <- function(pit,
                           species     = "Steelhead",
                           mig_year    = NULL,
                           rear_type   = "Wild Fish or Natural Production",
                           life_stage  = NULL,
                           event_types = c("Mark", "Observation")) {

  dat <- pit

  # ---- filters -------------------------------------------------------------
  if (!is.null(species) && "Mark.Species.Name" %in% names(dat))
    dat <- dat[!is.na(dat$Mark.Species.Name) &
               dat$Mark.Species.Name == species, ]

  if (!is.null(mig_year) && "Migration.Year.YYYY" %in% names(dat))
    dat <- dat[!is.na(dat$Migration.Year.YYYY) &
               as.character(dat$Migration.Year.YYYY) == as.character(mig_year), ]

  if (!is.null(rear_type) && "Mark.Rear.Type.Name" %in% names(dat))
    dat <- dat[!is.na(dat$Mark.Rear.Type.Name) &
               dat$Mark.Rear.Type.Name == rear_type, ]

  if (!is.null(life_stage) && "Event.Life.Stage.Value" %in% names(dat)) {
    # apply life stage filter only to Mark events; keep all Observations
    is_mark <- dat$Event.Type.Name == "Mark"
    keep    <- !is_mark | (!is.na(dat$Event.Life.Stage.Value) &
                            dat$Event.Life.Stage.Value == life_stage)
    dat <- dat[keep, ]
    # then drop tags whose Mark event was filtered out
    kept_tags <- unique(dat$Tag.Code[dat$Event.Type.Name == "Mark"])
    dat <- dat[dat$Tag.Code %in% kept_tags, ]
  }

  if (!is.null(event_types) && "Event.Type.Name" %in% names(dat))
    dat <- dat[dat$Event.Type.Name %in% event_types, ]

  if (nrow(dat) == 0)
    stop("No records remain after filtering. Check species, mig_year, ",
         "rear_type, and life_stage arguments.")

  # ---- extract site code ---------------------------------------------------
  # Event.Site.Name format: "CODE - Full Site Name"
  site_raw <- dat$Event.Site.Name
  dat$site <- trimws(sub(" -.*$", "", site_raw))

  # ---- parse date ----------------------------------------------------------
  dat$det_date <- as.Date(dat$Event.Date.MMDDYYYY, format = "%m/%d/%Y")

  # ---- extract mark RKM (passed through for use in prep_ge_data) ----------
  rkm_col  <- intersect(c("Mark.Site.RKM.Total", "Mark.Site.RKM.Value"),
                        names(dat))[1]
  mark_rkm <- if (!is.na(rkm_col))
    suppressWarnings(as.numeric(dat[[rkm_col]]))
  else
    rep(NA_real_, nrow(dat))

  # ---- select output -------------------------------------------------------
  data.frame(
    tag      = dat$Tag.Code,
    site     = dat$site,
    det_date = dat$det_date,
    mark_rkm = mark_rkm,
    stringsAsFactors = FALSE
  )
}
