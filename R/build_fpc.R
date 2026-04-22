#' @title Build passage FPC data frame by joining week numbers from trap data
#'
#' @description Joins calendar week numbers from the formatted trap data onto
#'   the filtered FPC passage records, producing a data frame suitable for use
#'   as \code{passageData} in \code{\link{SCRAPI2}} (after the user has added
#'   \code{GuidanceEfficiency} and \code{Collapse} columns separately).
#'
#'   Week numbers are pulled from the \code{WeekNumber} column of \code{trap}
#'   by matching \code{SampleEndDate} in \code{fpc} to \code{CollectionDate}
#'   in \code{trap}. Days in \code{fpc} with no matching trap date are dropped
#'   (\code{inner_join} behaviour); a message reports the count of dropped rows.
#'
#' @param fpc data frame returned by \code{\link{read_fpc}}.
#' @param trap data frame returned by \code{\link{lgr2SCRAPI}} (must contain
#'   \code{CollectionDate} and \code{WeekNumber} columns).
#' @param export_csv optional file path. When supplied the result is written as
#'   a CSV with dates formatted \code{"\%m/\%d/\%Y"} (the format expected by
#'   \code{\link{SCRAPI2}}'s \code{dateFormat} argument).
#'
#' @return A data frame with columns \code{SampleEndDate} (Date),
#'   \code{Week} (integer), and all remaining columns from \code{fpc}.
#'   \code{GuidanceEfficiency} and \code{Collapse} are intentionally absent —
#'   add \code{GuidanceEfficiency} from your Bayesian posterior mean and
#'   \code{Collapse} via the \code{strata=} argument of \code{\link{SCRAPI2}}.
#'
#' @seealso \code{\link{read_fpc}}, \code{\link{lgr2SCRAPI}},
#'   \code{\link{check_fpc}}, \code{\link{check_trap}},
#'   \code{\link{genstock_by_week}}
#'
#' @examples
#' \dontrun{
#' chnk_fpc <- build_fpc(chnk_fpc_raw, CHNK.trap,
#'                        export_csv = "MY2025CHNKsmolts_FPCdata.csv")
#' sthd_fpc <- build_fpc(sthd_fpc_raw, STHD.trap,
#'                        export_csv = "MY2025STHDsmolts_FPCdata.csv")
#' }
#'
#' @export
build_fpc <- function(fpc, trap, export_csv = NULL) {

  wk_lookup <- unique(data.frame(
    CollectionDate = as.Date(trap$CollectionDate),
    Week           = as.integer(trap$WeekNumber),
    stringsAsFactors = FALSE
  ))

  fpc2 <- fpc
  fpc2$SampleEndDate <- as.Date(fpc2$SampleEndDate)

  n_before <- nrow(fpc2)
  out <- merge(fpc2, wk_lookup,
               by.x = "SampleEndDate", by.y = "CollectionDate",
               all.x = FALSE)   # inner join — drop unmatched FPC dates

  n_dropped <- n_before - nrow(out)
  if (n_dropped > 0)
    message(n_dropped, " FPC row(s) dropped: no matching date in trap data.")

  # put Week immediately after SampleEndDate
  other_cols <- setdiff(names(out), c("SampleEndDate", "Week"))
  out <- out[, c("SampleEndDate", "Week", other_cols), drop = FALSE]

  if (!is.null(export_csv)) {
    out_write <- out
    out_write$SampleEndDate <- format(out_write$SampleEndDate, "%m/%d/%Y")
    write.csv(out_write, export_csv, row.names = FALSE, na = "NA")
    message("Written to: ", export_csv)
  }

  out
}


#' @title QC check for FPC passage data
#'
#' @description Checks critical columns of an FPC data frame for \code{NA}
#'   values and prints a summary. Designed to be run before \code{\link{SCRAPI2}}.
#'
#' @param df data frame returned by \code{\link{build_fpc}} (or similar).
#' @param label character label printed in the header (e.g. \code{"CHNK"}).
#' @param critical character vector of column names that must be non-NA.
#'
#' @return Invisible \code{NULL}. Results are printed to the console.
#' @export
check_fpc <- function(df, label = "",
                       critical = c("SampleEndDate", "Week",
                                    "SampleRate", "SampleCount",
                                    "CollectionCount")) {
  cat("===", label, "FPC ===  rows:", nrow(df), "\n")
  for (col in critical) {
    if (!col %in% names(df)) { cat("  MISSING column:", col, "\n"); next }
    n <- sum(is.na(df[[col]]))
    if (n > 0) cat("  FAIL", col, "NA:", n, "\n")
  }
  bad_cols <- intersect(critical, names(df))
  bad <- df[rowSums(is.na(df[, bad_cols, drop = FALSE])) > 0, ]
  if (nrow(bad)) {
    cat("  rows with NA in critical columns:\n")
    print(bad)
  } else {
    cat("  all critical columns clean\n")
  }
  invisible(NULL)
}


#' @title QC check for formatted trap data
#'
#' @description Checks critical columns of a formatted trap data frame for
#'   \code{NA} values and reports \code{GenStock} and (optionally)
#'   \code{BioScaleFinalAge} missingness rates.
#'
#' @param df data frame returned by \code{\link{lgr2SCRAPI}}.
#' @param label character label printed in the header (e.g. \code{"STHD"}).
#' @param check_age logical. When \code{TRUE} (steelhead) also reports
#'   unreadable \code{BioScaleFinalAge} values (missing, \code{"N:A"},
#'   or containing \code{"?"}).
#'
#' @return Invisible \code{NULL}. Results are printed to the console.
#' @export
check_trap <- function(df, label = "", check_age = FALSE) {
  cat("===", label, "trap ===  rows:", nrow(df), "\n")
  for (col in c("CollectionDate", "WeekNumber", "Rear")) {
    if (!col %in% names(df)) { cat("  MISSING column:", col, "\n"); next }
    n <- sum(is.na(df[[col]]))
    if (n > 0) cat("  FAIL", col, "NA:", n, "\n")
  }
  if ("GenStock" %in% names(df)) {
    n_gs <- sum(is.na(df$GenStock) | df$GenStock %in% c("NG", "N:A", ""))
    cat("  GenStock missing/NG:", n_gs,
        "(", round(100 * n_gs / nrow(df), 1), "%)\n")
  }
  if (check_age && "BioScaleFinalAge" %in% names(df)) {
    age <- as.character(df$BioScaleFinalAge)
    n_age <- sum(is.na(age) | age %in% c("N:A", "") | grepl("\\?", age))
    cat("  BioScaleFinalAge missing/unreadable:", n_age,
        "(", round(100 * n_age / nrow(df), 1), "%)\n")
  }
  invisible(NULL)
}


#' @title Check date alignment between FPC and trap data
#'
#' @description Compares the unique dates in \code{fpc} (\code{SampleEndDate})
#'   and \code{trap} (\code{CollectionDate}) and reports mismatches. Use this
#'   before running \code{\link{SCRAPI2}} to confirm the two data sources cover
#'   the same trapping days.
#'
#' @param fpc data frame with a \code{SampleEndDate} column (Date or character).
#' @param trap data frame with a \code{CollectionDate} column (Date or character).
#' @param label character label printed in the header.
#'
#' @return Invisible \code{NULL}. Results are printed to the console.
#' @export
check_date_alignment <- function(fpc, trap, label = "") {
  f <- unique(as.Date(fpc$SampleEndDate))
  t <- unique(as.Date(trap$CollectionDate))
  cat("===", label, "date alignment ===\n")
  miss_fpc <- setdiff(t, f)
  miss_trap <- setdiff(f, t)
  cat("  trap dates missing from FPC:", length(miss_fpc), "\n")
  if (length(miss_fpc)) cat("   ", format(sort(miss_fpc)), "\n")
  cat("  FPC dates missing from trap:", length(miss_trap), "\n")
  if (length(miss_trap)) cat("   ", format(sort(miss_trap)), "\n")
  invisible(NULL)
}


#' @title Tabulate GenStock by week with a usable-fish total column
#'
#' @description Returns a wide table of fish counts by \code{WeekNumber} and
#'   \code{GenStock}, with a \code{Total} column showing the number of fish
#'   usable for SCRAPI2 composition estimation in each week. For steelhead
#'   that means GenStock \strong{and} \code{fwAge} are both present; for
#'   Chinook it means GenStock is present.
#'
#' @param trap data frame returned by \code{\link{lgr2SCRAPI}}.
#' @param species one of \code{"chnk"} or \code{"sthd"}.
#' @param rear_only character. Filter to this rear type before tabulating.
#'   Default \code{"W"}. Set to \code{NULL} to include all fish.
#'
#' @return A data frame with one row per week. Columns are \code{Week},
#'   one column per GenStock code, and \code{Total} (usable fish per week).
#' @export
genstock_by_week <- function(trap, species, rear_only = "W") {

  species <- match.arg(species, c("chnk", "sthd"))

  df <- trap
  if (!is.null(rear_only) && "Rear" %in% names(df))
    df <- df[!is.na(df$Rear) & df$Rear == rear_only, ]

  counts <- table(df$WeekNumber, df$GenStock)
  out    <- as.data.frame.matrix(counts)
  out    <- cbind(Week = as.integer(rownames(out)), out)
  rownames(out) <- NULL

  # Total = fish usable for composition: GenStock present (+ fwAge for sthd)
  has_gs  <- !is.na(df$GenStock) & df$GenStock != "NA"
  if (species == "sthd") {
    usable <- has_gs & !is.na(df$fwAge)
  } else {
    usable <- has_gs
  }
  totals_per_week <- tapply(usable, df$WeekNumber, sum)
  out$Total <- as.integer(totals_per_week[as.character(out$Week)])

  out
}
