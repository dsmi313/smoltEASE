#' Three-row GRS detection histories for a day/night GE fit
#'
#' Counts each eligible spilled fish once on each row, using detections within
#' ten minutes of its first GRS observation. Uses the same eligible tags and
#' day/night classification as the GE data. Antennas within a row are pooled.
#' @param dat_up PIT events with tag, site, antenna and det_time columns.
#' @param ge_dn Named list returned by prep_ge_data_daynight(), one per rear.
#' @param passage_minutes Maximum elapsed minutes for a single spill passage.
#' @return List with counts (rear x week x period x seven nonzero histories),
#'   histories, rear_levels, weeks and patterns. Pass to fit_ge_rear_daynight().
#' @export
prep_grs_detection <- function(dat_up, ge_dn, passage_minutes = 10) {
  if (!all(c("tag", "site", "antenna", "det_time") %in% names(dat_up)))
    stop("dat_up needs tag, site, antenna and det_time.", call. = FALSE)
  if (!is.list(ge_dn) || is.null(names(ge_dn)) || anyDuplicated(names(ge_dn)))
    stop("ge_dn must be a named list with unique rear names.", call. = FALSE)
  if (length(passage_minutes) != 1L || !is.finite(passage_minutes) || passage_minutes <= 0)
    stop("passage_minutes must be positive.", call. = FALSE)
  weeks <- ge_dn[[1]]$base$weeks
  if (is.null(weeks)) weeks <- as.integer(sub("^wk", "", rownames(ge_dn[[1]]$n_day)))
  patterns <- as.matrix(expand.grid(U = 0:1, M = 0:1, D = 0:1))[-1, , drop = FALSE]
  keys <- apply(patterns, 1, paste0, collapse = "")
  counts <- array(0L, c(length(ge_dn), length(weeks), 2L, 7L),
                  dimnames = list(names(ge_dn), as.character(weeks), c("day", "night"), keys))
  spill <- dat_up[dat_up$site == "GRS" & !is.na(dat_up$det_time), , drop = FALSE]
  by_tag <- split(spill, spill$tag)
  output <- list()
  for (r in seq_along(ge_dn)) {
    tags <- ge_dn[[r]]$tags
    tags <- tags[tags$cell %in% c(3L, 4L), , drop = FALSE]
    for (i in seq_len(nrow(tags))) {
      x <- by_tag[[as.character(tags$tag[i])]]
      if (is.null(x) || !nrow(x)) stop("Eligible spill fish has no GRS reads: ", tags$tag[i], call. = FALSE)
      x <- x[as.numeric(difftime(x$det_time, min(x$det_time), units = "mins")) <= passage_minutes, , drop = FALSE]
      antenna <- toupper(trimws(as.character(x$antenna)))
      if (any(!antenna %in% c("01", "02", "03", "04", "05", "06", "07", "08", "09", "0A", "0B")))
        stop("Unknown GRS antenna; normalize antenna IDs before calling.", call. = FALSE)
      h <- paste0(as.integer(any(antenna %in% c("01", "02", "03", "04"))),
                  as.integer(any(antenna %in% c("05", "06", "07"))),
                  as.integer(any(antenna %in% c("08", "09", "0A", "0B"))))
      s <- match(tags$week[i], weeks); t <- match(tags$period[i], c("day", "night"))
      if (is.na(s) || is.na(t)) stop("Invalid week or period in eligible fish.", call. = FALSE)
      counts[r, s, t, match(h, keys)] <- counts[r, s, t, match(h, keys)] + 1L
      output[[length(output) + 1L]] <- data.frame(tag = tags$tag[i], rear = names(ge_dn)[r],
        week = tags$week[i], period = tags$period[i], history = h)
    }
  }
  list(counts = counts, histories = if (length(output)) do.call(rbind, output) else data.frame(),
       rear_levels = names(ge_dn), weeks = weeks, patterns = patterns)
}
