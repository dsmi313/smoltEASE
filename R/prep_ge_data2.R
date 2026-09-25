#' Construct six-cell data for the guidance-efficiency model
#'
#' @description Converts the standardized event table returned by
#'   \code{prep_pit_data2()} into the six mutually exclusive detection-history
#'   cells required by \code{fit_ge_model2()}. It also calculates weekly LGR
#'   and LGS percent spill, weekly mean LGR outflow, and the coarser parent
#'   groups used by the nested model.
#'
#' @param dat_up Data frame returned by \code{prep_pit_data2()}.
#' @param spill_data Hourly spill data containing \code{Site},
#'   \code{DateTime}, \code{HourlySpill}, and \code{HourlyFlow}. Multiple
#'   seasonal CSVs may be combined with \code{rbind()} before calling.
#' @param weeks Increasing ISO week numbers to fit.
#' @param adult_sites Six-character release-site codes treated as adult
#'   releases and excluded.
#' @param route_order Named numeric vector describing downstream route order.
#' @param transport_antennas Antenna IDs identifying transported fish at GRJ.
#' @param parent_floor Adjacent weeks are pooled until both \code{c4+c5} and
#'   \code{c3+c4} exceed this value. Default 10.
#' @param alpha_phi_mean,alpha_phi_sd,beta_phi_mean,beta_phi_sd Phi-prior
#'   values passed through to \code{fit_ge_model2()}.
#' @param tz Time zone used to parse spill timestamps. Default \code{"UTC"}.
#'
#' @details The six cells are: c1 GRJ only, c2 GRJ and GOJ, c3 GRS only,
#'   c4 GRS and GOJ, c5 GOJ only, and c6 transported. Supplemental transport
#'   records can flag a tag as transported but do not create a detection
#'   history for a tag absent from the main PIT data.
#'
#' @return A six-cell list accepted directly by \code{fit_ge_model2()}.
#'
#' @export
prep_ge_data2 <- function(
    dat_up,
    spill_data,
    weeks,
    adult_sites = c("LGRLDR", "PRDLD1", "SHERFT", "BONAFF"),
    route_order = c(GRJ = 1, GRS = 1, GOJ = 2, LMJ = 3, ICH = 4),
    transport_antennas = c("61", "62"),
    parent_floor = 10L,
    alpha_phi_mean = 0,
    alpha_phi_sd = 2,
    beta_phi_mean = 0,
    beta_phi_sd = 1,
    tz = "UTC") {
  required_pit <- c("tag", "site", "det_time", "antenna",
                    "release_site", "source")
  missing_pit <- setdiff(required_pit, names(dat_up))
  if (length(missing_pit)) {
    stop("dat_up is missing: ", paste(missing_pit, collapse = ", "),
         ". Use prep_pit_data2() first.", call. = FALSE)
  }
  required_spill <- c("Site", "DateTime", "HourlySpill", "HourlyFlow")
  missing_spill <- setdiff(required_spill, names(spill_data))
  if (length(missing_spill)) {
    stop("spill_data is missing: ", paste(missing_spill, collapse = ", "),
         call. = FALSE)
  }
  if (!is.numeric(weeks) || !length(weeks) || any(!is.finite(weeks)) ||
      any(weeks != floor(weeks)) || is.unsorted(weeks, strictly = TRUE)) {
    stop("weeks must be increasing integer ISO week numbers.", call. = FALSE)
  }
  if (length(parent_floor) != 1L || !is.finite(parent_floor) ||
      parent_floor < 0 || parent_floor != floor(parent_floor)) {
    stop("parent_floor must be one nonnegative integer.", call. = FALSE)
  }

  main <- dat_up[dat_up$source == "main", , drop = FALSE]
  route <- main[main$site %in% names(route_order), , drop = FALSE]
  route$route_number <- unname(route_order[route$site])
  route <- route[order(route$tag, route$det_time), , drop = FALSE]
  upstream_tags <- unique(route$tag[
    route$route_number < ave(route$route_number, route$tag, FUN = cummax)
  ])
  adult_release_tags <- unique(main$tag[main$release_site %in% adult_sites])
  excluded_tags <- union(upstream_tags, adult_release_tags)

  main <- main[!(main$tag %in% excluded_tags), , drop = FALSE]
  transported_tags <- unique(dat_up$tag[
    !(dat_up$tag %in% excluded_tags) & dat_up$site == "GRJ" &
      dat_up$antenna %in% as.character(transport_antennas)
  ])

  histories <- main[main$site %in% c("GRJ", "GRS", "GOJ"), , drop = FALSE]
  tags <- sort(unique(histories$tag))
  if (!length(tags)) stop("No usable route histories remain.", call. = FALSE)

  first_grj <- tapply(histories$det_time[histories$site == "GRJ"],
                      histories$tag[histories$site == "GRJ"], min)
  first_grs <- tapply(histories$det_time[histories$site == "GRS"],
                      histories$tag[histories$site == "GRS"], min)
  first_goj <- tapply(histories$det_time[histories$site == "GOJ"],
                      histories$tag[histories$site == "GOJ"], min)
  GRJ <- as.POSIXct(unname(first_grj[match(tags, names(first_grj))]),
                    origin = "1970-01-01", tz = tz)
  GRS <- as.POSIXct(unname(first_grs[match(tags, names(first_grs))]),
                    origin = "1970-01-01", tz = tz)
  GOJ <- as.POSIXct(unname(first_goj[match(tags, names(first_goj))]),
                    origin = "1970-01-01", tz = tz)

  LGR <- pmin(GRJ, GRS, na.rm = TRUE)
  paired <- !is.na(LGR) & !is.na(GOJ)
  goj_lag_days <- median(as.numeric(difftime(GOJ[paired], LGR[paired],
                                             units = "days")))
  fish_week <- as.integer(format(as.Date(LGR), "%V"))
  goj_only <- is.na(LGR) & !is.na(GOJ)
  if (any(goj_only) && !is.finite(goj_lag_days)) {
    stop("Cannot assign weeks to GOJ-only fish.", call. = FALSE)
  }
  fish_week[goj_only] <- as.integer(format(
    as.Date(GOJ[goj_only] - goj_lag_days * 86400), "%V"
  ))

  has_grj <- !is.na(GRJ)
  has_grs <- !is.na(GRS)
  has_goj <- !is.na(GOJ)
  has_grs[has_grj & has_grs] <- FALSE
  cell <- rep(NA_integer_, length(tags))
  cell[ has_grj & !has_goj] <- 1L
  cell[ has_grj &  has_goj] <- 2L
  cell[!has_grj &  has_grs & !has_goj] <- 3L
  cell[!has_grj &  has_grs &  has_goj] <- 4L
  cell[!has_grj & !has_grs &  has_goj] <- 5L
  cell[tags %in% transported_tags] <- 6L

  keep <- !is.na(cell) & !is.na(fish_week) & fish_week %in% weeks
  n <- matrix(as.integer(table(
    factor(fish_week[keep], levels = weeks),
    factor(cell[keep], levels = 1:6)
  )), nrow = length(weeks),
  dimnames = list(paste0("wk", weeks), paste0("c", 1:6)))

  identify_p <- n[, 4L] + n[, 5L]
  identify_phi_s <- n[, 3L] + n[, 4L]
  parent <- integer(nrow(n))
  parent_id <- 1L
  p_total <- 0L
  phi_total <- 0L
  open_rows <- integer(0)
  for (i in seq_len(nrow(n))) {
    p_total <- p_total + identify_p[i]
    phi_total <- phi_total + identify_phi_s[i]
    open_rows <- c(open_rows, i)
    if (p_total > parent_floor && phi_total > parent_floor) {
      parent[open_rows] <- parent_id
      parent_id <- parent_id + 1L
      p_total <- 0L
      phi_total <- 0L
      open_rows <- integer(0)
    }
  }
  if (length(open_rows)) {
    parent[open_rows] <- if (parent_id > 1L) parent_id - 1L else 1L
  }
  parent <- as.integer(factor(parent))

  raw_spill_time <- trimws(as.character(spill_data$DateTime))
  spill_time <- as.POSIXct(rep(NA_character_, length(raw_spill_time)), tz = tz)
  for (fmt in c("%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M",
                "%Y-%m-%d %H:%M:%S")) {
    unresolved <- which(is.na(spill_time) & !is.na(raw_spill_time) &
                          nzchar(raw_spill_time))
    if (!length(unresolved)) break
    spill_time[unresolved] <- as.POSIXct(raw_spill_time[unresolved],
                                        format = fmt, tz = tz)
  }
  spill <- data.frame(
    Site = trimws(as.character(spill_data$Site)), datetime = spill_time,
    HourlySpill = suppressWarnings(as.numeric(spill_data$HourlySpill)),
    HourlyFlow = suppressWarnings(as.numeric(spill_data$HourlyFlow))
  )
  spill <- spill[spill$Site %in% c("LGR", "LGS") &
                   !is.na(spill$datetime) & is.finite(spill$HourlySpill) &
                   is.finite(spill$HourlyFlow), , drop = FALSE]
  spill <- spill[!duplicated(spill[, c("Site", "datetime")]), , drop = FALSE]
  spill$date <- as.Date(spill$datetime)
  daily <- aggregate(cbind(HourlySpill, HourlyFlow) ~ Site + date,
                     data = spill, FUN = sum)
  daily_flow <- aggregate(HourlyFlow ~ Site + date,
                          data = spill, FUN = mean)
  names(daily_flow)[3] <- "outflow"
  daily$outflow <- daily_flow$outflow[
    match(paste(daily$Site, daily$date),
          paste(daily_flow$Site, daily_flow$date))]
  daily$spill_pct <- 100 * daily$HourlySpill / daily$HourlyFlow
  daily$week <- as.integer(format(daily$date, "%V"))
  weekly <- aggregate(spill_pct ~ Site + week, data = daily, FUN = mean)
  weekly_outflow <- aggregate(outflow ~ Site + week, data = daily, FUN = mean)
  lgr_spill <- weekly$spill_pct[
    match(paste("LGR", weeks), paste(weekly$Site, weekly$week))]
  lgs_spill <- weekly$spill_pct[
    match(paste("LGS", weeks), paste(weekly$Site, weekly$week))]
  lgr_outflow <- weekly_outflow$outflow[
    match(paste("LGR", weeks),
          paste(weekly_outflow$Site, weekly_outflow$week))]
  if (anyNA(lgr_spill) || anyNA(lgs_spill) || anyNA(lgr_outflow)) {
    stop("LGR/LGS spill or LGR outflow is missing for at least one requested week.",
         call. = FALSE)
  }

  n_seen <- as.integer(rowSums(n))
  message("Six-cell data: ", sum(n_seen), " histories; ",
          length(excluded_tags), " adult/upstream tags excluded; ",
          sum(tags %in% transported_tags), " transported.")
  list(
    weeks = as.integer(weeks), S = length(weeks), n = n, N_seen = n_seen,
    S_lik = sum(n_seen > 0), lik_idx = which(n_seen > 0), parent = parent,
    n_strat = length(unique(parent)),
    lgr_spill_std = as.numeric(scale(lgr_spill)),
    lgs_spill_std = as.numeric(scale(lgs_spill)),
    alpha_phi_mean = alpha_phi_mean, alpha_phi_sd = alpha_phi_sd,
    beta_phi_mean = beta_phi_mean, beta_phi_sd = beta_phi_sd,
    lgr_spill_pct = lgr_spill, lgs_spill_pct = lgs_spill,
    lgr_outflow = lgr_outflow
  )
}
