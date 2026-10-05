#' Classify times as day or night at Lower Granite Dam
#'
#' Uses sunrise and sunset (or civil dawn and dusk) from suncalc.
#' Detection times are treated as local clock times, whatever time zone label
#' they carry. Pass PTAGIS times parsed as clock strings (for example UTC).
#' @param time POSIXct detection times (local clock).
#' @param lat,lon Site coordinates; defaults are GRJ from PTAGIS.
#' @param tz_local Local time zone for sun times, with daylight saving.
#' @param twilight Use sunrise/sunset or civil dawn/dusk to delimit day.
#' @return Character vector day or night, NA where time is NA.
#' @export
classify_daynight <- function(time, lat = 46.660394, lon = -117.436261,
                             tz_local = "America/Los_Angeles",
                             twilight = c("sunrise_sunset", "civil")) {
  twilight <- match.arg(twilight)
  if (!requireNamespace("suncalc", quietly = TRUE)) {
    stop("Install the suncalc package.", call. = FALSE)
  }
  out <- rep(NA_character_, length(time))
  ok <- !is.na(time)
  if (!any(ok)) return(out)
  keep <- if (twilight == "civil") c("dawn", "dusk") else c("sunrise", "sunset")
  lt <- as.POSIXlt(time[ok])
  day_chr <- format(time[ok], "%Y-%m-%d")
  dates <- as.Date(unique(day_chr))
  sun <- suncalc::getSunlightTimes(date = dates, lat = lat, lon = lon,
                                  keep = keep, tz = tz_local)
  clock_min <- function(x) {
    y <- as.POSIXlt(x)
    y$hour * 60 + y$min + y$sec / 60
  }
  i <- match(as.Date(day_chr), sun$date)
  start <- clock_min(sun[[keep[1]]])[i]
  end <- clock_min(sun[[keep[2]]])[i]
  m <- lt$hour * 60 + lt$min + lt$sec / 60
  out[ok] <- ifelse(m >= start & m < end, "day", "night")
  out
}

#' Six-cell data split by day and night of LGR passage
#'
#' Splits histories by first LGR detection (GRJ or GRS), using the same cell
#' and week rules as [prep_ge_data2()]. GOJ-only fish stay pooled.
#' @inheritParams prep_ge_data2
#' @param lat,lon,tz_local,twilight Passed to [classify_daynight()].
#' @return List with base, n_day and n_night (S x 6 matrices with column 5
#'   zero), n_c5 (GOJ-only counts), twilight, and tags (tag, cell, week,
#'   lgr_time, period).
#' @export
prep_ge_data_daynight <- function(
    dat_up, spill_data, weeks,
    lat = 46.660394, lon = -117.436261,
    tz_local = "America/Los_Angeles", twilight = c("sunrise_sunset", "civil"),
    adult_sites = c("LGRLDR", "PRDLD1", "SHERFT", "BONAFF"),
    route_order = c(GRJ = 1, GRS = 1, GOJ = 2, LMJ = 3, ICH = 4),
    transport_antennas = c("61", "62"), parent_floor = 10L,
    alpha_phi_mean = 0, alpha_phi_sd = 2,
    beta_phi_mean = 0, beta_phi_sd = 1, tz = "UTC") {
  twilight <- match.arg(twilight)
  base <- prep_ge_data2(
    dat_up, spill_data = spill_data, weeks = weeks, adult_sites = adult_sites,
    route_order = route_order, transport_antennas = transport_antennas,
    parent_floor = parent_floor, alpha_phi_mean = alpha_phi_mean,
    alpha_phi_sd = alpha_phi_sd, beta_phi_mean = beta_phi_mean,
    beta_phi_sd = beta_phi_sd, tz = tz)
  main <- dat_up[dat_up$source == "main", , drop = FALSE]
  route <- main[main$site %in% names(route_order), , drop = FALSE]
  route$route_number <- unname(route_order[route$site])
  route <- route[order(route$tag, route$det_time), , drop = FALSE]
  upstream_tags <- unique(route$tag[
    route$route_number < ave(route$route_number, route$tag, FUN = cummax)])
  adult_release_tags <- unique(main$tag[main$release_site %in% adult_sites])
  excluded_tags <- union(upstream_tags, adult_release_tags)
  main <- main[!(main$tag %in% excluded_tags), , drop = FALSE]
  transported_tags <- unique(dat_up$tag[
    !(dat_up$tag %in% excluded_tags) & dat_up$site == "GRJ" &
      dat_up$antenna %in% as.character(transport_antennas)])
  histories <- main[main$site %in% c("GRJ", "GRS", "GOJ"), , drop = FALSE]
  tags <- sort(unique(histories$tag))
  first_time <- function(site) {
    x <- tapply(histories$det_time[histories$site == site],
                histories$tag[histories$site == site], min)
    as.POSIXct(as.numeric(x)[match(tags, names(x))], origin = "1970-01-01", tz = tz)
  }
  GRJ <- first_time("GRJ")
  GRS <- first_time("GRS")
  GOJ <- first_time("GOJ")
  LGR <- pmin(GRJ, GRS, na.rm = TRUE)
  paired <- !is.na(LGR) & !is.na(GOJ)
  goj_lag_days <- median(as.numeric(difftime(GOJ[paired], LGR[paired], units = "days")))
  fish_week <- as.integer(format(as.Date(LGR), "%V"))
  goj_only <- is.na(LGR) & !is.na(GOJ)
  fish_week[goj_only] <- as.integer(format(
    as.Date(GOJ[goj_only] - goj_lag_days * 86400), "%V"))
  has_grj <- !is.na(GRJ)
  has_grs <- !is.na(GRS)
  has_goj <- !is.na(GOJ)
  has_grs[has_grj & has_grs] <- FALSE
  cell <- rep(NA_integer_, length(tags))
  cell[has_grj & !has_goj] <- 1L
  cell[has_grj & has_goj] <- 2L
  cell[!has_grj & has_grs & !has_goj] <- 3L
  cell[!has_grj & has_grs & has_goj] <- 4L
  cell[!has_grj & !has_grs & has_goj] <- 5L
  cell[tags %in% transported_tags] <- 6L
  keep <- !is.na(cell) & !is.na(fish_week) & fish_week %in% weeks
  period <- classify_daynight(LGR, lat = lat, lon = lon,
                              tz_local = tz_local, twilight = twilight)
  if (any(keep & cell != 5L & is.na(period))) {
    stop("Some fish with an LGR detection could not be classified day or night.",
         call. = FALSE)
  }
  count <- function(sel) {
    matrix(as.integer(table(factor(fish_week[sel], levels = weeks),
                            factor(cell[sel], levels = 1:6))), nrow = length(weeks),
           dimnames = list(paste0("wk", weeks), paste0("c", 1:6)))
  }
  n_day <- count(keep & cell != 5L & period == "day")
  n_night <- count(keep & cell != 5L & period == "night")
  n_c5 <- as.integer(table(factor(fish_week[keep & cell == 5L], levels = weeks)))
  rebuilt <- n_day + n_night
  rebuilt[, 5L] <- n_c5
  if (!identical(unname(rebuilt), unname(base$n))) {
    stop("Day/night counts do not add up to the prep_ge_data2() six-cell counts.",
         call. = FALSE)
  }
  list(base = base, n_day = n_day, n_night = n_night, n_c5 = n_c5,
       twilight = twilight,
       tags = data.frame(tag = tags[keep], cell = cell[keep], week = fish_week[keep],
                         lgr_time = LGR[keep], period = period[keep],
                         stringsAsFactors = FALSE))
}
