#' @title Prepare stratum-level data for Bayesian GE model
#'
#' @description Wrangle raw PIT tag detection data into stratum-level summaries
#'   for the JAGS guidance efficiency model. Two distinct fish pools are created:
#'   (1) a psi estimation pool of upstream-tagged fish classified by their first
#'   route through LGR (GRS vs. UND), and (2) direct LGR counts of all fish
#'   detected at GRJ or GRS during each stratum period. A spill covariate is
#'   attached by joining \code{spill_weekly} to strata by overlapping date ranges.
#'
#' @param dat_up data frame of PIT tag detections as returned by
#'   \code{\link{prep_pit_data}}. Required columns: \code{tag} (character),
#'   \code{site} (character site code), \code{det_date} (Date). Optional column
#'   \code{mark_rkm} (numeric) is used with \code{min_mark_rkm} to restrict
#'   Pool A to upstream-tagged fish without affecting Pool B.
#' @param min_mark_rkm numeric. Restrict both pools to fish tagged strictly
#'   above this RKM (e.g. \code{695} for LGR). Downstream-tagged fish are
#'   excluded from both Pool A and Pool B because their LGR detections may
#'   reflect adult upstream movement rather than smolt downstream migration.
#'   Default \code{695}.
#' @param strat_assign data frame mapping weeks to strata. Accepts two formats:
#'   \itemize{
#'     \item \strong{Week/Collapse format} (preferred): columns \code{Week}
#'       (integer ISO week number) and \code{Collapse} (integer stratum ID).
#'       Dates are derived from ISO week numbers of detections in \code{dat_up}.
#'     \item \strong{Date format}: columns \code{date} (Date), \code{stratum},
#'       and \code{stratum_idx} (integer, 1-based, monotonically increasing).
#'   }
#' @param spill_data data frame of daily spill values. Required columns:
#'   \code{Date} (Date or character) and \code{spill.per} (spill in raw units,
#'   e.g. kcfs). Rows are joined to strata by date and averaged within each stratum.
#' @param species one of \code{"chnk"} or \code{"sthd"}.
#' @param downstream_sites character vector of site codes at downstream detection
#'   sites used to identify fish that have passed through LGR.
#'
#' @return A data frame with one row per stratum and columns:
#'   \code{stratum}, \code{stratum_idx}, \code{n_GRS_pool}, \code{n_UND},
#'   \code{n_pool}, \code{n_GRJ_obs}, \code{n_GRS_obs}, \code{spill_val}.
#'
#' @importFrom dplyr filter arrange group_by slice ungroup transmute left_join
#'   mutate case_when full_join summarise distinct count rename across
#' @importFrom tidyr pivot_wider replace_na
#' @export
prep_ge_data <- function(dat_up,
                         strat_assign,
                         spill_data,
                         species,
                         min_mark_rkm     = 695,
                         downstream_sites = c("GOJ","LMJ","MCJ","JDJ",
                                              "B2J","BCC","TWX",
                                              "PD5","PD6","PD7","PD8","PDW",
                                              "ICH","PDO","ESANIS","TTOWER",
                                              "ASMEBR","PIER3","MLRSNI",
                                              "LMILIS","FOUNDI","CRESIS")) {

  species <- match.arg(species, c("chnk", "sthd"))

  # Restrict to upstream-tagged fish for both pools.
  # Upstream fish are the representative sample for both psi estimation (Pool A)
  # and LGR passage counts (Pool B) — downstream-tagged fish may include adults
  # moving upstream through the adult fishway, not smolts migrating downstream.
  if (!is.null(min_mark_rkm) && "mark_rkm" %in% names(dat_up)) {
    n_before <- length(unique(dat_up$tag))
    dat_up   <- dat_up[!is.na(dat_up$mark_rkm) & dat_up$mark_rkm > min_mark_rkm, ]
    message("RKM filter (mark_rkm > ", min_mark_rkm, "): ",
            n_before - length(unique(dat_up$tag)), " tags removed, ",
            length(unique(dat_up$tag)), " retained.")
  }

  # --- Normalise strat_assign to date/stratum/stratum_idx format ---
  # Accept Week/Collapse tibble (e.g. from the user's strata table) and expand
  # to one row per calendar date using ISO week numbers from dat_up + spill_data.
  if (all(c("Week", "Collapse") %in% names(strat_assign)) &&
      !("date" %in% names(strat_assign))) {

    all_dates <- sort(unique(c(
      as.Date(dat_up$det_date),
      as.Date(spill_data$Date)
    )))
    wk_num <- as.integer(format(all_dates, "%V"))   # ISO week

    week_to_strat <- strat_assign
    names(week_to_strat)[names(week_to_strat) == "Collapse"] <- "stratum"
    week_to_strat$stratum_idx <- as.integer(
      factor(week_to_strat$stratum, levels = sort(unique(week_to_strat$stratum)))
    )

    strat_assign <- merge(
      data.frame(date = all_dates, Week = wk_num),
      week_to_strat,
      by = "Week", all.x = FALSE
    )[, c("date", "stratum", "stratum_idx")]
  }

  # --- Pool A: psi estimation pool ---
  # Each upstream-tagged fish is classified by its first LGR route (GRS or UND).
  # GRJ fish are excluded — they never entered the spillway, so they cannot inform
  # psi (P(detected at GRS | passed through spillway)).
  down_first <- dat_up %>%
    filter(site %in% downstream_sites) %>%
    arrange(tag, det_date) %>%
    group_by(tag) %>% slice(1) %>% ungroup() %>%
    transmute(tag, down_date = det_date)

  lgr_first <- dat_up %>%
    filter(site %in% c("GRJ", "GRS")) %>%
    arrange(tag, det_date) %>%
    group_by(tag) %>% slice(1) %>% ungroup() %>%
    transmute(tag, lgr_route = site)

  pool_a <- down_first %>%
    left_join(lgr_first, by = "tag") %>%
    mutate(lgr_class = case_when(
      lgr_route == "GRS" ~ "GRS",
      is.na(lgr_route)   ~ "UND",
      TRUE               ~ "GRJ"
    )) %>%
    filter(lgr_class %in% c("GRS", "UND")) %>%
    left_join(strat_assign, by = c("down_date" = "date")) %>%
    filter(!is.na(stratum))

  psi_pool <- pool_a %>%
    group_by(stratum, stratum_idx) %>%
    summarise(n_GRS_pool = sum(lgr_class == "GRS"),
              n_UND      = sum(lgr_class == "UND"),
              n_pool     = n_GRS_pool + n_UND,
              .groups = "drop")

  # --- Pool B: direct LGR counts (one row per fish per stratum) ---
  # distinct(tag, site, stratum) collapses multiple detections of the same fish
  # at the same site within a stratum to a single count.
  lgr_counts <- dat_up %>%
    filter(site %in% c("GRJ", "GRS")) %>%
    distinct(tag, site, det_date) %>%
    left_join(strat_assign, by = c("det_date" = "date")) %>%
    filter(!is.na(stratum)) %>%
    distinct(tag, site, stratum) %>%
    count(stratum, site) %>%
    pivot_wider(names_from = site, values_from = n, values_fill = 0)
  if (!"GRJ" %in% names(lgr_counts)) lgr_counts$GRJ <- 0L
  if (!"GRS" %in% names(lgr_counts)) lgr_counts$GRS <- 0L
  lgr_counts <- rename(lgr_counts, n_GRJ_obs = GRJ, n_GRS_obs = GRS)

  # --- Spill covariate: mean daily spill.per within each stratum ---
  # spill_data must have a Date column (Date or character) and spill.per
  spill_df <- spill_data
  spill_df$Date <- as.Date(spill_df$Date)

  spill_strat <- strat_assign %>%
    left_join(spill_df[, c("Date", "spill.per")], by = c("date" = "Date")) %>%
    group_by(stratum) %>%
    summarise(spill_val = mean(spill.per, na.rm = TRUE), .groups = "drop") %>%
    mutate(spill_val = replace_na(spill_val, 0))

  # --- Combine all stratum summaries ---
  full_join(psi_pool, lgr_counts, by = "stratum") %>%
    left_join(spill_strat, by = "stratum") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    mutate(psi         = ifelse(n_pool > 0, n_GRS_pool / n_pool, NA_real_),
           n_total_GRS = n_GRS_obs / psi) %>%
    arrange(stratum_idx)
}
