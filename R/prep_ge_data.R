#' @title Prepare stratum-level data for Bayesian GE model
#'
#' @description Wrangle raw PIT tag detection data into stratum-level summaries
#'   for the JAGS guidance efficiency model. Two distinct fish pools are created:
#'   (1) a psi estimation pool of upstream-tagged fish classified by their first
#'   route through LGR (GRS vs. UND), and (2) direct counts of all fish
#'   detected at GRJ, GRS, or GOJ during each stratum period. A spill covariate
#'   is attached by joining \code{spill_weekly} to strata by overlapping date
#'   ranges.
#'
#'   Five detection history counts are also computed for the multistate likelihood
#'   in \code{\link{fit_ge_model}}:
#'   \enumerate{
#'     \item \code{n_h1}: GRJ detected, not at GOJ
#'     \item \code{n_h2}: GRJ detected and GOJ detected
#'     \item \code{n_h3}: GRS detected, not at GOJ
#'     \item \code{n_h4}: GRS detected and GOJ detected
#'     \item \code{n_h5}: not detected at LGR, detected at GOJ
#'   }
#'
#' @param dat_up data frame of PIT tag detections as returned by
#'   \code{\link{prep_pit_data}}. Required columns: \code{tag} (character),
#'   \code{site} (character site code), \code{det_date} (Date). Optional column
#'   \code{mark_rkm} (numeric) is used with \code{min_mark_rkm} to restrict
#'   Pool A to upstream-tagged fish without affecting Pool B.
#' @param min_mark_rkm numeric. Restrict both pools to fish tagged strictly
#'   above this RKM (e.g. \code{695} for LGR). Default \code{695}.
#' @param strat_assign data frame mapping weeks to strata. Accepts two formats:
#'   \itemize{
#'     \item \strong{Week/Collapse format} (preferred): columns \code{Week}
#'       (integer ISO week number) and \code{Collapse} (integer stratum ID).
#'     \item \strong{Date format}: columns \code{date} (Date), \code{stratum},
#'       and \code{stratum_idx} (integer, 1-based, monotonically increasing).
#'   }
#' @param parent_strata optional two-column data frame giving a coarser parent
#'   grouping for a nested (two-level) GE fit. The first column must match the
#'   stratum ids produced here (i.e. the \code{Collapse}/stratum value used in
#'   \code{strat_assign}; for a weekly run that is the week number), the second
#'   column is the parent group id (e.g. the collapse-stratum). When supplied, a
#'   \code{parent} column is added to the output and \code{\link{fit_ge_model}}
#'   shrinks each stratum's psi toward its parent group's mean. Pass your
#'   Week/Collapse table here when running weekly. Default \code{NULL}
#'   (single-level fit).
#' @param spill_data data frame of daily LGR spill values. Required columns:
#'   \code{Date} (Date or character) and \code{spill.per} (numeric, percentage).
#' @param lgs_spill_data data frame of daily LGS spill values. Required columns:
#'   \code{Date} (Date or character) and \code{spill.per} (numeric, percentage).
#'   Used to fit the phi (GOJ detection) regression. If \code{NULL}, LGS spill
#'   is assumed equal to LGR spill (not recommended under post-ROD operations).
#' @param species one of \code{"chnk"} or \code{"sthd"}.
#' @param goj_site character. Site code for the Little Goose bypass array used
#'   to compute detection history counts. Default \code{"GOJ"}.
#' @param downstream_sites character vector of site codes used to identify
#'   fish that have passed through LGR for Pool A classification.
#'
#' @return A data frame with one row per stratum and columns:
#'   \code{stratum}, \code{stratum_idx}, \code{n_GRS_pool}, \code{n_UND},
#'   \code{n_pool}, \code{n_GRJ_obs}, \code{n_GRS_obs}, \code{n_GOJ_obs},
#'   \code{spill_val}, \code{lgs_spill_val}, \code{n_h1}, \code{n_h2},
#'   \code{n_h3}, \code{n_h4}, \code{n_h5}.
#'
#' @importFrom dplyr filter arrange group_by slice ungroup transmute left_join
#'   mutate case_when full_join summarise distinct count rename across
#' @importFrom tidyr pivot_wider replace_na complete
#' @export
prep_ge_data <- function(dat_up,
                         strat_assign,
                         spill_data,
                         species,
                         lgs_spill_data   = NULL,
                         goj_site         = "GOJ",
                         min_mark_rkm     = 695,
                         downstream_sites = c("GOJ","LMJ","MCJ","JDJ",
                                              "B2J","BCC","TWX",
                                              "PD5","PD6","PD7","PD8","PDW",
                                              "ICH","PDO","ESANIS","TTOWER",
                                              "ASMEBR","PIER3","MLRSNI",
                                              "LMILIS","FOUNDI","CRESIS"),
                         parent_strata    = NULL) {

  species <- match.arg(species, c("chnk", "sthd"))

  # Restrict to upstream-tagged fish
  if (!is.null(min_mark_rkm) && "mark_rkm" %in% names(dat_up)) {
    n_before <- length(unique(dat_up$tag))
    dat_up   <- dat_up[!is.na(dat_up$mark_rkm) & dat_up$mark_rkm > min_mark_rkm, ]
    message("RKM filter (mark_rkm > ", min_mark_rkm, "): ",
            n_before - length(unique(dat_up$tag)), " tags removed, ",
            length(unique(dat_up$tag)), " retained.")
  }

  # Normalise strat_assign to date/stratum/stratum_idx format
  if (all(c("Week", "Collapse") %in% names(strat_assign)) &&
      !("date" %in% names(strat_assign))) {

    all_dates <- sort(unique(c(
      as.Date(dat_up$det_date),
      as.Date(spill_data$Date)
    )))
    wk_num <- as.integer(format(all_dates, "%V"))

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

  # Pool A: psi estimation pool
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
              .groups    = "drop")

  # Pool B: direct LGR and GOJ counts
  # n_GOJ_obs is a raw diagnostic count of GOJ detections per stratum. It is
  # not used inside the JAGS likelihood (which runs on n_h1-n_h5), but lets
  # you quickly check GOJ detection volume per stratum, useful for the
  # identifiability discussion around p and phi at high spill.
  lgr_counts <- dat_up %>%
    filter(site %in% c("GRJ", "GRS", goj_site)) %>%
    distinct(tag, site, det_date) %>%
    left_join(strat_assign, by = c("det_date" = "date")) %>%
    filter(!is.na(stratum)) %>%
    distinct(tag, site, stratum) %>%
    count(stratum, site) %>%
    pivot_wider(names_from = site, values_from = n, values_fill = 0)
  if (!"GRJ" %in% names(lgr_counts)) lgr_counts$GRJ <- 0L
  if (!"GRS" %in% names(lgr_counts)) lgr_counts$GRS <- 0L
  if (!goj_site %in% names(lgr_counts)) lgr_counts[[goj_site]] <- 0L
  lgr_counts <- rename(lgr_counts, n_GRJ_obs = GRJ, n_GRS_obs = GRS,
                       n_GOJ_obs = !!goj_site)

  # Five detection history counts for multistate likelihood
  # First detection per tag per site (GRJ, GRS, GOJ)
  lgr_sites <- c("GRJ", "GRS", goj_site)

  first_det <- dat_up %>%
    filter(site %in% lgr_sites) %>%
    arrange(tag, det_date) %>%
    group_by(tag, site) %>% slice(1) %>% ungroup() %>%
    select(tag, site, det_date)

  tag_wide <- first_det %>%
    pivot_wider(names_from = site, values_from = det_date,
                values_fn  = min, values_fill = NA)

  # Ensure all three columns exist
  for (col in c("GRJ", "GRS", goj_site)) {
    if (!col %in% names(tag_wide)) tag_wide[[col]] <- as.Date(NA)
  }
  names(tag_wide)[names(tag_wide) == goj_site] <- "GOJ"

  fish <- tag_wide %>%
    mutate(
      a = case_when(
        !is.na(GRJ)              ~ 1L,
        !is.na(GRS) & is.na(GRJ) ~ 2L,
        TRUE                      ~ 0L
      ),
      b        = if_else(!is.na(GOJ), 1L, 0L),
      lgr_date = case_when(
        !is.na(GRJ) ~ GRJ,
        !is.na(GRS) ~ GRS,
        !is.na(GOJ) ~ GOJ
      ),
      h = case_when(
        a == 1 & b == 0 ~ 1L,
        a == 1 & b == 1 ~ 2L,
        a == 2 & b == 0 ~ 3L,
        a == 2 & b == 1 ~ 4L,
        a == 0 & b == 1 ~ 5L
      )
    ) %>%
    filter(!is.na(lgr_date), !is.na(h)) %>%
    left_join(strat_assign, by = c("lgr_date" = "date")) %>%
    filter(!is.na(stratum))

  hist_counts <- fish %>%
    count(stratum, h) %>%
    complete(stratum = unique(strat_assign$stratum), h = 1:5,
             fill = list(n = 0L)) %>%
    pivot_wider(names_from = h, values_from = n, values_fill = 0L,
                names_prefix = "n_h")

  # LGR spill covariate
  spill_df       <- spill_data
  spill_df$Date  <- as.Date(spill_df$Date)

  spill_strat <- strat_assign %>%
    left_join(spill_df[, c("Date", "spill.per")], by = c("date" = "Date")) %>%
    group_by(stratum) %>%
    summarise(spill_val = mean(spill.per, na.rm = TRUE), .groups = "drop") %>%
    mutate(spill_val = replace_na(spill_val, 0))

  # LGS spill covariate for phi regression
  if (!is.null(lgs_spill_data)) {
    lgs_df       <- lgs_spill_data
    lgs_df$Date  <- as.Date(lgs_df$Date)
    lgs_strat <- strat_assign %>%
      left_join(lgs_df[, c("Date", "spill.per")], by = c("date" = "Date")) %>%
      group_by(stratum) %>%
      summarise(lgs_spill_val = mean(spill.per, na.rm = TRUE), .groups = "drop") %>%
      mutate(lgs_spill_val = replace_na(lgs_spill_val, 0))
  } else {
    message("lgs_spill_data not provided; using LGR spill for phi regression.",
            " Provide LGS spill data for more accurate GOJ detection modelling.")
    lgs_strat <- spill_strat %>%
      rename(lgs_spill_val = spill_val) %>%
      select(stratum, lgs_spill_val)
  }

  # Combine all stratum summaries
  ge_out <- full_join(psi_pool, lgr_counts, by = "stratum") %>%
    left_join(hist_counts,  by = "stratum") %>%
    left_join(spill_strat,  by = "stratum") %>%
    left_join(lgs_strat,    by = "stratum") %>%
    mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
    arrange(stratum_idx)

  # Attach parent grouping for a nested fit. The parent table's first column
  # matches the stratum id; the second is the parent group. Strata with no
  # match are placed in their own singleton parent (so they are never pooled
  # into an unrelated group).
  if (!is.null(parent_strata)) {
    pmap   <- as.data.frame(parent_strata)
    praw   <- pmap[[2]][match(ge_out$stratum, pmap[[1]])]
    if (anyNA(praw)) {
      warning(sum(is.na(praw)), " stratum/strata have no row in 'parent_strata';",
              " each is placed in its own parent group.")
      praw[is.na(praw)] <- paste0("_solo_", ge_out$stratum[is.na(praw)])
    }
    ge_out$parent <- as.integer(factor(praw, levels = unique(praw)))
  }

  ge_out
}
