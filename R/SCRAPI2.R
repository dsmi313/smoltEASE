# SCRAPI2 - SCRAPI with posterior GSI and guidance efficiency uncertainty
# Extends SCRAPI v2.2 (Steinhorst & Ackerman) by propagating two additional
# uncertainty sources into the bootstrap confidence intervals:
#   1. Per-fish posterior draws of genetic stock assignment (GSI)
#   2. Uncertainty in the juvenile bypass guidance efficiency (GE) estimate

#' @title SCRAPI2: SCRAPI with posterior GSI and guidance efficiency uncertainty
#'
#' @description Extends \code{\link{SCRAPI}} by propagating two additional sources
#'   of uncertainty into the bootstrap confidence intervals: (1) per-fish posterior
#'   draws of genetic stock identification (GSI) assignments, and (2) uncertainty
#'   in the juvenile bypass system guidance efficiency (GE) estimate. When neither
#'   \code{gsiDraws} nor \code{geDraws} is supplied, GE and stock assignment are
#'   fixed at their point values in every row, matching \code{SCRAPI}'s inputs --
#'   but the reported \code{Estimate} is still the bootstrap mean over uniform
#'   resampling (see Details), not \code{SCRAPI}'s row-1 raw-data value, so the
#'   two will agree closely but not bit-for-bit.
#'
#' @inheritParams SCRAPI
#' @param gsiDraws a data frame where the first column contains individual fish IDs
#'   (matching \code{fishID} in \code{smoltData}) and each remaining column contains
#'   the genetic stock assignment for one posterior draw. Column names are not
#'   otherwise used. When supplied, each bootstrap iteration draws a fresh set
#'   of per-fish stock assignments rather than treating the observed
#'   \code{Primary} column as fixed. If the number of draw columns is less than
#'   \code{B}, GSI
#'   columns are sampled with replacement (useful when the genetics lab
#'   provides fewer posterior draws than the bootstrap size).
#' @param fishID column name in \code{smoltData} that contains individual fish IDs
#'   for matching to \code{gsiDraws}. Required when \code{gsiDraws} is supplied.
#'   Defaults to \code{"MasterID"}.
#' @param n_point deprecated, no longer used. Every bootstrap row now draws its
#'   own GSI column, so there is no longer a separate averaging step for the
#'   point estimate to control. Kept only so existing calls do not break; a
#'   message is emitted if it is set to anything other than its default.
#' @param geDraws a data frame of pre-computed GE posterior draws returned by
#'   \code{\link{generate_ge_draws}}. Must have a \code{SampleEndDate} column
#'   matching dates in \code{passageData} plus one or more numeric draw columns.
#'   In each bootstrap iteration the corresponding draw
#'   column replaces the fixed \code{GuidanceEfficiency} values. If the number
#'   of draw columns is less than \code{B}, columns are sampled with replacement
#'   (same policy as \code{gsiDraws}). Set to \code{NULL} (default) to treat GE
#'   as fixed.
#' @param pointEst always \code{"mean"}, the mean of the B bootstrap rows. Kept
#'   as an argument only for call compatibility with earlier versions; any
#'   other value stops with an error rather than silently falling back.
#'   \code{"median"} is no longer offered: every row's composition proportions
#'   sum to 1 by construction, so the mean preserves additivity exactly (a
#'   stock's estimate plus every other stock's equals WildSmolts), and the
#'   auxiliary outputs (\code{Rear.csv}, \code{Prime.csv}, \code{PxS.csv}) are
#'   themselves bootstrap means -- offering a median main-table estimate would
#'   make those disagree with \code{answer} for no real benefit.
#'
#' @details
#' \strong{No row is the raw data evaluated once.} Every one of the B rows is
#' a complete bootstrap/Monte Carlo replicate with its own GE draw and GSI
#' draw when those posterior inputs are supplied, its own binomial resample
#' of passage counts, and its own uniform within-stratum resample of rearing
#' and composition fish. \code{Estimate} and \code{LCI}/\code{UCI} are therefore two
#' summaries (mean, and quantiles) of the same distribution, not two
#' different distributions -- an earlier version held row 1 fixed at the
#' observed data with no GE or resampling variation, which is no longer the
#' case.
#'
#' \strong{Guidance efficiency basis.} Abundance expands as
#' \code{count / (rate * GE)}. Because the abundance estimator is linear in
#' \code{1/GE}, the \code{"mean"} point estimate is exactly the Monte Carlo average of the
#' estimator over the GE posterior, the same construction used for the
#' \code{gsiDraws} averaging -- not a separate harmonic-mean formula computed
#' outside the bootstrap.
#'
#' \strong{Point estimates outside the interval can still happen.} \code{mean}
#' and \code{quantile} are different functionals; for a cell whose bootstrap
#' distribution is heavily right-skewed (very few informative fish, GE
#' estimated near zero on a high-volume day) the mean can still sit outside a
#' percentile interval computed from the same rows. This is no longer an
#' artifact of mixing a fixed and a resampled process. It indicates a strongly
#' right-skewed or heavy-tailed bootstrap distribution whose estimate is
#' sensitive to rare, large replicates and should be inspected. A message
#' lists any rows where this happens.
#'
#' \strong{GE is excluded from the genetic-composition sampling rate SR.}
#' Within each replicate, the selected GE draw enters both the aggregate
#' passage expansion and the rearing inclusion rate \code{True}.
#' Consequently, uncertainty in GE can affect total passage and the
#' estimated wild proportion, but it does not directly reweight genetic
#' stock composition. \code{SR} is \code{trap_rate * genotype_rate}, calculated
#' from the observed trapping and genotyping data and then held fixed across
#' bootstrap replicates. Earlier versions folded
#' GE into \code{SR} as well, which meant a fish genotyped on a near-zero-GE
#' day was inflated the same way abundance is -- importing the same
#' right-skew instability into stock proportions with no data-driven reason
#' to think stock mix correlates with day-to-day guidance efficiency. That
#' argument does not extend to rear type: whether wild/hatchery run timing
#' correlates with GE is a real empirical question, not one
#' inverse-probability weighting gets to assume away, so \code{True} keeps
#' GE and is rebuilt from that replicate's GE draw.
#'
#' \strong{Resampling is uniform regardless.} \code{SR} is fixed and
#' \code{True} varies by replicate, but weighting the resample by them while
#' also applying their inverse inside the estimator would tend to cancel the
#' intended inverse-probability correction and shift the bootstrap toward the
#' unweighted sample composition -- a mismatch inherited from the original
#' \code{SCRAPI}. Resampling fish
#' and rearing units with equal probability within each stratum, while
#' keeping \code{1/SR} and \code{1/True} in the estimator itself, avoids that
#' cancellation.
#'
#' \strong{Every reported output shares the same B replicates.} \code{ProWild},
#' the stock proportions and abundances written to \code{Prime.csv} and
#' \code{PxS.csv}, and the weekly passage totals in \code{Rear.csv} and the
#' printed summary are all accumulated across the same bootstrap rows that
#' build \code{Estimate}/\code{LCI}/\code{UCI}, not computed separately from
#' the observed data at the point GE. They will reconcile with the main table
#' rather than needing independent verification.
#'
#' @return A list (returned invisibly) with two elements:
#'   \item{CI}{matrix of point estimates and bootstrap confidence intervals,
#'     identical in structure to the table printed by \code{SCRAPI}}
#'   \item{bootstrap}{full \code{B x p} matrix of bootstrap draws}
#'   CSV output files are also written using the same naming scheme as
#'   \code{SCRAPI}.
#'
#' @author Original SCRAPI by Kirk Steinhorst and Mike Ackerman.
#'   SCRAPI2 extensions by Thomas Delomas.
#'
#' @examples
#' \dontrun{
#' SCRAPI2(smoltData   = sthdScrapiInput,
#'         passageData = sthdSmoltPassData,
#'         Primary     = "GenStock",
#'         Secondary   = "fwAge",
#'         gsiDraws    = myGsiDraws,   # MasterID + boot_1 ... boot_2000
#'         fishID      = "MasterID",
#'         Run         = "sthdSmolt2_",
#'         RTYPE       = "W",
#'         alph        = 0.1,
#'         B           = 2000)
#' }
#'
#' @param strata optional data frame supplying the week-to-stratum collapse
#'   mapping directly, so the \code{Collapse} column does not need to be
#'   pre-baked into \code{passageData}. Must have exactly two columns: the first
#'   matching the \code{strat} column (e.g. \code{Week}), the second giving the
#'   stratum integer (e.g. \code{Collapse}). When supplied this overrides any
#'   existing \code{Collapse} column in \code{passageData}. Default \code{NULL}.
#' @param seed optional integer. When supplied, \code{set.seed(seed)} is called
#'   before any random draws (GSI/GE column sampling, the bootstrap resampling
#'   loop), making a run reproducible. Default \code{NULL} (no seeding).
#'
#' @importFrom stats rbinom quantile plogis
#' @importFrom Hmisc mApply
#' @export

SCRAPI2 <- function(smoltData = NULL, Dat = "CollectionDate", Rr = "Rear",
                    Primary = "GenStock", Secondary = NA, passageData = NULL,
                    strat = "Week", dat = "SampleEndDate", tally = "SampleCount",
                    samrate = "SampleRate", guidance = "GuidanceEfficiency",
                    collaps = "Collapse", Run = "output", RTYPE = "W",
                    REARSTRAT = TRUE, alph = 0.1, B = 5000,
                    dateFormat = "%m/%d/%Y",
                    gsiDraws = NULL, fishID = "MasterID", n_point = 100L,
                    geDraws = NULL, strata = NULL, seed = NULL,
                    pointEst = "mean")
{
  if (length(B) != 1L || !is.finite(B) || B < 1 || B %% 1 != 0)
    stop("B must be a positive integer.", call. = FALSE)
  B <- as.integer(B)

  if (length(alph) != 1L || !is.finite(alph) || alph <= 0 || alph >= 1)
    stop("alph must be strictly between 0 and 1.", call. = FALSE)

  if (!identical(pointEst, "mean"))
    stop("SCRAPI2 uses the bootstrap mean to preserve additive abundance ",
         "totals; pointEst = '", pointEst, "' is no longer supported.",
         call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(gsiDraws) && n_point != 100L)
    message("n_point is no longer used: every bootstrap row now draws its own ",
            "GSI column, so there is nothing separate for n_point to control.")

  # ---- import data -------------------------------------------------------
  if(is.character(smoltData))   { All  <- read.csv(smoltData,  header = TRUE) } else { All  <- smoltData  }
  if(is.character(passageData)) { pass <- read.csv(passageData, header = TRUE) } else { pass <- passageData }

  # Normalise date columns to Date class up front. POSIXct input (common when
  # loaded from DB drivers) prints with timestamp suffixes under table() and
  # then fails to re-parse under as.Date(..., format=dateFormat).
  All[[Dat]]  <- as.Date(All[[Dat]],  tryFormats = c(dateFormat, "%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y"))
  pass[[dat]] <- as.Date(pass[[dat]], tryFormats = c(dateFormat, "%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y"))

  if (anyNA(All[[Dat]]))
    stop("smoltData contains missing or unparseable dates in '", Dat, "'.",
         call. = FALSE)
  if (anyNA(pass[[dat]]))
    stop("passageData contains missing or unparseable dates in '", dat, "'.",
         call. = FALSE)

  if(!is.null(geDraws)) {
    if(!"SampleEndDate" %in% names(geDraws))
      stop("geDraws must have a 'SampleEndDate' column", call. = FALSE)
    geDraws$SampleEndDate <- as.Date(
      geDraws$SampleEndDate,
      tryFormats = c(dateFormat, "%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y")
    )
    if(anyNA(geDraws$SampleEndDate))
      stop("geDraws contains missing or unparseable SampleEndDate values.",
           call. = FALSE)
    if(anyDuplicated(geDraws$SampleEndDate))
      stop("geDraws contains duplicate SampleEndDate values.", call. = FALSE)
  }

  # ---- apply strata mapping if supplied ----------------------------------
  if (!is.null(strata)) {
    if (ncol(strata) != 2)
      stop("'strata' must be a two-column data frame: Week (or equivalent) and Collapse")
    pass[[collaps]] <- strata[[2L]][match(pass[[strat]], strata[[1L]])]
    if (anyNA(pass[[collaps]]))
      warning("Some weeks in passageData have no matching row in 'strata'; ",
              "those rows will have NA stratum and be excluded.")
  }

  # ---- validate new parameters -------------------------------------------
  if(!is.null(gsiDraws)) {
    if(ncol(gsiDraws) < 2L)
      stop("gsiDraws must contain an ID column and at least one draw column.",
           call. = FALSE)
    if(!fishID %in% names(All))
      stop("fishID column '", fishID, "' not found in smoltData")
    if(anyNA(gsiDraws[[1]]) || anyDuplicated(gsiDraws[[1]]))
      stop("The first gsiDraws column must contain unique, non-missing fish IDs.",
           call. = FALSE)
    if(!all(All[[fishID]] %in% gsiDraws[[1]]))
      stop("Not all fishID values in smoltData are present in gsiDraws")
    n_gsi_available <- ncol(gsiDraws) - 1L
    if(n_gsi_available < 1L)
      stop("gsiDraws must have at least one draw column in addition to the ID column")
    if(n_gsi_available < B)
      message("gsiDraws has ", n_gsi_available, " draw column(s) but B = ", B,
              "; sampling GSI draws with replacement.")
    # Always an independent random draw, never seq_len(B): even when there
    # are exactly B or more columns, taking them in raw column order would
    # pair column b of GSI with column b of GE below by shared position,
    # which could accidentally correlate two unrelated MCMC chains that
    # happen to share iteration count or thinning.
    gsi_idx_boot <- sample.int(n_gsi_available, B, replace = n_gsi_available < B)
  }
  if(!is.null(geDraws)) {
    n_ge_available <- ncol(geDraws) - 1L
    if(n_ge_available < 1L)
      stop("geDraws must have at least one draw column in addition to SampleEndDate")
    if(!all(vapply(geDraws[, -1, drop = FALSE], is.numeric, logical(1))))
      stop("All geDraws draw columns must be numeric.", call. = FALSE)
    ge_chk <- as.matrix(geDraws[, -1, drop = FALSE])
    if(any(!is.finite(ge_chk)) || any(ge_chk <= 0) || any(ge_chk > 1))
      stop("geDraws values must all be finite and in (0, 1].", call. = FALSE)
    missing_ge_dates <- unique(pass[[dat]][
      !pass[[dat]] %in% geDraws$SampleEndDate
    ])
    if(length(missing_ge_dates) > 0L)
      stop("No GE posterior draws were supplied for passage date(s): ",
           paste(missing_ge_dates, collapse = ", "), call. = FALSE)
    if(n_ge_available < B)
      message("geDraws has ", n_ge_available, " draw column(s) but B = ", B,
              "; sampling GE draws with replacement.")
    ge_idx_boot <- sample.int(n_ge_available, B, replace = n_ge_available < B)
  }

  # ---- header ------------------------------------------------------------
  cat("\nStart time:", date(), "\n")
  cat("\nThis is a SCRAPI2 run of", Run, "\n")
  cat("\nFocus is on fish of type", RTYPE, "\n")
  cat("\nPrimary composition variable:", Primary,
      "| Secondary:", Secondary, "\n")
  cat("\nWild adjustment by week:", REARSTRAT, "\n")
  cat("\nBootstrap iterations: B =", B, "| Alpha =", alph, "\n")
  cat("\nPoint estimate:", pointEst, "\n")
  if(!is.null(gsiDraws))
    cat("\nGSI posterior uncertainty: ON  (each of B =", B, "bootstrap rows draws its own GSI column)\n")
  if(!is.null(geDraws))
    cat("\nGuidance efficiency uncertainty: ON  (geDraws supplied, B =", B, "draws)\n")

  # ---- inner: average secondary prop fallback ----------------------------
  getAvgProp <- function(Fh) {
    keep_secondary <- !is.na(Fh$SGrp) & Fh$SGrp != "NA" &
      nzchar(trimws(as.character(Fh$SGrp)))
    Fh <- droplevels(Fh[keep_secondary, , drop = FALSE])
    if(nrow(Fh) == 0L)
      return(matrix(0, nrow = nPgrps, ncol = nSgrps,
                    dimnames = list(as.character(Pgrps),
                                    as.character(Sgrps))))
    Freqs <- Hmisc::mApply(1/Fh$SR,
                           list(factor(Fh$SGrp, levels = Sgrps),
                                factor(Fh$PGrp, levels = Pgrps)),
                           sum)
    Freqs[is.na(Freqs)] <- 0
    prop.table(Freqs, margin = 1)
  }

  # ---- inner: core estimator ---------------------------------------------
  thetahat <- function(passage, RearDat, Fish) {
    strat_key <- as.character(strats)

    dailypass <- passage$Tally / passage$Ptrue
    bystrata  <- tapply(dailypass, passage$Stratum, sum)
    bystrata  <- bystrata[strat_key]

    HNCWstrat <- Hmisc::mApply(1/RearDat$True,
                               list(RearDat$Stratum, RearDat$Rear), sum)
    HNCWstrat[is.na(HNCWstrat)] <- 0

    if(ncol(as.data.frame(HNCWstrat)) == 1) {
      PWild <- 1
    } else {
      HNCWprop <- prop.table(HNCWstrat, margin = 2)
      if(!REARSTRAT) {
        ColTotals   <- apply(HNCWstrat, 1, sum)
        Proportions <- ColTotals / sum(ColTotals)
        PWild       <- Proportions[2]
      } else {
        PWild <- HNCWprop[2, ]
      }
    }

    if (length(PWild) == 1L) {
      PWild <- rep(PWild, nstrats)
    } else {
      PWild <- PWild[strat_key]
    }

    if (anyNA(bystrata) || anyNA(PWild))
      stop("Could not align passage and rearing estimates by stratum.",
           call. = FALSE)

    WildStrata <- setNames(as.numeric(PWild) * as.numeric(bystrata), strat_key)
    TotalWild  <- sum(WildStrata)

    # Pin grouping to the full level sets so a bootstrap resample missing a
    # primary group (or stratum) still yields a full [nPgrps x nstrats] matrix;
    # otherwise Primaryproportions/Primaryests shrink and the assembled theta
    # vector no longer matches p. Mirrors the Freqs/getAvgProp handling below.
    Primarystrata <- Hmisc::mApply(1/Fish$SR,
                                   list(factor(Fish$Strat, levels = strats),
                                        factor(Fish$PGrp,  levels = Pgrps)),
                                   sum)
    Primarystrata[is.na(Primarystrata)] <- 0
    Primaryproportions <- prop.table(Primarystrata, margin = 2)

    if (is.null(colnames(Primaryproportions)) ||
        !all(strat_key %in% colnames(Primaryproportions)))
      stop("Primary composition strata could not be aligned.", call. = FALSE)

    Primaryproportions <- Primaryproportions[, strat_key, drop = FALSE]
    Primaryests        <- Primaryproportions %*% WildStrata[strat_key]

    if(!is.na(Secondary)) {
      SecondAbund <- array(0, dim = c(nPgrps, nSgrps, nstrats))
      keep_secondary <- !is.na(Fish$SGrp) & Fish$SGrp != "NA" &
        nzchar(trimws(as.character(Fish$SGrp)))
      Fish <- droplevels(Fish[keep_secondary, , drop = FALSE])
      if(nrow(Fish) == 0L) {
        PrimeBySecond <- matrix(0, nrow = nPgrps, ncol = nSgrps)
        Second <- numeric(nSgrps)
        return(list(TotalWild, WildStrata, Primaryproportions,
                    Primaryests, Second, PrimeBySecond, PWild))
      }
      AvgProp_b <- getAvgProp(Fish)
      AvgProp_b[!is.finite(AvgProp_b)] <- 0
      Freqs <- Hmisc::mApply(1/Fish$SR,
                             list(factor(Fish$PGrp, levels = Pgrps),
                                  factor(Fish$SGrp, levels = Sgrps),
                                  factor(Fish$Strat, levels = strats)),
                             sum)
      Freqs[is.na(Freqs)] <- 0
      Props <- prop.table(Freqs, margin = c(1, 3))
      if(any(is.nan(Props))) {
        for(h in 1:nstrats)
          for(i in 1:nPgrps)
            if(any(is.nan(Props[i, , h]))) Props[i, , h] <- AvgProp_b[i, ]
      }
      for(h in 1:nstrats) {
        ThisPrime          <- WildStrata[h] * Primaryproportions[, h]
        SecondAbund[, , h] <- as.vector(t(c(ThisPrime * Props[, , h])))
      }
      PrimeBySecond <- apply(SecondAbund, c(1, 2), sum)
      PrimeBySecond[is.na(PrimeBySecond)] <- 0
      Second <- apply(PrimeBySecond, 2, sum)
      return(list(TotalWild, WildStrata, Primaryproportions,
                  Primaryests, Second, PrimeBySecond, PWild))
    } else {
      return(list(TotalWild, WildStrata, Primaryproportions,
                  Primaryests, PWild))
    }
  }

  # ---- column indices ----------------------------------------------------
  FISHdate    <- which(Dat      == names(All))
  FISHrear    <- which(Rr       == names(All))
  FISHpndx    <- which(Primary  == names(All))
  FISHidx     <- if(!is.null(gsiDraws)) which(fishID == names(All)) else NULL
  if(!is.na(Secondary)) FISHsndx <- which(Secondary == names(All))

  PASSstrat   <- which(strat    == names(pass))
  PASSdate    <- which(dat      == names(pass))
  PASSrate    <- which(samrate  == names(pass))
  PASScounts  <- which(tally    == names(pass))
  PASSguideff <- which(guidance == names(pass))
  PASScollaps <- which(collaps  == names(pass))

  # Required columns other than GuidanceEfficiency (handled below) must exist.
  # Missing Collapse is common when users forget to pass 'strata'; surface the
  # fix rather than letting the failure cascade into a cryptic dimnames error.
  req <- list(strat = PASSstrat, dat = PASSdate, samrate = PASSrate,
              tally = PASScounts, collaps = PASScollaps)
  missing_nm <- names(req)[vapply(req, length, integer(1)) == 0]
  if (length(missing_nm) > 0) {
    lookup <- c(strat = strat, dat = dat, samrate = samrate,
                tally = tally, collaps = collaps)
    hint <- ""
    if ("collaps" %in% missing_nm)
      hint <- paste0(" (pass strata = <Week/Collapse data frame> to have ",
                     "SCRAPI2 build the '", collaps, "' column for you)")
    stop("passageData is missing required column(s): ",
         paste(shQuote(lookup[missing_nm]), collapse = ", "), hint)
  }

  # Baseline GE derived from the draws. This is used only for input setup and
  # validation; the reported Estimate is calculated from the B complete
  # bootstrap/Monte Carlo replicates below.
  ge_point <- function(gd, dates) {
    m   <- as.matrix(gd[, -1, drop = FALSE])
    idx <- match(dates, gd$SampleEndDate)
    day <- 1 / rowMeans(1 / m)                          # per-day harmonic mean
    day[idx]
  }

  # When geDraws is supplied the bootstrap uses per-draw GE, but the initial
  # setup still needs a baseline GE column. If GuidanceEfficiency is absent,
  # derive that setup value from the draws.
  if (length(PASSguideff) == 0) {
    if (is.null(geDraws))
      stop("passageData is missing the '", guidance, "' column and no geDraws supplied.")
    pass[[guidance]] <- ge_point(geDraws, pass[[dat]])
    PASSguideff <- which(guidance == names(pass))
  }

  ndays <- nrow(pass)

  # ---- passage setup -----------------------------------------------------
  Cpattern <- unique(cbind(pass[, PASSstrat], pass[, PASScollaps]))
  cat("\nStrata collapsed according to:\n")
  temp <- t(Cpattern); rownames(temp) <- c("Week","Strata")
  colnames(temp) <- rep("", ncol(temp)); print(temp)

  # The initial setup needs a valid baseline inclusion rate even though the
  # reported Estimate is calculated entirely from the B replicates below.
  if (!is.null(geDraws)) {
    ge_pe     <- ge_point(geDraws, pass[[dat]])
    pass$true <- pass[, PASSrate] * ge_pe
  } else {
    pass$true <- pass[, PASSrate] * pass[, PASSguideff]
  }

  bad_pass <- is.na(pass[, PASScollaps]) |
              !is.finite(pass[, PASSrate]) | pass[, PASSrate] <= 0 |
              pass[, PASSrate] > 1 |
              !is.finite(pass$true) | pass$true <= 0 |
              pass$true > 1 |
              !is.finite(pass[, PASScounts]) | pass[, PASScounts] < 0
  if (any(bad_pass))
    stop("Invalid passage rows on: ",
         paste(pass[[dat]][bad_pass], collapse = ", "), call. = FALSE)

  passdata <- data.frame(Stratum = pass[, PASScollaps],
                         Tally   = pass[, PASScounts],
                         Ptrue   = pass$true)

  # ---- pre-compute geDraws daily matrix (n_days x B) ---------------------
  if(!is.null(geDraws)) {
    ge_day_idx <- match(pass[[dat]], geDraws$SampleEndDate)
    ge_mat_raw <- as.matrix(geDraws[, -1, drop = FALSE])   # n_gedays x n_ge_available
    ge_day_mat <- ge_mat_raw[ge_day_idx, ge_idx_boot, drop = FALSE]
  }

  # ---- assign stratum and true rate to each fish -------------------------
  # Pre-parse dates once; use which() so NA rows don't end up as subscripts.
  # Previous approach compared raw character dates which breaks when the two
  # sources use different string formats (e.g. "2025-04-08" vs "4/08/2025").
  pass_dates_d  <- as.Date(pass[, PASSdate], format = dateFormat)
  all_dates_d   <- as.Date(All[,  FISHdate], format = dateFormat)

  if (all(is.na(pass_dates_d)))
    stop("Could not parse any '", dat, "' values in passageData with dateFormat='",
         dateFormat, "'. Example value: ", shQuote(pass[1, PASSdate]))
  if (all(is.na(all_dates_d)))
    stop("Could not parse any '", Dat, "' values in smoltData with dateFormat='",
         dateFormat, "'. Example value: ", shQuote(All[1, FISHdate]))

  All$Collaps <- NA_integer_
  for(d in unique(all_dates_d)) {
    if (is.na(d)) next
    pidx  <- which(!is.na(pass_dates_d) & pass_dates_d == d)
    fidx  <- which(!is.na(all_dates_d)  & all_dates_d  == d)
    if (length(pidx) > 0 && length(fidx) > 0)
      All$Collaps[fidx] <- pass[pidx[1], PASScollaps]
  }

  # intersect() would strip the Date class; use %in% to preserve it.
  all_valid  <- all_dates_d[!is.na(all_dates_d) & all_dates_d %in% pass_dates_d]
  set        <- sort(unique(all_valid))
  ndates     <- length(set)
  # True = trap_rate * GE, same as it always was. Unlike SR, this is not
  # dropping GE by design: whether wild/hatchery run timing correlates with
  # GE is an empirical question, not one inverse-probability weighting gets
  # to assume away, so GE stays in and gets redrawn per bootstrap iteration
  # from that iteration's GE draw (true_day carries the day index needed to
  # look it up; the point-GE value here is used only where geDraws is NULL
  # and for the checkRates validation below).
  All$true     <- NA_real_
  All$true_day <- NA_integer_
  for(nn in seq_len(ndates)) {
    set_d    <- set[nn]
    prow_idx <- which(!is.na(pass_dates_d) & pass_dates_d == set_d)
    ptmp     <- pass[prow_idx, , drop = FALSE]
    fish_idx <- which(!is.na(all_dates_d) & all_dates_d == set_d)
    if (nrow(ptmp) > 0 && length(fish_idx) > 0) {
      All$true[fish_idx]     <- ptmp$true[1]
      All$true_day[fish_idx] <- prow_idx[1]
    }
  }
  rear_ok <- !is.na(All$Collaps) & !is.na(All$true) & !is.na(All$true_day)
  if (any(!rear_ok))
    warning(sum(!rear_ok),
            " fish excluded because they lack a matching passage date.",
            call. = FALSE)

  RearData <- data.frame(Rear     = All[rear_ok, FISHrear],
                         Stratum  = All$Collaps[rear_ok],
                         True     = All$true[rear_ok],
                         true_day = All$true_day[rear_ok])
  nAll <- nrow(RearData)

  # ---- filter to RTYPE and non-NA Primary --------------------------------
  AllRTYPE   <- droplevels(All[All[, FISHrear] == RTYPE, ])
  AllPrimary <- droplevels(AllRTYPE[AllRTYPE[, FISHpndx] != "NA", ])

  primary_ok <- !is.na(AllPrimary$Collaps) & !is.na(AllPrimary$true_day)
  AllPrimary <- droplevels(AllPrimary[primary_ok, , drop = FALSE])
  nFISH      <- nrow(AllPrimary)

  # Pre-parse dates so NA rows don't leak into subscripts/comparisons.
  allp_dates_d <- as.Date(AllPrimary[, FISHdate], format = dateFormat)
  all_valid    <- allp_dates_d[!is.na(allp_dates_d) & allp_dates_d %in% pass_dates_d]
  set          <- sort(unique(all_valid))
  ndates       <- length(set)

  # Per-date primary sampling rate. Use the Date-native count of primary fish on
  # each date (length(fidx)); a prior table()/as.Date(names(...)) re-parse broke
  # once date columns became Date class, leaving SR all-zero (see git 6b1e1d8).
  #
  # SR = trap_rate * subrate. Both are calculated from the observed data and
  # then held fixed across bootstrap replicates:
  # trap_rate is how much of the day the trap operated, subrate is
  # genotyped/trapped. GE is deliberately excluded. A genotyped fish's 1/SR
  # weight answers "how many fish does this one represent within what was
  # trapped" -- that's a question about trap and lab subsampling, not about
  # guidance efficiency. Folding GE in would inflate the same few genotyped
  # fish on a near-zero-GE day the same way it inflates abundance, importing
  # the 1/psi tail instability into composition with no data-driven reason to
  # think stock mix correlates with day-to-day GE. GE enters exactly once,
  # in the aggregate daily passage expansion (Ptrue_b) inside the bootstrap
  # loop, and nowhere else.
  AllPrimary$subrate  <- numeric(nFISH)
  AllPrimary$trap_rate <- numeric(nFISH)
  for(nn in seq_len(ndates)) {
    set_d <- set[nn]
    ptmp  <- pass[which(!is.na(pass_dates_d) & pass_dates_d == set_d), , drop = FALSE]
    fidx  <- which(!is.na(allp_dates_d) & allp_dates_d == set_d)
    if (nrow(ptmp) > 0 && length(fidx) > 0) {
      AllPrimary$subrate[fidx]   <- length(fidx) / ptmp[1, PASScounts]
      AllPrimary$trap_rate[fidx] <- ptmp[1, PASSrate]
    }
  }
  AllPrimary$SR <- AllPrimary$trap_rate * AllPrimary$subrate

  # ---- strata / group metadata -------------------------------------------
  observed_groups <- as.character(AllPrimary[, FISHpndx])
  observed_groups <- observed_groups[
    !is.na(observed_groups) & observed_groups != "NA" &
      nzchar(trimws(observed_groups))
  ]

  if (!is.null(gsiDraws)) {
    gsi_rows <- match(AllPrimary[[fishID]], gsiDraws[[1]])
    posterior_groups <- unique(unlist(
      lapply(gsiDraws[gsi_rows, -1, drop = FALSE], as.character),
      use.names = FALSE
    ))
    posterior_groups <- posterior_groups[
      !is.na(posterior_groups) & posterior_groups != "NA" &
        nzchar(trimws(posterior_groups))
    ]
    Pgrps <- sort(unique(c(observed_groups, posterior_groups)))
  } else {
    Pgrps <- sort(unique(observed_groups))
  }

  nPgrps  <- length(Pgrps)
  if(nPgrps == 0L)
    stop("No usable levels were found in Primary = '", Primary, "'.",
         call. = FALSE)
  p       <- 1 + nPgrps
  strats  <- unique(Cpattern[, 2])
  nstrats <- length(strats)

  # ---- validate SR and True before the bootstrap runs --------------------
  # A non-finite or non-positive rate is a data problem (a date with zero
  # trapped or genotyped fish, a formatting error in the passage file), not
  # something to paper over. Stop rather than invent an inclusion probability
  # via stratum-mean imputation, which could silently mask exactly this kind
  # of error. Must happen before AllPrime below snapshots AllPrimary$SR.
  checkRates <- function(x, strat_of_x, label) {
    for(h in strats) {
      idx <- which(strat_of_x == h)
      if (length(idx) == 0L)
        stop(label, " has no usable observations in stratum ", h, ".",
             call. = FALSE)
      if(any(!is.finite(x[idx]) | x[idx] <= 0 | x[idx] > 1))
        stop("Invalid ", label, " in stratum ", h,
             " (values must be finite and in (0, 1]).")
    }
  }
  checkRates(RearData$True, RearData$Stratum, "True")
  checkRates(AllPrimary$SR, AllPrimary$Collaps, "SR")

  if(!is.na(Secondary)) {
    AllPrime <- data.frame(Strat = AllPrimary$Collaps,
                           PGrp  = AllPrimary[, FISHpndx],
                           SGrp  = AllPrimary[, FISHsndx],
                           SR    = AllPrimary$SR)
    secondary_groups <- as.character(AllPrime$SGrp)
    Sgrps <- sort(unique(secondary_groups[
      !is.na(secondary_groups) & secondary_groups != "NA" &
        nzchar(trimws(secondary_groups))
    ]))
    nSgrps  <- length(Sgrps)
    if(nSgrps == 0L)
      stop("No usable levels were found in Secondary = '", Secondary, "'.",
           call. = FALSE)
    p       <- p + nSgrps + nPgrps * nSgrps
  } else {
    AllPrime <- data.frame(Strat = AllPrimary$Collaps,
                           PGrp  = AllPrimary[, FISHpndx],
                           SR    = AllPrimary$SR)
  }

  # ---- helper: build AllPrime from AllPrimary ----------------------------
  makeFishDat <- function(ap) {
    if(!is.na(Secondary))
      data.frame(Strat = ap$Collaps, PGrp = ap[, FISHpndx],
                 SGrp  = ap[, FISHsndx], SR = ap$SR)
    else
      data.frame(Strat = ap$Collaps, PGrp = ap[, FISHpndx], SR = ap$SR)
  }

  nsFish <- if(!is.na(Secondary)) {
    secondary_value <- as.character(AllPrime$SGrp)
    sum(!is.na(secondary_value) & secondary_value != "NA" &
          nzchar(trimws(secondary_value)))
  } else NA

  # ---- bootstrap ---------------------------------------------------------
  # Every one of the B rows is a complete bootstrap/Monte Carlo replicate: its
  # own GE draw and GSI draw when supplied, its own binomial passage resample,
  # and its own uniform resample of composition and rearing fish. No row is
  # the raw data evaluated once. Resampling is uniform even though True varies with GE_b
  # and SR is fixed: in both cases the 1/SR and 1/True inverse-probability
  # weighting happens inside thetahat, not duplicated into the resampling
  # probabilities as well (that duplication is what cancelled the estimator's
  # own weights in the inherited SCRAPI design).
  #
  # ProWild, WildCollaps, pPropTable, sAbunTable, and the weekly passage
  # totals are all accumulated here across the same B rows that build
  # theta.b, so Rear.csv, Prime.csv, PxS.csv, and the printed weekly totals
  # reconcile with the Estimate column in `answer` -- none of them come from
  # a separate single-pass calculation at the point GE.
  theta.b <- matrix(0, nrow = B, ncol = p)

  strat_key <- as.character(strats)
  WildStrata_sum <- setNames(numeric(nstrats), strat_key)
  # pAbun_sum accumulates each iteration's actual stock-by-stratum abundance
  # (proportion x that iteration's WildStrata), not proportion and total
  # separately: E[P] * E[N] != E[P*N] in general, and averaging them apart
  # would make Prime.csv disagree with the main table's per-stock Estimate,
  # which is the mean of est_b[[4]] (Primaryproportions %*% WildStrata) --
  # exactly this product, just already summed across strata there.
  pAbun_sum <- matrix(0, nrow = nstrats, ncol = nPgrps,
                      dimnames = list(strat_key, as.character(Pgrps)))
  if(!is.na(Secondary)) sAbun_sum <- matrix(0, nPgrps, nSgrps)
  bystrata_sum <- numeric(nstrats); names(bystrata_sum) <- strat_key
  byweek_levels <- sort(unique(pass[, PASSstrat]))
  byweek_sum    <- numeric(length(byweek_levels)); names(byweek_sum) <- as.character(byweek_levels)

  for(b in 1:B) {
    # GE draw for this iteration, decided before anything that depends on it.
    # This replicate's GE draw enters daily passage here and is reused
    # below to construct the replicate-specific rearing weights, True_b.
    # GE does not enter genetic-composition SR.
    if(!is.null(geDraws)) {
      ptrue_b <- pass[, PASSrate] * ge_day_mat[, b]
    } else {
      ptrue_b <- passdata$Ptrue
    }

    # binomial resample of daily counts, on the same GE basis as this draw
    est_daily_b <- round(passdata$Tally / ptrue_b)
    cntstar     <- numeric(ndays)
    for(i in 1:ndays)
      if(est_daily_b[i] > 0) cntstar[i] <- rbinom(1, est_daily_b[i], ptrue_b[i])

    dailyStar <- data.frame(Stratum = passdata$Stratum,
                            Tally   = cntstar,
                            Ptrue   = ptrue_b)

    bystrata_sum <- bystrata_sum + tapply(dailyStar$Tally / dailyStar$Ptrue,
                                          dailyStar$Stratum, sum)[as.character(strats)]
    byweek_sum   <- byweek_sum + tapply(dailyStar$Tally / dailyStar$Ptrue,
                                        pass[, PASSstrat], sum)[as.character(byweek_levels)]

    # uniform bootstrap of rearing data by stratum. True is rebuilt from this
    # iteration's ptrue_b (trap_rate * GE_b), restoring GE to the wild/
    # hatchery weight -- unlike SR, whether rear-type timing correlates with
    # GE is an open empirical question, not something to assume away.
    # Resampling itself stays uniform; only the 1/True weight inside
    # thetahat should carry GE, not the resampling probability.
    WHstar <- NULL
    for(h in strats) {
      jw <- RearData[which(RearData$Stratum == h), ]
      if(nrow(jw) == 0) next
      True_b     <- ptrue_b[jw$true_day]
      idx        <- sample.int(nrow(jw), replace = TRUE)
      jw_star    <- jw[idx, ]
      jw_star$True <- True_b[idx]
      WHstar <- rbind(WHstar, jw_star)
    }
    RearStar <- WHstar

    # uniform bootstrap of fish by stratum; SR is fixed (trap rate x
    # genotype rate), so there is nothing to weight the resample by either
    ap_boot <- NULL
    for(h in strats) {
      jw  <- AllPrimary[which(AllPrimary$Collaps == h), ]
      if(nrow(jw) == 0) next
      idx <- sample.int(nrow(jw), replace = TRUE)
      ap_boot <- rbind(ap_boot, jw[idx, ])
    }
    ap_b <- ap_boot

    # GSI draw for this iteration
    if(!is.null(gsiDraws)) {
      im <- match(ap_b[[fishID]], gsiDraws[[1]])
      gsi_b <- as.character(gsiDraws[[gsi_idx_boot[b] + 1L]][im])

      if(anyNA(gsi_b) || any(gsi_b == "NA") ||
         any(!nzchar(trimws(gsi_b))))
        stop("Missing GSI assignment in bootstrap iteration ", b,
             call. = FALSE)

      # Replacing the complete column avoids factor-level coercion.
      ap_b[[FISHpndx]] <- gsi_b
    }

    fd_b  <- makeFishDat(ap_b)
    est_b <- thetahat(dailyStar, RearStar, fd_b)

    wild_b <- est_b[[2]][strat_key]
    prop_b <- est_b[[3]][, strat_key, drop = FALSE]

    if (anyNA(wild_b) ||
        !identical(names(wild_b), strat_key) ||
        !identical(colnames(prop_b), strat_key))
      stop("Stratum alignment failed in bootstrap iteration ", b, call. = FALSE)

    WildStrata_sum <- WildStrata_sum + unname(wild_b)
    pAbun_sum      <- pAbun_sum + sweep(t(prop_b), 1, unname(wild_b), `*`)
    if(!is.na(Secondary)) sAbun_sum <- sAbun_sum + est_b[[6]]

    theta.b[b, ] <- if(!is.na(Secondary))
      c(est_b[[1]], t(est_b[[4]]), est_b[[5]], as.vector(t(est_b[[6]])))
    else
      c(est_b[[1]], t(est_b[[4]]))
  }

  WildCollaps <- WildStrata_sum / B
  pAbunRaw    <- pAbun_sum / B
  ProWild     <- WildCollaps / (bystrata_sum / B)
  pPropTable  <- sweep(pAbunRaw, 1, WildCollaps, `/`)
  pPropTable[!is.finite(pPropTable)] <- NA_real_
  if(!is.na(Secondary)) sAbunTable <- sAbun_sum / B
  rpasscollaps <- round(bystrata_sum / B)
  byweek_boot  <- round(byweek_sum / B)

  cat("\nTotal smolts by week:\n"); print(byweek_boot)
  cat("\nTotal smolts by statistical week:\n"); print(rpasscollaps)
  cat("\nTotal smolts:", round(sum(bystrata_sum / B)), "\n")

  # ---- confidence intervals ----------------------------------------------

  # Mean, not median: within every single bootstrap row the composition
  # proportions sum to 1 by construction (prop.table), so stock/age/sex cells
  # already sum exactly to their parent total in that row. Averaging preserves
  # an identity that holds row-wise; taking the median across rows does not.
  pointVec <- colMeans(theta.b)

  CI <- matrix(0, nrow = p, ncol = 3)
  for(j in 1:p) {
    cij      <- quantile(theta.b[, j], c(alph/2, 1 - alph/2))
    CI[j, ]  <- c(pointVec[j], cij)
  }
  CI <- round(CI)

  answer <- cbind(CI,
                  round((((CI[, 3] - CI[, 2]) / CI[, 1]) / 2) * 100, 1))
  colnames(answer) <- c("Estimate", "LCI", "UCI", "P1")

  if(!is.na(Secondary)) {
    grpnams <- ""
    for(prim in Pgrps)
      for(nam in Sgrps) grpnams <- c(grpnams, paste0(nam, prim))
    rownames(answer) <- c("WildSmolts", as.character(Pgrps), Sgrps, grpnams[-1])
  } else {
    rownames(answer) <- c("WildSmolts", as.character(Pgrps))
  }

  # Flag rather than hide any row where the point estimate leaves its interval.
  outside <- which(answer[, "Estimate"] < answer[, "LCI"] |
                   answer[, "Estimate"] > answer[, "UCI"])
  if(length(outside) > 0)
    message("Estimate falls outside its interval for ", length(outside), " row(s): ",
            paste(rownames(answer)[outside], collapse = ", "),
            ". Mean and quantile are different summaries of the same bootstrap ",
            "distribution; this flags a cell whose distribution is heavily ",
            "skewed, not a coding artifact. See ?SCRAPI2 Details.")

  cat("\n"); print(answer)

  # ---- write outputs (same layout as SCRAPI) -----------------------------
  tSmolts <- t(rbind(rpasscollaps, round(ProWild, 4), round(WildCollaps, 0)))
  colnames(tSmolts) <- c("TotalSmolts", "p(Wild)", "WildSmolts")
  write.csv(tSmolts, file = paste0(Run, "Rear.csv"))

  ciFile <- paste0(Run, "CIs.csv")
  header <- if(is.na(Secondary)) paste(RTYPE, "-", Primary) else
    paste(RTYPE, "-", Primary, "-", Secondary)
  write.table(header, file = ciFile, append = FALSE,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  suppressWarnings(write.table(answer, file = ciFile,
                               col.names = NA, sep = ",", append = TRUE))
  write.table(paste("RearSampleSize =",  nAll),  file = ciFile, append = TRUE,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(paste("PrimeSampleSize =", nFISH), file = ciFile, append = TRUE,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  if(!is.na(Secondary))
    write.table(paste("SecondSampleSize =", nsFish), file = ciFile,
                append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

  primeFile   <- paste0(Run, "Prime.csv")
  pAbunTable  <- round(pAbunRaw, 0)
  PrimeTotals <- round(colSums(pAbunRaw), 0)
  write.table(table(AllPrimary$Collaps, AllPrimary[, Primary]),
              file = primeFile, col.names = NA, sep = ",", append = FALSE)
  write.table(pPropTable,   file = primeFile, row.names = TRUE,
              col.names = FALSE, append = TRUE, sep = ",")
  write.table(pAbunTable,   file = primeFile, row.names = TRUE,
              col.names = FALSE, append = TRUE, sep = ",")
  write.table(t(PrimeTotals), file = primeFile, row.names = "PrimeTotals",
              col.names = FALSE, append = TRUE, sep = ",")

  if(!is.na(Secondary)) {
    sAbunTable <- cbind(sAbunTable, apply(sAbunTable, 1, sum))
    sAbunTable <- rbind(sAbunTable, apply(sAbunTable, 2, sum))
    rownames(sAbunTable) <- c(as.character(Pgrps), "sTotals")
    colnames(sAbunTable) <- c(Sgrps, "pTotals")
    write.table(round(sAbunTable, 0), file = paste0(Run, "PxS.csv"),
                col.names = NA, append = FALSE, sep = ",")
  }

  cat("\nEnd time:", date(), "\n")
  invisible(list(CI = answer, bootstrap = theta.b))
}
