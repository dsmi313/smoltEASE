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
#'   \code{gsiDraws} nor \code{geDraws} is supplied the function behaves identically
#'   to \code{SCRAPI}.
#'
#' @inheritParams SCRAPI
#' @param gsiDraws a data frame where the first column contains individual fish IDs
#'   (matching \code{fishID} in \code{smoltData}) and each remaining column contains
#'   the genetic stock assignment for one posterior draw, named \code{boot_1},
#'   \code{boot_2}, etc. When supplied, each bootstrap iteration draws a fresh set
#'   of per-fish stock assignments rather than treating the observed \code{Primary}
#'   column as fixed. If the number of draw columns is less than
#'   \code{max(B, n_point)}, GSI columns are sampled with replacement (useful when
#'   the genetics lab provides fewer posterior draws than the bootstrap size).
#' @param fishID column name in \code{smoltData} that contains individual fish IDs
#'   for matching to \code{gsiDraws}. Required when \code{gsiDraws} is supplied.
#'   Defaults to \code{"MasterID"}.
#' @param n_point number of GSI posterior draws used to compute the point estimate
#'   when \code{gsiDraws} is supplied. The point estimate is the mean across these
#'   draws. Default is 100.
#' @param geDraws a data frame of pre-computed GE posterior draws returned by
#'   \code{\link{generate_ge_draws}}. Must have a \code{SampleEndDate} column
#'   matching dates in \code{passageData} plus draw columns named \code{boot_1},
#'   \code{boot_2}, etc. In each bootstrap iteration the corresponding draw
#'   column replaces the fixed \code{GuidanceEfficiency} values. If the number
#'   of draw columns is less than \code{B}, columns are sampled with replacement
#'   (same policy as \code{gsiDraws}). Set to \code{NULL} (default) to treat GE
#'   as fixed.
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
#'         n_point     = 100,
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
                    gsiDraws = NULL, fishID = "MasterID", n_point = 100,
                    geDraws = NULL, strata = NULL)
{
  # ---- import data -------------------------------------------------------
  if(is.character(smoltData))   { All  <- read.csv(smoltData,  header = TRUE) } else { All  <- smoltData  }
  if(is.character(passageData)) { pass <- read.csv(passageData, header = TRUE) } else { pass <- passageData }

  # Normalise date columns to Date class up front. POSIXct input (common when
  # loaded from DB drivers) prints with timestamp suffixes under table() and
  # then fails to re-parse under as.Date(..., format=dateFormat).
  All[[Dat]]  <- as.Date(All[[Dat]])
  pass[[dat]] <- as.Date(pass[[dat]])
  if(!is.null(geDraws))
    geDraws$SampleEndDate <- as.Date(geDraws$SampleEndDate)

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
    if(!fishID %in% names(All))
      stop("fishID column '", fishID, "' not found in smoltData")
    if(!all(All[[fishID]] %in% gsiDraws[[1]]))
      stop("Not all fishID values in smoltData are present in gsiDraws")
    n_gsi_available <- ncol(gsiDraws) - 1L
    if(n_gsi_available < 1L)
      stop("gsiDraws must have at least one draw column in addition to the ID column")
    n_needed <- max(B, n_point)
    if(n_gsi_available < n_needed) {
      message("gsiDraws has ", n_gsi_available, " draw column(s) but max(B, n_point) = ",
              n_needed, "; sampling GSI draws with replacement.")
      gsi_idx_point <- sample.int(n_gsi_available, n_point, replace = TRUE)
      gsi_idx_boot  <- sample.int(n_gsi_available, B,       replace = TRUE)
    } else {
      gsi_idx_point <- seq_len(n_point)
      gsi_idx_boot  <- seq_len(B)
    }
  }
  if(!is.null(geDraws)) {
    if(!"SampleEndDate" %in% names(geDraws))
      stop("geDraws must have a 'SampleEndDate' column")
    n_ge_available <- ncol(geDraws) - 1L
    if(n_ge_available < 1L)
      stop("geDraws must have at least one draw column in addition to SampleEndDate")
    if(n_ge_available < B) {
      message("geDraws has ", n_ge_available, " draw column(s) but B = ", B,
              "; sampling GE draws with replacement.")
      ge_idx_boot <- sample.int(n_ge_available, B, replace = TRUE)
    } else {
      ge_idx_boot <- seq_len(B)
    }
  }

  # ---- header ------------------------------------------------------------
  cat("\nStart time:", date(), "\n")
  cat("\nThis is a SCRAPI2 run of", Run, "\n")
  cat("\nFocus is on fish of type", RTYPE, "\n")
  cat("\nPrimary composition variable:", Primary,
      "| Secondary:", Secondary, "\n")
  cat("\nWild adjustment by week:", REARSTRAT, "\n")
  cat("\nBootstrap iterations: B =", B, "| Alpha =", alph, "\n")
  if(!is.null(gsiDraws))
    cat("\nGSI posterior uncertainty: ON  (n_point =", n_point, "draws for point estimate)\n")
  if(!is.null(geDraws))
    cat("\nGuidance efficiency uncertainty: ON  (geDraws supplied, B =", B, "draws)\n")

  # ---- inner: average secondary prop fallback ----------------------------
  getAvgProp <- function(Fh) {
    Fh    <- droplevels(Fh[Fh$SGrp != "NA", ])
    Freqs <- Hmisc::mApply(1/Fh$SR,
                           list(factor(Fh$SGrp, levels = Sgrps),
                                factor(Fh$PGrp, levels = Pgrps)),
                           sum)
    Freqs[is.na(Freqs)] <- 0
    prop.table(Freqs, margin = 1)
  }

  # ---- inner: core estimator ---------------------------------------------
  thetahat <- function(passage, RearDat, Fish) {
    dailypass <- passage$Tally / passage$Ptrue
    bystrata  <- tapply(dailypass, passage$Stratum, sum)

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
    WildStrata <- PWild * bystrata
    TotalWild  <- sum(WildStrata)

    Primarystrata <- Hmisc::mApply(1/Fish$SR,
                                   list(Fish$Strat, Fish$PGrp), sum)
    Primarystrata[is.na(Primarystrata)] <- 0
    Primaryproportions <- prop.table(Primarystrata, margin = 2)
    Primaryests        <- Primaryproportions %*% WildStrata

    if(!is.na(Secondary)) {
      SecondAbund <- array(0, dim = c(nPgrps, nSgrps, nstrats))
      Fish <- droplevels(Fish[Fish$SGrp != "NA", ])
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
            if(any(is.nan(Props[i, , h]))) Props[i, , h] <- AvgProp[i, ]
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

  # When geDraws is supplied the bootstrap uses per-draw GE, but the initial
  # setup still needs a point-estimate GE column. If GuidanceEfficiency is
  # absent, derive it from the row means of geDraws.
  if (length(PASSguideff) == 0) {
    if (is.null(geDraws))
      stop("passageData is missing the '", guidance, "' column and no geDraws supplied.")
    ge_date_idx <- match(as.Date(pass[, PASSdate], format = dateFormat),
                         as.Date(geDraws$SampleEndDate, format = dateFormat))
    ge_means    <- rowMeans(as.matrix(geDraws[, -1, drop = FALSE]), na.rm = TRUE)
    pass[[guidance]] <- ifelse(is.na(ge_date_idx), mean(ge_means, na.rm = TRUE),
                               ge_means[ge_date_idx])
    PASSguideff <- which(guidance == names(pass))
  }

  ndays <- nrow(pass)

  # ---- passage setup -----------------------------------------------------
  Cpattern <- unique(cbind(pass[, PASSstrat], pass[, PASScollaps]))
  cat("\nStrata collapsed according to:\n")
  temp <- t(Cpattern); rownames(temp) <- c("Week","Strata")
  colnames(temp) <- rep("", ncol(temp)); print(temp)

  pass$true      <- pass[, PASSrate] * pass[, PASSguideff]
  pass$estimated <- pass[, PASScounts] / pass$true

  passdata <- data.frame(Stratum = pass[, PASScollaps],
                         Tally   = pass[, PASScounts],
                         Ptrue   = pass$true)

  # ---- pre-compute geDraws daily matrix (n_days x B) ---------------------
  if(!is.null(geDraws)) {
    ge_day_idx <- match(as.Date(pass[, dat], format = dateFormat),
                        as.Date(geDraws$SampleEndDate, format = dateFormat))
    ge_mat_raw <- as.matrix(geDraws[, -1, drop = FALSE])   # n_gedays x n_ge_available
    ge_season  <- colMeans(ge_mat_raw, na.rm = TRUE)[ge_idx_boot]  # season fallback per boot
    ge_day_mat <- matrix(ge_season, nrow = ndays, ncol = B, byrow = TRUE)
    valid_ge   <- !is.na(ge_day_idx)
    ge_day_mat[valid_ge, ] <- ge_mat_raw[ge_day_idx[valid_ge], ge_idx_boot, drop = FALSE]
  }

  passcollaps  <- tapply(pass$estimated, pass[, PASScollaps], sum)
  rpasscollaps <- round(passcollaps)
  cat("\nTotal smolts by week:\n"); print(round(tapply(pass$estimated, pass[, PASSstrat], sum)))
  cat("\nTotal smolts by statistical week:\n"); print(rpasscollaps)
  cat("\nTotal smolts:", round(sum(passcollaps)), "\n")

  # ---- assign stratum and true rate to each fish -------------------------
  # Pre-parse dates once; use which() so NA rows don't end up as subscripts.
  # Previous approach compared raw character dates which breaks when the two
  # sources use different string formats (e.g. "2025-04-08" vs "4/08/2025").
  nAll          <- nrow(All)
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
  All$true <- NA_real_
  for(nn in seq_len(ndates)) {
    set_d    <- set[nn]
    ptmp     <- pass[which(!is.na(pass_dates_d) & pass_dates_d == set_d), , drop = FALSE]
    fish_idx <- which(!is.na(all_dates_d) & all_dates_d == set_d)
    if (nrow(ptmp) > 0 && length(fish_idx) > 0)
      All$true[fish_idx] <- ptmp$true[1]
  }
  n_unmatched <- sum(is.na(All$Collaps) | is.na(All$true))
  if (n_unmatched > 0)
    message(n_unmatched, " fish in smoltData have no matching date in passageData and will be excluded.")

  RearData <- data.frame(Rear    = All[, FISHrear],
                         Stratum = All$Collaps,
                         True    = All$true)

  # ---- filter to RTYPE and non-NA Primary --------------------------------
  AllRTYPE   <- droplevels(All[All[, FISHrear] == RTYPE, ])
  AllPrimary <- droplevels(AllRTYPE[AllRTYPE[, FISHpndx] != "NA", ])
  nFISH      <- nrow(AllPrimary)

  # Pre-parse dates so NA rows don't leak into subscripts/comparisons.
  allp_dates_d <- as.Date(AllPrimary[, FISHdate], format = dateFormat)
  all_valid    <- allp_dates_d[!is.na(allp_dates_d) & allp_dates_d %in% pass_dates_d]
  set          <- sort(unique(all_valid))
  ndates       <- length(set)

  tabl      <- table(AllPrimary[, FISHdate], AllPrimary[, FISHpndx])
  nPrime    <- apply(tabl, 1, sum)
  nPrime_d  <- as.Date(names(nPrime), format = dateFormat)

  AllPrimary$SR <- numeric(nFISH)
  for(nn in seq_len(ndates)) {
    set_d <- set[nn]
    ptmp  <- pass[which(!is.na(pass_dates_d) & pass_dates_d == set_d), , drop = FALSE]
    pc    <- nPrime[which(!is.na(nPrime_d) & nPrime_d == set_d)]
    fidx  <- which(!is.na(allp_dates_d) & allp_dates_d == set_d)
    if (nrow(ptmp) > 0 && length(fidx) > 0 && length(pc) > 0)
      AllPrimary$SR[fidx] <- ptmp$true[1] * pc[1] / ptmp[1, PASScounts]
  }

  # ---- strata / group metadata -------------------------------------------
  Pgrps   <- sort(unique(AllPrimary[, FISHpndx]))
  nPgrps  <- length(Pgrps)
  p       <- 1 + nPgrps
  strats  <- unique(Cpattern[, 2])
  nstrats <- length(strats)

  if(!is.na(Secondary)) {
    AllPrime <- data.frame(Strat = AllPrimary$Collaps,
                           PGrp  = AllPrimary[, FISHpndx],
                           SGrp  = AllPrimary[, FISHsndx],
                           SR    = AllPrimary$SR)
    Sgrps   <- sort(unique(as.character(AllPrime$SGrp[AllPrime$SGrp != "NA"])))
    nSgrps  <- length(Sgrps)
    p       <- p + nSgrps + nPgrps * nSgrps
    AvgProp <- getAvgProp(AllPrime)
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

  # ---- point estimate ----------------------------------------------------
  if(!is.null(gsiDraws)) {
    ind_matches <- match(AllPrimary[[fishID]], gsiDraws[[1]])
    accumPE <- NULL
    for(i in 1:n_point) {
      ap_i <- AllPrimary
      ap_i[, FISHpndx] <- as.character(gsiDraws[[gsi_idx_point[i] + 1L]][ind_matches])
      fd_i  <- makeFishDat(ap_i)
      est_i <- thetahat(passdata, RearData, fd_i)
      # accumulate theta vector
      ncomp <- if(!is.na(Secondary)) 4L else 3L  # components before PWild
      tvec  <- if(!is.na(Secondary))
        c(est_i[[1]], t(est_i[[4]]), est_i[[5]], as.vector(t(est_i[[6]])))
      else
        c(est_i[[1]], t(est_i[[4]]))
      accumPE <- if(is.null(accumPE)) tvec else accumPE + tvec
    }
    pe_vec   <- accumPE / n_point
    # use last thetahat result for ProWild / WildCollaps / pPropTable
    ests_pe  <- thetahat(passdata, RearData, makeFishDat(AllPrimary))
  } else {
    ests_pe <- thetahat(passdata, RearData, AllPrime)
  }

  pwild_idx <- if(!is.na(Secondary)) 7L else 5L
  ProWild    <- ests_pe[[pwild_idx]]
  WildCollaps <- ests_pe[[2]]
  pPropTable  <- t(ests_pe[[3]])
  if(!is.na(Secondary)) sAbunTable <- ests_pe[[6]]
  nsFish <- if(!is.na(Secondary)) sum(AllPrime$SGrp != "NA") else NA

  # ---- bootstrap ---------------------------------------------------------
  theta.b <- matrix(0, nrow = B, ncol = p)

  for(b in 1:B) {
    # --- passage counts ---
    if(b == 1) {
      dailyStar <- passdata
      RearStar  <- RearData
      ap_b      <- AllPrimary
    } else {
      # binomial resample of daily counts
      est_daily <- round(passdata$Tally / passdata$Ptrue)
      cntstar   <- numeric(ndays)
      for(i in 1:ndays)
        if(est_daily[i] > 0) cntstar[i] <- rbinom(1, est_daily[i], passdata$Ptrue[i])

      # GE uncertainty: use pre-computed geDraws column for this iteration
      if(!is.null(geDraws)) {
        ptrue_b <- pass[, PASSrate] * ge_day_mat[, b]
      } else {
        ptrue_b <- passdata$Ptrue
      }
      dailyStar <- data.frame(Stratum = passdata$Stratum,
                              Tally   = cntstar,
                              Ptrue   = ptrue_b)

      # weighted bootstrap of rearing data by stratum
      WHstar <- NULL
      for(h in strats) {
        jw  <- RearData[which(RearData$Stratum == h), ]
        idx <- sample.int(nrow(jw), replace = TRUE, prob = jw$True)
        WHstar <- rbind(WHstar, jw[idx, ])
      }
      RearStar <- WHstar

      # weighted bootstrap of fish by stratum
      ap_boot <- NULL
      for(h in strats) {
        jw  <- AllPrimary[which(AllPrimary$Collaps == h), ]
        idx <- sample.int(nrow(jw), replace = TRUE, prob = jw$SR)
        ap_boot <- rbind(ap_boot, jw[idx, ])
      }
      ap_b <- ap_boot
    }

    # --- GSI draw for this iteration ---
    if(!is.null(gsiDraws)) {
      im    <- match(ap_b[[fishID]], gsiDraws[[1]])
      ap_b[, FISHpndx] <- as.character(gsiDraws[[gsi_idx_boot[b] + 1L]][im])
    }

    fd_b    <- makeFishDat(ap_b)
    est_b   <- thetahat(dailyStar, RearStar, fd_b)

    theta.b[b, ] <- if(!is.na(Secondary))
      c(est_b[[1]], t(est_b[[4]]), est_b[[5]], as.vector(t(est_b[[6]])))
    else
      c(est_b[[1]], t(est_b[[4]]))
  }

  # ---- if gsiDraws, replace row 1 with averaged point estimate ----------
  if(!is.null(gsiDraws)) theta.b[1, ] <- pe_vec

  # ---- confidence intervals ----------------------------------------------
  CI <- matrix(0, nrow = p, ncol = 3)
  for(j in 1:p) {
    cij      <- quantile(theta.b[, j], c(alph/2, 1 - alph/2))
    CI[j, ]  <- c(theta.b[1, j], cij)
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
  pAbunTable  <- round(sweep(pPropTable, 1, WildCollaps, `*`), 0)
  PrimeTotals <- round(apply(pAbunTable, 2, sum), 0)
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
    rownames(sAbunTable) <- c(levels(Pgrps), "sTotals")
    colnames(sAbunTable) <- c(Sgrps, "pTotals")
    write.table(round(sAbunTable, 0), file = paste0(Run, "PxS.csv"),
                col.names = NA, append = FALSE, sep = ",")
  }

  cat("\nEnd time:", date(), "\n")
  invisible(list(CI = answer, bootstrap = theta.b))
}
