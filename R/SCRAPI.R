# SCRAPI v2.2 - Salmonid Composition and Run Analyses for Pacific salmon Indices
# Original authors: Kirk Steinhorst and Mike Ackerman (mackerman44/SCOBI)
# Ported into smoltEASE for use alongside EASE uncertainty framework.
# mApply (plyr) replaced with base R tapply; behavior is identical.

#' @title SCRAPI v2.2: The smolt companion to SCOBI
#'
#' @description Perform compositional analyses of smolts at Lower Granite Dam.
#'
#' @param smoltData the data.frame (or path to .csv) containing biological data for
#'   smolts to be analyzed. Should contain all smolts trapped in a given migratory year
#'   and species (sthd, ch0, or ch1). Use \code{lgr2SCRAPI()} to format raw data from
#'   LGTrappingDB.
#' @param Dat column in \code{smoltData} containing sample date for each smolt
#' @param Rr column in \code{smoltData} containing rear type (W or HNC)
#' @param Primary primary category in \code{smoltData} to estimate (e.g. GenStock)
#' @param Secondary secondary category to estimate nested within Primary; use
#'   \code{NA} if no secondary decomposition is desired
#' @param passageData data.frame (or path to .csv) containing smolt passage data with
#'   calendar week, sampling rate, daily smolt count, guidance efficiency, and
#'   collapsing scheme
#' @param strat column in \code{passageData} containing calendar weeks
#' @param dat column in \code{passageData} containing sample dates
#' @param tally column in \code{passageData} containing daily smolt count in the trap
#' @param samrate column in \code{passageData} containing daily trap sample rate
#' @param guidance column in \code{passageData} containing daily guidance efficiency
#' @param collaps column in \code{passageData} containing the collapsing scheme
#' @param Run label used as prefix for all output files
#' @param RTYPE rear type to analyse: "W" (wild) or "HNC" (hatchery no-clip)
#' @param REARSTRAT if \code{TRUE}, proportion wild is calculated on a
#'   time-stratified basis; if \code{FALSE} it is pooled across the emigration
#' @param alph alpha for confidence intervals (e.g. 0.10 gives 90\% CIs)
#' @param B number of bootstrap iterations
#' @param dateFormat date format string used in fish and passage data
#'   (default \code{"\%m/\%d/\%Y"})
#'
#' @return NULL (results written to CSV files and printed to console)
#'
#' @author Kirk Steinhorst and Mike Ackerman
#'
#' @examples
#' \dontrun{
#' SCRAPI(smoltData = sthdScrapiInput, Primary = "GenStock", Secondary = "fwAge",
#'        passageData = sthdSmoltPassData, Run = "sthdSmoltTest",
#'        RTYPE = "W", alph = 0.1, B = 200)
#' }
#'
#' @importFrom stats rbinom quantile
#' @export

SCRAPI <- function(smoltData = NULL, Dat = "CollectionDate", Rr = "Rear",
                   Primary = "GenStock", Secondary = NA, passageData = NULL,
                   strat = "Week", dat = "SampleEndDate", tally = "SampleCount",
                   samrate = "SampleRate", guidance = "GuidanceEfficiency",
                   collaps = "Collapse", Run = "output", RTYPE = "W",
                   REARSTRAT = TRUE, alph = 0.1, B = 5000,
                   dateFormat = "%m/%d/%Y")
{
  # Import data
  if(is.character(smoltData))   { All  <- read.csv(file = smoltData,  header = TRUE) } else { All  <- smoltData  }
  if(is.character(passageData)) { pass <- read.csv(file = passageData, header = TRUE) } else { pass <- passageData }

  cat("\nStart time: ", date(), "\n")
  cat("\nThis is a run of", Run, "\n")
  cat("\nFocus is on fish of type", RTYPE, "\n")
  cat("\nPrimary composition variable is", Primary, "and the secondary variable is", Secondary, "\n")
  cat("\nWild adjustment by week is", REARSTRAT, "\n")
  cat("\nNumber of bootstrap iterations: B =", B, "\n")
  cat("\nAlpha =", alph, "\n")

  # Average secondary proportion across all strata (fallback when a stratum has no fish)
  getAvgProp <- function(Fh) {
    Fh    <- droplevels(Fh[which(Fh$SGrp != "NA"), ])
    Freqs <- tapply(1/Fh$SR,
                    list(factor(Fh$SGrp, levels = Sgrps),
                         factor(Fh$PGrp, levels = Pgrps)),
                    sum)
    Freqs[is.na(Freqs)] <- 0
    prop.table(Freqs, margin = 1)
  }

  # Core estimator: total wild smolt passage and composition by category
  thetahat <- function(passage, RearDat, Fish) {
    dailypass <- passage$Tally / passage$Ptrue
    bystrata  <- tapply(dailypass, passage$Stratum, sum)

    HNCWstrat <- tapply(1/RearDat$True,
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

    if(assignflg) {
      assign("ProWild",    PWild,      envir = .GlobalEnv)
      assign("WildCollaps", WildStrata, envir = .GlobalEnv)
    }

    TotalWild          <- sum(WildStrata)
    Primarystrata      <- tapply(1/Fish$SR,
                                 list(Fish$Strat, Fish$PGrp), sum)
    Primarystrata[is.na(Primarystrata)] <- 0
    Primaryproportions <- prop.table(Primarystrata, margin = 2)
    Primaryests        <- Primaryproportions %*% WildStrata

    if(!is.na(Secondary)) {
      SecondAbund <- array(numeric(nPgrps * nSgrps * nstrats),
                           dim = c(nPgrps, nSgrps, nstrats))
      Fish <- droplevels(Fish[which(Fish$SGrp != "NA"), ])
      if(assignflg) assign("nsFish", nrow(Fish), envir = .GlobalEnv)
      Freqs <- tapply(1/Fish$SR,
                      list(factor(Fish$PGrp, levels = Pgrps),
                           factor(Fish$SGrp, levels = Sgrps),
                           factor(Fish$Strat, levels = strats)),
                      sum)
      Freqs[is.na(Freqs)] <- 0
      Props <- prop.table(Freqs, margin = c(1, 3))
      if(any(is.nan(Props))) {
        for(h in 1:nstrats) {
          for(i in 1:nPgrps) {
            if(any(is.nan(Props[i, , h]))) Props[i, , h] <- AvgProp[i, ]
          }
        }
      }
      for(h in 1:nstrats) {
        ThisPrime          <- WildStrata[h] * Primaryproportions[, h]
        SecondAbund[, , h] <- as.vector(t(c(ThisPrime * Props[, , h])))
      }
      PrimeBySecond <- apply(SecondAbund, c(1, 2), sum)
      PrimeBySecond[is.na(PrimeBySecond)] <- 0
      if(assignflg) assign("sTable", PrimeBySecond, envir = .GlobalEnv)
      Second <- apply(PrimeBySecond, 2, sum)
      return(list(TotalWild, WildStrata, Primaryproportions, Primaryests,
                  Second, PrimeBySecond))
    } else {
      return(list(TotalWild, WildStrata, Primaryproportions, Primaryests))
    }
  }

  # Bootstrap function
  bootsmolt <- function(FishWH, FishDat, LGDdaily) {
    theta.b <- matrix(numeric(p * B), ncol = p)
    for(b in 1:B) {
      if(b == 1) {
        dailyStar <- passdata
        RearStar  <- RearData
        indivStar <- FishDat
      } else {
        dailypass <- round(LGDdaily$Tally / LGDdaily$Ptrue)
        cntstar   <- numeric(ndays)
        for(i in 1:ndays) {
          if(dailypass[i] != 0)
            cntstar[i] <- rbinom(1, dailypass[i], LGDdaily$Ptrue[i])
        }
        dailyStar <- data.frame(Stratum = LGDdaily$Stratum,
                                Tally   = cntstar,
                                Ptrue   = LGDdaily$Ptrue)

        # Weighted bootstrap of rearing data by stratum
        WHstar <- NULL
        for(h in strats) {
          justwk <- FishWH[FishWH$Stratum == h, ]
          nwk    <- nrow(justwk)
          idx    <- sample.int(nwk, replace = TRUE, prob = unlist(justwk$True))
          wkstar <- justwk[idx, ]
          WHstar <- if(is.null(WHstar)) wkstar else rbind(WHstar, wkstar)
        }
        RearStar <- WHstar

        # Weighted bootstrap of composition fish by stratum
        indivStar <- NULL
        for(h in strats) {
          justwk <- FishDat[FishDat$Strat == h, ]
          nwk    <- nrow(justwk)
          idx    <- sample.int(nwk, replace = TRUE, prob = unlist(justwk$SR))
          wkstar <- justwk[idx, ]
          indivStar <- if(is.null(indivStar)) wkstar else rbind(indivStar, wkstar)
        }
      }

      eststar <- thetahat(dailyStar, RearStar, indivStar)
      if(!is.na(Secondary))
        theta.b[b, ] <- c(eststar[[1]], t(eststar[[4]]), eststar[[5]],
                          as.vector(t(eststar[[6]])))
      else
        theta.b[b, ] <- c(eststar[[1]], t(eststar[[4]]))
    }

    CI <- matrix(numeric(ncol(theta.b) * 3), ncol = 3)
    for(j in 1:p) {
      CIj      <- quantile(theta.b[, j], c(alph/2, 1 - alph/2))
      CI[j, ]  <- c(theta.b[1, j], CIj)
    }
    round(CI)
  }

  # ------------ MAIN ----------------------------------------------------

  assignflg <- TRUE  # thetahat will set global ProWild / WildCollaps on first call

  FISHdate <- which(Dat     == names(All))
  FISHrear <- which(Rr      == names(All))
  FISHpndx <- which(Primary == names(All))
  if(!is.na(Secondary)) {
    FISHsndx  <- which(Secondary == names(All))
    FISHsgrps <- unique(All[, FISHsndx])
    nSgrps    <- length(FISHsgrps)
  }

  PASSstrat   <- which(strat    == names(pass))
  PASSdate    <- which(dat      == names(pass))
  PASSrate    <- which(samrate  == names(pass))
  PASScounts  <- which(tally    == names(pass))
  PASSguideff <- which(guidance == names(pass))
  PASScollaps <- which(collaps  == names(pass))

  ndays <- nrow(pass)

  Cpattern <- unique(cbind(pass[, PASSstrat], pass[, PASScollaps]))
  cat("\nStrata are collapsed according to:\n")
  temp <- t(Cpattern)
  rownames(temp) <- c("Week", "Strata")
  colnames(temp) <- rep("", ncol(temp))
  print(temp)

  pass$true      <- pass[, PASSrate] * pass[, PASSguideff]
  pass$estimated <- pass[, PASScounts] / pass$true

  passdata <- data.frame(Stratum = pass[, PASScollaps],
                         Tally   = pass[, PASScounts],
                         Ptrue   = pass$true)

  passweeks    <- tapply(pass$estimated, pass[, PASSstrat],   sum)
  rpassweeks   <- round(passweeks)
  cat("\nEstimate of total smolts by week:\n")
  print(rpassweeks)

  passcollaps  <- tapply(pass$estimated, pass[, PASScollaps], sum)
  rpasscollaps <- round(passcollaps)
  cat("\nEstimate of total smolts by statistical week:\n")
  print(rpasscollaps)

  totalpassage <- round(sum(passcollaps))
  cat("\nEstimate of total smolts:", totalpassage, "\n")

  # Assign collapsed stratum to each fish
  nAll        <- nrow(All)
  All$Collaps <- numeric(nAll)
  AllDates    <- unique(All[, FISHdate])
  for(d in AllDates) {
    CollStrat                          <- pass[which(pass[, PASSdate] == d), PASScollaps]
    All$Collaps[which(All[, FISHdate] == d)] <- CollStrat
  }

  # Assign true sampling rate to each fish
  set    <- intersect(All[, FISHdate], pass[, PASSdate])
  ndates <- length(set)
  All$true <- numeric(nAll)
  for(nn in 1:ndates) {
    passtemp <- pass[as.Date(pass[, PASSdate],  format = dateFormat) ==
                       as.Date(set[nn], origin = "1970-01-01", format = dateFormat), ]
    All$true[as.Date(All[, FISHdate], format = dateFormat) ==
               as.Date(set[nn], origin = "1970-01-01", format = dateFormat)] <- passtemp$true
  }

  RearData <- data.frame(Rear    = All[, FISHrear],
                         Stratum = All$Collaps,
                         True    = All$true)

  AllRTYPE   <- droplevels(All[which(All[, FISHrear] == RTYPE), ])
  AllPrimary <- droplevels(AllRTYPE[which(AllRTYPE[, FISHpndx] != "NA"), ])
  nFISH      <- nrow(AllPrimary)

  set    <- intersect(AllPrimary[, FISHdate], pass[, PASSdate])
  ndates <- length(set)

  Pgrps   <- sort(unique(AllPrimary[, FISHpndx]))
  nPgrps  <- length(Pgrps)
  p       <- 1 + nPgrps
  strats  <- unique(Cpattern[, 2])
  nstrats <- length(strats)

  tabl   <- table(AllPrimary[, FISHdate], AllPrimary[, FISHpndx])
  nPrime <- apply(tabl, 1, sum)

  AllPrimary$SR <- numeric(nFISH)
  for(nn in 1:ndates) {
    passtemp   <- pass[as.Date(pass[, PASSdate],  format = dateFormat) ==
                         as.Date(set[nn], origin = "1970-01-01", format = dateFormat), ]
    Primecount <- nPrime[which(as.Date(names(nPrime), format = dateFormat) ==
                                 as.Date(set[nn], origin = "1970-01-01", format = dateFormat))]
    AllPrimary$SR[as.Date(AllPrimary[, FISHdate], format = dateFormat) ==
                    as.Date(set[nn], origin = "1970-01-01", format = dateFormat)] <-
      passtemp$true * Primecount / passtemp[, PASScounts]
  }

  if(!is.na(Secondary)) {
    AllPrime <- data.frame(Strat = AllPrimary$Collaps,
                           PGrp  = AllPrimary[, FISHpndx],
                           SGrp  = AllPrimary[, FISHsndx],
                           SR    = AllPrimary$SR)
    Sgrps   <- unique(as.character(AllPrime$SGrp))
    Sgrps   <- sort(Sgrps[Sgrps != "NA"])
    nSgrps  <- length(Sgrps)
    p       <- p + nSgrps + nPgrps * nSgrps
    AvgProp <- getAvgProp(AllPrime)
  } else {
    AllPrime <- data.frame(Strat = AllPrimary$Collaps,
                           PGrp  = AllPrimary[, FISHpndx],
                           SR    = AllPrimary$SR)
  }

  ests      <- thetahat(passdata, RearData, AllPrime)
  assignflg <- FALSE
  pPropTable <- t(ests[[3]])
  if(!is.na(Secondary)) sAbunTable <- ests[[6]]

  if(B > 0) {
    answer <- bootsmolt(RearData, AllPrime, passdata)
    answer <- cbind(answer,
                    round((((answer[, 3] - answer[, 2]) / answer[, 1]) / 2) * 100, 1))
    colnames(answer) <- c("Estimate", "LCI", "UCI", "P1")

    if(!is.na(Secondary)) {
      grpnams <- ""
      for(prim in Pgrps)
        for(nam in Sgrps) grpnams <- c(grpnams, paste0(nam, prim))
      grpnams <- grpnams[-1]
      rownames(answer) <- c("WildSmolts", as.character(Pgrps), Sgrps, grpnams)
    } else {
      rownames(answer) <- c("WildSmolts", as.character(Pgrps))
    }
    cat("\n")
    print(answer)

    # Write strata / rearing summary
    tSmolts <- t(rbind(rpasscollaps, round(ProWild, 4), round(WildCollaps, 0)))
    colnames(tSmolts) <- c("TotalSmolts", "p(Wild)", "WildSmolts")
    write.csv(tSmolts, file = paste0(Run, "Rear.csv"))

    # Write CI output
    ciFile <- paste0(Run, "CIs.csv")
    header <- if(is.na(Secondary)) paste(RTYPE, "-", Primary) else
                                   paste(RTYPE, "-", Primary, "-", Secondary)
    write.table(header, file = ciFile, append = FALSE, row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    suppressWarnings(write.table(answer, file = ciFile, col.names = NA,
                                 sep = ",", append = TRUE))
    write.table(paste("RearSampleSize =",  nAll),  file = ciFile, append = TRUE,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(paste("PrimeSampleSize =", nFISH), file = ciFile, append = TRUE,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    if(!is.na(Secondary))
      write.table(paste("SecondSampleSize =", nsFish), file = ciFile, append = TRUE,
                  row.names = FALSE, col.names = FALSE, quote = FALSE)

    # Write primary composition table
    primeFile   <- paste0(Run, "Prime.csv")
    pFreqTable  <- table(AllPrimary$Collaps, AllPrimary[, Primary])
    pAbunTable  <- round(sweep(pPropTable, MARGIN = 1, WildCollaps, `*`), 0)
    PrimeTotals <- round(apply(pAbunTable, 2, sum), 0)
    write.table(pFreqTable,   file = primeFile, col.names = NA, sep = ",", append = FALSE)
    write.table(pPropTable,   file = primeFile, row.names = TRUE, col.names = FALSE,
                append = TRUE, sep = ",")
    write.table(pAbunTable,   file = primeFile, row.names = TRUE, col.names = FALSE,
                append = TRUE, sep = ",")
    write.table(t(PrimeTotals), file = primeFile, row.names = "PrimeTotals",
                col.names = FALSE, append = TRUE, sep = ",")

    # Write primary x secondary table
    if(!is.na(Secondary)) {
      sAbunTable <- cbind(sAbunTable, apply(sAbunTable, 1, sum))
      sAbunTable <- rbind(sAbunTable, apply(sAbunTable, 2, sum))
      rownames(sAbunTable) <- c(levels(Pgrps), "sTotals")
      colnames(sAbunTable) <- c(Sgrps, "pTotals")
      sAbunTable <- round(sAbunTable, 0)
      write.table(sAbunTable, file = paste0(Run, "PxS.csv"),
                  col.names = NA, append = FALSE, sep = ",")
    }
  } else {
    cat("\nNo bootstrap; just print estimates\n")
    print(ests)
  }

  cat("\nEnd time:", date(), "\n")
}
