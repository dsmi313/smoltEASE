test_that("GRS preparation pools antennas by row and excludes later visits", {
  tags <- data.frame(tag=c("x","y"),cell=c(3L,4L),week=16L,period=c("day","night"))
  gd <- list(H=list(base=list(weeks=16L),tags=tags),
             W=list(base=list(weeks=16L),tags=tags[FALSE,]))
  dat <- data.frame(tag=c("x","x","x","y"),site="GRS",
    antenna=c("01","02","08","06"),
    det_time=as.POSIXct(c("2025-04-14 12:00:00","2025-04-14 12:00:01",
                         "2025-04-15 12:00:00","2025-04-14 23:00:00"),tz="UTC"))
  x <- prep_grs_detection(dat,gd)
  expect_equal(x$histories$history,c("100","010"))
  expect_equal(sum(x$counts),2)
  expect_equal(x$counts[1,1,1,"100"],1L)
  dat$antenna[1] <- "FF"
  expect_error(prep_grs_detection(dat,gd),"Unknown GRS")
})

test_that("joint row likelihood samples effective p=q*d and rejects mismatches", {
  skip_if_not_installed("rjags")
  day <- rbind(c(8,2,12,2,0,0),c(7,2,10,2,0,3),c(6,1,8,1,0,4),c(7,1,10,1,0,0))
  c5 <- rep(2L,4); n <- day*2; n[,5] <- c5
  g <- list(base=list(n=n,weeks=13:16,parent=1:4,
    lgr_spill_pct=c(25,45,65,35),lgs_spill_pct=c(20,40,60,30),
    lgr_outflow=c(70,85,100,80),alpha_phi_mean=0,alpha_phi_sd=2,
    beta_phi_mean=0,beta_phi_sd=1),n_day=day,n_night=day,n_c5=c5)
  dn <- list(H=g,W=g)
  cnt <- array(0L,c(2,4,2,7))
  for (r in 1:2) for (s in 1:4) for (t in 1:2) {
    total <- sum(day[s,3:4]); cnt[r,s,t,] <- c(1,1,1,1,1,1,total-6)
  }
  grs <- list(counts=cnt,weeks=13:16,rear_levels=c("H","W"))
  bad <- grs; bad$counts[1,1,1,1] <- 0L
  expect_error(fit_ge_rear_daynight(dn,grs_detection=bad),"exactly")
  fit <- suppressWarnings(fit_ge_rear_daynight(dn,grs_detection=grs,day_ge="shared",
    n_adapt=100,n_burnin=100,n_iter=200,n_chains=2,n_thin=2,verbose=FALSE))
  m <- do.call(rbind,lapply(fit$samples,as.matrix))
  expect_true(fit$joint_grs_detection)
  expect_equal(m[,"p[1,1]"],m[,"q[1,1]"]*m[,"array_detection[1]"],tolerance=1e-10)
  expect_equal(m[,"p_night[1,1]"],m[,"q_night[1,1]"]*m[,"array_detection[2]"],tolerance=1e-10)
  d <- generate_rear_ge_draws_daynight(fit,pass_dates=as.Date("2025-03-24"),B=5,
    daily_spill=data.frame(Date=as.Date("2025-03-24"),spill.per=25,outflow=70),
    spill_counts=data.frame(Date=as.Date("2025-03-24"),S_day=5,S_night=5),seed=1)
  expect_true(all(is.finite(as.matrix(d[,-1]))))
})
