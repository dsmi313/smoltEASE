# Run AFTER your MY2025 GE comparison script in the same R session.
# Uses its cleaned events, day/night prep, covariates and trap/FPC data.
# Antenna rows are assumed independent conditional on passage period.
# Array detection pooled across rear/week; spillbay routing q is rear/week-specific.
# q_night_offset=0 assumes equal day/night spillbay routing; test sensitivity.
stopifnot(MY == 2025L, packageVersion("smoltEASE") >= "0.1.9")
DET_DIR <- "C:/Users/david.smith/OneDrive - State of Idaho/Desktop/SCRAPI/SCRAPI2/Multi-Year Model/Results no U/MY Tests/2025/Results/detection"
dir.create(DET_DIR, recursive=TRUE, showWarnings=FALSE)
grs <- smoltEASE::prep_grs_detection(events,ge_dn)
readr::write_csv(grs$histories,file.path(DET_DIR,"GRS_row_histories.csv"))
saveRDS(grs,file.path(DET_DIR,"GRS_detection_data.rds"))
# New decomposition changes the prior from effective p to spillbay routing q.
# The before/after comparison includes this change, not just added row data.
after <- fit_cached(file.path(DET_DIR,"fit_joint_GRS.rds"),
  list(grs=grs$counts,ge_dn=ge_dn,psi_spill=psi_spill,psi_outflow=psi_outflow,
       settings=dn_settings,trans_on=TRANS_ON,q_night_offset=0),
  function() do.call(smoltEASE::fit_ge_rear_daynight,c(list(ge_dn=ge_dn,
    weeks=WEEKS,target_rear="W",parent=ge_data$W$parent,psi_spill=psi_spill,
    psi_outflow=psi_outflow,trans_on=TRANS_ON,day_ge="shared",
    grs_detection=grs,q_night_offset=0),dn_settings)))
if (!isTRUE(after$diagnostics$core_ge_rhat_pass))
  warning("Joint detection fit has flagged core convergence; results provisional.")
readr::write_csv(after$summary,file.path(DET_DIR,"joint_posterior_summary.csv"))
readr::write_csv(after$summary[grepl("^(array_detection|row_detection|q|q_night|p|p_night)\\[",after$summary$parameter),],
                file.path(DET_DIR,"routing_and_detection_summary.csv"))
smolt_data <- read.csv(FILES$smolt,check.names=FALSE,stringsAsFactors=FALSE)
smolt_data$CollectionDate <- parse_date(smolt_data$CollectionDate)
comparison <- list(); stratum_comparison <- list()
for (label in c("before","after")) {
  fit <- if (label=="before") dn_fit else after
  draws <- smoltEASE::generate_rear_ge_draws_daynight(fit,rear_type="W",
    pass_dates=pass_dates,B=B,daily_spill=lagged,
    spill_counts=dplyr::filter(spill_counts,rear_type=="W"),shrink_k=SHRINK_K,seed=SEED+1L)
  saveRDS(draws,file.path(DET_DIR,paste0("daily_GE_",label,".rds")))
  pass_k <- passage
  pass_k$GuidanceEfficiency <- 1/rowMeans(1/as.matrix(draws[,-1]))
  prefix <- file.path(DET_DIR,paste0("MY2025_W_",label,"_"))
  res <- smoltEASE::SCRAPI2(smoltData=smolt_data,Dat="CollectionDate",Rr="Rear",
    Primary=PRIMARY,Secondary=SECONDARY,passageData=pass_k,strat="Week",dat="SampleEndDate",
    tally="SampleCount",samrate="SampleRate",guidance="GuidanceEfficiency",collaps="Collapse",
    Run=prefix,RTYPE="W",REARSTRAT=TRUE,alph=ALPHA,B=B,dateFormat="%Y-%m-%d",
    gsiDraws=NULL,fishID="MasterID",geDraws=draws,seed=SEED+2L,pointEst="mean")
  saveRDS(res,file.path(DET_DIR,paste0("SCRAPI2_",label,".rds")))
  comparison[[label]] <- ci_table(res,label)
  stratum_comparison[[label]] <- read_rear(paste0(prefix,"Rear.csv")) |>
    dplyr::mutate(run=label)
}
comparison <- dplyr::bind_rows(comparison)
ref <- comparison |> dplyr::filter(run=="before") |>
  dplyr::select(group,before_estimate=estimate)
comparison <- comparison |> dplyr::left_join(ref,by="group") |>
  dplyr::mutate(change=estimate-before_estimate,percent_change=100*change/before_estimate)
readr::write_csv(comparison,file.path(DET_DIR,"SCRAPI_before_after.csv"))
readr::write_csv(dplyr::bind_rows(stratum_comparison),file.path(DET_DIR,"SCRAPI_strata_before_after.csv"))
print(as.data.frame(comparison |> dplyr::filter(group %in% c("WildSmolts","GRROND","IMNAHA"))))
cat("Outputs:",DET_DIR,"\n")
