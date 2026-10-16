test_that("exported GE fits accept structure and route prior options", {
  for (name in c("fit_ge_rear_model2", "fit_ge_rear_daynight")) {
    args <- formals(getExportedValue("smoltEASE", name))
    expect_true(all(c("structure", "delta_prior") %in% names(args)))
    expect_identical(eval(args$structure), c("original", "consistent"))
    expect_identical(eval(args$delta_prior), c("tight", "loose"))
  }
})
