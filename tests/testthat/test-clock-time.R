test_that("prep_pit_data2 reads PTAGIS AM/PM times on a 24-hour clock", {
  pit <- data.frame(
    `Tag Code` = c("A", "B", "C", "D", "E"),
    `Species Name` = "Steelhead", `Rear Type Code` = "W",
    `Site Name` = "GRJ - Lower Granite Dam Juvenile",
    `Obs Time Value` = c("5/1/2020 5:05:06 PM", "5/1/2020 5:05:06 AM",
                         "5/1/2020 12:10:00 AM", "5/1/2020 12:10:00 PM",
                         "05/01/2020 17:05:06"),
    `Antenna ID` = "01", `Release Site Name` = "JOHTRP - Johnson Creek Trap",
    check.names = FALSE, stringsAsFactors = FALSE)
  out <- prep_pit_data2(pit)
  expect_equal(format(out$det_time, "%H:%M:%S"),
               c("17:05:06", "05:05:06", "00:10:00", "12:10:00", "17:05:06"))
  expect_equal(out$det_date, rep(as.Date("2020-05-01"), 5))
})

test_that("unparseable PIT times still stop prep_pit_data2", {
  pit <- data.frame(
    `Tag Code` = "A", `Species Name` = "Steelhead", `Rear Type Code` = "W",
    `Site Name` = "GRJ - Lower Granite Dam Juvenile",
    `Obs Time Value` = "not a time", `Antenna ID` = "01",
    `Release Site Name` = "JOHTRP - Johnson Creek Trap",
    check.names = FALSE, stringsAsFactors = FALSE)
  expect_error(prep_pit_data2(pit), "unparseable")
})
