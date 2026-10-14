# Parse clock times written either on a 12-hour clock with AM/PM (as PTAGIS
# writes Obs Time Value, e.g. "5/1/2020 5:05:06 PM") or on a 24-hour clock.
# strptime() ignores trailing text, so a 24-hour format would read
# "5:05:06 PM" as 05:05:06. Times ending in AM or PM are therefore parsed
# only with the 12-hour formats, and all other times only with the 24-hour
# formats. Unparsed values are returned as NA for the caller to handle.
.parse_clock_time <- function(x, tz) {
  x <- trimws(as.character(x))
  out <- as.POSIXct(rep(NA_character_, length(x)), tz = tz)
  ok <- !is.na(x) & nzchar(x)
  ampm <- ok & grepl("[AaPp][Mm]$", x)
  groups <- list(
    list(rows = which(ampm),
         formats = c("%m/%d/%Y %I:%M:%S %p", "%m/%d/%Y %I:%M %p",
                     "%Y-%m-%d %I:%M:%S %p")),
    list(rows = which(ok & !ampm),
         formats = c("%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M",
                     "%Y-%m-%d %H:%M:%S")))
  for (g in groups) {
    for (fmt in g$formats) {
      todo <- g$rows[is.na(out[g$rows])]
      if (!length(todo)) break
      out[todo] <- as.POSIXct(x[todo], format = fmt, tz = tz)
    }
  }
  out
}
