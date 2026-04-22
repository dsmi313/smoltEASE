#' @title Query LGD trapping database
#'
#' @description Pulls smolt trap biological data from the LGDTrapping SQL
#'   Server database (\code{TempLGTrappingExcelPivotfromvwLGDMasterCombine}).
#'   Returns the result as a data frame and optionally writes a CSV.
#'
#' @param spawn_year migratory year string, e.g. \code{"MY2025"}.
#' @param srr_prefix first character(s) of the SRR code used to filter species.
#'   Use \code{"1"} for Chinook and \code{"3"} for steelhead.
#' @param life_stage LGDLifeStage code. Use \code{"JY"} for Chinook smolts and
#'   \code{"JV"} for steelhead smolts.
#' @param server SQL Server instance name. Default \code{"idfgnrsql2"}.
#' @param database database name. Default \code{"LGDTrapping"}.
#' @param driver ODBC driver name. Default \code{"ODBC Driver 17 for SQL Server"}.
#' @param export_csv optional file path. When supplied the result is written as
#'   a CSV (no row names) in addition to being returned.
#'
#' @return A data frame of trap records matching the query filters.
#'
#' @details Requires the \pkg{DBI} and \pkg{odbc} packages and Windows
#'   Integrated Authentication (\code{Trusted_Connection = yes}) to the
#'   IDFG network. The underlying view is
#'   \code{TempLGTrappingExcelPivotfromvwLGDMasterCombine}.
#'
#'   Before formatting with \code{lgr2SCRAPI()} check that \code{GenStock}
#'   uses the correct 6-letter codes (e.g. \code{LOSALM} not \code{LOWSALM}).
#'
#' @examples
#' \dontrun{
#' chnk_raw <- query_lgd_trap("MY2025", srr_prefix = "1", life_stage = "JY",
#'                             export_csv = "MY2025.CHNKsmolts_trapBioData_allVariables.csv")
#' sthd_raw <- query_lgd_trap("MY2025", srr_prefix = "3", life_stage = "JV",
#'                             export_csv = "MY2025.STHDsmolts_trapBioData_allVariables.csv")
#' }
#'
#' @export
query_lgd_trap <- function(spawn_year,
                            srr_prefix,
                            life_stage,
                            server     = "idfgnrsql2",
                            database   = "LGDTrapping",
                            driver     = "ODBC Driver 17 for SQL Server",
                            export_csv = NULL) {

  if (!requireNamespace("DBI",  quietly = TRUE) ||
      !requireNamespace("odbc", quietly = TRUE))
    stop("Packages 'DBI' and 'odbc' are required. Install with:\n",
         "  install.packages(c('DBI', 'odbc'))")

  con <- DBI::dbConnect(
    odbc::odbc(),
    Driver             = driver,
    Server             = server,
    Database           = database,
    Trusted_Connection = "yes"
  )
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  sql <- sprintf(
    "SELECT * FROM TempLGTrappingExcelPivotfromvwLGDMasterCombine
     WHERE SpawnYear = '%s' AND SRR LIKE '%s%%' AND LGDLifeStage = '%s'",
    spawn_year, srr_prefix, life_stage
  )

  out <- DBI::dbGetQuery(con, sql)

  if (!is.null(export_csv))
    write.csv(out, export_csv, row.names = FALSE)

  out
}
