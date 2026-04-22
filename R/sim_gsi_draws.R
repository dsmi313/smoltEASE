#' Simulate gsiDraws for testing SCRAPI2
#'
#' Creates a placeholder \code{gsiDraws} data frame for use with
#' \code{\link{SCRAPI2}} while waiting for real posterior draws from the
#' genetics lab.
#'
#' Two modes are available:
#' \describe{
#'   \item{\code{"fixed"}}{Every draw repeats the fish's observed stock
#'     assignment. GSI contributes zero additional uncertainty — useful for
#'     verifying the pipeline before real draws arrive.}
#'   \item{\code{"population"}}{Each draw independently samples every fish's
#'     assignment from the observed population-level frequency distribution.
#'     This is a maximum-uncertainty baseline: it ignores all individual-level
#'     information and treats every fish as equally likely to belong to any
#'     stock.}
#' }
#'
#' @param smoltData data frame of individual fish records, or a character path
#'   to a CSV file.
#' @param Primary column name in \code{smoltData} that holds the stock
#'   assignment. Default \code{"GenStock"}.
#' @param fishID column name in \code{smoltData} that holds individual fish IDs.
#'   Default \code{"MasterID"}.
#' @param B number of bootstrap draws (should match the \code{B} you will pass
#'   to \code{SCRAPI2}). Default \code{5000}.
#' @param n_point number of point-estimate draws (should match the
#'   \code{n_point} you will pass to \code{SCRAPI2}). Default \code{100}.
#' @param method \code{"fixed"} (default) or \code{"population"}. See Details.
#' @param seed optional integer passed to \code{set.seed()} for reproducibility
#'   when \code{method = "population"}.
#'
#' @return A data frame with \code{max(B, n_point) + 1} columns. Column 1 is
#'   named after \code{fishID}; the remaining columns are named
#'   \code{boot_1}, \code{boot_2}, \ldots{} as required by
#'   \code{\link{SCRAPI2}}.
#'
#' @examples
#' \dontrun{
#' # Pipeline test: fixed assignments, no GSI uncertainty
#' gsiDraws <- sim_gsi_draws(STHD.trap, Primary = "GenStock", B = 2000)
#' result   <- SCRAPI2(smoltData = STHD.trap, ..., gsiDraws = gsiDraws, B = 2000)
#'
#' # Maximum-uncertainty baseline
#' gsiDraws <- sim_gsi_draws(STHD.trap, method = "population", B = 2000, seed = 42)
#' }
#'
#' @export
sim_gsi_draws <- function(smoltData, Primary = "GenStock", fishID = "MasterID",
                          B = 5000, n_point = 100,
                          method = c("fixed", "population"),
                          seed = NULL) {
  method <- match.arg(method)
  if (is.character(smoltData)) smoltData <- read.csv(smoltData, header = TRUE)

  if (!fishID %in% names(smoltData))
    stop("fishID column '", fishID, "' not found in smoltData")
  if (!Primary %in% names(smoltData))
    stop("Primary column '", Primary, "' not found in smoltData")

  ids    <- smoltData[[fishID]]
  stocks <- as.character(smoltData[[Primary]])
  n_fish <- length(ids)
  n_draws <- max(B, n_point)

  valid <- stocks[!is.na(stocks) & stocks != "NA"]
  if (length(valid) == 0)
    stop("No non-NA stock assignments found in '", Primary, "' column")
  stock_levels <- sort(unique(valid))

  if (!is.null(seed)) set.seed(seed)

  draws_mat <- if (method == "fixed") {
    matrix(stocks, nrow = n_fish, ncol = n_draws)
  } else {
    stock_probs <- tabulate(factor(valid, levels = stock_levels)) / length(valid)
    replicate(n_draws,
              sample(stock_levels, size = n_fish, replace = TRUE, prob = stock_probs))
  }

  result <- as.data.frame(draws_mat, stringsAsFactors = FALSE)
  colnames(result) <- paste0("boot_", seq_len(n_draws))
  result <- cbind(setNames(data.frame(ids, stringsAsFactors = FALSE), fishID),
                  result)
  result
}
