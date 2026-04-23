#' Simulate gsiDraws for testing SCRAPI2
#'
#' Creates a placeholder \code{gsiDraws} data frame for use with
#' \code{\link{SCRAPI2}} while waiting for real posterior draws from the
#' genetics lab.
#'
#' Three modes are available:
#' \describe{
#'   \item{\code{"fixed"}}{Every draw repeats the fish's observed stock
#'     assignment. GSI contributes zero additional uncertainty — useful for
#'     verifying the pipeline before real draws arrive.}
#'   \item{\code{"population"}}{Each draw independently samples every fish's
#'     assignment from the observed population-level frequency distribution.
#'     This is a maximum-uncertainty baseline: it ignores all individual-level
#'     information and treats every fish as equally likely to belong to any
#'     stock.}
#'   \item{\code{"dirichlet"}}{Each draw samples a stock assignment for each
#'     fish from a Dirichlet-multinomial distribution. The fish's observed stock
#'     gets concentration \code{alpha_self}; all other stocks get
#'     \code{alpha_other}. Fish with NA stock retain NA in all draws. Lower
#'     \code{alpha_self} relative to \code{alpha_other} injects more
#'     uncertainty; higher values approach \code{"fixed"} behaviour.}
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
#' @param n_point number of point-estimate draws. Default \code{100}.
#' @param method \code{"fixed"} (default), \code{"population"}, or
#'   \code{"dirichlet"}. See Details.
#' @param alpha_self Dirichlet concentration on the observed stock
#'   (\code{method = "dirichlet"} only). Higher values = less uncertainty.
#'   Default \code{10}.
#' @param alpha_other Dirichlet concentration on all other stocks
#'   (\code{method = "dirichlet"} only). Default \code{1}.
#' @param seed optional integer passed to \code{set.seed()} for reproducibility.
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
#'
#' # Dirichlet uncertainty: moderate (alpha_self = 10, alpha_other = 1)
#' gsiDraws <- sim_gsi_draws(STHD.trap, method = "dirichlet",
#'                            alpha_self = 10, alpha_other = 1, B = 2000, seed = 42)
#'
#' # Maximum-uncertainty baseline
#' gsiDraws <- sim_gsi_draws(STHD.trap, method = "population", B = 2000, seed = 42)
#' }
#'
#' @export
sim_gsi_draws <- function(smoltData, Primary = "GenStock", fishID = "MasterID",
                          B = 5000, n_point = 100,
                          method = c("fixed", "dirichlet", "population"),
                          alpha_self = 10, alpha_other = 1,
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
  n_stocks <- length(stock_levels)

  if (!is.null(seed)) set.seed(seed)

  draws_mat <- if (method == "fixed") {

    matrix(stocks, nrow = n_fish, ncol = n_draws)

  } else if (method == "population") {

    stock_probs <- tabulate(factor(valid, levels = stock_levels)) / length(valid)
    replicate(n_draws,
              ifelse(is.na(stocks) | stocks == "NA", NA_character_,
                     sample(stock_levels, size = n_fish, replace = TRUE,
                            prob = stock_probs)))

  } else {
    # dirichlet: per-fish Dirichlet-multinomial draws
    # For each fish with a known stock, build a Dirichlet concentration vector:
    # alpha_self on the observed stock, alpha_other on all others.
    # Sample one stock label per draw from the resulting probabilities.

    rdirichlet_one <- function(alpha) {
      x <- rgamma(length(alpha), shape = alpha, rate = 1)
      x / sum(x)
    }

    draws_mat <- matrix(NA_character_, nrow = n_fish, ncol = n_draws)

    for (i in seq_len(n_fish)) {
      s <- stocks[i]
      if (is.na(s) || s == "NA") next
      alpha_vec <- rep(alpha_other, n_stocks)
      names(alpha_vec) <- stock_levels
      alpha_vec[s] <- alpha_self
      # Draw n_draws probability vectors, sample one stock each
      probs_mat <- replicate(n_draws, rdirichlet_one(alpha_vec))  # n_stocks x n_draws
      draws_mat[i, ] <- stock_levels[apply(probs_mat, 2,
                                           function(p) sample.int(n_stocks, 1, prob = p))]
    }

    draws_mat
  }

  result <- as.data.frame(draws_mat, stringsAsFactors = FALSE)
  colnames(result) <- paste0("boot_", seq_len(n_draws))
  result <- cbind(setNames(data.frame(ids, stringsAsFactors = FALSE), fishID),
                  result)
  result
}
