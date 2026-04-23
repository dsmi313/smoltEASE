#' Simulate gsiDraws for testing SCRAPI2
#'
#' @param smoltData data frame or CSV path.
#' @param Primary stock assignment column. Default \code{"GenStock"}.
#' @param fishID fish ID column. Default \code{"MasterID"}.
#' @param B number of bootstrap draws. Default \code{5000}.
#' @param n_point number of point-estimate draws. Default \code{100}.
#' @param method \code{"fixed"}, \code{"dirichlet"} (default), or \code{"population"}.
#' @param alpha_self Dirichlet concentration on observed stock. Default \code{10}.
#' @param alpha_other Dirichlet concentration on all other stocks. Default \code{1}.
#' @param seed optional integer for reproducibility.
#'
#' @export
sim_gsi_draws <- function(smoltData, Primary = "GenStock", fishID = "MasterID",
                          B = 5000, n_point = 100,
                          method = c("fixed", "dirichlet", "population"),
                          alpha_self = 10, alpha_other = 1,
                          seed = NULL) {

  method <- match.arg(method)
  if (is.character(smoltData)) smoltData <- read.csv(smoltData, header = TRUE)
  if (!fishID  %in% names(smoltData)) stop("fishID column '",  fishID,  "' not found")
  if (!Primary %in% names(smoltData)) stop("Primary column '", Primary, "' not found")

  ids     <- smoltData[[fishID]]
  stocks  <- as.character(smoltData[[Primary]])
  n_fish  <- length(ids)
  n_draws <- max(B, n_point)

  # NA guard: treat both R NA and the string "NA" as missing
  is_na <- is.na(stocks) | stocks == "NA"

  valid <- stocks[!is_na]
  if (length(valid) == 0) stop("No non-NA stock assignments found in '", Primary, "'")
  stock_levels <- sort(unique(valid))
  n_stocks     <- length(stock_levels)

  if (!is.null(seed)) set.seed(seed)

  draws_mat <- if (method == "fixed") {

    m <- matrix(stocks, nrow = n_fish, ncol = n_draws)
    m[is_na, ] <- NA_character_
    m

  } else if (method == "population") {

    stock_probs <- tabulate(factor(valid, levels = stock_levels)) / length(valid)
    m <- matrix(
      sample(stock_levels, size = n_fish * n_draws, replace = TRUE, prob = stock_probs),
      nrow = n_fish, ncol = n_draws
    )
    m[is_na, ] <- NA_character_
    m

  } else {
    # Dirichlet-multinomial, vectorized
    # alpha matrix: n_fish x n_stocks
    alpha_mat <- matrix(alpha_other, nrow = n_fish, ncol = n_stocks,
                        dimnames = list(NULL, stock_levels))
    fish_with_stock <- which(!is_na)
    obs_idx <- match(stocks[fish_with_stock], stock_levels)
    alpha_mat[cbind(fish_with_stock, obs_idx)] <- alpha_self

    # Draw gammas: need n_fish x n_stocks x n_draws
    # Do it draw by draw to keep memory manageable
    m <- matrix(NA_character_, nrow = n_fish, ncol = n_draws)
    for (d in seq_len(n_draws)) {
      g     <- matrix(rgamma(n_fish * n_stocks, shape = as.vector(alpha_mat)),
                      nrow = n_fish, ncol = n_stocks)
      probs <- g / rowSums(g)
      m[fish_with_stock, d] <- stock_levels[
        apply(probs[fish_with_stock, , drop = FALSE], 1,
              function(p) sample.int(n_stocks, 1L, prob = p))
      ]
    }
    m
  }

  result           <- as.data.frame(draws_mat, stringsAsFactors = FALSE)
  colnames(result) <- paste0("boot_", seq_len(n_draws))
  cbind(setNames(data.frame(ids, stringsAsFactors = FALSE), fishID), result)
}
