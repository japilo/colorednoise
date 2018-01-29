#' Estimate the Temporal Autocorrelation of a Numeric Vector
#'
#' A wrapper for the \code{\link[stats]{acf}} function from the stats package that
#' extracts only the temporal autocorrelation at a lag of one timestep (which is
#' the type of temporal autocorrelation that this package simulates).
#'
#' @param x A numeric vector.
#' @param na.action The way the function will handle NAs in the vector. Set to \code{\link[stats]{na.pass}} by default.
#' @return A single numeric value: the estimate of the temporal autocorrelation with a lag of 1.
#' @examples
#' rednoise <- raw_noise(timesteps = 50, mu = 0.5, sigma = 0.2, phi = 0.3)
#' autocorrelation(rednoise)
#' @export
autocorrelation <- function(x, na.action = na.omit) {
    acf(x, plot = F, na.action = na.action)[[1]][2]
}
