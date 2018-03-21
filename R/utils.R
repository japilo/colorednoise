#' Estimate the Temporal Autocorrelation of a Numeric Vector
#'
#' A wrapper for the \code{\link[stats]{acf}} function from the stats package that
#' extracts only the temporal autocorrelation at a lag of one timestep (which is
#' the type of temporal autocorrelation that this package simulates).
#'
#' @importFrom stats acf
#' @param x A numeric vector.
#' @param biasCorrection Autocorrelation estimates are biased for short time series. The function can
#' correct for this bias in the manner proposed by Quenouille (1949). Set to TRUE by default.
#' @param na.action The way the function will handle NAs in the vector. Set to \code{\link[stats]{na.pass}} by default.
#' @return A single numeric value: the estimate of the temporal autocorrelation with a lag of 1.
#' @examples
#' rednoise <- colored_noise(timesteps = 50, mu = 0.5, sigma = 0.2, phi = 0.3)
#' autocorrelation(rednoise)
#' @export
autocorrelation <- function(x, biasCorrection = TRUE, na.action = na.omit) {
  if (biasCorrection == TRUE) {
    r <- acf(x, plot = F, na.action = na.action)[[1]][2]
    r1 <- acf(x[1:round(length(x)/2)], plot = F, na.action = na.action)[[1]][2]
    r2 <- acf(x[round(length(x)/2)+1:length(x)], plot = F, na.action = na.action)[[1]][2]
    2*r - 0.5*(r1 + r2)
  } else if (biasCorrection == FALSE) {
    acf(x, plot = F, na.action = na.action)[[1]][2]}
}
