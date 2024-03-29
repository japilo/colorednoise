#' @keywords internal
"_PACKAGE"

#' \code{colorednoise} package
#'
#' Simulate Temporally Autocorrelated Populations
#'
#' See the README on \href{https://github.com/japilo/colorednoise#readme}{GitHub}
#'
#' @name colorednoise-package
#' @useDynLib colorednoise, .registration = TRUE
#' @import purrr
#' @importFrom Rcpp evalCpp
#' @rawNamespace import(data.table, except = transpose)
#' @importFrom stats sd acf na.omit plogis
NULL

## quiets concerns of R CMD check re: the .'s that appear in
## pipelines
if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "mean.trans",
    "sd.trans", "noise", "timestep", "dist", "zero", "ref"))
