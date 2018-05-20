#' \code{colorednoise} package
#'
#' Simulate Temporally Autocorrelated Populations
#'
#' See the README on
#' \href{https://cran.r-project.org/package=colorednoise/README.html}{CRAN}
#' or \href{https://github.com/japilo/colorednoise#readme}{GitHub}
#'
#' @docType package
#' @name colorednoise
#' @useDynLib colorednoise, .registration = TRUE
#' @import purrr
#' @import dplyr
#' @importFrom stats sd acf na.omit plogis
#' @importFrom purrrlyr by_row
NULL

## quiets concerns of R CMD check re: the .'s that appear in
## pipelines
if (getRversion() >= "2.15.1") utils::globalVariables(c(".", "mean.trans",
    "sd.trans", "noise", "timestep", "dist", "zero"))
