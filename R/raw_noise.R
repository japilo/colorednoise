#' Generate Autocorrelated Noise
#'
#' This function generates temporally autocorrelated random numbers with a mean,
#' standard deviation, and autocorrelation you specify.
#' @param timesteps The number of temporally autocorrelated random numbers (one
#'   per timestep) you want.
#' @param mu The mean of the temporally autocorrelated random numbers.
#' @param sigma The standard deviation of the temporally autocorrelated random
#'   numbers.
#' @param phi The temporal autocorrelation. 0 is white noise (uncorrelated),
#'   positive values are red noise (directly correlated) and negative values are blue
#'   noise (inversely correlated).
#' @return A vector of temporally autocorrelated random numbers.
#' @examples
#' rednoise <- raw_noise(timesteps = 30, mu = 0.5, sigma = 0.2, phi = 0.3)
#' rednoise
#' @export
raw_noise <- function(timesteps, mu, sigma, phi) {
    delta <- mu * (1 - phi)
    variance <- sigma^2 * (1 - phi^2)
    noise <- vector(mode = "double", length = timesteps)
    noise[1] <- c(rnorm(1, mu, sigma))
    for (i in (1:(timesteps - 1))) {
        noise[i + 1] <- delta + phi * noise[i] + rnorm(1, 0, sqrt(variance))
    }
    return(noise)
}
#' Estimate Mean, SD, and Autocorrelation of Sample Noise.
#'
#' This function estimates the temporal autocorrelation of a vector of random
#' numburs, as well as the sample mean and standard deviation. Try feeding it the output of \code{\link{raw_noise}}.
#' @param noise The vector of random numbers.
#' @return A labeled vector with the sample mean, sample SD, and sample
#'   autocorrelation.
#' @examples
#' rednoise <- raw_noise(timesteps = 30, mu = 0.5, sigma = 0.2, phi = 0.3)
#' raw_estim(rednoise)
#' @export
raw_estim <- function(noise) {
    ac <- autocorrelation(noise)
    xbar <- mean(noise)
    stdev <- sd(noise)
    output <- c(xbar, stdev, ac)
    names(output) <- c("Mean", "SD", "Autocorrelation")
    return(output)
}
#' Generate Autocorrelated Noise for Every Combination of the Given Parameters.
#'
#' This function generates sets of temporally autocorrelated random numbers for
#' every possible combination of parameter values you specify. Essentially a
#' loop of \code{\link{raw_noise}} that outputs a list. All parameters can be
#' given as single values or vectors of values.
#' @param timesteps How many timesteps you want in each set. Can be scalar or
#'   vector.
#' @param mu The mean of the temporally autocorrelated random numbers. Can be
#'   scalar or vector.
#' @param sigma The standard deviation of the temporally autocorrelated random
#'   numbers. Can be scalar or vector.
#' @param phi The temporal autocorrelation. 0 is white noise (uncorrelated),
#'   positive values are red noise (correlated) and negative values are blue
#'   noise (inversely correlated). Can be scalar or vector.
#' @param replicates How many replicates you would like of each possible
#'   combination of parameters.
#' @return A list of vectors of temporally autocorrelated random numbers. Each
#'   element in the list is named with the combination of parameters that
#'   generated it.
#' @examples
#' loop <- raw_noise_loop(timesteps=c(5:10), mu=c(0.2, 0.5), sigma=c(0.2, 0.5),
#'                        phi=c(0, 0.1), replicates=10)
#' loop[[1]]
#' @export
raw_noise_loop <- function(timesteps, mu, sigma, phi, replicates) {
    . <- NULL
    # Simulates raw noise for each combination of parameters -----------
    raw_sims <- cross_df(list(timesteps = timesteps, phi = phi, mu = mu, sigma = sigma)) %>%
        rerun(.n = replicates, pmap(., raw_noise))
    # Labels each element of the list with the parameters that generated it -----
    for (j in 1:length(raw_sims)) {
        for (i in 1:length(raw_sims[[j]][[2]])) {
            names(raw_sims[[j]][[2]])[i] <- paste(c(rbind(names(raw_sims[[j]][[1]]),
                as.character(raw_sims[[j]][[1]][i, ]))), collapse = " ")
        }
    }
    # Unnests the list ------------------
    map(raw_sims, 2) %>% flatten()
}
#' Simulation and Estimation of Colored Noise
#'
#' This function simulates generates sets of temporally autocorrelated random
#' numbers for every possible combination of parameter values you specify, then
#' estimates the mean, standard deviation, and autocorrelation of each set of
#' random numbers so generated. Internally, the function does the same thing as
#' \code{\link{raw_noise_loop}}, but instead of outputting the random numbers, it
#' outputs measures of each set of random numbers.
#' @importFrom dplyr bind_rows
#' @inheritParams raw_noise_loop
#' @return A data frame with one row for each set of random numbers. The
#'   variables are mean, SD, autocorrelation, and the four parameters used to
#'   generate the random numbers.
#' @examples
#' estimates <- raw_estim_loop(timesteps=c(5:10), mu=c(0.2, 0.5), sigma=c(0.2, 0.5),
#'                             phi=c(0, 0.1), replicates=10)
#' head(estimates)
#' @export
raw_estim_loop <- function(timesteps, mu, sigma, phi, replicates) {
    . <- NULL
    raw_sims <- cross_df(list(timesteps = timesteps, phi = phi, mu = mu, sigma = sigma)) %>%
        rerun(.n = replicates, pmap(., raw_noise))
    # Output a data frame where each row is the mean, SD, and AC of a noise vector
    labels <- raw_sims %>% map(1) %>% bind_rows()
    noise <- raw_sims %>% map(2) %>% flatten()
    estimates <- data.frame(mean = map_dbl(noise, mean), sd = map_dbl(noise, sd), autocorrelation = map_dbl(noise,
        autocorrelation))
    cbind(labels, estimates)
}
