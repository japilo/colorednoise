#' This function generates temporally autocorrelated random numbers with a mean, standard deviation,
#' and autocorrelation you specify.
#' @param timesteps The number of temporally autocorrelated random numbers (one per timestep) you want.
#' @param mu The mean of the temporally autocorrelated random numbers.
#' @param sigma The standard deviation of the temporally autocorrelated random numbers.
#' @param phi The temporal autocorrelation. 0 is white noise (uncorrelated), positive values are red noise (correlated) and negative values are blue noise (inversely correlated).
#' @return A vector of temporally autocorrelated random numbers.
#' @example
#' raw_noise(timesteps = 30, mu = 0.5, sigma = 0.2, phi = 0.3)
raw_noise <- function(timesteps, mu, sigma, phi) {
    delta <- mu * (1 - phi)
    variance <- sigma^2 * (1 - phi^2)
    noise <- c(rnorm(1, mu, sigma))
    for (i in (1:(timesteps - 1))) {
      noise[i + 1] <- delta + phi * noise[i] + rnorm(1, 0, sqrt(variance))
    }
    return(noise)
}
raw_estim <- function(noise) {
    ac <- acf(x, lag = 1, plot = F)[[1]][2]
    xbar <- mean(x)
    stdev <- sd(x)
    output <- c(xbar, stdev, ac)
    names(output) <- c("Mean", "SD", "Autocorrelation")
    return(output)
}
raw_noise_loop <- function(years, mu, sigma, phi, replicates) {
  # Simulates raw noise for each combination of parameters -----------
    raw_sims <- cross_df(list(years = years, phi = phi, mu = mu, sigma = sigma)) %>%
        rerun(.n = replicates, pmap(., raw))
  # Labels each simulation with the parameters that generated it ------
    labeled_sims <- vector(mode = "list", replicates)
    for (i in 1:length(raw_sims)) {
        labeled_sims[[i]] <- vector(mode = "list", nrow(raw_sims[[1]][[1]]))
        for (j in 1:length(raw_sims[[i]][[2]])) {
            labeled_sims[[i]][[j]] <- cbind(t(raw_sims[[i]][[2]][[j]]), raw_sims[[i]][[1]][j,
                ])
        }
    }
  # Unnests the list and binds into one data frame ---------------------
    labeled_sims %>% flatten() %>% bind_rows()
}
