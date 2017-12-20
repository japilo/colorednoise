
#' Simulate a Time Series of a Temporally Autocorrelated Population
#'
#' The backbone of this package, this function simulates a population with
#' temporally autocorrelated vital rates (survival and fertility). In other
#' words, this function will show you the dynamics over time of a population
#' whose survival and fertility is stochastic, but also correlated to the
#' survival and fertility in the previous year, respectively. The assumptions of the
#' simulation are that the population is asexually reproducing or female-only,
#' individuals are able to reproduce in their first timestep, survival and
#' fertility are the same at all ages, and that individuals continue to be reproductively capable
#' until they die.
#'
#' Be advised that not all combinations of values will work. Certain combinations
#' of mean and SD for survival and fertility are very unrealistic for real
#' populations, and you'll end up with mostly NAs in your output. Use common sense
#' as a demogropher / population biologist.
#' @import stats
#' @param start The starting population size.
#' @param timesteps The number of timesteps you want to simulate. Individuals
#'   are added and killed off every timestep according to the survival and
#'   fertility rates. In ecological applications, timesteps are usually years,
#'   but theoretically they can be any length of time.
#' @param survPhi The temporal autocorrelation of survival. 0 is white noise (uncorrelated),
#'   positive values are red noise (directly correlated) and negative values are
#'   blue noise (inversely correlated).
#' @param fecundPhi The temporal autocorrelation of fecundity. As above.
#' @param survMean The mean survival from timestep to timestep. Must be a value
#'   between 0 (all individuals die) and 1 (all individuals live).
#' @param survSd The standard deviation of the survival from timestep to
#'   timestep. Must be a value between 0 and 1.
#' @param fecundMean The mean fertility: mean offspring produced by each individual per timestep.
#' @param fecundSd The standard deviation of the fertility.
#' @return A data frame with five variables: timestep, newborns (new individuals
#'   added this timestep), survivors (individuals alive last year who survived
#'   this timestep), population (total individuals alive), and growth (the
#'   increase or decrease in population size from last year).
#' @examples
#' series1 <- timeseries(start = 20, timesteps = 10, survPhi = 0.7, fecundPhi = -0.1, survMean = 0.6,
#' survSd = 0.52, fecundMean = 1.2, fecundSd = 0.7)
#' head(series1)
#' @export
timeseries <- function(start, timesteps, survPhi, fecundPhi, survMean, survSd, fecundMean, fecundSd) {
    # These lines generate the temporally autocorrelated random numbers ----
    delta <- qlogis(survMean) * (1 - survPhi)
    variance <- qlogis(survSd)^2 * (1 - survPhi^2)
    x <- c(rnorm(1, qlogis(survMean), qlogis(survSd)))
    for (i in (1:(timesteps - 1))) {
        x[i + 1] <- delta + survPhi * x[i] + rnorm(1, 0, sqrt(variance))
    }
    delta <- log(1/fecundMean) * (1 - fecundPhi)
    variance <- abs(log(fecundSd))^2 * (1 - fecundPhi^2)
    y <- c(rnorm(1, log(1/fecundMean), abs(log(fecundSd))))
    for (i in (1:(timesteps - 1))) {
        y[i + 1] <- delta + fecundPhi * y[i] + rnorm(1, 0, sqrt(variance))
    }
    # Transforming the random numbers into probabilities of fecundity and survival
    St <- plogis(x)
    Ft <- exp(y)
    # This while loop kills some individuals according to St probabilities, creates
    # new ones according to Ft probabilities, calculates population size and growth,
    # and steps time until 'timesteps' has been reached.
    count <- 0
    population <- vector(mode = "integer", length = timesteps + 1)
    population[1] <- start
    survivors <- vector(mode = "integer", length = timesteps)
    newborns <- vector(mode = "integer", length = timesteps)
    growth <- vector(mode = "integer", length = timesteps)
    while (count < timesteps) {
        count <- count + 1
        survivors[count] <- rbinom(n = 1, size = population[count], prob = St[count])
        newborns[count] <- sum(rpois(n = survivors[count], lambda = Ft[count]))
        growth[count] <- -(population[count] - survivors[count]) + newborns[count]
        population[count + 1] <- population[count] + growth[count]
    }
    # Getting rid of the starting population value
    population <- population[2:length(population)]
    timestep <- 1:timesteps
    data <- data.frame(timestep, newborns, survivors, population, growth)
    return(data)
}
#' Simulate Temporally Autocorrelated Populations for Every Combination of
#' Parameters
#'
#' Essentially a loop of \code{\link{timeseries}}, this function simulates a
#' population with temporally autocorrelated vital rates for every combination
#' of parameters you specify, with as many replicates as desired. It also
#' estimates the sample mean survival and fertility for each simulated
#' population. Please be advised that this function can be very computationally
#' intensive if you provide many possible parameter values and/or ask for many
#' replicates.
#' @import purrr
#' @importFrom dplyr mutate
#' @param timesteps The number of timesteps you want to simulate. Individuals
#'   are added and killed off every timestep according to the survival and
#'   fertility rates. Can be a scalar or a vector of values to loop over.
#' @param start The starting population size. Can be a scalar or vector.
#' @param survPhi The temporal autocorrelation of survival. 0 is white noise (uncorrelated),
#'   positive values are red noise (directly correlated) and negative values are
#'   blue noise (inversely correlated). Can be a scalar or a vector.
#' @param fecundPhi The temporal autocorrelation of fecundity. As above.
#' @param survMean The mean survival from timestep to timestep. Must be a value
#'   between 0 (all individuals die) and 1 (all individuals live). Can be a scalar
#'   or a vector.
#' @param survSd The standard deviation of the survival from timestep to
#'   timestep. Must be a value between 0 and 1. Can be a scalar or a vector.
#' @param fecundMean The mean fertility: mean offspring produced by each individual per timestep. Can be a scalar or a vector.
#' @param fecundSd The standard deviation of the fertility. Can be a scalar or a vector of values.
#' @param replicates How many replicates you would like of each possible
#'   combination of parameters.
#' @return A list of data frames, each with fourteen variables: timestep,
#'   newborns (new individuals added this timestep), survivors (individuals
#'   alive last year who survived this timestep), population (total individuals
#'   alive), growth (the increase or decrease in population size from last
#'   year), estimated survival in the timestep, estimated fecundity in the
#'   timestep, and the seven parameters used to generate the simulation.
#' @examples
#' survival_range <- autocorr_sim(timesteps = 30, start = 200, survPhi = 0.3, fecundPhi = 0.1,
#'                                survMean = c(0.2, 0.3, 0.4, 0.5, 0.6), survSd = 0.5,
#'                                fecundMean = 1.6, fecundSd = 0.5, replicates = 50)
#' head(survival_range[[1]])
#' @export
autocorr_sim <- function(timesteps, start, survPhi, fecundPhi, survMean, survSd, fecundMean, fecundSd,
    replicates) {
    . <- NULL
    # Simulates a population for each combination of parameters -----------
    raw_sims <- cross_df(list(start = start, timesteps = timesteps, survPhi = survPhi, fecundPhi = fecundPhi, survMean = survMean,
        survSd = survSd, fecundMean = fecundMean, fecundSd = fecundSd)) %>% rerun(.n = replicates,
        pmap(., timeseries))
    # Labels each simulation with the parameters that generated it ------
    labeled_sims <- vector(mode = "list", replicates)
    for (i in 1:length(raw_sims)) {
        labeled_sims[[i]] <- vector(mode = "list", nrow(raw_sims[[1]][[1]]))
        for (j in 1:length(raw_sims[[i]][[2]])) {
            labeled_sims[[i]][[j]] <- cbind(raw_sims[[i]][[2]][[j]], raw_sims[[i]][[1]][j,])
        }
    }
    # Unnests the list and adds estimates of survival and fertility
    sims <- labeled_sims %>% flatten() %>%
      map(~mutate(., est_surv = survivors/(population - growth),
                  est_fecund = newborns/survivors))
    return(sims)
}
