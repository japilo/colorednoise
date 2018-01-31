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
#'                                fecundMean = 1.1, fecundSd = 0.5, replicates = 50)
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
      map(~mutate(., est_surv = survivors/population,
                  est_fecund = newborns/survivors))
    return(sims)
}

#' Temporally Autocorrelated Matrix Models

matrix_model <- function(mean.mat, sd.mat, phi.mat, pop.init, timesteps){
  mean.mat.trans <- matrix(c(log(mean.mat[1,]), qlogis(mean.mat[2:nrow(mean.mat),])), byrow=T, nrow = nrow(mean.mat))
  sd.mat.trans <- matrix(c(map2_dbl(mean.mat[1,], sd.mat[1,], variancefix, "log"),
                           map2_dbl(mean.mat[2:nrow(mean.mat),], sd.mat[2:nrow(sd.mat),], variancefix, "logis")),
                         byrow=T, nrow = nrow(sd.mat))
  noise.mat = list()
  for (i in 1:length(phi.mat)) {
    noise.mat[[i]] = raw_noise(timesteps, mu = mean.mat.trans[i], sigma = sd.mat.trans[i], phi = phi.mat[i])
  }
  fecundity.index <- seq(1, length(noise.mat), nrow(mean.mat))
  fecundity <- noise.mat[fecundity.index] %>% map(exp)
  transitions <- noise.mat[-fecundity.index] %>% map(plogis)
  yearly.mat <- list()
  for(i in 1:timesteps) {
    yearly.mat[[i]] <- matrix(c(map_dbl(1:length(fecundity), function(x) {fecundity[[x]][i]}),
                                map_dbl(seq(1, length(transitions), nrow(mean.mat)-1), function(x) {transitions[[x]][i]})),
                              byrow = T, nrow = nrow(mean.mat))
  }
  population <- list(matrix(pop.init, ncol = ncol(mean.mat)))
  for (i in 1:(timesteps-1)){
    population[[i + 1]] <- population[[i]] %*% yearly.mat[[i]]
  }
  population %>% map(as.data.frame) %>% bind_rows(.id = "timestep")
}
