#' Simulate Temporally Autocorrelated Populations for Every Combination of
#' Parameters
#'
#' Essentially a loop of \code{\link{unstructured_pop}}, this function simulates a
#' population with temporally autocorrelated vital rates for every combination
#' of parameters you specify, with as many replicates as desired. It also
#' estimates the sample mean survival and fertility for each simulated
#' population. Please be advised that this function can be very computationally
#' intensive if you provide many possible parameter values and/or ask for many
#' replicates.
#'
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
autocorr_sim <- function(timesteps, start, survPhi, fecundPhi, survMean,
    survSd, fecundMean, fecundSd, replicates) {
    # Simulates a population for each combination of parameters
    # -----------
    raw_sims <- cross_df(list(start = start, timesteps = timesteps,
        survPhi = survPhi, fecundPhi = fecundPhi, survMean = survMean,
        survSd = survSd, fecundMean = fecundMean, fecundSd = fecundSd)) %>%
        rerun(.n = replicates, pmap(., unstructured_pop))
    # Labels each simulation with the parameters that generated it
    # ------
    labeled_sims <- vector(mode = "list", replicates)
    for (i in 1:length(raw_sims)) {
        labeled_sims[[i]] <- vector(mode = "list", nrow(raw_sims[[1]][[1]]))
        for (j in 1:length(raw_sims[[i]][[2]])) {
            labeled_sims[[i]][[j]] <- cbind(raw_sims[[i]][[2]][[j]],
                raw_sims[[i]][[1]][j, ])
        }
    }
    # Unnests the list and adds estimates of survival and fertility
    sims <- labeled_sims %>% flatten() %>% map(~mutate(., est_surv = survivors/population,
        est_fecund = newborns/survivors))
    return(sims)
}

#' Temporally Autocorrelated Matrix Population Models
#'
#' Simulate a structured population with temporal autocorrelation using standard Leslie matrices.
#' Each element in the Leslie matrix has a specified mean, variance, and temporal autocorrelation value.
#' The matrix can have arbitrary dimensions and can have transitions besides linear survival. This model
#' includes environmental stochasticity with colored noise as well as demographic stochasticity. Density
#' dependence is not currently supported.
#'
#' @param data The input data can be one of two formats: a list of three matrices, or a data frame
#' with three columns. \cr
#' If it is a list of three matrices, they must be standard Leslie matrices: the first
#' a matrix of mean values for each matrix element, the second a matrix of standard deviations, and the third
#' a matrix of temporal autocorrelations. \cr
#' If it is a data frame, there must be three columns, one for mean vital rates, one for standard deviations, and one labeled 'autocorrelation.' \cr
#' If the population has n stages, the first n rows of the data frame must be the matrix elements for the first stage,
#' and the next n*(1-n) rows must be the transition probabilities, each row of the matrix from first to last transposed vertically. \cr
#' If you want to run a matrix population model without temporal autocorrelation, simply set all autocorrelation values to zero.
#' @param initialPop An initial population vector. The length must be the same as the number of classes in the matrices.
#' @param timesteps The number of timesteps you would like to simulate the population.
#' @param corrMatrix Optional: add a correlation matrix describing within-year correlations between vital rates. The vital rates must be
#' in the same order as they are in the data frame format above: a Leslie matrix turned into a vector row-wise.
#' @param colNames Optional: If the mean, sd, and autocorrelation columns of your data frame input are not
#' named 'mean', 'sd', and 'autocorrelation', provide their names here in a character vector, e.g.,
#' c(mean = 'Mean', sd = 'Standard Deviation', autocorrelation = 'phi')
#' @param matrixStructure By default, the function assumes that the first row of the matrix gives fecundities while
#' the rest of the matrix gives transition or survival probabilities. However, these assumptions do not apply to
#' many plant matrices. If your matrix has transition probabilities in the first row or fecundities beyond the first row
#' (e.g., clonal reproduction), provide a character matrix here with the same dimensions as your matrix that gives in
#' strings whether each element is 'fecundity' or 'transition'.
#' @param demoStochasticity By default, the matrix model is run without demographic stochasticity - that is,
#' it simply matrix-multiplies the population vector by a stochastic matrix of vital rates. If you would prefer a
#' stochastic projection with demographic stochasticity, that is, the given vital rates are used as rates for
#' probability distributions of survival and fertility, set this argument to TRUE. Please note that the function becomes
#' much more computationally intensive if you turn on demographic stochasticity. It is also important to note that the
#' deomgraphic stochasticity helper function cannot handle population sizes larger than the maximum for 32-bit integers,
#' and will therefore crash if asked to model large populations with high population growth for long time periods.
#' @return A data frame with n + 2 columns, where n is the number of stages in the matrix. One column indicates the timestep,
#' there is one column with the population size for each stage, and one column for total population size.
#' @examples
#' meanMat <- matrix(c(0.55, 0.6, 0.24, 0.4), byrow = TRUE, ncol = 2)
#' sdMat <- matrix(c(0.3, 0.35, 0.05, 0.1), byrow = TRUE, ncol = 2)
#' phiMat <- matrix(c(-0.2, -0.2, 0, 0), byrow = TRUE, ncol = 2)
#' initialPop <- c(100, 100)
#' sim <- matrix_model(list(meanMat, sdMat, phiMat), initialPop, 50)
#' head(sim)
#' @export
matrix_model <- function(data, initialPop, timesteps, corrMatrix = NULL,
    colNames = NULL, matrixStructure = NULL, demoStochasticity = FALSE) {
    stages <- length(initialPop)
    if (is.data.frame(data) == T) {
        if (is.null(colNames) == F) {
            data <- data[, colNames] %>% rename(!!!colNames)
        }
        if (all(names(data) == c("mean", "sd", "autocorrelation")) == F) {
            stop("Please name data frame columns correctly")
        }
        if (stages^2 != nrow(data)) {
            stop("Number of stages in initial population vector does not match number of rows in data. Do you have a row for every matrix element?")
        }
        if (all(data$mean > 0) == F) {
            stop("Invalid values in mean column")
        }
        if (all(data$sd > 0) == F) {
            stop("Invalid values in SD column")
        }
        dat <- data %>% as_tibble()
    } else if (is.list(data) == T) {
        if (length(data) > 3) {
            stop("List data should only have 3 elements")
        }
        if (all(data[[1]] >= 0) == F) {
            stop("Invalid values in mean matrix")
        }
        if (all(data[[2]] >= 0) == F) {
            stop("Invalid values in SD matrix")
        }
        dat <- tibble(mean = as.vector(t(data[[1]])), sd = as.vector(t(data[[2]])),
            autocorrelation = as.vector(t(data[[3]])))
    } else {
        stop("Invalid data type. Must be a list of three matrices or a data frame with three columns.")
    }
    if (is.null(matrixStructure) == T) {
        df <- dat %>% cbind(dist = c(rep("log", stages), rep("qlogis",
            nrow(dat) - stages))) %>% mutate(dist = as.character(dist))
    } else if (all(dim(matrixStructure) == c(stages, stages)) == T) {
        stopifnot(matrixStructure == c("fecundity", "transition"))
        dists <- ifelse(as.vector(t(matrixStructure)) == "fecundity",
            "log", "qlogis")
        df <- dat %>% cbind(dist = dists) %>% mutate(dist = as.character(dist))
    } else {
        stop("Either your initial population vector has an invalid length or your matrixStructure is invalid")
    }
    if (is.null(corrMatrix) == T) {
        elements <- df %>% rowwise() %>% mutate(mean.trans = ifelse(mean ==
            0, 0, invoke(dist, list(mean))), sd.trans = ifelse(sd ==
            0, 0, variancefix(mean, sd, dist)), noise = list(colored_noise(timesteps,
            mean.trans, sd.trans, autocorrelation)), natural.noise = ifelse(all(noise ==
            0) == T, list(noise), ifelse(dist == "log", list(exp(noise)),
            list(plogis(noise)))))
    } else if (is.matrix(corrMatrix) == T) {
        elements <- df %>% rowwise() %>% mutate(mean.trans = invoke(dist,
            list(mean)), sd.trans = variancefix(mean, sd, dist))
        elements$noise <- colored_multi_rnorm(100, elements$mean.trans,
            elements$sd.trans, elements$autocorrelation, corrMatrix) %>%
            split(rep(1:ncol(.), each = nrow(.)))
        elements <- elements %>% mutate(natural.noise = ifelse(dist ==
            "log", list(exp(noise)), list(plogis(noise))))
    } else {
        stop("Correlation matrix must be in matrix format")
    }
    if (demoStochasticity == T) {
        population <- demo_stochasticity(initialPop, elements$natural.noise)
        population %>% as_tibble() %>%
          by_row(., sum, .collate = "cols", .to = "total") %>%
          mutate(timestep = 1:n()) %>% select(timestep, everything()) %>%
          set_names(c("timestep", paste0("stage", 1:stages), "total"))
    } else if (demoStochasticity == F) {
        population <- projection(initialPop, elements$natural.noise)
        population %>% map(as_tibble) %>% bind_rows() %>%
          by_row(., sum, .collate = "cols", .to = "total") %>%
          mutate(timestep = 1:n()) %>% select(timestep, everything()) %>%
          set_names(c("timestep", paste0("stage", 1:stages), "total"))
    }
}
