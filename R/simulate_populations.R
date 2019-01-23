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
#' includes environmental stochasticity with colored noise. Density dependence and demographic stochasticity not currently supported.
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
#' @param covMatrix Optional: Add a covariance matrix describing within-year covariances between matrix elements. The matrix elements must be
#' in the same order as they are in the data frame format above: a Leslie matrix turned into a vector row-wise. There should be as many
#' columns as matrix elements, excluding repeat elements (see below) or structural zeros.
#' @param colNames Optional: If the mean, sd, and autocorrelation columns of your data frame input are not
#' named 'mean', 'sd', and 'autocorrelation', provide their names here in a character vector, e.g.,
#' `c(mean = 'Mean', sd = 'Standard Deviation', autocorrelation = 'phi')`
#' @param matrixStructure Optional: By default, the function assumes that the first row of the matrix gives fecundities while
#' the rest of the matrix gives transition or survival probabilities. However, these assumptions do not apply to
#' many plant matrices. If your matrix has transition probabilities in the first row or fecundities beyond the first row
#' (e.g., clonal reproduction), provide a character matrix here with the same dimensions as your matrix that gives in
#' strings whether each element is 'fecundity' or 'transition'.
#' @param repeatElements Optional: Sometimes not all matrix elements can be measured, and some transitions or fertilities
#' are generalized across classes. If you have any matrix elements that are copies of other matrix elements (e.g., stage 3
#' is assumed to have the same fertility as stage 4) indicate them here with a matrix of \emph{rowwise} (not column-wise)
#' indices that show which elements are repeats and which are unique. For example in a 2x2 matrix where both classes are
#' assumed to have the same fertility, input `matrix(c(1, 1, 3, 4), byrow = T, ncol = 2)`. If you indicate repeat elements
#' and you include a covariance matrix, the covariance matrix must only have as many columns as \emph{unique matrix elements}.
#' Structural zeros should \emph{not} be included here as repeats, as they are automatically detected in the function.
#' @param survivalOverflow If the survival for a stage is very high or very variable, the function may sometimes generate
#' projection matrices with survival that exceeds 1 for that stage. The function has two methods of dealing with this problem:
#' either discard all projection matrices and generate new ones until the survival falls within acceptable bounds ("redraw") or
#' divide all the non-fertility matrix elements for that stage by the survival such that they add to 1 ("scale"). The default
#' is "scale".
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
matrix_model <- function(data, initialPop, timesteps, covMatrix = NULL,
    colNames = NULL, matrixStructure = NULL, repeatElements = NULL,
    survivalOverflow = "scale") {
    stages <- length(initialPop)
    # Regularize all valid data inputs to the same format
    if (is.data.frame(data) == T) {
        if (is.null(colNames) == F) {
            data <- data[, colNames] %>% rename(!(!(!colNames)))
        }
        if (all(names(data) == c("mean", "sd", "autocorrelation")) ==
            F) {
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
        if (length(unique(map_int(data, length)))>1) {
          stop("Matrices are not equal dimensions")
        }
        dat <- tibble(mean = as.vector(t(data[[1]])), sd = as.vector(t(data[[2]])),
            autocorrelation = as.vector(t(data[[3]])))
    } else {
        stop("Invalid data type. Must be a list of three matrices or a data frame with three columns.")
    }
    # Generate neutral / null versions of matrixStructure, repeatElements, and
    # covMatrix if not specified in function input
    if (is.null(matrixStructure) == T) {
      matrixStructure <- matrix(
        c(rep("fecundity", stages), rep("transition", stages^2-stages)),
        byrow = T, ncol = stages)
    }
    dists <- ifelse(as.vector(t(matrixStructure)) == "fecundity",
                    "log", "qlogis")
    if (is.null(repeatElements) == T) {
      repeatElements <- matrix(seq(1:stages^2), ncol = stages, byrow = T)
    }
    dat <- dat %>% mutate(dist = dists, zero = mean==0&sd==0,
                          ref = as.vector(t(repeatElements)))
    repeats <- repeatElements == matrix(seq(1:stages^2), ncol = stages, byrow = T)
    # Create version of data that can be used to generate colored noise
    inputs <- dat %>% slice(which(t(repeats))) %>% filter(zero == F) %>%
      rowwise() %>% mutate(
        mean.trans = ifelse(mean == 0, 0, invoke(dist, list(mean))),
        sd.trans = ifelse(sd == 0, 0, variancefix(mean, sd, dist))
      )
    if (is.null(covMatrix) == T) {
      covMatrix <- cor2cov(inputs$sd, diag(nrow(inputs)))
    }
    # Create colored noise, discard if invalid matrix
    if(survivalOverflow == "redraw") {
      repeat {
        inputs$noise <- colored_multi_rnorm(timesteps, inputs$mean.trans, inputs$sd.trans,
                                            inputs$autocorrelation, covMatrix) %>% split(col(.))
        result <- left_join(dat, inputs, by = c("mean", "sd", "autocorrelation", "dist", "zero", "ref"))
        result$noise[dat$zero==T] <- rep(list(rep.int(0, timesteps)), sum(dat$zero==T))
        # checking for >1 probability
        result <- result %>% rowwise() %>% mutate(
          natural.noise = ifelse(zero == T, list(noise),
                                 ifelse(dist == "log", list(exp(noise)),
                                        list(plogis(noise))))
        )
        matrices <- map(1:timesteps, function(x) {
          matrix(map_dbl(result$natural.noise, x), byrow = T, ncol = stages)
        })
        matrices %>% map(replace, which(matrixStructure=="fecundity"), 0) %>%
          map(colSums) %>% map_lgl(function(x) all(x <= 1)) -> check
        if (all(check == T)) {
          break
        }
      }
    } else if(survivalOverflow == "scale") {
      inputs$noise <- colored_multi_rnorm(timesteps, inputs$mean.trans, inputs$sd.trans,
                                          inputs$autocorrelation, covMatrix) %>% split(col(.))
      result <- left_join(dat, inputs, by = c("mean", "sd", "autocorrelation", "dist", "zero", "ref"))
      result$noise[dat$zero==T] <- rep(list(rep.int(0, timesteps)), sum(dat$zero==T))
      # checking for >1 probability
      result <- result %>% rowwise() %>% mutate(
        natural.noise = ifelse(zero == T, list(noise),
                               ifelse(dist == "log", list(exp(noise)),
                                      list(plogis(noise)))))
      matrices <- map(1:timesteps, function(x) {
        matrix(map_dbl(result$natural.noise, x), byrow = T, ncol = stages)
      })
      matrices %>% map(replace, which(matrixStructure=="fecundity"), 0) %>%
        map(colSums) %>% modify_depth(1, function(x) x <= 1) -> check
      if (all(map_lgl(check, all)) == F) {
        matrices <- map(1:length(matrices), function(x) {
          matrices[[x]] %>% replace(which(matrixStructure=="fecundity"), 0) %>%
            t() %>% split(seq(nrow(.))) %>%
            map_if(!check[[x]], function(y) {y/sum(y)}) %>%
            reduce(rbind) %>% t() %>%
            replace(which(matrixStructure=="fecundity"), matrices[[x]][matrixStructure=="fecundity"])
        })
      }
    } else {stop("survivalOverflow must be set to 'redraw' or 'scale'")}
    pop <- projection(initialPop, matrices)
    pop %>% map(as_tibble, .name_repair = ~ c(paste0("stage", 1:stages))) %>%
      bind_rows() %>% group_by(timestep = row_number()) %>% nest() %>%
      mutate(total = map(data, sum)) %>% unnest()
}
