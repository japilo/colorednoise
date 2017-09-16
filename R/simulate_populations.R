timeseries <- function(start, years, phi, surv.mean, surv.sd, fecund.mean, fecund.sd) {
    require(stats)
    require(boot)
    # These lines generate the temporally autocorrelated random numbers ----
    delta <- surv.mean * (1 - phi)
    variance <- surv.sd^2 * (1 - phi^2)
    x <- c(rnorm(1, surv.mean, surv.sd))
    for (i in (1:(years - 1))) {
        x[i + 1] <- delta + phi * x[i] + rnorm(1, 0, sqrt(variance))
    }
    delta <- fecund.mean * (1 - phi)
    variance <- fecund.sd^2 * (1 - phi^2)
    y <- c(rnorm(1, fecund.mean, fecund.sd))
    for (i in (1:(years - 1))) {
        y[i + 1] <- delta + phi * y[i] + rnorm(1, 0, sqrt(variance))
    }
    # Transforming the random numbers into probabilities of fecundity and survival
    St <- inv.logit(x)
    Ft <- inv.logit(y)
    # This while loop kills some individuals according to St probabilities, creates new
    # ones according to Ft probabilities, calculates population size and growth, and
    # steps time until 'years' has been reached.
    count <- 0
    population <- start
    survivors <- vector()
    fecundity <- vector()
    growth <- vector()
    while (count < years) {
        count <- count + 1
        survivors[count] <- rbinom(n = 1, size = population[count], prob = St[count])
        fecundity[count] <- rbinom(n = 1, size = survivors[count], prob = Ft[count])
        growth[count] <- -(population[count] - survivors[count]) + fecundity[count]
        population[count + 1] <- population[count] + growth[count]
    }
    # Getting rid of the starting population value
    population <- population[2:length(population)]
    year <- 1:years
    data <- data.frame(year, fecundity, survivors, population, growth)
    return(data)
}
autocorr_sim <- function(length, start, phi, surv.mean, surv.sd, fecund.mean, fecund.sd,
    replicates) {
  # Simulates a population for each combination of parameters -----------
    raw_sims <- cross_df(list(start = start, years = length, phi = phi, surv.mean = surv.mean,
        surv.sd = surv.sd, fecund.mean = fecund.mean, fecund.sd = fecund.sd)) %>%
        rerun(.n = replicates, pmap(., timeseries))
  # Labels each simulation with the parameters that generated it ------
    labeled_sims <- vector(mode = "list", replicates)
    for (i in 1:length(raw_sims)) {
        labeled_sims[[i]] <- vector(mode = "list", nrow(raw_sims[[1]][[1]]))
        for (j in 1:length(raw_sims[[i]][[2]])) {
            labeled_sims[[i]][[j]] <- cbind(raw_sims[[i]][[2]][[j]], raw_sims[[i]][[1]][j,
                ])
        }
    }
    # Unnests the list and adds estimates of survival and fecundity
    sims <- labeled_sims %>% flatten() %>% map(~mutate(., est.surv = survivors/(growth -
        population), est.fecund = fecundity/survivors))
    return(sims)
}
autocorr_estim <- function(listdf, groupvars) {
    listdf %>% map(~group_by_at(., groupvars)) %>% map(~summarize(., acf.surv = acf(est.surv,
        plot = F, na.action = na.pass)[[1]][2], acf.fecund = acf(est.fecund, plot = F,
        na.action = na.pass)[[1]][2])) %>% bind_rows()
}
