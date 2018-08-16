library(purrr)
library(dplyr)
context("Consistency in population simulations")

test_that("unstructured_pop can produce blue noise populations", {
    test_blue <- rerun(.n = 1000, unstructured_pop(start = 5000, timesteps = 50,
        survPhi = -0.5, survMean = 0.4, survSd = 0.05, fecundPhi = 0,
        fecundMean = 1.5, fecundSd = 0.2)) %>% map(~mutate(., est_surv = survivors/population,
        est_fecund = newborns/survivors)) %>% map_dbl(~autocorrelation(.$est_surv))
    expect_true(mean(test_blue) < -0.2)
})

test_that("unstructured_pop can produce red noise populations", {
    test_red <- rerun(.n = 1000, unstructured_pop(start = 5000, timesteps = 50,
        survPhi = 0.5, survMean = 0.4, survSd = 0.05, fecundPhi = 0,
        fecundMean = 1.5, fecundSd = 0.2)) %>% map(~mutate(., est_surv = survivors/population,
        est_fecund = newborns/survivors)) %>% map_dbl(~autocorrelation(.$est_surv))
    expect_true(mean(test_red) > 0.2)
})

test_that("matrix_model can produce cross-correlated autocorrelated populations without demographic stochasticity", {
  meanMat <- matrix(c(0.6087097, 0.2480645, 1.6687097, 0.4335484), ncol=2)
  sdMat <- matrix(c(0.0442929, 0.03251947, 0.34437133, 0.10898160), ncol=2)
  phiMat <- matrix(c(-0.20349906,  0.05242292, -0.20349906,  0.02614703), ncol=2)
  covMatrix <- matrix(c(1.0000000, -0.4251527,  0.5000000, -0.4412049,
                        -0.4251527,  1.0000000, -0.4251527, -0.4000000,
                        0.5000000, -0.4251527,  1.0000000, -0.4412049,
                        -0.4412049, -0.4000000, -0.4412049,  1.0000000), byrow=T, ncol = 4) %>%
    cor2cov(sigma = as.vector(sdMat))
  matrixStructure <- matrix(c("transition", "transition", "fecundity", "transition"), ncol = 2)
  test <- matrix_model(list(meanMat, sdMat, phiMat), c(100, 100), 50, covMatrix,
                       matrixStructure = matrixStructure)
  expect_true(sum(test[nrow(test), 2:3]) > sum(test[1, 2:3]))
})

test_that("output of variancefix transforms back to the natural scale", {
  mu <- 0.5
  sigma <- 0.2
  mu.log <- log(mu)
  sigma.log <- variancefix(mu, sigma, "log")
  rand <- rnorm(1000, mu.log, sigma.log) %>% exp()
  expect_true(sd(rand)>=sigma && sd(rand)<0.25)
})
