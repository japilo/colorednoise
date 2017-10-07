library(colorednoise)
library(purrr)
library(dplyr)
context("Autocorrelation of timeseries output")

test_that("timeseries can produce blue noise populations", {
  rerun(.n = 1000,
        timeseries(start = 20, timesteps = 10, phi = -0.4, survMean = 0.6,
        survSd = 0.52, fecundMean = 0.3, fecundSd = 0.7)) %>%
    map(~mutate(., est_surv = survivors/(growth - population),
                est_fecund = newborns/survivors)) %>%
    map_dbl(~autocorrelation(.$est_surv)) -> test_blue
  expect_true(mean(test_blue) < -0.1)
})

test_that("timeseries can produce red noise populations", {
  rerun(.n = 1000,
        timeseries(start = 20, timesteps = 10, phi = 0.4, survMean = 0.6,
                   survSd = 0.52, fecundMean = 0.3, fecundSd = 0.7)) %>%
    map(~mutate(., est_surv = survivors/(growth - population),
                est_fecund = newborns/survivors)) %>%
    map_dbl(~autocorrelation(.$est_surv)) -> test_red
  expect_true(mean(test_red) > 0.1)
})
