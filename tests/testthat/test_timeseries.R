library(purrr)
library(dplyr)
context("Autocorrelation of timeseries output")

test_that("timeseries can produce blue noise populations", {
  rerun(.n = 1000,
        timeseries(start = 200, timesteps = 50, survPhi = -0.5, survMean = 0.4,
        survSd = 0.05, fecundPhi = 0, fecundMean = 1.5, fecundSd = 0.2)) %>%
    map(~mutate(., est_surv = survivors/population,
                est_fecund = newborns/survivors)) %>%
    map_dbl(~autocorrelation(.$est_surv)) -> test_blue
  expect_true(mean(test_blue) < -0.2)
})

test_that("timeseries can produce red noise populations", {
  rerun(.n = 1000,
        timeseries(start = 200, timesteps = 50, survPhi = 0.5, survMean = 0.4,
                   survSd = 0.05, fecundPhi = 0, fecundMean = 1.5, fecundSd = 0.2)) %>%
    map(~mutate(., est_surv = survivors/population,
                est_fecund = newborns/survivors)) %>%
    map_dbl(~autocorrelation(.$est_surv)) -> test_red
  expect_true(mean(test_red) > 0.2)
})
