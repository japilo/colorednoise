library(colorednoise)
library(purrr)
library(dplyr)
context("Autocorrelation of colored noise")

test_that("colored noise can produce blue noise", {
    test_blue <- rerun(.n = 1000, colored_noise(timesteps = 100, mu = 0.5,
        sigma = 0.2, phi = -0.5)) %>% map_dbl(autocorrelation) %>%
        mean()
    expect_true(-0.55 < test_blue && test_blue < -0.45)
})

test_that("colored noise can produce red noise", {
    test_red <- rerun(.n = 1000, colored_noise(timesteps = 100, mu = 0.5,
        sigma = 0.2, phi = 0.5)) %>% map_dbl(autocorrelation) %>% mean()
    expect_true(0.55 > test_red && test_red > 0.45)
})

test_that("cor2cov output is correct", {
  corr <- matrix(c(1, 0.53, 0.73, 0.53, 1, 0.44, 0.73, 0.44, 1), nrow = 3)
  sigmas <- c(2, 0.3, 1.2)
  covar <- cor2cov(sigmas, corr)
  expect_true(all.equal(cov2cor(covar), corr))
})
