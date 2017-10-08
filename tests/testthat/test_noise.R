library(colorednoise)
library(purrr)
library(dplyr)
context("Autocorrelation of raw noise")

test_that("raw noise can produce blue noise", {
  rerun(.n = 1000,
    raw_noise(timesteps = 100, mu = 0.5, sigma = 0.2, phi = -0.5)) %>%
    map_dbl(autocorrelation) %>% mean() -> test_blue
  expect_true(-0.55 < test_blue && test_blue < -0.45)
})

test_that("raw noise can produce red noise", {
  rerun(.n = 1000,
        raw_noise(timesteps = 100, mu = 0.5, sigma = 0.2, phi = 0.5)) %>%
    map_dbl(autocorrelation) %>% mean() -> test_red
  expect_true(0.55 > test_red && test_red > 0.45)
})
