library(colorednoise)
library(purrr)
library(dplyr)
context("Autocorrelation of colored noise")

test_that("colored noise can produce blue noise", {
    test_blue <- rerun(.n = 1000, colored_noise(timesteps = 100, mean = 0.5,
        sd = 0.2, phi = -0.5)) %>% map_dbl(autocorrelation) %>%
        mean()
    expect_true(-0.55 < test_blue && test_blue < -0.45)
})

test_that("colored noise can produce red noise", {
    test_red <- rerun(.n = 1000, colored_noise(timesteps = 100, mean = 0.5,
        sd = 0.2, phi = 0.5)) %>% map_dbl(autocorrelation) %>% mean()
    expect_true(0.55 > test_red && test_red > 0.45)
})

test_that("cor2cov output is correct", {
  corr <- matrix(c(1, 0.53, 0.73, 0.53, 1, 0.44, 0.73, 0.44, 1), nrow = 3)
  sigmas <- c(2, 0.3, 1.2)
  covar <- cor2cov(sigmas, corr)
  expect_true(all.equal(cov2cor(covar), corr))
})

test_that("colored_multi_rnorm can produce red noise", {
  set.seed(989)
  corr <- matrix(c(1, 0.53, 0.73, 0.53, 1, 0.44, 0.73, 0.44, 1), nrow = 3)
  test <- colored_multi_rnorm(100, c(0, 3, 5), c(1, 0.5, 1), c(0.5, 0.5, 0.5), corr) %>%
    as_tibble() %>% summarise_all(autocorrelation) %>% as.numeric()
  expect_true(all(test > 0.5)==T)
})

test_that("colored_multi_rnorm can produce blue noise", {
  set.seed(19045)
  corr <- matrix(c(1, 0.53, 0.73, 0.53, 1, 0.44, 0.73, 0.44, 1), nrow = 3)
  test <- colored_multi_rnorm(100, c(0, 3, 5), c(1, 0.5, 1), c(-0.5, -0.5, -0.5), corr) %>%
    as_tibble() %>% summarise_all(autocorrelation) %>% as.numeric()
  expect_true(all(test < -0.4)==T)
})

test_that("multi_rnorm produces correlated variables", {
  set.seed(20)
  mus <- c(0, 3, 5)
  sigmas <- matrix(c(1, 0.265, 2.19, 0.265, 0.25, 0.66, 2.19, 0.66, 9), ncol = 3)
  mat <- multi_rnorm(1000, mus, sigmas)
  expect_true(all.equal(cov(mat), matrix(
    c(0.99893199634091711658, 0.26444969493511205627, 2.19288232024640672435,
      0.26444969493511205627, 0.25975361085119447191, 0.75318406596432485589,
      2.19288232024640672435, 0.75318406596432485589, 8.99581067172942994148),
  byrow = T, ncol=3)))
})
