test_that("colored noise can produce blue noise", {
    test_blue <- map(1:1000, ~colored_noise(timesteps = 100, mean = 0.5,
        sd = 0.2, phi = -0.5)) %>% map_dbl(autocorrelation) %>%
        mean()
    expect_true(-0.55 < test_blue && test_blue < -0.45)
})

test_that("colored noise can produce red noise", {
    test_red <- map(1:1000, ~colored_noise(timesteps = 100, mean = 0.5,
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
    as.data.table() %>%
    .[, .(V1_autocorr = autocorrelation(V1), V2_autocorr = autocorrelation(V2), V3_autocorr = autocorrelation(V3))] %>%
    as.numeric()
  expect_true(all(test > 0.5)==T)
})

test_that("colored_multi_rnorm can produce blue noise", {
  set.seed(19045)
  corr <- matrix(c(1, 0.53, 0.73, 0.53, 1, 0.44, 0.73, 0.44, 1), nrow = 3)
  test <- colored_multi_rnorm(100, c(0, 3, 5), c(1, 0.5, 1), c(-0.5, -0.5, -0.5), corr) %>%
    as.data.table() %>%
    .[, .(V1_autocorr = autocorrelation(V1), V2_autocorr = autocorrelation(V2), V3_autocorr = autocorrelation(V3))] %>%
    as.numeric()
  expect_true(all(test < -0.4)==T)
})

test_that("Bias correction increases estimated autocorrelation", {
  set.seed(180)
  test_bias <- map(1:1000, ~colored_noise(timesteps = 100, mean = 0.5,
        sd = 0.2, phi = 0.5)) %>% map_dbl(autocorrelation) %>% mean()
  test_nobias <- map(1:1000, ~colored_noise(timesteps = 100, mean = 0.5,
        sd = 0.2, phi = 0.5)) %>% map_dbl(autocorrelation, biasCorrection = F) %>% 
        mean()
  expect_true(test_bias > test_nobias)
})