## colorednoise 1.0.0

* `colorednoise` can now execute matrix population models with the `matrix_model` function, along with its helper functions `colored_multi_rnorm`, `multi_rnorm`, and `cor2cov`.
* Renamed `raw_noise` to `colored_noise` and `timeseries` to `unstructured_pop`.
* The package now has a more extensive testing suite.

## colorednoise 0.0.2

* Recoded the `raw_noise` and `timeseries` functions in Rcpp.
* Changed `timeseries` to take separate autocorrelation values for survival and fertility.
* The output variance of survival probabilities in `timeseries` now reliably matches the input variance.
