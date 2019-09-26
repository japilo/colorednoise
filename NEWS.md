## colorednoise 1.0.5

* Updated to be compatible with tidyr v1.0.0

## colorednoise 1.0.4

* Updated to be compatible with tibble v2.0.0

## colorednoise 1.0.3

* `matrix_model` now has more options for dealing with survival values erroneously set to >1.

## colorednoise 1.0.2

* The Makevars have been tweaked so the package runs on more operating systems and compilers.
* `colorednoise` no longer depends on the deprecated package `purrrlyr`.

## colorednoise 1.0.0

* `colorednoise` can now execute matrix population models with the `matrix_model` function, along with its helper functions `colored_multi_rnorm`, `multi_rnorm`, and `cor2cov`.
* Renamed `raw_noise` to `colored_noise` and `timeseries` to `unstructured_pop`.
* The package now has a more extensive testing suite.

## colorednoise 0.0.2

* Recoded the `raw_noise` and `timeseries` functions in Rcpp.
* Changed `timeseries` to take separate autocorrelation values for survival and fertility.
* The output variance of survival probabilities in `timeseries` now reliably matches the input variance.
