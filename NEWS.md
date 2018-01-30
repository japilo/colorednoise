## colorednoise 0.0.2

* Recoded the `raw_noise` and `timeseries` functions in Rcpp.
* Changed `timeseries` to take separate autocorrelation values for survival and fertility.
* The output variance of survival probabilities in `timeseries` now reliably matches the input variance.
