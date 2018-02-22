#ifndef NOISE_H
#define NOISE_H

#include <RcppArmadillo.h>
Rcpp::NumericVector colored_noise(int timesteps, double mu, double sigma, double phi);
arma::mat multi_rnorm(int n, Rcpp::NumericVector mu, Rcpp::NumericMatrix sigma);
arma::mat cor2cov(Rcpp::NumericVector sigma, Rcpp::NumericMatrix corrMatrix);
Rcpp::NumericMatrix colored_multi_rnorm(int timesteps, Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector phi, Rcpp::NumericMatrix corrMatrix);

#endif
