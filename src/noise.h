#ifndef NOISE_H
#define NOISE_H

#include <Rcpp.h>
Rcpp::NumericVector raw_noise(int timesteps, double mu, double sigma, double phi);

#endif
