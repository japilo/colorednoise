#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' Generate Autocorrelated Noise
//'
//' This function generates temporally autocorrelated random numbers with a mean,
//' standard deviation, and autocorrelation you specify.
//' @param timesteps The number of temporally autocorrelated random numbers (one
//'   per timestep) you want.
//' @param mu The mean of the temporally autocorrelated random numbers.
//' @param sigma The standard deviation of the temporally autocorrelated random
//'   numbers.
//' @param phi The temporal autocorrelation. 0 is white noise (uncorrelated),
//'   positive values are red noise (directly correlated) and negative values are blue
//'   noise (inversely correlated).
//' @return A vector of temporally autocorrelated random numbers.
//' @examples
//' rednoise <- raw_noise(timesteps = 30, mu = 0.5, sigma = 0.2, phi = 0.3)
//' rednoise
//' @export
// [[Rcpp::export]]
NumericVector raw_noise(int timesteps, double mu, double sigma, double phi) {
  double delta = mu * (1 - phi);
  double variance = pow(sigma, 2.0) * (1 - pow(phi, 2.0));
  NumericVector noise(timesteps);
  noise[0] = R::rnorm(mu, sigma);
  for(int i = 0; i < timesteps-1; ++i) {
    noise[i+1] = delta + phi*noise[i] + R::rnorm(0, sqrt(variance));
  }
  return noise;
}

// Generate Correlated Normal Random Numbers
arma::mat mvrnorm(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// Convert from Correlation Matrix to Covariance Matrix
arma::mat cor2cov(arma::vec sigmas, arma::vec corrMatrix) {
  arma::mat m1 = diagmat(sigmas);
  return m1 * corrMatrix * m1;
}

// Generate Multiple Cross-Correlated & Autocorrelated Variables
NumericMatrix colored_mvrnorm(int timesteps, NumericVector mu, NumericVector sigma, NumericMatrix corrMatrix, NumericVector phi) {
  NumericVector sigma2(sigma.length());
  // NumericVector delta(mu.length());
  // for (int i = 0; i < mu.length(); ++i) {
    // delta[i] = mu[i] * (1 - phi[i]);
  //}
  for (int i = 0; i < sigma.length(); ++i) {
    sigma2[i] = sqrt(pow(sigma[i], 2.0) * (1 - pow(phi[i], 2.0)));
  }
  NumericMatrix cov = cor2cov(sigma2, corrMatrix);
  NumericVector zeroes(mu.length());
  NumericMatrix draws = mvrnorm(timesteps, zeroes, cov);
  return draws;
  // After verifying output, add autocorrelation pulling from draws
  // NumericMatrix noise(timesteps, sigma.length());
  // for (int i = 0; i < sigma.length(); ++i) {
  //   noise(0,i) = draws(0,i)
  //}
  // for (int i = 0; i < noise.ncols(); ++i) {
  //   for (int j = 0; j < noise.nrows(); ++j) {
  //     noise(j+1,i) = delta[i] + phi[i]*noise(j,i) + draws(j,i)
  // }
  //}
  // return noise;
}

// Test code

/*** R
raw_noise(timesteps = 30, mu = 0.5, sigma = 0.2, phi = 0.3)
*/
