#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' Generate Autocorrelated Noise
//'
//' Generates temporally autocorrelated random numbers with a mean,
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
//' rednoise <- colored_noise(timesteps = 30, mu = 0.5, sigma = 0.2, phi = 0.3)
//' rednoise
//' @export
// [[Rcpp::export]]
NumericVector colored_noise(int timesteps, double mu, double sigma, double phi) {
  double delta = mu * (1 - phi);
  double variance = pow(sigma, 2.0) * (1 - pow(phi, 2.0));
  NumericVector noise(timesteps);
  noise[0] = R::rnorm(mu, sigma);
  for(int i = 0; i < timesteps-1; ++i) {
    noise[i+1] = delta + phi*noise[i] + R::rnorm(0, sqrt(variance));
  }
  return noise;
}

//' Generate Correlated Normal Random Numbers
//'
//' Generate random numbers from a multivariate normal distribution.
//' It can be used to create correlated random numbers.
//' @param n The number of samples desired for each variable.
//' @param mu A vector giving the mean of each variable.
//' @param sigma A valid covariance matrix.
//' @return A matrix with n rows and as many columns as mu values.
//' @examples
//' mus <- c(0, 3, 5)
//' sigmas <- matrix(c(1, 0.265, 2.19, 0.265, 0.25, 0.66, 2.19, 0.66, 9), ncol = 3)
//' mat <- multi_rnorm(100, mus, sigmas)
//' var(mat)
//' @export
// [[Rcpp::export]]
arma::mat multi_rnorm(int n, NumericVector mu, NumericMatrix sigma) {
  arma::vec mu2 = as<arma::vec>(mu);
  arma::mat sigma2 = as<arma::mat>(sigma);
  int ncols = sigma2.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu2, 1, n).t() + Y * arma::chol(sigma2);
}

//' Convert from Correlation Matrix to Covariance Matrix
//'
//' Convert a correlation matrix to a covariance matrix.
//' @param sigma A vector of standard deviations for the variables you're describing. Length must be the same as the number of rows/columns of CorrMatrix.
//' @param corrMatrix A valid correlation matrix.
//' @return A covariance matrix with the same dimensions as corrMatrix.
//' @examples
//' corr <- matrix(c(1, 0.53, 0.73, 0.53, 1, 0.44, 0.73, 0.44, 1), nrow = 3)
//' sigmas <- c(2, 0.3, 1.2)
//' covar <- cor2cov(sigmas, corr)
//' cov2cor(covar)
//' @export
// [[Rcpp::export]]
arma::mat cor2cov(NumericVector sigma, NumericMatrix corrMatrix) {
  arma::vec sigmas2 = as<arma::vec>(sigma);
  arma::mat corrs = as<arma::mat>(corrMatrix);
  arma::mat m1 = diagmat(sigmas2);
  return m1 * corrs * m1;
}

//' Generate Multiple Cross-Correlated & Autocorrelated Variables
//'
//' Generates random variables that are correlated to each other and temporally autocorrelated.
//' @param timesteps The number of temporally autocorrelated random numbers (one
//'   per timestep) you want.
//' @param mu A vector giving the mean of each variable.
//' @param sigma A vector giving the standard deviation of each variable.
//' @param phi A vector giving the temporal autocorrelation of each variable.
//' @param corrMatrix A valid correlation matrix. The number of rows/columns must match the length of the mu, sigma, and phi vectors.
//' @return A matrix with as many rows as timesteps and as many columns as mu/sigma/phi values.
//' @examples
//' corr <- matrix(c(1, 0.53, 0.73, 0.53, 1, 0.44, 0.73, 0.44, 1), nrow = 3)
//' test <- colored_multi_rnorm(100, c(0, 3, 5), c(1, 0.5, 1), corr, c(0.5, -0.3, 0))
//' var(test)
//' test %>% as.data.frame() %>% summarize_all(.funs = c("mean", "sd", "autocorrelation"))
//' @export
// [[Rcpp::export]]
NumericMatrix colored_multi_rnorm(int timesteps, NumericVector mu, NumericVector sigma, NumericVector phi, NumericMatrix corrMatrix) {
  // Convert mu and sigma to delta and modified SD
  NumericVector sigma2(sigma.length());
  NumericVector delta(mu.length());
  for (int i = 0; i < mu.length(); ++i) {
    delta[i] = mu[i] * (1 - phi[i]);
  }
  for (int i = 0; i < sigma.length(); ++i) {
    sigma2[i] = sqrt(pow(sigma[i], 2.0) * (1 - pow(phi[i], 2.0)));
  }
  // Generate cross-correlated noise with modified SD around zero
  NumericMatrix cov = wrap(cor2cov(sigma2, corrMatrix));
  NumericVector zeroes(mu.length());
  NumericMatrix draws = wrap(multi_rnorm(timesteps, zeroes, cov));
  // initialize colored noise vectors with corresponding mu and sigma
  NumericMatrix noise(timesteps, sigma.length());
  for (int i = 0; i < sigma.length(); ++i) {
    noise(0,i) = R::rnorm(mu[i], sigma[i]);
  }
  // Generate colored noise using draws from correlated noise
  int ncols = noise.ncol();
  int nrows = noise.nrow();
  for (int i = 0; i < ncols; ++i) {
     for (int j = 0; j < nrows; ++j) {
       noise(j+1,i) = delta[i] + phi[i]*noise(j,i) + draws(j,i);
   }
  }
  return noise;
}

// Test code

/*** R
colored_noise(timesteps = 30, mu = 0.5, sigma = 0.2, phi = 0.3)
*/
