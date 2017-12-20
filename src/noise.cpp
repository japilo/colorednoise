#include <Rcpp.h>
using namespace Rcpp;

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

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
raw_noise(timesteps = 30, mu = 0.5, sigma = 0.2, phi = 0.3)
*/
