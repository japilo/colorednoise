#include "noise.h"
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
DataFrame timeseries(int start, int timesteps, double survPhi, double fecundPhi, double survMean, double survSd, double fecundMean, double fecundSd) {
  // These lines generate the temporally autocorrelated random numbers
  double survmu = R::qlogis(survMean, 0, 1, true, false);
  double survsigma = R::qlogis(survSd, 0, 1, true, false);
  NumericVector surv = raw_noise(timesteps, survmu, survsigma, survPhi);
  double fecundmu = log(fecundMean);
  double fecundsigma = abs(log(fecundSd));
  NumericVector fecund = raw_noise(timesteps, fecundmu, fecundsigma, fecundPhi);
  // Transforming the random numbers into fecundity and survival on the response scale
  NumericVector St = plogis(surv);
  NumericVector Ft = exp(fecund);
  // This loop kills some individuals according to St probabilities, creates
  // new ones according to Ft counts, calculates population size and growth,
  // and steps time until 'timesteps' has been reached.
  IntegerVector population(timesteps);
  population[0] = start;
  IntegerVector survivors(timesteps);
  IntegerVector growth(timesteps);
  IntegerVector newborns(timesteps);
  for(int i = 0; i < timesteps; ++i) {
    if (i < timesteps - 1) {
    survivors[i] = R::rbinom(population[i], St[i]);
    newborns[i] = sum(Rcpp::rpois(survivors[i], Ft[i]));
    growth[i] = -(population[i] - survivors[i]) + newborns[i];
    population[i + 1] = population[i] + growth[i];
    } else {
      survivors[i] = R::rbinom(population[i], St[i]);
      newborns[i] = sum(Rcpp::rpois(survivors[i], Ft[i]));
      growth[i] = -(population[i] - survivors[i]) + newborns[i];
    }
  }
 IntegerVector timestep = seq_len(timesteps);
  DataFrame output = DataFrame::create(
    Named("timestep") = timestep,
    _("newborns") = newborns,
    _("survivors") = survivors,
    _("population") = population,
    Named("growth") = growth
  );
  return output;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timeseries(start = 20, timesteps = 10, survPhi = 0.7, fecundPhi = -0.1, survMean = 0.6,
           survSd = 0.52, fecundMean = 1.2, fecundSd = 0.7)
*/
