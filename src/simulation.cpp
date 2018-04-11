#include <RcppArmadillo.h>
#include "noise.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Variance Fix
//
// This function changes the variance so that once it is back-transformed
// from the logit scale, the original variance is recovered.
// [[Rcpp::export]]
double variancefix(double mu, double sigma, std::string dist){
  if (sigma == 0) {
    return 0;
  }
  if (dist == "qlogis") {
  double mn = R::qlogis(mu, 0, 1, true, false);
  double out = (sigma*pow(exp(mn) + 1, 2.0))/exp(mn);
  return out;
  } else if (dist == "log") {
    double mn = log(mu);
    double out = sigma*(1/exp(mn));
    return out;
  } else {
    stop("Inadmissible value");
  }
}

//' Simulated Time Series of an Unstructured Temporally Autocorrelated Population
//'
//' This function simulates an unstructured population with
//' temporally autocorrelated vital rates (survival and fertility). In other
//' words, this function will show you the dynamics over time of a population
//' whose survival and fertility is stochastic, but also correlated to the
//' survival and fertility in the previous year, respectively. The assumptions of the
//' simulation are that the population is asexually reproducing or female-only,
//' survival and fertility are the same at all ages / stages,
//' and that individuals continue to be reproductively capable until they die.
//'
//' Be advised that not all combinations of values will work. If you set survival and
//' fertility unrealistically high, the population size will tend toward infinity and
//' the simulation will fail because the numbers are too large to handle. Use your
//' common sense as a demographer / population biologist.
//' @param start The starting population size.
//' @param timesteps The number of timesteps you want to simulate. Individuals
//'   are added and killed off every timestep according to the survival and
//'   fertility rates. In ecological applications, timesteps are usually years,
//'   but theoretically they can be any length of time.
//' @param survPhi The temporal autocorrelation of survival. 0 is white noise (uncorrelated),
//'   positive values are red noise (directly correlated) and negative values are
//'   blue noise (inversely correlated).
//' @param fecundPhi The temporal autocorrelation of fecundity. As above.
//' @param survMean The mean survival from timestep to timestep. Must be a value
//'   between 0 (all individuals die) and 1 (all individuals live).
//' @param survSd The standard deviation of the survival from timestep to
//'   timestep. Must be a value between 0 and 1.
//' @param fecundMean The mean fertility: mean offspring produced by each individual per timestep.
//' @param fecundSd The standard deviation of the fertility.
//' @return A data frame with four variables: timestep, population (total individuals
//'   alive at the start of the timestep), newborns (new individuals
//'   born this timestep), and survivors (individuals who survive this timestep).
//' @examples
//' series1 <- unstructured_pop(start = 20, timesteps = 10, survPhi = 0.7, fecundPhi = -0.1,
//' survMean = 0.6, survSd = 0.52, fecundMean = 1.2, fecundSd = 0.7)
//' head(series1)
//' @export
// [[Rcpp::export]]
DataFrame unstructured_pop(int start, int timesteps, double survPhi, double fecundPhi, double survMean, double survSd, double fecundMean, double fecundSd) {
  // These lines generate the temporally autocorrelated random numbers
  double survmu = R::qlogis(survMean, 0, 1, true, false);
  double survsigma = variancefix(survMean, survSd, "qlogis");
  Rcpp::NumericVector survnoise = colored_noise(timesteps, survmu, survsigma, survPhi);
  Rcpp::NumericVector St = plogis(survnoise);
  double fecundmu = log(fecundMean);
  double fecundsigma = variancefix(fecundMean, fecundSd, "log");
  Rcpp::NumericVector fecundnoise = colored_noise(timesteps, fecundmu, fecundsigma, fecundPhi);
  Rcpp::NumericVector Ft = exp(fecundnoise);
  // This loop kills some individuals according to St probabilities, creates
  // new ones according to Ft counts, calculates population size and growth,
  // and steps time until 'timesteps' has been reached.
  std::vector<long> population(timesteps);
  population[0] = start;
  std::vector<long> survivors(timesteps);
  std::vector<long> newborns(timesteps);
  for(int i = 0; i < timesteps; ++i) {
    if (i < timesteps - 1) {
    survivors[i] = R::rbinom(population[i], St[i]);
    newborns[i] = sum(Rcpp::rpois(survivors[i], Ft[i]));
    population[i + 1] = survivors[i] + newborns[i];
    } else {
      survivors[i] = R::rbinom(population[i], St[i]);
      newborns[i] = sum(Rcpp::rpois(survivors[i], Ft[i]));
    }
  }
 IntegerVector timestep = seq_len(timesteps);
  DataFrame output = DataFrame::create(
    Named("timestep") = timestep,
    _("newborns") = wrap(newborns),
    _("survivors") = wrap(survivors),
    Named("population") = wrap(population)
  );
  return output;
}

// Matrix Projection Function without Demographic Stochasticity

// Feed in an initial population vector and a list of values for each matrix vital rate
// in each timestep. Get out a list of matrices with the population for each timestep
// [[Rcpp::export]]
Rcpp::List projection(arma::vec initialPop, List noise) {
  List vectors(noise.length());
  for (int i = 0; i < noise.length(); ++i) {
    vectors[i] = as<vec>(noise[i]);
  }
  vec sample = vectors[1];
  int timesteps = sample.size();
  List matrices(timesteps);
  for (int i = 0; i < timesteps; ++i) {
    arma::vec timestep(vectors.length());
    for (int j = 0; j < vectors.length(); ++j) {
      vec v = vectors[j];
      timestep(j) = v[i];
    }
    mat projection = mat(timestep);
    projection.reshape(initialPop.size(), initialPop.size());
    matrices[i] = projection.t();
  }
  List population(timesteps);
  population[0] = initialPop;
  for (int i = 0; i < timesteps-1; ++i) {
    rowvec pop = population[i];
    mat proj = matrices[i];
    population[i + 1] = pop*proj;
  }
  population[0] = initialPop.t();
  return population;
}

// Matrix Projection Function with Demographic Stochasticity

// Feed in an initial population vector and a matrix of vital rates for each year
// Get out a matrix with stage-specific population each year
Rcpp::NumericMatrix demo_stochasticity(arma::Row<long> initialPop, List noise) {
  Rcpp::NumericVector sample = noise[1];
  int timesteps = sample.length();
  int stages = initialPop.size();
  arma::mat rates(timesteps, stages*stages);
  for (int i = 0; i < noise.length(); ++i) {
    arma::vec newcol = noise[i];
    rates.col(i) = newcol;
  }
  arma::Mat<long> population(timesteps, stages);
  arma::Mat<long> elements(timesteps, stages*stages);
  population.row(0) = initialPop;
  for(int i = 0; i < timesteps-1; ++i) {
    for(int j = 0; j < stages; ++j) {
      if (j == 0) {
        for(int k = 0; k < stages; ++k) {
          elements(i, k) = sum(Rcpp::rpois(population(i, k), rates(i, k)));
        }
        arma::Row<long> offspring = elements.row(i);
        population(i+1, j) = sum(offspring.subvec(0, stages-1));
      } else {
        for (int k = 0; k < stages; ++k) {
          elements(i, j+1+k) = R::rbinom(population(i, k), rates(i, stages+k));
        }
        arma::Row<long> recruits = elements.row(i);
        population(i+1,j) = sum(recruits.subvec(stages, stages*stages-1));
      }
    }
  }
  return wrap(population);
}
