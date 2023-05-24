#include <Rcpp.h>
using namespace Rcpp;

// Define the Ricker function in C++
  // [[Rcpp::export]]
DataFrame rickerCpp(
  double K, double N0, int timesteps, double distP = 0.05, double env = 0, double env_trend = 0,
  int spinup = 0, int sample_interval = 1,
  double b0_m = 0.01, double b1_m = 0, double b2_m = 0,
  double b0_g = 0, double b1_g = 0, double b2_g = 0,
  double b0_r = 0.2, double b1_r = 0, double b2_r = 0
) {
  NumericVector N(timesteps);  // Create a vector to store population sizes
  
  N[0] = N0;  // Set initial population size
  
  IntegerVector t_vec, N_vec, env_vec, r_vec, actualr_vec, m_vec, g_vec, R_vec, disturbance_vec;
  
  for (int t = 1; t < timesteps; t++) {
    env = R::rnorm(env + env_trend * t, 0.1);
    double g = b0_g + b1_g * env + b2_g * N[t - 1];
    double m = b0_m * R::plogis(b1_m * env + b2_m * N[t - 1]);
    double r = b0_r * exp(b1_r * env + b2_r * N[t - 1] / K);
    
    double G = exp(g * (1 - N[t - 1] / K));
    double M = m * N[t - 1];
    int R = R::rpois(r);
    int disturbance = R::rbinom(1, 1, 1 - distP);
    
    N[t] = (N[t - 1] * G - M) * disturbance + R;
    N[t] = std::max(N[t], 0.0);
    
    t_vec.push_back(t);
    N_vec.push_back(N[t]);
    env_vec.push_back(env);
    r_vec.push_back(r);
    actualr_vec.push_back(N[t - 1] / N[t]);
    m_vec.push_back(m);
    g_vec.push_back(g);
    R_vec.push_back(R);
    disturbance_vec.push_back(disturbance);
  }
  
  int num_samples = (timesteps - 1 - spinup) / sample_interval + 1;
  IntegerVector sample_vec = seq(spinup, timesteps - 1, sample_interval);
  
  DataFrame out_df = DataFrame::create(
    Named("t") = t_vec[sample_vec - 1],
    Named("N") = N_vec[sample_vec - 1],
    Named("env") = env_vec[sample_vec - 1],
    Named("r") = r_vec[sample_vec - 1],
    Named("actualr") = actualr_vec[sample_vec - 1],
    Named("m") = m_vec[sample_vec - 1],
    Named("g") = g_vec[sample_vec - 1],
    Named("R") = R_vec[sample_vec - 1],
    Named("disturbance") = disturbance_vec[sample_vec - 1]
  );
  
  return out_df;
}
