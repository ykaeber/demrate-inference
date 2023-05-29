#include <Rcpp.h>
using namespace Rcpp;

// Function to simulate Beverton-Holt model
// [[Rcpp::export]]
DataFrame beverton_holt(
    double N0, int timesteps, int spinup = 0, int sample_interval = 1,
    double b0_e = 0.0, double b1_e = 0.0,
    double b0_k = 100, double b1_k = 0.3,
    double b1_g = 0.0, double b2_g = 0,
    double b0_r = 1, double b1_recr = 0.1, double b2_recr = -0.001, int opt_x1 = 1,
    double distP = 0.01, int Nrep = 1
    ) {
  
  std::vector<int> rep_vec, t_vec, disturbance_vec;
  std::vector<double> N_vec, x1_vec, K_vec, g_vec, recr_vec, recr_mean_vec, actualR_vec;
  
  for(int rep_i = 1; rep_i <= Nrep; rep_i++){
    
    int N_size = timesteps + 1;
    std::vector<double> N(N_size);
    
    // Set initial N size
    N[0] = N0;
    double recr = 1.0;
    
    for (int t = 1; t <= timesteps; t++) {
      double x1 = b0_e + R::rnorm(0, 0.02);
      x1 = std::min(1.0, std::max(0.0, x1));
      
      double K = b0_k + b1_k * x1;
      double g = R::plogis(b1_g * x1 + b2_g * N[t-1], 0, 1, TRUE, FALSE) * 2.0;
      
      int disturbance = R::rbinom(1, 1 - distP);
      
      N[t] = (g * N[t-1]) / (1 + (N[t-1] / K)) * disturbance;
      
      recr = std::max(recr, 1.0);
      double K_r = 10.0;
      double g_recr = R::plogis(b0_r + b1_recr * x1 + b2_recr * N[t], 0, 1, TRUE, FALSE) * 2.0;
      if (g_recr < 0.0001) g_recr = 0.0;
      double recr_mean = (g_recr * recr) / (1 + (recr / K_r));
      
      recr = R::rpois(recr_mean);
      
      N[t] += recr;
      
      rep_vec.push_back(rep_i);
      t_vec.push_back(t);
      N_vec.push_back(N[t]);
      x1_vec.push_back(x1);
      K_vec.push_back(K);
      g_vec.push_back(g);
      recr_mean_vec.push_back(recr_mean);
      recr_vec.push_back(recr);
      disturbance_vec.push_back(disturbance);
      actualR_vec.push_back(N[t-1] / N[t]);
    }
  }
  
  int num_samples = (timesteps - spinup) / sample_interval + 1;
  std::vector<int> sample_vec(num_samples);
  for (int i = 0; i < num_samples; i++) {
    sample_vec[i] = spinup + i * sample_interval;
  }
  
  DataFrame out_df = DataFrame::create(
    Named("rep") = rep_vec,
    Named("t") = t_vec,
    Named("N") = N_vec,
    Named("x1") = x1_vec,
    Named("K") = K_vec,
    Named("g") = g_vec,
    Named("recr_mean") = recr_mean_vec,
    Named("recr") = recr_vec,
    Named("disturbance") = disturbance_vec,
    Named("actualR") = actualR_vec
  );
  
  return out_df;
}
