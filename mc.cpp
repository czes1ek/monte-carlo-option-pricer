//Working: Anithetic variates

#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <vector>

double normal_cdf(double x) {
      return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

//Monte Carlo Function
double monte_carlo_call(double S0, double K, double r, double sigma, double T, int N,
                        unsigned int seed, bool antithetic) {
      std::mt19937 gen(seed);
      std::normal_distribution<double> norm(0.0, 1.0);

      double sum = 0.0;

      if (antithetic) {
            for (int i = 0; i < N/2; i++) {
                  double Z = norm(gen);

                  double S_T_up     = S0 * std::exp((r - 0.5*sigma*sigma)*T + sigma*std::sqrt(T)*Z);
                  double S_T_down   = S0 * std::exp((r - 0.5*sigma*sigma)*T +sigma*std::sqrt(T)*(-Z));

                  sum += std::max(S_T_up - K, 0.0);
                  sum += std::max(S_T_down - K, 0.0);
            }
      } else {
            for (int i = 0; i < N; i++) {
                  double S_T;
                  double Z = norm(gen);

                  S_T = S0 * std::exp((r - 0.5*sigma*sigma)*T + sigma*std::sqrt(T)*Z);
                  sum += std::max(S_T - K, 0.0);    
            }
      }

      sum = std::exp(-r*T) * (sum/N);

      return sum;
}

//Black Scholes Function
double black_scholes_call(double S0, double K, double r, double sigma, double T) {
      double C;
      double d1;
      double d2;

      d1 = ( std::log(S0/K) + (r + 0.5*sigma*sigma)*T ) / ( sigma*std::sqrt(T) );
      d2 = d1 - sigma*std::sqrt(T);

      C = S0 * normal_cdf(d1) - K * std::exp(-r*T) * normal_cdf(d2);

      return C;
}

void run_test () {
            double S0 = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
            int N = 1000000;

            double mc = monte_carlo_call(S0,  K,  r,  sigma,  T,  N, 12345, true);
            double bs = black_scholes_call(S0, K, r, sigma, T);

            double tolerance = 0.05;

            if (std::abs(mc - bs) < tolerance){
                  std::cout << "[Pass] Monte Carlo matches Black-Scholes... " << mc << " - " << bs << " <= 0.05" << std::endl;
            } else {
                  std:: cout << "[Fail] Monte Carlo does not match Black-Scholes... " << mc << " - " << bs << " > 0.05"<< std::endl;
            }
      }

double measure_std(double S0, double K, double r, double sigma, double T, int N,
                  int num_trials, bool antithetic) {

                  std::vector<double> prices;

                  //collect
                  for (int i = 0; i < num_trials; i++) {
                        prices.push_back(monte_carlo_call(S0, K, r, sigma, T, N, i, antithetic));
                        }

                        double sum = 0.0;
                        for (double p : prices) {
                              sum += p;
                        }

                        double mean = sum / num_trials;

                        double variance = 0.0;
                        for (double p : prices) {
                              variance += (p - mean) * (p - mean);
                        }
                        
                        variance = variance / num_trials;
                  
                  return std::sqrt(variance);

                  }

int main()
{     
      double S0=100, K=100, r=0.05, sigma=0.20, T=1.0;
      int N = 100000;
      int trials = 50;

      run_test();

      double std_plain = measure_std(S0, K, r, sigma, T, N, trials, false);
      double std_anti  = measure_std(S0, K, r, sigma, T, N, trials, true);

      std::cout << "Std dev without antithetic: " << std_plain << std::endl;
      std::cout << "Std dev with antithetic:    " << std_anti  << std::endl;
      
      return 0;
}