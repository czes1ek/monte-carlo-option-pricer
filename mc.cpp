//Working: Anithetic variates

#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <cstring>

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
                        prices.push_back(monte_carlo_call(S0, K, r, sigma, T, N, static_cast<unsigned int>(i), antithetic));
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

struct Estimate {
      double price;
      double std_error;
};

//Monte Carlo with standard error. Uses the same draws as monte_carlo_call, so its price
//agrees with it to rounding. The standard error comes from independent samples:
//one discounted payoff each, or the average of an antithetic pair.
Estimate monte_carlo_call_with_error(double S0, double K, double r, double sigma, double T, int N,
                                     unsigned int seed, bool antithetic) {
      std::mt19937 gen(seed);
      std::normal_distribution<double> norm(0.0, 1.0);

      const int samples = antithetic ? N/2 : N;
      double sum = 0.0;          // payoffs, accumulated as in monte_carlo_call
      double sample_sum = 0.0;
      double sample_sq_sum = 0.0;

      for (int i = 0; i < samples; i++) {
            double Z = norm(gen);
            double up = std::max(S0 * std::exp((r - 0.5*sigma*sigma)*T + sigma*std::sqrt(T)*Z) - K, 0.0);
            double sample = up;
            sum += up;
            if (antithetic) {
                  double down = std::max(S0 * std::exp((r - 0.5*sigma*sigma)*T +sigma*std::sqrt(T)*(-Z)) - K, 0.0);
                  sum += down;
                  sample = 0.5 * (up + down);
            }
            sample_sum += sample;
            sample_sq_sum += sample * sample;
      }

      double price = std::exp(-r*T) * (sum/N);
      double mean = sample_sum / samples;
      double variance = (sample_sq_sum - samples * mean * mean) / (samples - 1);
      double std_error = std::exp(-r*T) * std::sqrt(variance / samples);

      return {price, std_error};
}

//Prints CSV: price (from monte_carlo_call), its standard error and the absolute error
//against Black-Scholes for 1e3..1e6 paths, plain and antithetic, all with one fixed seed.
void print_convergence() {
      double S0 = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
      unsigned int seed = 12345;
      double bs = black_scholes_call(S0, K, r, sigma, T);

      std::cout << std::setprecision(10);
      std::cout << "variant,paths,seed,S0,K,r,sigma,T,mc_price,std_error,bs_price,abs_error\n";
      for (bool antithetic : {false, true}) {
            for (int N : {1000, 10000, 100000, 1000000}) {
                  double price = monte_carlo_call(S0, K, r, sigma, T, N, seed, antithetic);
                  Estimate e = monte_carlo_call_with_error(S0, K, r, sigma, T, N, seed, antithetic);
                  std::cout << (antithetic ? "antithetic" : "plain") << ',' << N << ',' << seed << ','
                            << S0 << ',' << K << ',' << r << ',' << sigma << ',' << T << ','
                            << price << ',' << e.std_error << ',' << bs << ',' << std::abs(price - bs) << '\n';
            }
      }
}

//Runs the original test, then checks the convergence mode. Returns 0 if everything passes.
int run_tests() {
      run_test();

      double S0 = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
      int N = 1000000;
      unsigned int seed = 12345;
      double tolerance = 0.05;
      double bs = black_scholes_call(S0, K, r, sigma, T);
      int failures = 0;

      for (bool antithetic : {false, true}) {
            const char* name = antithetic ? "antithetic" : "plain";
            double price = monte_carlo_call(S0, K, r, sigma, T, N, seed, antithetic);
            double error = std::abs(price - bs);
            bool ok = error < tolerance;
            std::cout << (ok ? "[Pass] " : "[Fail] ") << "Convergence " << name << " at 1e6 paths: |"
                      << price << " - " << bs << "| = " << error << (ok ? " < " : " >= ") << tolerance << std::endl;
            failures += ok ? 0 : 1;

            //The standard error must describe the same samples as the reported price.
            Estimate e = monte_carlo_call_with_error(S0, K, r, sigma, T, N, seed, antithetic);
            bool same = std::abs(e.price - price) < 1e-9;
            std::cout << (same ? "[Pass] " : "[Fail] ") << "Convergence " << name
                      << " standard error uses the same paths as the price" << std::endl;
            failures += same ? 0 : 1;
      }
      return failures == 0 ? 0 : 1;
}

int main(int argc, char** argv)
{     
      // ./mc                prints the original test and variance-reduction result
      // ./mc --convergence  prints the convergence table as CSV
      // ./mc --test         runs the tests; exits non-zero on failure
      if (argc > 1 && std::strcmp(argv[1], "--convergence") == 0) {
            print_convergence();
            return 0;
      }
      if (argc > 1 && std::strcmp(argv[1], "--test") == 0) {
            return run_tests();
      }

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