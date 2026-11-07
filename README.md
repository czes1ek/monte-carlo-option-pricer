# Monte Carlo Option Pricer

A European call option pricer written in C++, built from scratch. It prices
options via Monte Carlo simulation, validates the result against the
closed-form Black-Scholes solution, and implements antithetic variates as a
variance reduction technique.

## What it does

- Prices a European call option by simulating asset price paths under
  geometric Brownian motion and averaging the discounted payoffs.
- Validates the simulation against the analytical Black-Scholes price with an
  automated tolerance-based test.
- Implements antithetic variates and measures the resulting variance reduction
  across repeated trials.

## Result

Across 50 trials of 100,000 paths each, antithetic variates reduced the
estimator's standard deviation from 0.0447 to 0.0339, a reduction of roughly
24% with no increase in path count.

## Build and run
g++ -std=c++17 -O2 mc.cpp -o mc

./mc

./mc --test           # runs the tests; exits non-zero on failure

./mc --convergence    # prints a CSV convergence table

`--convergence` prints, for 1e3, 1e4, 1e5 and 1e6 paths (seed 12345, S0=100,
K=100, r=0.05, sigma=0.2, T=1), the Monte Carlo price, its standard error, the
Black-Scholes price and the absolute error, with and without antithetic
variates. All rows share the seed, so each run extends the previous one's draws.

## Background

A European call gives the right to buy an asset at a fixed strike price K at
expiration T. Its value is the expected discounted payoff under risk-neutral
dynamics. Monte Carlo estimates this by simulating many possible terminal
prices and averaging the payoffs. Black-Scholes gives the exact value for this
case, which makes it a useful check on the simulation's correctness.