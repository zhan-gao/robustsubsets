


# robustsubsets

## Overview

An R implementation of [robust subset
selection](https://arxiv.org/abs/2005.08217).

Robust subset selection is a robust adaption of the classic best subset
selection estimator, and is defined by the constrained least squares
problem:

![](man/figures/README-equation.png)<!-- -->

Robust subsets seeks out the best subset of predictors and observations
and performs a least squares fit on this subset. The number of
predictors used in the fit is controlled by the parameter `k` and the
observations by the parameter `h`.

## Installation

You should install Gurobi and the associated R package gurobi before
installing robustsubsets. Gurobi is available for free under academic
license at <https://www.gurobi.com/>.

To install `robustsubsets` from GitHub, run the following code:

``` r
devtools::install_github('ryan-thompson/robustsubsets')
```

## Usage

The `rss()` function fits a robust subset regression model for a grid of
`k` and `h`. The `cv.rss()` function provides a convenient way to
automatically cross-validate these parameters.

``` r
library(robustsubsets)

# Generate training data with contaminated predictor matrix
set.seed(0)
n <- 100 # Number of observations
p <- 10 # Number of predictors
p0 <- 5 # Number of relevant predictors
ncontam <- 10 # Number of contaminated observations
beta <- c(rep(1, p0), rep(0, p - p0))
x <- matrix(rnorm(n * p), n, p)
e <- rnorm(n, c(rep(10, ncontam), rep(0, n - ncontam)))
y <- x %*% beta + e

# Fit using robust subset selection
fit <- rss(x, y)
coef(fit, k = p0, h = n - ncontam)
```

    ##  [1] -0.1424147  0.7800366  1.0083840  1.0259560  1.0409005  0.9616482
    ##  [7]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000

``` r
# Cross-validate using robust subset selection
cl <- parallel::makeCluster(2)
fit <- cv.rss(x, y, cluster = cl)
parallel::stopCluster(cl)
coef(fit)
```

    ##  [1] -0.1424147  0.7800366  1.0083840  1.0259560  1.0409005  0.9616482
    ##  [7]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000

## Documentation

See the package [reference manual](robustsubsets_1.1.0.pdf).
