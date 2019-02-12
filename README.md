<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.com/btaschler/mcap.svg?branch=master)](https://travis-ci.com/btaschler/mcap)

## MCAP

`mcap` provides a model-based clustering approach in very high
dimensions (especially when `p` is much larger than `n`) via adaptive
projections. Clustering is based on full variances Gaussian mixture
modelling in a lower dimensional (projected) space. The projection
dimension is set adaptively in a data-driven manner based on a cluster
stability criterion. Available projection variants (so far) include PCA
and random Projections (Gaussian as well as sparse methods).

### Resources

See our paper: *currently under review*

preprint: …

### Getting Started

Clone or download the [code](https://github.com/btaschler/mcap) from
github.

Alternatively, you can install `mcap` directly from github with:

``` r
# install.packages("devtools")
devtools::install_github("btaschler/mcap")
```

### Prerequisites

Dependencies on other packages:

  - for parallelisation: `foreach`, `doParallel`, `parallel`

  - for clustering: `pcaMethods`, `nethet`, `mclust`, `RandPro`,
    `kernlab`

  - misc: `iterators`, `magrittr`, `stats`, `dplyr`, `tidyverse`,
    `utils, methods`, `data.table`, `RevoUtilsMath`

### Quick demo

This is a basic example showing how to use `mcap` to cluster two (known)
groups:

``` r
library(mcap)

### basic example code
K <- 2       #number of clusters (groups)
n_k <- 200   #number of samples per group
p <- 1000    #number of features (dimension)
A <- matrix(rnorm(n_k*p), n_k, p)            #data for group 1
B <- matrix(rnorm(n_k*p, mean = 1), n_k, p)  #data for group 2
X <- rbind(A, B)                   #input matrix
Y <- c(rep(0, n_k), rep(1, n_k))   #known labels
           
## using PCA projections
model_fit <- MCAPfit(X, k = K, projection = 'PCA', centering_per_group = FALSE,
                     true_labels = Y, parallel = TRUE)

## sparse random projection
model_fit <- MCAPfit(X, k = K, projection = 'li', centering_per_group = FALSE,
                     true_labels = Y, parallel = TRUE)

## adjusted Rand index
print(model_fit$fit_gmm$aRI)

## display assigned cluster labels for each sample
print(model_fit$fit_gmm$model_fit$comp)

## show optimised projection dimension
print(model_fit$fit_q_opt$q_opt)
```

### Versioning

For all available versions, see
[releases](https://github.com/btaschler/mcap/releases). We use [Semantic
Versioning](http://semver.org/).

### Authors

  - **[Bernd Taschler](https://github.com/btaschler), Sach Mukherjee**

List of
[contributors](https://github.com/btaschler/mcap/graphs/contributors).

### License

This project is licensed under the GNU General Public License – see the
[GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) for details.

### Acknowledgments

  - [Konstantinos Perrakis](https://github.com/kperrakis) for valuable
    discussions.

  - The coffee machine for mental and physical support.
