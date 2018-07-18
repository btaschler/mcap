<!-- README.md is generated from README.Rmd. Please edit that file -->
mcap
----

The goal of `mcap` is to provide a model-based clustering approach in very high dimensions (especially when *p* ≫ *n*) via adaptive projections. Clustering is based on Gaussian mixture models in a lower dimensional (projected) space. Projection dimension is set adaptively based on a cluster stability criterion. Available projection variants (so far) include PCA and Random Projections.

### Resources

See our paper: ...

### Getting Started

Clone or download the [code](https://github.com/btaschler/MCAP) from github.

Alternatively, you can install `mcap` directly from github with:

``` r
# install.packages("devtools")
devtools::install_github("btaschler/MCAP")
```

### Prerequisites

Dependencies on other packages: - `huge` - `pcaMethods` - `nethet` - ...

    [ ... examples ... ]

### Quick demo

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```

Using the Gram matrix to compute PCA components: Input matrix *X* ∼ MVN(0, 1), output matrix *Z* consists of the first *n**p**c* = 2 principal components.

``` r
library(mcap)
X <- matrix(rnorm(2000),20,100)
Z <- GramPCA(X, npc = 2)
print(Z)
#>               [,1]         [,2]
#>  [1,] -0.112476751  0.135196472
#>  [2,]  0.021942693  0.127673443
#>  [3,]  0.177794912 -0.034702578
#>  [4,] -0.043324750 -0.168622561
#>  [5,] -0.345966418 -0.348263476
#>  [6,] -0.392370608 -0.009594412
#>  [7,]  0.137982361  0.422597804
#>  [8,]  0.343726442 -0.156445340
#>  [9,]  0.165096801  0.023553258
#> [10,]  0.230089922 -0.005554391
#> [11,] -0.299886794  0.116836876
#> [12,]  0.136500259 -0.115854485
#> [13,] -0.115429222 -0.084472274
#> [14,] -0.364219954 -0.105020620
#> [15,]  0.048155690  0.240377440
#> [16,]  0.166415701 -0.153754406
#> [17,]  0.315929543  0.231336978
#> [18,]  0.000618085 -0.030285042
#> [19,] -0.229565212  0.419043866
#> [20,]  0.158987298 -0.504046552
```

### Running the tests

``` r
## basic example code
```

### Versioning

We use [Semantic Versioning](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/btaschler/MCAP/releases).

### Authors

-   **Bernd Taschler, Sach Mukherjee** -- *Initial work* -- [BT](https://github.com/btaschler)

See also the list of [contributors](https://github.com/btaschler/MCAP/graphs/contributors) who participated in this project.

### License

This project is licensed under the GNU General Public License -- see the [LICENSE.md](LICENSE.md) file for details.

### Acknowledgments

-   Konstantinos Perrakis for valuable discussions.
