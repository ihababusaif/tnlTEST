
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tnlTEST

<!-- badges: start -->
<!-- badges: end -->

The goal of tnlTEST is to provide functions to perform the hypothesis
tests for the two sample problem based on order statistics and power
comparisons.

## Installation

You can install the released version of tnlTEST from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("tnlTEST")
```

## Details

A non-parametric two-sample test is performed for testing null
hypothesis *H*<sub>0</sub> : *F* = *G* against the alternative
hypothesis *H*<sub>1</sub> : *F* ≠ *G*. The assumptions of the
*T*<sub>*n*</sub><sup>(ℓ)</sup> test are that both samples should come
from a continuous distribution and the samples should have the same
sample size.<br /> Missing values are silently omitted from x and
y.<br /> Exact and simulated p-values are available for the
*T*<sub>*n*</sub><sup>(ℓ)</sup> test. If exact =“NULL” (the default) the
p-value is computed based on exact distribution when the sample size is
less than 11. Otherwise, p-value is computed based on a Monte Carlo
simulation. If exact =“TRUE”, an exact p-value is computed. If
exact=“FALSE”, a Monte Carlo simulation is performed to compute the
p-value. It is recommended to calculate the p-value by a Monte Carlo
simulation (use exact=“FALSE”), as it takes too long to calculate the
exact p-value when the sample size is greater than 10. <br /> The
probability mass function (pmf), cumulative density function (cdf) and
quantile function of *T*<sub>*n*</sub><sup>(ℓ)</sup> are also available
in this package, and the above-mentioned conditions about exact =“NULL”,
exact =“TRUE” and exact=“FALSE” is also valid for these functions.<br />
Exact distribution of *T*<sub>*n*</sub><sup>(ℓ)</sup> test is also
computed under Lehman alternative.<br /> Random number generator of
*T*<sub>*n*</sub><sup>(ℓ)</sup> test statistic are provided under null
hypothesis in the library.

## Examples

`tnl.test` function performs a nonparametric test for two sample test on
vectors of data.

``` r
library(tnlTEST)
require(stats)
 x=rnorm(7,2,0.5)
 y=rnorm(7,0,1)
 tnl.test(x,y,l=2)
#> $statistic
#> [1] 2
#> 
#> $p.value
#> [1] 0.02447552
```

`ptnl` gives the distribution function of
*T*<sub>*n*</sub><sup>(ℓ)</sup> against the specified quantiles.

``` r
library(tnlTEST)
 ptnl(q=2,n=6,l=2,exact="NULL",trial = 100000)
#> $method
#> [1] "exact"
#> 
#> $cdf
#> [1] 0.03030303
```

`dtnl` gives the density of *T*<sub>*n*</sub><sup>(ℓ)</sup> against the
specified quantiles.

``` r
library(tnlTEST)
 dtnl(k=3,n=7,l=2,exact="TRUE")
#> $method
#> [1] "exact"
#> 
#> $pmf
#> [1] 0.05710956
```

`qtnl` gives the quantile function of *T*<sub>*n*</sub><sup>(ℓ)</sup>
against the specified probabilities.

``` r
library(tnlTEST)
 qtnl(.3,8,1,exact="FALSE",trial = 100000)
#> $method
#> [1] "Monte Carlo simulation"
#> 
#> $quantile
#> [1] 3
```

`rtnl` generates random values from *T*<sub>*n*</sub><sup>(ℓ)</sup>.

``` r
library(tnlTEST)
 rtnl(N=15,n=7,l=2)
#>  [1] 7 6 5 5 3 5 6 6 6 7 6 7 6 6 7
```

`dtnl.lehmann` gives the density of *T*<sub>*n*</sub><sup>(ℓ)</sup>
under Lehmann alternatives.

``` r
library(tnlTEST)
 dtnl.lehmann(l=2,n=6,gamma=0.8)
#> [1] 0.00000000 0.03757203 0.08230829 0.14229514 0.23276690 0.50505764
```

`ptnl.lehmann` gives the distribution function of
*T*<sub>*n*</sub><sup>(ℓ)</sup> under Lehmann alternatives.

``` r
library(tnlTEST)
 ptnl.lehmann(l=2,5,gamma=1.2)
#> [1] 0.0000000 0.0444164 0.1529147 0.3694389 1.0000000
```
