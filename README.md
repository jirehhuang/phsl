phsl: Partitioned Hybrid Structure Learning
================
Jireh Huang
(<jirehhuang@ucla.edu>)

This package implements the Bayesian network structure learning
algorithms developed in [Huang and Zhou
(2022)](https://doi.org/10.1007/s10994-022-06145-4), interfacing with
the R package [bnlearn](https://www.bnlearn.com/). In particular, this
package features the partitioned PC (pPC), *p*-value adjacency
thresholding (PATH), and hybrid greedy initialization (HGI) algorithms,
which culminate in the partitioned hybrid greedy search (pHGS)
algorithm.

# Installation

Run the following code to install the package from GitHub.

``` r
devtools::install_github("jirehhuang/phsl")
```

# Usage

The pPC, PATH, and HGI algorithms are integrated into the `bnsl()`
function, which is a mix-and-match Bayesian network structure learning
function modeled after the `rsmax2()` function from the
[bnlearn](https://www.bnlearn.com/) package. The pPC algorithm is
available as a constraint-based algorithm in the `restrict` argument,
PATH may be specified by the `path` argument when `restrict = "ppc"`,
and HGI may be activated by the `hgi` argument for constraint-based and
hybrid approaches. See `help(bnsl)` for more details and examples.
Wrapper functions `ppc()` and `phgs()` contain presets to implement pPC
with PATH and pHGS, respectively.

# Reference

Please cite the following paper when using any part of this package,
modified or as is.

[Huang, J., & Zhou, Q. (2022). Partitioned hybrid learning of Bayesian
network structures. *Machine Learning*.
https://doi.org/10.1007/s10994-022-06145-4](https://doi.org/10.1007/s10994-022-06145-4)
