
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The disjointCycles package

<!-- badges: start -->
<!-- badges: end -->

This R package contains source code and replication files for the paper
\`\`Causal Discovery for Linear Non-Gaussian Models with Disjoint
Cycles’’ by Drton, Garrote-Lopez, Nikov, Robeva and Wang. Specifically,
the project seeks to estimate a causal graph from observational data
when the data follow a linear structural equation model. Most notably,
the procedure allows for cycles in the graph. In contrast to previous
work, we test whether certain higher order moments are zero or non-zero
instead of employing an ICA based approach. As a result, our procedure
scales to problems where the number of variables is prohibitively large
for an ICA based approach. As a tradeoff, we require the cycles in the
true graph to be disjoint.

## Installation

You can install the R package from it’s [GitHub
repository](https://github.com/ysamwang/disjointCycles) with the
following code:

``` r
# install.packages("pak")
pak::pak("ysamwang/disjointCycles")
```

## Replicating the numerical results in the paper

The simulations provided in the paper are quite extensive and would
require a high performance cluster to replicate the exact settings used
in the paper. Nonetheless, we provide the code used to generate
numerical results and under settings where the problem size and number
of replications can run on most personal machines.

- To replicate the simulations used to create Figure 4 and 5 use the
  script large_sim.R
- To replicate the simulations in Section 5 regarding the ICA based
  procedure use the script ica_cluster_sim.R. You will also need
  pytetrad which can be accessed at this [github
  repository](https://github.com/cmu-phil/py-tetrad).
