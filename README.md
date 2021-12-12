R/StageWise
================
Jeffrey Endelman

Plant breeders typically use multi-environment datasets that are
unbalanced with respect to individuals over time and often within year
as well. Experimental designs and the best models for
micro-environmental variation may also vary with environment, defined as
a location x year combination. These complexities can be handled
efficiently with a two-stage analysis, in which best linear unbiased
estimates (BLUEs) of genotypic value are computed for each environment
in Stage 1 and then used as the response variable in Stage 2. To fully
utilize the data, however, the variance-covariance matrix of the
estimates from Stage 1 must be included in Stage 2 ([Piepho et
al. 2012](https://doi.org/10.1002/bimj.201100219); [Damesa et
al. 2017](https://doi.org/10.2134/agronj2016.07.0395)), which is not
possible with most software packages for genomics-assisted breeding.

R package StageWise was developed to provide a simple interface for
genomic selection and GWAS based on a fully efficient, two-stage
analysis. The software assumes genetic effects from a single population,
which is suitable for many crops but not recommended when there are
multiple heterotic groups. The example datasets are for tetraploid
potato, but the software works for diploids and higher polyploids too.
The ASReml-R package is used for variance component estimation, which
requires a license from [VSN
International](https://www.vsni.co.uk/software/asreml-r).

Funding for software development has come from USDA Hatch Project
1013047 and NIFA Award 2020-51181-32156, and the vignette datasets were
generated with support from NIFA Awards 2016-34141-25707 and
2019-34141-30284, Potatoes USA, the Wisconsin Potato and Vegetable
Growers Association, and the University of Wisconsin-Madison.

To install and load the package:

``` r
install.packages("devtools")
devtools::install_github("jendelman/StageWise", build_vignettes=FALSE)
library(StageWise)
```

There are three vignettes in the package:

-   [Vignette 1](https://jendelman.github.io/StageWise/Vignette1.html)
    illustrates analysis of a single trait under a compound symmetry
    model for GxE, which is appropriate when the genetic correlation
    between environments is similar.

-   [Vignette 2](https://jendelman.github.io/StageWise/Vignette2.html)
    illustrates analysis of a single trait under a heterogeneous
    covariance structure, which is often needed with diverse locations.

-   [Vignette 3](https://jendelman.github.io/StageWise/Vignette3.html)
    illustrates the analysis of correlated traits under a compound
    symmetry model for GxE.

For a complete specification of package functions, consult the
[reference manual.](https://jendelman.github.io/StageWise/manual.pdf)
