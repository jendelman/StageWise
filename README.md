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
estimates from Stage 1 should be included in Stage 2 ([Piepho et
al. 2012](https://doi.org/10.1002/bimj.201100219); [Damesa et
al. 2017](https://doi.org/10.2134/agronj2016.07.0395)), which is not
possible with most software packages for genomics-assisted breeding.

R package StageWise was developed to provide a simple interface for
genomic selection and GWAS based on a fully efficient, two-stage
analysis [(Endelman 2023)](https://doi.org/10.1007/s00122-023-04298-x).
Each phenotypic record is associated with a single genotype id, which is
suitable for clonal and inbred lines but not hybrid crops. The ASReml-R
package is used for variance component estimation, which requires a
license from [VSN
International](https://www.vsni.co.uk/software/asreml-r).

Software development has been supported by USDA Hatch Project 1013047
and the USDA National Institute of Food and Agriculture (NIFA) Award
2020-51181-32156. The potato datasets were generated with support from
NIFA Awards 2016-34141-25707 and 2019-34141-30284, Potatoes USA, the
Wisconsin Potato and Vegetable Growers Association, and the University
of Wisconsin-Madison.

To install and load the package:

``` r
install.packages("devtools")
devtools::install_github("jendelman/StageWise", build_vignettes=FALSE)
library(StageWise)
```

There are three vignettes in the package:

-   [Vignette 1](https://jendelman.github.io/StageWise/Vignette1.html)
    illustrates analysis of a single trait with a homogeneous GxE model
    (same genetic correlation between all environments).

-   [Vignette 2](https://jendelman.github.io/StageWise/Vignette2.html)
    illustrates analysis of a single trait with heterogeneous GxE
    covariance, which is often needed with diverse locations.

-   [Vignette 3](https://jendelman.github.io/StageWise/Vignette3.html)
    illustrates the analysis of correlated traits under a homogeneous
    GxE model.

For a complete specification of package functions, consult the
[reference manual.](https://jendelman.github.io/StageWise/manual.pdf)
