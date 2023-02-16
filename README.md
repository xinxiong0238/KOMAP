
<!-- README.md is generated from README.Rmd. Please edit that file -->

# KOMAP

<!-- badges: start -->
<!-- badges: end -->

The goal of KOMAP is to fit an unsupervised phenotyping system which
does not require manual chart review to acquire gold standard labels for
algorithm training and validation. The online algorithm only requires
summary-level covariance matrix for coefficient training. To get a quick
model evaluation without individual-level data, users can upload
conditional mean vector and covariance matrix to obtain a simulated
model AUC. KOMAP can deal with cases when both codified and NLP features
are available, or when only one feature type is at hand (just leave
`target.code` or `target.cui` equal to `NULL` if not available).

## Installation

You can install the development version of KOMAP from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xinxiong0238/KOMAP")
```

## Toy Example for Input data

To illustrate the requirement for input data format, users can refer to
the toy examples stored in the package.

``` r
## Input covariance matrix
?cov_RA
#> No documentation for 'cov_RA' in specified packages and libraries:
#> you could try '??cov_RA'

## Input dictionary (can be NULL)
?dict_RA
#> No documentation for 'dict_RA' in specified packages and libraries:
#> you could try '??dict_RA'

## Input individual data to predict disease score (can be NULL)
?dat_part
#> No documentation for 'dat_part' in specified packages and libraries:
#> you could try '??dat_part'

## Input gold label to generate real AUC (can be NULL)
?gold_label
#> No documentation for 'gold_label' in specified packages and libraries:
#> you could try '??gold_label'
```

## Phenotyping Example

This is a basic example of how to use the main function in the KOMAP
package to perform model training, score prediction and model evaluation
with built-in toy rheumatoid arthritis data. To do phenotyping on other
diseases, you should specify the related codified features and/or NLP
features either with prior knowledge or with ONCE searching results
(<https://shiny.parse-health.org/CUIsearch-dev/>). The csv file you
downloaded from ONCE app has the same format as `codify_RA` (or
`cui_RA`) stored in the KOMAP pacakge.

``` r
library(KOMAP)
codify.feature <- codify_RA$Variable[codify_RA$select == 1]
cuisearch.feature <- cui_RA$cui[cui_RA$select == 1]
input.cov <- cov_RA
target.code <- 'PheCode:714.1'
target.cui <- 'C0003873'
nm.utl <- 'utl'
nm.pi <- 'pi'
nm.id <- 'patient_num'
nm.y <- 'Y'
dat.part <- dat_part
gold.label <- gold_label

## Only fit the model without any validation and without feature screening
out_0 <- KOMAP(input.cov, target.code, target.cui, nm.utl, nm.multi = NULL, dict_RA,
             pred = FALSE, eval.real = FALSE, eval.sim = FALSE)
#> 
#> Finish estimating coefficients.

## Only fit the model without any validation
out_1 <- KOMAP(input.cov, target.code, target.cui, nm.utl, nm.multi = NULL, dict_RA,
             codify.feature, cuisearch.feature,               
             pred = FALSE, eval.real = FALSE, eval.sim = FALSE)
#> Check feature format in `input.cov`, `codify.feature` and/or `cuisearch.feature`...
#> Num of total feat: 53
#> Num of selected codify feat: 19
#> Num of selected codify feat after intersection: 16
#> Num of selected NLP feat: 20
#> Num of selected NLP feat after intersection: 16
#> 
#> More detailed info...
#> Num of PheCode in input.cov: 8
#> Num of PheCode in codify.feature: 7
#> Num of CCS in input.cov: 0
#> Num of CCS in codify.feature: 0
#> Num of LOINC in input.cov: 2
#> Num of LOINC in codify.feature: 0
#> Num of RXNORM in input.cov: 17
#> Num of RXNORM in codify.feature: 12
#> Num of CUI in input.cov: 25
#> Num of CUI in cuisearch.feature: 20
#> 
#> Finish estimating coefficients.

## Fit the model and calculate simulated AUC
out_2 <- KOMAP(input.cov, target.code, target.cui, nm.utl, nm.multi = NULL, dict_RA,
             codify.feature, cuisearch.feature,               
             pred = FALSE, eval.real = FALSE, eval.sim = TRUE,
             mu0, mu1, var0, var1, prev_Y, B = 10000)
#> Check feature format in `input.cov`, `codify.feature` and/or `cuisearch.feature`...
#> Num of total feat: 53
#> Num of selected codify feat: 19
#> Num of selected codify feat after intersection: 16
#> Num of selected NLP feat: 20
#> Num of selected NLP feat after intersection: 16
#> 
#> More detailed info...
#> Num of PheCode in input.cov: 8
#> Num of PheCode in codify.feature: 7
#> Num of CCS in input.cov: 0
#> Num of CCS in codify.feature: 0
#> Num of LOINC in input.cov: 2
#> Num of LOINC in codify.feature: 0
#> Num of RXNORM in input.cov: 17
#> Num of RXNORM in codify.feature: 12
#> Num of CUI in input.cov: 25
#> Num of CUI in cuisearch.feature: 20
#> 
#> Finish estimating coefficients.
#> Finish estimating AUC.

## If individual data is provided, KOMAP can perform disease score prediction
out_3 <- KOMAP(input.cov, target.code, target.cui, nm.utl, nm.multi = NULL, dict_RA,
             codify.feature, cuisearch.feature,               
             pred = TRUE, eval.real = FALSE, eval.sim = FALSE,
             dat.part = dat.part, nm.id = nm.id)
#> Check feature format in `input.cov`, `codify.feature` and/or `cuisearch.feature`...
#> Num of total feat: 53
#> Num of selected codify feat: 19
#> Num of selected codify feat after intersection: 16
#> Num of selected NLP feat: 20
#> Num of selected NLP feat after intersection: 16
#> 
#> More detailed info...
#> Num of PheCode in input.cov: 8
#> Num of PheCode in codify.feature: 7
#> Num of CCS in input.cov: 0
#> Num of CCS in codify.feature: 0
#> Num of LOINC in input.cov: 2
#> Num of LOINC in codify.feature: 0
#> Num of RXNORM in input.cov: 17
#> Num of RXNORM in codify.feature: 12
#> Num of CUI in input.cov: 25
#> Num of CUI in cuisearch.feature: 20
#> 
#> Finish estimating coefficients.
#> Finish predicting scores.

## If individual data and gold label are provided, KOMAP can perform disease score prediction and calculate the true AUC
out_4 <- KOMAP(input.cov, target.code, target.cui, nm.utl, nm.multi = NULL, dict_RA,
             codify.feature, cuisearch.feature,               
             pred = TRUE, eval.real = TRUE, eval.sim = FALSE,
             dat.part = dat.part, nm.id = nm.id, 
             gold.label = gold.label, nm.pi = nm.pi, nm.y = nm.y)
#> Check feature format in `input.cov`, `codify.feature` and/or `cuisearch.feature`...
#> Num of total feat: 53
#> Num of selected codify feat: 19
#> Num of selected codify feat after intersection: 16
#> Num of selected NLP feat: 20
#> Num of selected NLP feat after intersection: 16
#> 
#> More detailed info...
#> Num of PheCode in input.cov: 8
#> Num of PheCode in codify.feature: 7
#> Num of CCS in input.cov: 0
#> Num of CCS in codify.feature: 0
#> Num of LOINC in input.cov: 2
#> Num of LOINC in codify.feature: 0
#> Num of RXNORM in input.cov: 17
#> Num of RXNORM in codify.feature: 12
#> Num of CUI in input.cov: 25
#> Num of CUI in cuisearch.feature: 20
#> 
#> Finish estimating coefficients.
#> Finish predicting scores.
#> Finish evaluating model prediction.
```
