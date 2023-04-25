
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

## Dara preprocessing

### Input covariance matrix (wide format)

One of the KOMAP algorithm inputs is the covariance matrix of the
clinical features (i.e., log(count+1)) derived from your target
population. To obtain the covariance matrix, we start from the raw EHR
records. Note in the stored fake ehr data, we only include 5 patients
and 2 different ehr codes - lab code and diagnostic ICD code
(DIAG-ICD9).

``` r
library(KOMAP)
head(fake_ehr)
#>   patient_num days_since_admission concept_type concept_code
#> 1          10                    1    DIAG-ICD9          285
#> 2          10                    1    DIAG-ICD9          599
#> 3          10                    1    LAB-LOINC       6690-2
#> 4          10                    1    LAB-LOINC       6690-2
#> 5          10                    1    LAB-LOINC       2532-0
#> 6          10                    1    LAB-LOINC       2160-0
```

We need to first calculate the raw concept count once per day for each
patient(i.e., for those raw codes that show up more than one time in a
day, treat as one time):

``` r
library(dplyr)
fake_ehr_raw_count <- fake_ehr %>% 
  group_by(patient_num, concept_type, concept_code) %>%
  summarise(count_once_per_day = length(unique(days_since_admission)))
head(fake_ehr_raw_count)
#> # A tibble: 6 × 4
#> # Groups:   patient_num, concept_type [1]
#>   patient_num concept_type concept_code count_once_per_day
#>         <int> <chr>        <chr>                     <int>
#> 1          10 DIAG-ICD9    285                           1
#> 2          10 DIAG-ICD9    296                           1
#> 3          10 DIAG-ICD9    394                           1
#> 4          10 DIAG-ICD9    402                           1
#> 5          10 DIAG-ICD9    440                           1
#> 6          10 DIAG-ICD9    480                           1
```

Then for diagnostic ICD code, we need to roll them up to phecode level.
An toy example of the mapping file is stored in the package:

``` r
head(phecode_map)
#>   ICD910 concept_code       Phecode                    Description
#> 1   ICD9          079   PheCode:079                Viral infection
#> 2   ICD9          112   PheCode:112                    Candidiasis
#> 3   ICD9          153   PheCode:153              Colorectal cancer
#> 4   ICD9          153 PheCode:153.2                   Colon cancer
#> 5   ICD9          250   PheCode:250              Diabetes mellitus
#> 6   ICD9          272   PheCode:272 Disorders of lipoid metabolism
```

For each patient, the phecode count is equal to the sum of all its
descendant ICD code once-per-day count:

Then, we combine the rolled up phecode count and raw loinc count
together:

``` r
fake_ehr_count_loinc <- fake_ehr_raw_count %>% filter(concept_type == 'LAB-LOINC') %>%
  mutate(count = count_once_per_day,
         concept_code = paste0('LOINC:', concept_code)) %>%
  select(!count_once_per_day)
fake_ehr_count <- rbind(fake_ehr_count_loinc, fake_ehr_rollup_count_phecode)
head(fake_ehr_count)
#> # A tibble: 6 × 4
#> # Groups:   patient_num, concept_type [1]
#>   patient_num concept_type concept_code count
#>         <int> <chr>        <chr>        <int>
#> 1          10 LAB-LOINC    LOINC:1742-6     1
#> 2          10 LAB-LOINC    LOINC:1920-8     4
#> 3          10 LAB-LOINC    LOINC:1975-2     2
#> 4          10 LAB-LOINC    LOINC:1988-5     3
#> 5          10 LAB-LOINC    LOINC:2160-0     7
#> 6          10 LAB-LOINC    LOINC:2276-4     1
```

Then we transform the data to the wide format, conduct logarithm
transformation, and calculate the covariance matrix:

``` r
fake_ehr_count_wide <- tidyr::pivot_wider(fake_ehr_count, id_cols = c('patient_num'),
                                         names_from = 'concept_code', values_from = 'count')
fake_ehr_count_wide[is.na(fake_ehr_count_wide)] <- 0
fake_ehr_logcount_wide <- fake_ehr_count_wide
fake_ehr_logcount_wide[,-1] <- log(fake_ehr_logcount_wide[, -1] + 1)
fake_ehr_cov <- cov(fake_ehr_logcount_wide[, -1])
fake_ehr_cov[1:3,1:3]
#>              LOINC:1742-6 LOINC:1920-8 LOINC:1975-2
#> LOINC:1742-6   0.41764272  -0.04620953  -0.15587476
#> LOINC:1920-8  -0.04620953   0.11930231   0.08004193
#> LOINC:1975-2  -0.15587476   0.08004193   0.34109280
```

Note that it may be beneficial to screen out patients who never have the
diagnostic code of your target disease before calculating the covariance
matrix. For example, if you are interesting in rheumatoid arthritis
phenotyping task, it is suggested to focus on patients who had at least
1 count of `PheCode:714.1` in their EHR records.

### Input covariance matrix (long format)

You can also input the covariance matrix in a long format. The data must
contain three columns with the first two indicating the node pair names.
Note that when KOMAP tranforms the input to a wide format, it will
assign 0 to missing pairwise covariance value. See an example below:

``` r
head(cov_RA_long)
#> # A tibble: 6 × 3
#>   from  to                   cov
#>   <chr> <chr>              <dbl>
#> 1 utl   utl             1.02    
#> 2 utl   PheCode:714.1   0.319   
#> 3 utl   C0003873        0.359   
#> 4 utl   C0549206        0.122   
#> 5 utl   RXNORM:32624    0.0258  
#> 6 utl   PheCode:286.11 -0.000248
```

### Input feature filters (optional)

To further filter out unrelated codes, you can specify a vector of
codified features as well as a vector of NLP features (e.g., CUIs in the
stored example). You can refer to
[ONCE](https://shiny.parse-health.org/ONCE/), an online searching app
that returns a dictionary of related codified features as well as NLP
features based on knowledge extracted from online article corpus. The
csv file you downloaded from ONCE app has the same format as `codify_RA`
(or `cui_RA`) stored in the KOMAP pacakge.

``` r
codify.feature <- codify_RA$Variable[codify_RA$high_confidence_level == 1]
head(codify.feature)
#> [1] "PheCode:714.1" "PheCode:714"   "RXNORM:614391" "RXNORM:214555"
#> [5] "RXNORM:6851"   "RXNORM:27169"
nlp.feature <- cui_RA$cui[cui_RA$high_confidence_level == 1]
head(nlp.feature)
#> [1] "C0003873" "C0035450" "C0242708" "C0063041" "C0409651" "C0717758"
```

If `codify.feature` or `nlp.feature` argument is not equal to NULL, then
KOMAP will only take a subset of the covariance matrix corresponding to
your input. You can also do the filtering on the input covariance matrix
by yourself.

### Input conditional covariance matrix, conditional mean vector and prevalance (optional)

If you have some labeled data at hand and you want to quickly check
KOMAP performance, this package provides a simulated AUC that only needs
conditional covariance matrix, conditional mean vector and prevalence
without any individual data. To generate such AUC, we start from
`fake_ehr_label_logcount_wide` which was generated by the aforementioned
process using another set of fake patients. The format is very similar
to `fake_ehr_logcount_wide`, except that `fake_ehr_label_logcount_wide`
has an additional column indicating the true label for each patient.

``` r
fake_ehr_label_logcount_wide[1:3,1:5]
#>   Y patient_num LOINC:1742-6 LOINC:1751-7 LOINC:1920-8
#> 1 1          19     1.098612    0.6931472     1.791759
#> 2 0         161     1.609438    0.6931472     1.098612
#> 3 0         179     1.791759    0.0000000     1.098612
```

The conditional statistics are derived based on this label column.

``` r
fake_ehr_var0 = cov(fake_ehr_label_logcount_wide %>% filter(Y == 0) %>% select(-c(Y, patient_num)))
fake_ehr_var1 = cov(fake_ehr_label_logcount_wide %>% filter(Y == 1) %>% select(-c(Y, patient_num)))
fake_ehr_mu0 = colMeans(fake_ehr_label_logcount_wide %>% filter(Y == 0) %>% select(-c(Y, patient_num)))
fake_ehr_mu1 = colMeans(fake_ehr_label_logcount_wide %>% filter(Y == 1) %>% select(-c(Y, patient_num)))
fake_ehr_prev = mean(fake_ehr_label_logcount_wide$Y)
```

Note that the two conditional matrices must have colnames and rownames,
and the conditional mean vectors must have names.

### Input gold label individual data (optional)

If you want to obtain the accurate AUC value by providing individual
count data, the input gold label data should have the format as
`fake_ehr_label_logcount_wide` with one column indicating patient labels
(the column name should go to the `nm.y` argument). Additionally, you
can have another column suggesting the sample probability (the column
name should go to the `nm.pi` argument) if each patient in the label
data was sampled with unequal probability.

``` r
fake_ehr_label_logcount_wide_pi[1:3,1:6]
#>    pi Y patient_num LOINC:1742-6 LOINC:1751-7 LOINC:1920-8
#> 1 0.3 1          19     1.098612    0.6931472     1.791759
#> 2 0.7 0         161     1.609438    0.6931472     1.098612
#> 3 0.4 0         179     1.791759    0.0000000     1.098612
```

### Input individual data for disease score prediction (optional)

If you have some unlabeled data and want to use KOMAP to perform disease
score prediction, the input data should be similar to
`fake_ehr_logcount_wide`. If only part of your individual data has
labels (i.e., some patients have NA value on the `Y` column), KOMAP can
still perform prediction by ignoring the `Y` column. You also need to
assign the unique patient id to a column (i.e, `nm.id` argument) for
your input individual data, so that the output prediction data frame
will also have that column in case of mismatching issue.

## Toy example for how use stored input data

### Stored input data for phenotyping

To illustrate the requirement for input data format, users can refer to
the toy examples stored in the package.

``` r
library(KOMAP)
## Input covariance matrix
?cov_RA

## Input dictionary (can be NULL)
?dict_RA

## Input individual data to predict disease score and/or evaluate model performance (can be NULL)
?dat_part
```

## Phenotyping Example

This is a basic example of how to use the main function in the KOMAP
package to perform model training, score prediction and model evaluation
with built-in toy rheumatoid arthritis data.

### Only output with regression coefficients

``` r
codify.feature <- codify_RA$Variable[codify_RA$high_confidence_level == 1]
nlp.feature <- cui_RA$cui[cui_RA$high_confidence_level == 1]
input.cov <- cov_RA
target.code <- 'PheCode:714.1'
target.cui <- 'C0003873'
nm.utl <- 'utl'
nm.pi <- 'pi'
nm.id <- 'patient_num'
nm.y <- 'Y'
dat.part <- dat_part

## When the input is in a long format:
out_input_long <- KOMAP(cov_RA_long, is.wide = FALSE, target.code, target.cui, nm.utl,
                        nm.multi = NULL, dict_RA, pred = FALSE, eval.real = FALSE, eval.sim = FALSE)
#> 
#> Input long format data, transformed to wide format covariance matrix (42 unique nodes).
#> 
#> Finish estimating coefficients.
```

If you do not specify `pred`, `eval.real` or `eval.sim` equal to be
`TRUE`, then KOMAP only output the regression coefficients which can be
used for disease scoring calculation with additional individual data.
Inside the `est` output, the `long_df` data frame in the output pulls
coefficients from different model settings together if any (e.g., only
use main ICD as surrogate or use both main ICD and main NLP as
surrogates), while `lst` organizes the regression output by model
settings.

``` r
head(out_input_long$est$long_df)
#>          disease               method        target          feat
#> 1: PheCode:714.1 mainICD + allfeature PheCode:714.1 PheCode:714.1
#> 2: PheCode:714.1 mainICD + allfeature PheCode:714.1           utl
#> 3: PheCode:714.1 mainICD + allfeature PheCode:714.1      C0003873
#> 4: PheCode:714.1 mainICD + allfeature PheCode:714.1      C1959609
#> 5: PheCode:714.1 mainICD + allfeature PheCode:714.1      C0919386
#> 6: PheCode:714.1 mainICD + allfeature PheCode:714.1      C0409651
#>                                 desc       coeff
#> 1:              rheumatoid arthritis  0.67831962
#> 2:                Healthcare Utility -0.21875372
#> 3:              Rheumatoid Arthritis  0.15097257
#> 4:                    Erosion lesion  0.04129051
#> 5:               Pathology procedure -0.03343481
#> 6: Seropositive rheumatoid arthritis  0.03003230
```

Similarly, you can input the covariance matrix in wide format:

``` r
out_0 <- KOMAP(input.cov, is.wide = TRUE, target.code, target.cui, nm.utl, nm.multi = NULL, dict_RA,
             pred = FALSE, eval.real = FALSE, eval.sim = FALSE)
#> 
#> Finish estimating coefficients.
```

You can also specify a list of codified/narrative features that are
relevant to the target disease. KOMAP will perform feature screening
before regression:

``` r
out_1 <- KOMAP(input.cov, is.wide = TRUE, target.code, target.cui, nm.utl, nm.multi = NULL, dict_RA,
             codify.feature, nlp.feature,               
             pred = FALSE, eval.real = FALSE, eval.sim = FALSE)
#> Check feature format in `input.cov`, `codify.feature` and/or `cuisearch.feature`...
#> Num of total feat: 184
#> Num of selected codify feat: 48
#> Num of selected codify feat after intersection: 37
#> Num of selected NLP feat: 71
#> Num of selected NLP feat after intersection: 48
#> 
#> More detailed info...
#> Num of PheCode in input.cov: 41
#> Num of PheCode in codify.feature: 19
#> Num of CCS in input.cov: 4
#> Num of CCS in codify.feature: 1
#> Num of LOINC in input.cov: 17
#> Num of LOINC in codify.feature: 6
#> Num of RXNORM in input.cov: 40
#> Num of RXNORM in codify.feature: 18
#> Num of CUI in input.cov: 81
#> Num of CUI in cuisearch.feature: 71
#> 
#> Finish estimating coefficients.
```

### Output regression coefficients and simulated AUC

When specifying `eval.sim=TRUE` and input conditional mean vectors,
conditional covaraince matrices and the prevalence value, KOMAP will
return not only the regression coefficients, but also a simulated AUC
that helps a grasp of model performance without any individual labeled
data.

``` r
out_2 <- KOMAP(input.cov, is.wide = TRUE, target.code, target.cui, nm.utl, nm.multi = NULL, dict_RA,
             codify.feature, nlp.feature,                            
             pred = FALSE, eval.real = FALSE, eval.sim = TRUE,
             mu0, mu1, var0, var1, prev_Y, B = 10000)
#> Check feature format in `input.cov`, `codify.feature` and/or `cuisearch.feature`...
#> Num of total feat: 184
#> Num of selected codify feat: 48
#> Num of selected codify feat after intersection: 37
#> Num of selected NLP feat: 71
#> Num of selected NLP feat after intersection: 48
#> 
#> More detailed info...
#> Num of PheCode in input.cov: 41
#> Num of PheCode in codify.feature: 19
#> Num of CCS in input.cov: 4
#> Num of CCS in codify.feature: 1
#> Num of LOINC in input.cov: 17
#> Num of LOINC in codify.feature: 6
#> Num of RXNORM in input.cov: 40
#> Num of RXNORM in codify.feature: 18
#> Num of CUI in input.cov: 81
#> Num of CUI in cuisearch.feature: 71
#> 
#> Finish estimating coefficients.
#> Finish estimating AUC.
out_2$sim_eval
#>                      method       auc
#> 1          mainICD + codify 0.9499877
#> 2 mainICDNLP + codify & NLP 0.9690431
```

### Output regression coefficients and predicted disease scores

When specifying `pred=TRUE` and input a group of individual data (i.e.,
log count of EHR features), KOMAP will return both the regression
coefficients, and predicted disease scores for the input patients.
Remember to specify the `nm.id` argument (unique patient id) to avoid
score matcing issue.

``` r
out_3 <- KOMAP(input.cov, is.wide = TRUE, target.code, target.cui, nm.utl, nm.multi = NULL, dict_RA,
             codify.feature, nlp.feature,                          
             pred = TRUE, eval.real = FALSE, eval.sim = FALSE,
             dat.part = dat.part, nm.id = nm.id)
#> Check feature format in `input.cov`, `codify.feature` and/or `cuisearch.feature`...
#> Num of total feat: 184
#> Num of selected codify feat: 48
#> Num of selected codify feat after intersection: 37
#> Num of selected NLP feat: 71
#> Num of selected NLP feat after intersection: 48
#> 
#> More detailed info...
#> Num of PheCode in input.cov: 41
#> Num of PheCode in codify.feature: 19
#> Num of CCS in input.cov: 4
#> Num of CCS in codify.feature: 1
#> Num of LOINC in input.cov: 17
#> Num of LOINC in codify.feature: 6
#> Num of RXNORM in input.cov: 40
#> Num of RXNORM in codify.feature: 18
#> Num of CUI in input.cov: 81
#> Num of CUI in cuisearch.feature: 71
#> 
#> Finish estimating coefficients.
#> Finish predicting scores.
head(out_3$pred_prob)
#>    patient_num mainICD + codify mainICDNLP + codify & NLP
#> 3           s1        0.4132152                -1.1677991
#> 5           s2       -0.2087015                 0.4529423
#> 6           s3       -1.2699338                -2.6168660
#> 10          s4       -0.6606499                -1.4898302
#> 13          s5        0.2164402                 0.2210345
#> 20          s6       -0.7194731                -0.3074381
```

### Output regression coefficients and true AUC

When specifying `eval.real=TRUE` and input a group of individual data
(i.e., log count of EHR features) as well as their labels, KOMAP will
return both the regression coefficients, and real AUC evaluated based on
your gold label data. Remember to specify the `nm.y` argument (gold
label for each patient). It will also aumatically return the predicted
disease score for each patient.

``` r
out_4 <- KOMAP(input.cov, is.wide = TRUE, target.code, target.cui, nm.utl, nm.multi = NULL, dict_RA,
             codify.feature, nlp.feature,                           
             pred = FALSE, eval.real = TRUE, eval.sim = FALSE,
             dat.part = dat.part, nm.id = nm.id, nm.pi = nm.pi, nm.y = nm.y)
#> Check feature format in `input.cov`, `codify.feature` and/or `cuisearch.feature`...
#> Num of total feat: 184
#> Num of selected codify feat: 48
#> Num of selected codify feat after intersection: 37
#> Num of selected NLP feat: 71
#> Num of selected NLP feat after intersection: 48
#> 
#> More detailed info...
#> Num of PheCode in input.cov: 41
#> Num of PheCode in codify.feature: 19
#> Num of CCS in input.cov: 4
#> Num of CCS in codify.feature: 1
#> Num of LOINC in input.cov: 17
#> Num of LOINC in codify.feature: 6
#> Num of RXNORM in input.cov: 40
#> Num of RXNORM in codify.feature: 18
#> Num of CUI in input.cov: 81
#> Num of CUI in cuisearch.feature: 71
#> 
#> Finish estimating coefficients.
#> Finish evaluating model prediction.
out_4$real_eval
#>                      method       auc F_score_max
#> 1          mainICD + codify 0.9370693   0.8762717
#> 2 mainICDNLP + codify & NLP 0.9588487   0.8904735
```
