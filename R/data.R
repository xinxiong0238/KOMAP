#' A sample covariance matrix for rheumatoid arthritis
#'
#' A covariance matrix of rheumatoid arthritis (RA) related feature counts (after log transformation). Must have colnames and rownames.
#'
#' @docType data
#' @keywords datasets
#' @name cov_RA
#' @usage data(cov_RA)
#' @format A numeric matrix with 53 rows and 53 columns
NULL

#' A long-format sample covariance table for rheumatoid arthritis
#'
#' A long-format covariance data frame of rheumatoid arthritis (RA) related feature counts (after log transformation). Must have three columns with the first two indicating the node name and the third one indicating covariance value.
#'
#' @docType data
#' @keywords datasets
#' @name cov_RA_long
#' @usage data(cov_RA_long)
#' @format A numeric matrix with 1681 rows and 3 columns
NULL


#' Codified features related to rheumatoid arthritis
#'
#' A dataframe containing codified features related to rheumatoid disease (PheCode:714.1) generated by ONCE. The variables are as follows:
#'
#' \itemize{
#'   \item Variable. codified feature id (possible categories: PheCode/RXNORM/CCS/LOINC/Other lab)
#'   \item Description. full feature description
#'   \item relatedness_to_target_disease. Cosine similarity between the feature and PheCode:714.1.
#'   \item high_confidence_level. Selected features with further screening strategy.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name codify_RA
#' @usage data(codify_RA)
#' @format A data frame with 66 rows and 4 variables
NULL



#' NLP features related to rheumatoid arthritis
#'
#' A dataframe containing NLP features related to rheumatoid disease (C0003873) generated by ONCE. The variables are as follows:
#'
#' \itemize{
#'   \item cui. NLP feature id (possible categories: CUI)
#'   \item term. full feature description
#'   \item relatedness_to_target_disease. cosine similarity between the feature and C0003873
#'   \item importance_score. Weight combining cosine similarity, article count and term frequency that shows relatedness between features and C0003873
#'   \item high_confidence_level. Selected features with further screening strategy.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name cui_RA
#' @usage data(cui_RA)
#' @format A data frame with 897 rows and 5 variables
NULL



#' Dictionary mapping feature id to feature description
#'
#' A dataframe containing id and description of features appearing in \code{input.cov}, \code{codify.feature} and \code{cuisearch.feature}. The variables are as follows:
#'
#' \itemize{
#'   \item Variable. feature id (possible categories: PheCode/CCS/LOINC/RXNORM/CUI)
#'   \item Description. full feature description
#' }
#'
#' @docType data
#' @keywords datasets
#' @name dict_RA
#' @usage data(dict_RA)
#' @format A data frame with 183 rows and 2 variables
NULL



#' Dictionary mapping from ICD9 to Phecode
#'
#' A dataframe containing ICD9 and the mapped phecodes for features appearing in \code{fake_ehr}.
#' @docType data
#' @keywords datasets
#' @name phecode_map
#' @usage data(phecode_map)
#' @format A data frame with 48 rows and 4 variables
NULL


#' Pseudo-individual data for score prediction and model evaluation
#'
#' A dataframe containing pseudo patient id and codified and NLP feature counts (after log transformation). The variables are as follows:
#'
#' \itemize{
#'   \item patient_num unique patient id (pseudo)
#'   \item Y gold label (optional)
#'   \item pi sample probability (optional if uniformly sampled)
#'   \item utl health utility score
#'   \item ... Codified and NLP features
#' }
#'
#' @docType data
#' @keywords datasets
#' @name dat_part
#' @usage data(dat_part)
#' @format A data frame with 242 rows and 187 variables
NULL



#' Fake EHR data for preprocessing illustration
#'
#' A dataframe containing pseudo patient id and longitudinal EHR codes. The variables are as follows:
#'
#' \itemize{
#'   \item patient_num unique patient id (fake)
#'   \item days_since_admission time indicator
#'   \item concept_type two different EHR code types (`DIAG-ICD9` and `LAB-LOINC`) in this fake data
#'   \item concept_code EHR code name
#'   \item value lab values. -999 if the record belongs to `DIAG-ICD9`.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name fake_ehr
#' @usage data(fake_ehr)
#' @format A data frame with 346 rows and 5 variables
NULL



#' Fake EHR data for preprocessing illustration
#'
#' A dataframe containing pseudo patient id and longitudinal EHR codes. The variables are as follows:
#'
#' \itemize{
#'   \item patient_num unique patient id (fake)
#'   \item days_since_admission time indicator
#'   \item concept_type two different EHR code types (`DIAG-ICD9` and `LAB-LOINC`) in this fake data
#'   \item concept_code EHR code name
#' }
#'
#' @docType data
#' @keywords datasets
#' @name fake_ehr
#' @usage data(fake_ehr)
#' @format A data frame with 346 rows and 5 variables
NULL


#' Fake EHR data in wide format with gold labels
#'
#' A dataframe containing pseudo patient id, gold lables and rollup counts for EHR codes. The variables are as follows:
#'
#' \itemize{
#'   \item Y gold labels
#'   \item patient_num unique patient id (fake)
#'   \item ... EHR code counts
#' }
#'
#' @docType data
#' @keywords datasets
#' @name fake_ehr_label_logcount_wide
#' @usage data(fake_ehr_label_logcount_wide)
#' @format A data frame with 8 rows and 50 variables
NULL


#' Fake EHR data in wide format with gold labels and sample probabilities
#'
#' A dataframe containing pseudo patient id, gold lables, sample probabilities and rollup counts for EHR codes. The variables are as follows:
#'
#' \itemize{
#'   \item pi sample probability of the patient
#'   \item Y gold labels
#'   \item patient_num unique patient id (fake)
#'   \item ... EHR code counts
#' }
#'
#' @docType data
#' @keywords datasets
#' @name fake_ehr_label_logcount_wide_pi
#' @usage data(fake_ehr_label_logcount_wide_pi)
#' @format A data frame with 8 rows and 51 variables
NULL


#' Fake EHR data for preprocessing illustration (main input of the function \code{gen_cov_input})
#'
#' A dataframe containing pseudo patient id, days since admission date and EHR codes. The variables are as follows:
#'
#' \itemize{
#'   \item patient_num unique patient id (fake)
#'   \item days_since_admission time indicator
#'   \item code EHR code name
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ehr_data
#' @usage data(ehr_data)
#' @format A data frame with 346 rows and 3 variables
NULL



#' Fake rollup dictionary for preprocessing illustration (main input of the function \code{gen_cov_input})
#'
#' A dataframe containing EHR codes and rollup group codes. The variables are as follows:
#'
#' \itemize{
#'   \item code EHR code name
#'   \item group Rolled-up group code
#' }
#'
#' @docType data
#' @keywords datasets
#' @name rollup_dict
#' @usage data(rollup_dict)
#' @format A data frame with 48 rows and 2 variables
NULL

#' Fake filters for preprocessing illustration (main input of the function \code{gen_cov_input})
#'
#' A dataframe containing EHR codes and the least frequency required for patient inclusion. The variables are as follows:
#'
#' \itemize{
#'   \item code EHR code name
#'   \item freq The least frequency that is required for patients to be included to calculate the covariance matrix
#' }
#'
#' @docType data
#' @keywords datasets
#' @name filter_df
#' @usage data(filter_df)
#' @format A data frame with 1 rows and 2 variables
NULL

#' Conditional feature mean vector (Y=0)
#'
#' A named numeric vector containing the mean of log(feature count) given negative labels for RA disease(Y=0).
#'
#' @docType data
#' @keywords datasets
#' @name mu0
#' @usage data(mu0)
#' @format A named vector with 184 entries.
NULL


#' Conditional feature mean vector (Y=1)
#'
#' A named numeric vector containing the mean of log(feature count) given positive labels for RA disease(Y=1).
#'
#' @docType data
#' @keywords datasets
#' @name mu1
#' @usage data(mu1)
#' @format A named vector with 184 entries.
NULL


#' Conditional feature covariance matrix (Y=0)
#'
#' A numeric covariance matrix of log(feature count) given negative labels for RA disease(Y=0). Must have colnames and rownames.
#'
#' @docType data
#' @keywords datasets
#' @name var0
#' @usage data(var0)
#' @format A 184*184 covariance matrix.
NULL


#' Conditional feature covariance matrix (Y=1)
#'
#' A numeric covariance matrix of log(feature count) given positive labels for RA disease(Y=1). Must have colnames and rownames.
#'
#' @docType data
#' @keywords datasets
#' @name var1
#' @usage data(var1)
#' @format A 184*184 covariance matrix.
NULL


#' Prevalence of PheCode714.1
#'
#' A number indicating the prevalence of RA disease in the pseudo-population.
#'
#' @docType data
#' @keywords datasets
#' @name prev_Y
#' @usage data(prev_Y)
#' @format A number.
NULL
