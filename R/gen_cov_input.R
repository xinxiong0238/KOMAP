#' Generate covariance matrix from individual ehr data
#'
#' @param ehr_data Raw EHR data records with 3 colnames indicating the unique patient id, the number of days since admission and the concept code. Must be in order.
#' @param rollup_dict Dataframe with the first column indicating the raw code (e.g., ICD code) and the second the rolled-up group code (e.g., phecode). Can be \code{NULL} if no rollup is needed.
#' @param filter_df Dataframe with the first column indicating the concept code and the second indicating the least frequency that is required for patients to be included. Multiple filters are allowed by specifying different rows. Can be \code{NULL} if no filter is needed.
#' Note that the code count calculation for each patient includes two steps. First, we calculate the raw concept count once per day for each patient(i.e., for those raw codes that show up more than one time in a day, treat as one time). Second, the rolled-up group code count is equal to the sum of all its descendant code once-per-day count.
#' @returns A covariance matrix.
#' @examples
#' data(ehr_data)
#' data(rollup_dict)
#' data(filter_df)
#' input_cov <- gen_cov_input(ehr_data, rollup_dict, filter_df)
#' @export
#' @importFrom rlang .data
gen_cov_input <- function(ehr_data, rollup_dict, filter_df){
  colnames(ehr_data) <- c('patient_num', 'days_since_admission', 'code')
  colnames(filter_df) <- c('code', 'freq')

  if(!is.null(rollup_dict)){
    colnames(rollup_dict) <- c('code', 'group')
  }else{
    rollup_dict <- data.frame(`code` = NA, `group` = NA)
  }

  ehr_raw_count <- ehr_data %>%
    dplyr::group_by(.data$patient_num, .data$code) %>%
    dplyr::summarise(count_once_per_day = length(unique(.data$days_since_admission)))
  ehr_raw_count <- dplyr::left_join(ehr_raw_count, rollup_dict, multiple = "all")
  ehr_raw_count$group[is.na(ehr_raw_count$group)] <-
    ehr_raw_count$code[is.na(ehr_raw_count$group)]

  ehr_rollup_count <- ehr_raw_count %>%
    dplyr::group_by(.data$patient_num, .data$group) %>%
    dplyr::summarise(count = sum(.data$count_once_per_day))

  ehr_count_wide <- tidyr::pivot_wider(ehr_rollup_count, id_cols = c('patient_num'),
                                       names_from = 'group', values_from = 'count')
  ehr_count_wide[is.na(ehr_count_wide)] <- 0
  if(!is.null(filter_df)){
    filter_df <- dplyr::left_join(filter_df, rollup_dict)
    filter_df$group[is.na(filter_df$group)] = filter_df$code[is.na(filter_df$group)]
    if(nrow(filter_df) > 0){
      id_filter = ehr_rollup_count$patient_num
      for(i in 1:nrow(filter_df)){
        id = ehr_rollup_count %>%
          dplyr::filter(.data$group == filter_df$group[i] & .data$count >= filter_df$freq[i]) %>%
          dplyr::select(.data$patient_num)
        id_filter = intersect(id_filter, id$patient_num)
      }
      ehr_count_wide <- ehr_count_wide %>% dplyr::filter(.data$patient_num %in% id_filter)
    }
  }

  ehr_logcount_wide <- ehr_count_wide
  ehr_logcount_wide[,-1] <- log(ehr_logcount_wide[, -1] + 1)
  ehr_cov <- stats::cov(ehr_logcount_wide[, -1])
  return(ehr_cov)
}
