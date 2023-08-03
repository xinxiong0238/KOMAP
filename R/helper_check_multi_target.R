node_format_check_multi_target <- function(targets, input.cov, related.feature, nm.utl){
  message('Check feature format in `input.cov` and `related.feature`...')
  feat.cov = input.cov
  if(is.null(colnames(feat.cov)) | is.null(rownames(feat.cov))){
    stop('The input covariance matrix must have colnames and rownames!')
  }
  if(is.null(nm.utl)){
    stop('You must provide the variable name for the healthcare utility score!')
  }

  nm.feat.all = sort(rownames(feat.cov))
  if(!all(nm.utl %in% nm.feat.all)){
    stop('The utl variable does not exist in your input covariance matrix!')
  }
  select.feat = related.feature
  if(length(intersect(nm.feat.all, select.feat)) == 0){
    stop('No intersection between features in your input covariance matrix and selected features!')
  }

  message(paste0('Num of total feat: ', length(nm.feat.all)))
  message(paste0('Num of related feat: ', length(related.feature)))
  nm.code.filter = intersect(stats::na.omit(nm.feat.all), colnames(feat.cov))
  message(paste0('Num of related feat after intersection: ', length(nm.code.filter)))

  message('\nMore detailed info...')
  if(!is.null(related.feature)){
    ## check phecode
    allfeat.phecode = nm.feat.all[stringr::str_detect(nm.feat.all, '^PheCode\\:[0-9\\.]+$')]
    message(paste0('Num of PheCode in input.cov: ', length(allfeat.phecode)))
    selectfeat.phecode = related.feature[stringr::str_detect(related.feature, '^PheCode\\:[0-9\\.]+$')]
    message(paste0('Num of PheCode in related.feature: ', length(selectfeat.phecode)))


    ## check CCS
    allfeat.ccs = nm.feat.all[stringr::str_detect(nm.feat.all, '^CCS\\:[0-9]+$')]
    message(paste0('Num of CCS in input.cov: ', length(allfeat.ccs)))
    selectfeat.ccs = related.feature[stringr::str_detect(related.feature, '^CCS\\:[0-9]+$')]
    message(paste0('Num of CCS in related.feature: ', length(selectfeat.ccs)))

    ## check LOINC
    allfeat.loinc = nm.feat.all[stringr::str_detect(nm.feat.all, '^LOINC\\:[XC0-9\\-]+$')]
    message(paste0('Num of LOINC in input.cov: ', length(allfeat.loinc)))
    selectfeat.loinc = related.feature[stringr::str_detect(related.feature, '^LOINC\\:[0-9\\-]+$')]
    message(paste0('Num of LOINC in related.feature: ', length(selectfeat.loinc)))

    ## check RXNORM
    allfeat.rxnorm = nm.feat.all[stringr::str_detect(nm.feat.all, '^RXNORM\\:[0-9]+$')]
    message(paste0('Num of RXNORM in input.cov: ', length(allfeat.rxnorm)))
    selectfeat.rxnorm = related.feature[stringr::str_detect(related.feature, '^RXNORM\\:[0-9]+$')]
    message(paste0('Num of RXNORM in related.feature: ', length(selectfeat.rxnorm)))

    ## check CUI
    allfeat.cui = nm.feat.all[stringr::str_detect(nm.feat.all, '^C[0-9]+$')]
    message(paste0('Num of CUI in input.cov: ', length(allfeat.cui)))
    selectfeat.cui = related.feature[stringr::str_detect(related.feature, '^C[0-9]+$')]
    message(paste0('Num of CUI in cuisearch.feature: ', length(selectfeat.cui)))

  }else{
    allfeat.phecode = allfeat.ccs = allfeat.loinc = allfeat.rxnorm =
      selectfeat.phecode = selectfeat.ccs = selectfeat.loinc = selectfeat.rxnorm = NULL
  }


  ## check others
  allfeat.good = c(allfeat.phecode, allfeat.ccs,
                   allfeat.loinc, allfeat.rxnorm, allfeat.cui)
  selectfeat.good = c(selectfeat.phecode, selectfeat.ccs,
                      selectfeat.loinc, selectfeat.rxnorm, selectfeat.cui)
  allfeat.other = setdiff(nm.feat.all, allfeat.good)
  selectfeat.other = setdiff(c(related.feature), selectfeat.good)


  if(length(allfeat.other) > 0){
    warning('There exists some features in your input covariance matrix that do not belong to standard categories (PheCode/CCS/RXNORM/LOINC)!')
  }
  if(length(selectfeat.other) > 0){
    warning('There exists some selected features that do not belong to standard categories (PheCode/CCS/RXNORM/LOINC)!')
  }

}


node_format_check_part_multi_target <- function(targets, input.cov, nm.utl){
  message('Check feature format in `input.cov`...')
  feat.cov = input.cov
  if(is.null(colnames(feat.cov)) | is.null(rownames(feat.cov))){
    stop('The input covariance matrix must have colnames and rownames!')
  }
  if(is.null(nm.utl)){
    stop('You must provide the variable name for the healthcare utility score!')
  }

  nm.feat.all = sort(rownames(feat.cov))
  if(!all(nm.utl %in% nm.feat.all)){
    stop('The utl variable does not exist in your input covariance matrix!')
  }
  message(paste0('Num of total feat: ', length(nm.feat.all)))

}




KOMAP.est.check_multi_target <- function(input.cov, targets, nm.utl, related.feature){

  ### Check conditions for KOMAP.est
  if(any(!targets %in% rownames(input.cov))) stop('Your main surrogates does not exist in the covariance matrix!')

  ### Check input feature format
  node_format_check_multi_target(targets, input.cov, related.feature, nm.utl)
}

KOMAP.est.check.part_multi_target <- function(input.cov, targets, nm.utl){

  ### Check conditions for KOMAP.est
  if(any(!targets %in% rownames(input.cov))) stop('Your main surrogates does not exist in the covariance matrix!')

  ### Check input feature format
  node_format_check_part_multi_target(targets, input.cov, nm.utl)
}


