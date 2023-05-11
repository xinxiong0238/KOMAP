node_format_check <- function(target.code, target.cui, input.cov, codify.feature, cuisearch.feature, nm.utl){
  message('Check feature format in `input.cov`, `codify.feature` and/or `cuisearch.feature`...')
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
  select.feat = c(codify.feature, cuisearch.feature)
  if(length(intersect(nm.feat.all, select.feat)) == 0){
    stop('No intersection between features in your input covariance matrix and selected features!')
  }


  #### Both provide target code and target cui
  if(!is.null(target.code) & !is.null(target.cui)){
    nm.code.filter = codify.feature
    nm.cui.filter.final = cuisearch.feature

    message(paste0('Num of total feat: ', length(nm.feat.all)))
    message(paste0('Num of selected codify feat: ', length(nm.code.filter)))
    nm.code.filter = intersect(stats::na.omit(nm.code.filter), colnames(feat.cov))
    message(paste0('Num of selected codify feat after intersection: ', length(nm.code.filter)))

    message(paste0('Num of selected NLP feat: ', length(nm.cui.filter.final)))
    nm.cui.filter.final = intersect(stats::na.omit(nm.cui.filter.final), colnames(feat.cov))
    message(paste0('Num of selected NLP feat after intersection: ', length(nm.cui.filter.final)))
  }else{
    if(!is.null(target.code)){
      nm.code.filter = codify.feature

      message(paste0('Num of total feat: ', length(nm.feat.all)))
      message(paste0('Num of selected codify feat: ', length(nm.code.filter)))
      nm.code.filter = intersect(stats::na.omit(nm.code.filter), colnames(feat.cov))
      message(paste0('Num of selected codify feat after intersection: ', length(nm.code.filter)))

    }else{
      nm.cui.filter.final = cuisearch.feature

      message(paste0('Num of total feat: ', length(nm.feat.all)))
      message(paste0('Num of selected NLP feat: ', length(nm.cui.filter.final)))
      nm.cui.filter.final = intersect(stats::na.omit(nm.cui.filter.final), colnames(feat.cov))
      message(paste0('Num of selected NLP feat after intersection: ', length(nm.cui.filter.final)))
    }
  }

  message('\nMore detailed info...')
  if(!is.null(codify.feature)){
    ## check phecode
    allfeat.phecode = nm.feat.all[stringr::str_detect(nm.feat.all, '^PheCode\\:[0-9\\.]+$')]
    message(paste0('Num of PheCode in input.cov: ', length(allfeat.phecode)))
    selectfeat.phecode = codify.feature[stringr::str_detect(codify.feature, '^PheCode\\:[0-9\\.]+$')]
    message(paste0('Num of PheCode in codify.feature: ', length(selectfeat.phecode)))


    ## check CCS
    allfeat.ccs = nm.feat.all[stringr::str_detect(nm.feat.all, '^CCS\\:[0-9]+$')]
    message(paste0('Num of CCS in input.cov: ', length(allfeat.ccs)))
    selectfeat.ccs = codify.feature[stringr::str_detect(codify.feature, '^CCS\\:[0-9]+$')]
    message(paste0('Num of CCS in codify.feature: ', length(selectfeat.ccs)))

    ## check LOINC
    allfeat.loinc = nm.feat.all[stringr::str_detect(nm.feat.all, '^LOINC\\:[XC0-9\\-]+$')]
    message(paste0('Num of LOINC in input.cov: ', length(allfeat.loinc)))
    selectfeat.loinc = codify.feature[stringr::str_detect(codify.feature, '^LOINC\\:[0-9\\-]+$')]
    message(paste0('Num of LOINC in codify.feature: ', length(selectfeat.loinc)))

    ## check RXNORM
    allfeat.rxnorm = nm.feat.all[stringr::str_detect(nm.feat.all, '^RXNORM\\:[0-9]+$')]
    message(paste0('Num of RXNORM in input.cov: ', length(allfeat.rxnorm)))
    selectfeat.rxnorm = codify.feature[stringr::str_detect(codify.feature, '^RXNORM\\:[0-9]+$')]
    message(paste0('Num of RXNORM in codify.feature: ', length(selectfeat.rxnorm)))
  }else{
    allfeat.phecode = allfeat.ccs = allfeat.loinc = allfeat.rxnorm =
      selectfeat.phecode = selectfeat.ccs = selectfeat.loinc = selectfeat.rxnorm = NULL
  }

  if(!is.null(cuisearch.feature)){
    ## check CUI
    allfeat.cui = nm.feat.all[stringr::str_detect(nm.feat.all, '^C[0-9]+$')]
    message(paste0('Num of CUI in input.cov: ', length(allfeat.cui)))
    selectfeat.cui = cuisearch.feature[stringr::str_detect(cuisearch.feature, '^C[0-9]+$')]
    message(paste0('Num of CUI in cuisearch.feature: ', length(selectfeat.cui)))
  }else{
    allfeat.cui = selectfeat.cui = NULL
  }

  ## check others
  allfeat.good = c(allfeat.phecode, allfeat.ccs,
                   allfeat.loinc, allfeat.rxnorm, allfeat.cui)
  selectfeat.good = c(selectfeat.phecode, selectfeat.ccs,
                      selectfeat.loinc, selectfeat.rxnorm, selectfeat.cui)
  allfeat.other = setdiff(nm.feat.all, allfeat.good)
  selectfeat.other = setdiff(c(codify.feature, cuisearch.feature), selectfeat.good)


  if(length(allfeat.other) > 0){
    warning('There exists some features in your input covariance matrix that do not belong to standard categories (PheCode/CCS/RXNORM/LOINC)!')
  }
  if(length(selectfeat.other) > 0){
    warning('There exists some selected features that do not belong to standard categories (PheCode/CCS/RXNORM/LOINC)!')
  }

}


node_format_check_part <- function(target.code, target.cui, input.cov, nm.utl){
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




KOMAP.est.check <- function(input.cov, target.code, target.cui, nm.utl,
                        codify.feature, cuisearch.feature){

  ### Check conditions for KOMAP.est
  if(is.null(target.code)){
    stop('You must specify the main phecode for the target disease!')
  }else{
    if(!target.code %in% rownames(input.cov)) stop('Your main phecode does not exist in the covariance matrix!')
  }

  if(!is.null(target.cui)){
    if(!target.cui %in% rownames(input.cov)) stop('Your main CUI does not exist in the covariance matrix!')
  }

  ### Check input feature format
  node_format_check(target.code, target.cui, input.cov, codify.feature, cuisearch.feature, nm.utl)
}

KOMAP.est.check.part <- function(input.cov, target.code, target.cui, nm.utl){

  ### Check conditions for KOMAP.est
  if(is.null(target.code)){
    stop('You must specify the main phecode for the target disease!')
  }else{
    if(!target.code %in% rownames(input.cov)) stop('Your main phecode does not exist in the covariance matrix!')
  }

  if(!is.null(target.cui)){
    if(!target.cui %in% rownames(input.cov)) stop('Your main CUI does not exist in the covariance matrix!')
  }

  ### Check input feature format
  node_format_check_part(target.code, target.cui, input.cov, nm.utl)
}



KOMAP.sim_auc.check <- function(out, feat.out, mu0, mu1, var0, var1, prev_Y){
  if(any(is.null(mu0), is.null(mu1), is.null(var0), is.null(var1), is.null(prev_Y))){
    stop('You must provide conditional cov matrices, mean vectors and prevalence value to get simulated AUC!')
  }

  mu0_nm = names(mu0); mu1_nm = names(mu1); var0_nm = colnames(var0); var1_nm = colnames(var1)
  if(any(is.null(mu0_nm), is.null(mu1_nm))){
    stop('Conditional mean vectors must be assigned with feature names!')
  }
  if(any(is.null(var0_nm), is.null(var1_nm))){
    stop('Coditional mean vectors must have colnames corresponding to the feature names!')
  }
  if(any(!mu0_nm == mu1_nm)){
    stop('The name of mu0 and mu1 should be in the same order!')
  }
  if(any(!var0_nm == var1_nm)){
    stop('The name of var0 and var1 should be in the same order!')
  }
  if(any(!var0_nm == mu0_nm)){
    stop('The name of var0 and mu0 should be in the same order!')
  }
  if(any(!var1_nm == mu1_nm)){
    stop('The name of var1 and mu1 should be in the same order!')
  }

  if(length(intersect(feat.out, mu0_nm)) == 0){
    stop('Features in the conditional covariance matrix have no overlap with features in the trained model!')
  }
}


KOMAP.pred.check <- function(out, feat.out, dat.part, nm.utl, nm.id){
  ### Check conditions for KOMAP.pred
  if(any(is.null(dat.part), is.null(nm.id))){
    stop('You must provide individual data and its id column name (dat.part and nm.id argument) to get score prediction!')
  }
  if(!all(nm.utl %in% colnames(dat.part))){
    stop('The health utility variable does not exist in `dat.part`!')
  }
  if(!nm.id %in% colnames(dat.part)){
    stop('The patient id column does not exist in `dat.part`!')
  }
  if(any(duplicated(dat.part[,nm.id]))){
    stop('The id column you specify has duplicated entries!')
  }
  if(length(intersect(feat.out, colnames(dat.part))) == 0){
    stop('Features in the `dat.part` data have no overlap with features in the trained model!')
  }
}


KOMAP.eval.check <- function(pred.prob, gold.label, nm.pi, nm.y, nm.id, method_nm){
  ### Check conditions for KOMAP.eval
  if(!nm.id %in% colnames(gold.label)){
    stop('The patient id column does not exist in `gold.label`!')
  }
  if(!nm.y %in% colnames(gold.label)){
    stop('The nm.y column does not exist in `gold.label`!')
  }
  dat.merge = dplyr::left_join(gold.label, pred.prob, by = nm.id)
  if(nrow(dat.merge) == 0){
    stop('There is no patient in the intersection of `gold.label` and `dat.part`!')
  }
}
