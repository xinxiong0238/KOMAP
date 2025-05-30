#' Run komap main function
#'
#' @param input.cov A covariance matrix of the target disease with colnames and rownames equal to variable names.
#' @param is.wide If \code{TRUE}, \code{input.cov} is treated as a normal p by p covariance matrix. If \code{FALSE}, treat \code{input.cov} as in a long format with three columns (the first two indicating row and column names and the last one indicating covariance value). If only inputing the lower/upper triangle part, the function will automatically fill the other side assuming symmetric covariance matrix.
#' @param targets Main surrogates of the target disease.
#' @param related.feature A vector of features related to the target code. Can be generated by KESER or ONCE. If \code{NULL}, then use all features in \code{input.cov}.
#' @param nm.utl Variable name of health utility score.
#' @param dict Dictionary mapping feature id (the first column) to full description (the second column). Default is set to be NULL.
#' @param pred If \code{TRUE}, predict disease scores for individual data.
#' @param eval.real If \code{TRUE}, evaluate model performance based on prediction of individual data with gold labels.
#' @param eval.sim If \code{TRUE}, evaluate model performance based on conditional covariance matrices, conditional mean vectors and prevalence of the target disease.
#' @param mu0 Conditional mean vector of features given a negative label. Must be provided if \code{eval.sim}=\code{TRUE}.
#' @param mu1 Conditional mean vector of features given a positive label. Must be provided if \code{eval.sim}=\code{TRUE} and must have the same length and in the same feature order as \code{mu0}.
#' @param var0 Conditional covariance matrix of features given a negative label. Must be provided if \code{eval.sim}=\code{TRUE}.
#' @param var1 Conditional covariance matrix of features given a positive label. Must be provided if \code{eval.sim}=\code{TRUE} and must have the same dimension and in the same feature order as \code{var0}.
#' @param prev_Y Prevalence of the target disease. Must be provided if \code{eval.sim}=\code{TRUE}.
#' @param B The number of samples to calculate simulated AUC, given the conditional Gaussian assumption. Default is 10000.
#' @param dat.part Patient-level data for prediction. Must be provided if \code{pred}=\code{TRUE} or \code{eval.real}.
#' @param nm.id Name of the patient id column in \code{dat.part} and \code{gold.label}. Must be provided if \code{pred}=\code{TRUE} or \code{eval.real}.
#' @param nm.pi The column name of sample probability in \code{gold.label} if individuals in \code{gold.label} are generated with different probabilities. Default is \code{NULL}, meaning a uniform sample distribution.
#' @param nm.y The column name of labels in \code{gold.label}. Must be provided if \code{eval.real}=\code{TRUE}.
#' @returns A list containing estimation and/or prediction and/or evaluation results by running KOMAP.
#' @importFrom rlang .data
#' @export
KOMAP_corrupt_multi_target <- function(input.cov.train, input.cov.valid, is.wide = TRUE, targets, related.feature,
                                       nm.disease, nm.utl, nm.corrupts,
                                       dict = NULL, pred = FALSE, eval.real = FALSE, eval.sim = TRUE,
                                       mu0 = NULL, mu1 = NULL, var0 = NULL, var1 = NULL, prev_Y = NULL, B = 10000,
                                       dat.part = NULL, nm.id = NULL, nm.pi = NULL, nm.y = NULL){
  oldw <- getOption("warn")
  options(warn = -1)
  if(!is.wide){
    input.cov = as.data.frame(input.cov)
    unique.node = unique(c(input.cov[,1], input.cov[,2]))
    colnames(input.cov) = c('from', 'to', 'cov')

    input.cov.wide = stats::reshape(as.data.frame(input.cov), idvar = "from", timevar = "to",
                                    direction = "wide")
    rownames(input.cov.wide) = input.cov.wide$from; input.cov.wide$from = NULL
    colnames(input.cov.wide) = stringr::str_remove(colnames(input.cov.wide), '^cov\\.')
    miss.row = setdiff(unique.node, rownames(input.cov.wide))
    miss.row.matrix = matrix(NA, nrow = length(miss.row), ncol = ncol(input.cov.wide))
    rownames(miss.row.matrix) = miss.row; colnames(miss.row.matrix) = colnames(input.cov.wide)
    input.cov.wide = rbind(input.cov.wide, miss.row.matrix)

    miss.col = setdiff(unique.node, colnames(input.cov.wide))
    miss.col.matrix = matrix(NA, nrow = nrow(input.cov.wide), ncol = length(miss.col))
    colnames(miss.col.matrix) = miss.col; rownames(miss.col.matrix) = rownames(input.cov.wide)
    input.cov.wide = cbind(input.cov.wide, miss.col.matrix)

    input.cov = input.cov.wide
    input.cov.lower = input.cov[lower.tri(input.cov)]; input.cov.upper = input.cov[upper.tri(input.cov)]
    if(all(is.na(input.cov.lower))){
      input.cov[lower.tri(input.cov)] = input.cov[upper.tri(input.cov)]
    }else{
      if(all(is.na(input.cov.upper))){
        input.cov[upper.tri(input.cov)] = input.cov[lower.tri(input.cov)]
      }
    }
    input.cov[is.na(input.cov)] = 0
    message(paste0('\nInput long format data, transformed to wide format covariance matrix (',
                   length(unique.node),' unique nodes).'))
  }
  ### Check user's input
  if(!is.null(related.feature)){
    KOMAP.est.check_multi_target(input.cov.train, targets, nm.utl, related.feature)
  }else{
    KOMAP.est.check.part_multi_target(input.cov.train, targets, nm.utl)
  }
  out_main = KOMAP.est.corrupt.multi.target(input.cov.train, input.cov.valid, targets, related.feature,
                                            nm.disease, nm.utl, nm.corrupts, dict)

  message('\nFinish estimating coefficients.')
  nm.multi = NULL
  out = out_main$lst
  method_nm = names(out)
  feat.out = unique(unlist(sapply(out, function(x){x['feat']})))
  out_return = c(list(`est` = out_main))

  if(eval.sim){
    KOMAP.sim_auc.check(out, feat.out, mu0, mu1, var0, var1, prev_Y)
    out.sim = KOMAP.sim_auc(out, mu0, mu1, var0, var1, prev_Y, B = B)
    out_return = c(out_return, `sim_eval` = list(out.sim))
    message('Finish estimating AUC.')
  }

  if(pred){
    KOMAP.pred.check(out, feat.out, dat.part, nm.utl, nm.id)
    pred.prob = KOMAP.pred.corrupt.multi(out, targets, dat.part, nm.utl, nm.corrupts, nm.multi, nm.id)
    out_return = c(out_return, `pred_prob` = list(pred.prob))
    message('Finish predicting scores.')
  }

  if(eval.real){
    if(pred){
      gold.label = dat.part[,c(nm.id, nm.y, nm.pi)]
      gold.label = stats::na.omit(gold.label)
      KOMAP.eval.check(pred.prob, gold.label, nm.pi, nm.y, nm.id, method_nm)
      out.eval = KOMAP.eval(pred.prob, gold.label, nm.pi, nm.y, nm.id, method_nm)
      out_return = c(out_return,
                     `real_eval` = list(out.eval))
      message('Finish evaluating model prediction.')
    }else{
      KOMAP.pred.check(out, feat.out, dat.part, nm.utl, nm.id)
      pred.prob = KOMAP.pred.corrupt.multi(out, targets, dat.part, nm.utl, nm.corrupts, nm.multi, nm.id)
      gold.label = dat.part[,c(nm.id, nm.y, nm.pi)]
      gold.label = stats::na.omit(gold.label)
      KOMAP.eval.check(pred.prob, gold.label, nm.pi, nm.y, nm.id, method_nm)
      out.eval = KOMAP.eval(pred.prob, gold.label, nm.pi, nm.y, nm.id, method_nm)
      out_return = c(out_return, `pred_prob` = list(pred.prob),
                     `real_eval` = list(out.eval))
      message('Finish evaluating model prediction.')
    }
  }
  options(warn = oldw)
  return(out_return)
}


KOMAP.est.corrupt.multi.target <- function(input.cov.train, input.cov.valid, targets, related.feature,
                                           nm.disease, nm.utl, nm.corrupts, dict){
  out.all = c()
  for(target in targets){
    if(str_detect(target, '^PheCode\\:')){
      parent.code0 = stringr::str_extract(target, 'PheCode\\:[0-9]+')
      parent.code1 = stringr::str_extract(target, 'PheCode\\:[0-9]+\\.[0-9]{1}')
      parent.code2 = stringr::str_extract(target, 'PheCode\\:[0-9]+\\.[0-9]{2}')
      if(!is.na(parent.code2)){
        idd = which(!rownames(input.cov.train) %in% c(parent.code0, parent.code1))
        input.cov.train = input.cov.train[idd, idd]
        idd = which(!rownames(input.cov.valid) %in% c(parent.code0, parent.code1))
        input.cov.valid = input.cov.valid[idd, idd]
      }else{
        if(!is.na(parent.code1)){
          idd = which(!rownames(input.cov.train) %in% c(parent.code0))
          input.cov.train = input.cov.train[idd, idd]
          idd = which(!rownames(input.cov.valid) %in% c(parent.code0))
          input.cov.valid = input.cov.valid[idd, idd]
        }
      }
    }
  }


  if(is.null(related.feature)){
    related.feature = setdiff(colnames(input.cov.train), nm.utl)
  }
  related.feature = unique(c(nm.corrupts, related.feature))

  out = gen.KOMAP.est.table.corrupt.multi.target(input.cov.train, input.cov.valid, nm.disease, nm.utl,
                                                 nm.corrupts, targets, dict, related.feature)
  out.df = sapply(1:length(out), function(i){
    a = data.frame(`disease` = nm.disease,
                   `method` = names(out)[i],
                   `target` = paste0(out[[i]]$target, collapse = ', '),
                   `feat` = out[[i]]$beta$feat,
                   `coeff` = out[[i]]$beta$theta)
    return(a)
  }, simplify = FALSE)
  out.re = data.table::rbindlist(out.df)
  # out.re = out.re[out.re$method %in% str_remove(input$komap_feat,"combine\\_"), ]
  out.all = rbind(out.all, out.re)
  if(!is.null(dict)){
    out.all$desc = dict[match(out.all$feat, dict[,1]),2]
    out.all = out.all[,c('disease', 'method', 'target', 'feat', 'desc', 'coeff')]
    out.all$coeff = as.numeric(out.all$coeff)
    out.all = out.all[order(out.all$coeff, decreasing = TRUE), ]
    out.all$desc[out.all$feat == nm.utl] = 'Healthcare Utility'
  }else{
    out.all$desc = NA
    out.all = out.all[,c('disease', 'method', 'target', 'feat', 'desc', 'coeff')]
    out.all$coeff = as.numeric(out.all$coeff)
    out.all = out.all[order(out.all$coeff, decreasing = TRUE), ]
    out.all$desc[out.all$feat == nm.utl] = 'Healthcare Utility'
  }
  out.all = out.all %>% dplyr::arrange(.data$disease,.data$ method, -abs(.data$coeff))
  return(list(`long_df` = out.all,
              `lst` = out))
}


gen.KOMAP.est.table.corrupt.multi.target <- function(input.cov.train, input.cov.valid,
                                                     nm.disease, nm.utl, nm.corrupts, targets,
                                                     dict, related.feature){
  alpha.glmnent = 0.15
  colnames(input.cov.train)[stringr::str_detect(colnames(input.cov.train), '^CCS-PCS')] =
    stringr::str_replace(colnames(input.cov.train)[stringr::str_detect(colnames(input.cov.train), '^CCS-PCS')],
                         '^CCS-PCS', 'CCS')
  rownames(input.cov.train)[stringr::str_detect(rownames(input.cov.train), '^CCS-PCS')] =
    stringr::str_replace(rownames(input.cov.train)[stringr::str_detect(rownames(input.cov.train), '^CCS-PCS')],
                         '^CCS-PCS', 'CCS')

  colnames(input.cov.valid)[stringr::str_detect(colnames(input.cov.valid), '^CCS-PCS')] =
    stringr::str_replace(colnames(input.cov.valid)[stringr::str_detect(colnames(input.cov.valid), '^CCS-PCS')],
                         '^CCS-PCS', 'CCS')
  rownames(input.cov.valid)[stringr::str_detect(rownames(input.cov.valid), '^CCS-PCS')] =
    stringr::str_replace(rownames(input.cov.valid)[stringr::str_detect(rownames(input.cov.valid), '^CCS-PCS')],
                         '^CCS-PCS', 'CCS')

  nm.others = colnames(input.cov.train)
  nm.others = nm.others[!nm.others %in% c(targets, nm.utl)]
  all.nm = c(targets, nm.others, nm.utl)

  feat.cov.train = input.cov.train[all.nm, all.nm]
  feat.cov.valid = input.cov.valid[all.nm, all.nm]
  nm.feat.all = sort(colnames(feat.cov.train))
  nm.feat.filter = intersect(c(related.feature, nm.corrupts), stats::na.omit(nm.feat.all))
  out = c()
  for(jj in 1:length(targets)){
    method =  list(`feat` = unique(c(nm.feat.filter, targets)),
                   `target` = targets[jj],
                   `corrupt_target` = nm.corrupts[jj])
    out = c(out, list(method))
  }
  names(out) = paste0('Target_', targets, ' + related_feature')
  if(length(targets) > 1){
    method =  list(`feat` = unique(c(nm.feat.filter, targets)),
                   `target` = targets,
                   `corrupt_target` = nm.corrupts)
    out = c(out, `allTargets + related_feature` = list(method))
  }

  for(i in 1:length(out)){
    plot_data = c()
    method = out[[i]]
    feat = method$feat; out_targets = method$target; out_corrupts = method$corrupt_target
    n_komap = length(setdiff(feat, c(nm.corrupts, out_targets)))

    b.all = score = alpha = alpha_ma = c()

    all.nm =  c(out_targets, setdiff(feat, out_targets), nm.utl)
    feat.cov.part = feat.cov.train[all.nm, all.nm]
    feat.cov.part.valid = feat.cov.valid[all.nm, all.nm]
    U = svd(feat.cov.part)
    U = t(U$u %*% sqrt(diag(U$d)))
    U.valid = svd(feat.cov.part.valid)
    U.valid = t(U.valid$u %*% sqrt(diag(U.valid$d)))
    colnames(U) = rownames(U) = colnames(U.valid) = rownames(U.valid) = colnames(feat.cov.part)
    U_new = U
    for(j in 1:length(out_targets)){
      target.j = out_targets[j]
      reg1 = stats::lm(U[, target.j] ~ 0 + U[, nm.utl])
      U_new[, target.j] = reg1$residuals
      nm.corrupt = out_corrupts[j]

      U_new[, nm.corrupt] = U_new[, nm.corrupt] - reg1$coefficients * U_new[, nm.utl]
      cor(U_new[, target.j], U_new[, nm.corrupt])

      U.valid[, target.j] = U.valid[, target.j] - reg1$coefficients * U.valid[, nm.utl]
      U.valid[, nm.corrupt] = U.valid[, nm.corrupt] - reg1$coefficients * U.valid[, nm.utl]
      cor(U.valid[, target.j], U.valid[, nm.corrupt])

      alpha = cbind(alpha, reg1$coefficients)
      junk = rep(0, length(out_targets))
      junk[j] = 1
      alpha_ma = cbind(alpha_ma,
                       c(junk, rep(0, n_komap), -reg1$coefficients))

    }
    alpha_ma = cbind(alpha_ma,
                     rbind(matrix(0, length(out_targets), n_komap + length(nm.utl)),
                           diag(1, nrow = n_komap + length(nm.utl))))
    colnames(alpha_ma) = c(out_targets, setdiff(feat, c(nm.corrupts, out_targets)), nm.utl)

    for(j in 1:length(out_targets)){
      target.j = out_targets[j]
      corrupt.j = out_corrupts[j]
      Y = U_new[, target.j]
      id_target = which(colnames(U_new) == target.j)
      id_target_corrupt = which(colnames(U_new) == corrupt.j)
      id_elsetarget_corrupt = match(setdiff(nm.corrupts, corrupt.j), colnames(U_new))
      set.seed(10101)
      cv.fit <- glmnet::cv.glmnet(U_new[, -c(id_target_corrupt, id_elsetarget_corrupt)], Y,
                                  intercept = F, alpha = alpha.glmnent, standardize = F)
      lambda_vec = cv.fit$lambda
      mse_vec = sapply(lambda_vec, function(lambda){
        model.fit <- glmnet::glmnet(U_new[, -c(id_target_corrupt, id_elsetarget_corrupt)], Y, intercept = F, alpha = alpha.glmnent, standardize = F,
                                    lambda = lambda)
        junk = U.valid[, -c(id_target, id_elsetarget_corrupt)]
        colnames(junk)[colnames(junk) == corrupt.j] = target.j
        coef.beta = as.matrix(model.fit$beta)
        predict.valid.fit <- junk %*% coef.beta[colnames(junk), ]
        mse.fit = mean((U.valid[, id_target] - predict.valid.fit)^2)
        return(mse.fit)
      })
      # if(length(out_targets) == 1){
      #   plot(log(lambda_vec), mse_vec, type = 'l', main = paste(nm.disease, ' Corrupt single, ', target.j))
      #   abline(v = log(lambda_vec[which.min(mse_vec)]), col = 'red')
      # }else{
      #   plot(log(lambda_vec), mse_vec, type = 'l', main = paste(nm.disease, ' Corrupt multiple, ', target.j))
      #   abline(v = log(lambda_vec[which.min(mse_vec)]), col = 'red')
      # }
      plot_data = c(plot_data, list(lambda = lambda_vec, mse = mse_vec))
      model.fit <- glmnet::glmnet(U_new[, -c(id_target_corrupt, id_elsetarget_corrupt)], Y, intercept = F, alpha = alpha.glmnent, standardize = F,
                                  lambda = lambda_vec[which.min(mse_vec)])

      score0 = U_new[, -c(id_target_corrupt, id_elsetarget_corrupt)] %*% model.fit$beta
      rownames(model.fit$beta)[which(rownames(model.fit$beta) == corrupt.j)] = target.j
      b.all0 <- as.vector(model.fit$beta); names(b.all0) = rownames(model.fit$beta)
      # b.all0 = solve(XTX + W) %*% XTY
      b.all <- cbind(b.all, b.all0[match(names(b.all0), colnames(alpha_ma))])
      # score = cbind(score, score0[match(rownames(score0), colnames(alpha_ma))])
      score = cbind(score, score0)
    }
    if(length(out_targets) > 1){
      pca.score = stats::princomp(score)
      gamma = pca.score$loadings[,1]

      theta = alpha_ma %*% b.all %*% gamma
      beta.all <- data.frame(`feat` = colnames(alpha_ma), `theta` = theta)
    }else{
      theta = alpha_ma %*% b.all
      beta.all <- data.frame(`feat` = colnames(alpha_ma), `theta` = as.vector(theta))
    }
    method$beta = beta.all
    method$plot_data = plot_data
    out[[i]] = method
  }
  return(out)
}

#' @import mclust
KOMAP.pred.corrupt.multi <- function(out, targets, dat.part, nm.utl,
                                     nm.corrupts, nm.multi, nm.id = 'patient_num'){
  pred.cluster = pred.prob = pred.score = data.frame(`patient_num` = dat.part[,nm.id])
  for(i in 1:length(out)){
    method = out[[i]]
    feat = method$beta$feat
    # feat = setdiff(feat, c(nm.corrupt.code, nm.corrupt.cui))
    # print(setdiff(feat,  colnames(dat.part)))
    feat = intersect(feat, colnames(dat.part))
    if(length(feat) == 1){
      dat.part.filter = data.frame(dat.part[, feat])
      colnames(dat.part.filter) = feat
    }else{
      dat.part.filter = dat.part[, feat]
    }
    dat.part.filter = as.data.frame(cbind(dat.part[,nm.utl], dat.part.filter))
    colnames(dat.part.filter)[1] = nm.utl
    if(!is.null(nm.multi)){
      dat.part.filter = as.data.frame(cbind(dat.part[,nm.multi], dat.part.filter))
      colnames(dat.part.filter)[1] = nm.multi
    }
    b.all = data.frame(`coeff` = method$beta$theta, `feat` = method$beta$feat)
    b.all = b.all[b.all$feat %in% feat, ]
    S.norm <- as.matrix(dat.part.filter[,b.all$feat]) %*% matrix(b.all$coeff, ncol=1)

    ## Fit gaussian mixture model on S.norm, we only use the "length(nm.logS.ori) = 1" case:
    fit = S.norm
    pred.score = cbind(pred.score, fit)
    junk = mclust::Mclust(S.norm, G = 2, verbose = FALSE)
    cor1 = cor(dat.part[, targets[1]], junk$classification, method = 'kendall')
    cor2 = cor(dat.part[, targets[1]], 3-junk$classification, method = 'kendall')
    # cluster_i = ifelse(cor1 > cor2, junk$classification - 1, 3 - junk$classification)
    if(cor1 > cor2){
      cluster_i = factor(junk$classification, levels = c(1, 2), labels = c('no disease', 'disease'))
      prob_i = junk$z[, 2]
    }else{
      cluster_i = factor(junk$classification, levels = c(1, 2), labels = c('disease', 'no disease'))
      prob_i = junk$z[, 1]
    }
    pred.cluster = cbind(pred.cluster, cluster_i)
    pred.prob = cbind(pred.prob, prob_i)
  }
  colnames(pred.prob)[-1] = names(out)
  colnames(pred.cluster)[-1] = names(out)
  colnames(pred.score)[-1] = names(out)
  return(list(`pred.score` = stats::na.omit(pred.score),
              `pred.prob` = stats::na.omit(pred.prob),
              `pred.cluster` = stats::na.omit(pred.cluster)))
}



KOMAP.eval <- function(pred.prob, gold.label, nm.pi = NULL, nm.y = 'Y', nm.id, method_nm){
  dat.merge = dplyr::left_join(gold.label, pred.prob$pred.score, by = nm.id)
  dat.merge = stats::na.omit(dat.merge)
  Y <- dat.merge[, nm.y]
  auc_table <- c()
  F_table <- c()
  prev_vec <- c()
  label_num_vec <- c()
  # Since the sampling of gold label set is not completely at random, we need to include weights:
  if(!is.null(nm.pi)){
    w <- 1 / dat.merge[,nm.pi]
  }else{
    dat.merge$pi = w <- rep(1/nrow(dat.merge), nrow(dat.merge))
  }
  prev_vec <- c(prev_vec, mean(Y))
  label_num_vec <- c(label_num_vec, length(Y))
  num_method <- length(method_nm)
  auc_vec <- rep(0, num_method)
  F_vec <- rep(0, num_method)
  roc_lst <- vector('list', num_method)

  for (t in 1:num_method){
    # AUC
    auc_vec[t] <- AUC(Y, as.vector(dat.merge[,method_nm[t]]), wgt=w)
    # ROC table
    roc_lst[[t]] <- ROC(Y, as.vector(dat.merge[,method_nm[t]]), wgti=w, seq = seq(.01,.99,by=1e-5))[-1,]
    # Maximum F-score
    F_score_all <- 2 / (1 / roc_lst[[t]][,4] + 1 / roc_lst[[t]][,5])
    F_vec[t] <- max(F_score_all)
  }
  return(data.frame(`method` = method_nm,
                    `auc` = auc_vec,
                    `F_score_max` = F_vec)
         # `F_score_max` = F_vec
  )
}




