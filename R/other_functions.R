AUC = function(D,score,wgt=NULL){
  if(is.null(wgt)) {wgt=rep(1,length(D))}
  auc = sum(S.FUN(score,Yi=score,D*wgt,yes.smooth=F)*(1-D)*wgt)/sum((1-D)*wgt)
  return(auc)
}

ROC = function(D, score, wgti=NULL, seq = seq(.01,.99,by=.01)){
  roc = ROC.Est.FUN(D,score,0.5,seq,wgti = wgti)
  roc = matrix(roc[-1],ncol=6)
  colnames(roc) = c("cutoff","est.pos.rate","FPR","TPR","PPV","NPV")
  return(roc)
}

S.FUN <- function(yy,Yi,Di,yes.smooth=F){
  # if(yes.smooth){
  #   Y1i = Yi[Di==1]; n1 = sum(Di); bw = stats::bw.nrd(Y1i)/n1^0.6
  #   c(t(rep(1/n1,n1))%*%stats::pnorm((Y1i-VTM(yy,n1))/bw))
  # }else{
  #   return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  # }
  return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  ##sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
}

sum.I <- function(yy,FUN,Yi,Vi=NULL){
  ## sum_i I(yy FUN Yi)Vi
  # Vi weight
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of decending
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}


ROC.Est.FUN <- function(Di,yyi,yy0,fpr0=NULL,wgti=NULL,yes.smooth=F){
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
  if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));
  mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0
  for(k in 1:pp)
  {
    yy = yy0;
    if(!is.null(fpr0)){
      tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth);
      fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
      TPR = stats::approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y;
      TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR);
      yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }else{
      TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth);
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }
    out.yy = cbind(out.yy, yy)
    out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
    out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
    PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
    out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
    #AUC <- sum((sum.I(yyi[,k],"<=",Yi=yyi[,k],Vi=Di*wgti)+sum.I(yyi[,k],"<",Yi=yyi[,k],Vi=Di*wgti))*(1-Di)*wgti/2
    #             )/(sum((1-Di)*wgti)*sum(Di*wgti))
    AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
    out.AUC <- c(out.AUC, AUC)
  }
  out = c(out.AUC,out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV)
  out
}


Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F){
  yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth)
  return(stats::approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}
