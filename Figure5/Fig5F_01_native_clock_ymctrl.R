suppressMessages({library(glmnet);library(limma);library(data.table)})
set.seed(1)
ages<-c(Y=5,M=11,CTRL=22)   # drop public O; our CTRL is the sole OLD anchor
tissues<-c("Brain","Liver","Kidney","Spleen")
res<-list();pp<-list();tp<-list()
for(t in tissues){
  cnt<-as.matrix(read.csv(sprintf("merged/%s_counts.csv",t),row.names=1,check.names=FALSE))
  meta<-read.csv(sprintf("merged/%s_meta.csv",t),row.names=1)
  cpm<-log2(sweep(cnt,2,colSums(cnt),"/")*1e6+1)
  tr<-rownames(meta)[meta$group%in%c("Y","M","CTRL")]   # TRAIN = Y,M,CTRL (NO public O)
  age<-ages[meta[tr,"group"]];names(age)<-tr
  keep<-rowMeans(cpm[,tr]>1)>=0.5;X<-cpm[keep,,drop=FALSE]
  tt<-topTable(eBayes(lmFit(X[,tr],model.matrix(~age))),coef="age",number=Inf,sort.by="none")
  agenes<-rownames(tt)[tt$P.Value<0.05]
  xp<-t(X[agenes,tr]);best<-NULL
  for(a in seq(0.05,0.95,0.05)){cv<-cv.glmnet(xp,age,alpha=a,nfolds=length(age),grouped=FALSE,type.measure="mae");m<-min(cv$cvm);if(is.null(best)||m<best$m)best<-list(a=a,m=m,lam=cv$lambda.min)}
  loo<-setNames(numeric(length(tr)),tr)
  for(i in seq_along(tr)){s<-tr[-i];cvi<-cv.glmnet(t(X[agenes,s]),age[s],alpha=best$a,nfolds=length(s),grouped=FALSE,type.measure="mae");loo[tr[i]]<-as.numeric(predict(cvi,t(X[agenes,tr[i],drop=FALSE]),s="lambda.min"))}
  fit<-glmnet(xp,age,alpha=best$a,lambda=best$lam)
  te<-rownames(meta)[meta$group=="EGM"]
  pr<-setNames(as.numeric(predict(fit,t(X[agenes,te]))),te)
  ct_loo<-loo[meta[names(loo),"group"]=="CTRL"]
  p<-tryCatch(t.test(pr,ct_loo)$p.value,error=function(e)NA)
  res[[t]]<-data.table(tissue=t,n_agenes=length(agenes),alpha=best$a,LOO_MAE=round(mean(abs(loo-age)),2),LOO_R=round(cor(loo,age),3),clock_CTRL_LOO=round(mean(ct_loo),2),clock_EGM=round(mean(pr),2),EGM_minus_CTRL=round(mean(pr)-mean(ct_loo),2),p=round(p,4))
  pp[[t]]<-data.table(sample=tr,group=meta[tr,"group"],age=age,clock=round(loo,2),tissue=t)
  tp[[t]]<-data.table(sample=te,group="EGM",age=22,clock=round(pr,2),tissue=t)
  cat("done",t,"agenes",length(agenes),"CTRL_LOO",round(mean(ct_loo),2),"EGM",round(mean(pr),2),"EGM-CTRL",round(mean(pr)-mean(ct_loo),2),"p",round(p,4),"\n")
}
fwrite(rbindlist(res),"merged/native_clock_ymctrl_summary.csv")
fwrite(rbindlist(pp),"merged/native_clock_ymctrl_pubpred.csv")
fwrite(rbindlist(tp),"merged/native_clock_ymctrl_testpred.csv")
cat("\n===== train on public Y/M + our CTRL@22 (drop public O), project EGM =====\n");print(rbindlist(res))
