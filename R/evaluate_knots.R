#############################################################################
# This script has the r functions
#   - to create a ranking of kernel covariates per covariate and response
#############################################################################


# evaluate.knots: function to obtain the the best combination between 'nk' and 'bs' for each kernel covariate according to
#                 the AIC and the minimised generalized cross-validation (GCV) score
# - data: database with response and covariates
# - y: a vector indicating the response
# - x: a vector indicating the kernel covariate
# - bs: a two letter character string indicating the (penalized) smoothing basis to use
# - nk: the dimension of the basis used to represent the smooth term
# - family: family distribution ('bernoulli','gamma','poisson') [you can include more distributions by adding other if loops]
# - link: link function, 'log' by default, but it is not used in the case of 'bernoulli' family

evaluate.knots<-function(data,y.name,cov.kernel.name,bs,nk,family,link="log"){
  cov.NK<-data.frame()
  for(j in 1:length(bs)){
    for(i in 1:length(nk)){
      if(family=="bernoulli"){
        if(length(unique(data[,y.name]))==2){
          y=as.numeric(I(data[,y.name]>0))
          x=as.numeric(data[,cov.kernel.name])
          gam.model<-gam(y ~ s(x,bs=bs[j],k=nk[i],fx=TRUE), family = binomial)
        }
        if(length(unique(data[,y.name]))!=2){
          stop("y must have 2 unique levels")
        }
      }
      if(family=="gamma"){
        y=as.numeric(data[,y.name])
        x=as.numeric(data[,cov.kernel.name])
        gam.model<-gam(y ~ s(x,bs=bs[j],k=nk[i],fx=TRUE) ,family = Gamma(link=link))
      }
      if(family=="poisson"){
        y=as.numeric(data[,y.name])
        x=as.numeric(data[,cov.kernel.name])
        gam.model<-gam(y ~ s(x,bs=bs[j],k=nk[i],fx=TRUE) ,family = Poisson(link=link))
      }
      res<-data.frame(nk=nk[i],AIC=extractAIC(gam.model)[2],bs=bs[j],GCV=gam.model$gcv.ubre)
      cov.NK<-rbind(cov.NK,res)
    }
  }
  cov.NK$comp.GCV<-c(NA,cov.NK$GCV[1:(nrow(cov.NK)-1)])
  cov.NK$more<-cov.NK$GCV>cov.NK$comp.GCV
  cov.NK$more[which(cov.NK$nk==min(nk))]<-rep(NA,sum(cov.NK$nk==min(nk)))
  #determine the optimum number of knots and the effect
  res.df1<-data.frame()
  for (k in unique(cov.NK$bs)){
    first.TRUE<-(which(cov.NK[cov.NK$bs==k,"more"]==TRUE))[1]
    res.df1<-rbind(res.df1,cov.NK[cov.NK$bs==k,][first.TRUE-1,])
    res.df1<-res.df1[which.min(res.df1$AIC),]
  }
  res.df1[,c("nk","AIC","bs")]
}

# apply.evaluate.knots: allows to apply the evaluate.knots for all kernel covariates

apply.evaluate.knots=function(names.cov.kernel,data,y.name,bs,nk,family,link="log"){
  res.df=sapply(names.cov.kernel,function(x) evaluate.knots(data,y.name,x,bs,nk,family,link=link))
  res.df1=data.frame(cov=colnames(res.df),nk=as.numeric(res.df[1,]),AIC=as.numeric(res.df[2,]),bs=as.character(res.df[3,]))
  res.df1[order(res.df1$AIC),]
}
