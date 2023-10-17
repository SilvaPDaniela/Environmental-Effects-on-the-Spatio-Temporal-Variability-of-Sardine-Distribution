#############################################################################
# This script has the r functions
#   - to obtain kernel covariates for observed data and evaluate 
#############################################################################


# Libraries
library(kedd)

# Request: All covariates should consider the name of the measured covariate, then "_before_", followed by the number of lag in days and finally by the term "days"
# Example: If the the covariate name is X1, the covariates measured for the entire time lag interval should be named as "X1_before_0days", "X1_before_1days", ..., "X1_before_28days"

# Arguments of the function:
# - data: data.frame containing the covariates
# - common_col_name: 'X1_before_' (following the previous example)
# - max.l: maximum for the l parameter (Equation (2) from the paper)

apply.kernel<-function(data,common_col_name,max.l){
  max.time.lag<-sum(substr(names(data),1,nchar(common_col_name))==as.character(common_col_name))
  avail.time.lags=0:(max.time.lag-1)
  lag.from.center<-c(0:max.l)
  df.par<-data.frame(center=rep(avail.time.lags,each=length(lag.from.center)),lag.from.center=rep(lag.from.center,length(avail.time.lags)))
  df.par=df.par[-which(df.par[,1]-df.par[,2]<0),]
  df.par=df.par[-which(df.par[,1]+df.par[,2]>(max.time.lag-1)),]
  data.cov=data[,which(substr(names(data),1,nchar(common_col_name))==as.character(common_col_name))]
  kernel.cov<-data.frame(matrix(NA,ncol=nrow(df.par),nrow=nrow(data)))
  pos_und=unlist(gregexpr('_', common_col_name))[1]
  names(kernel.cov)<-paste(substr(names(data.cov),1,pos_und),"center_",df.par[,"center"],"_lag_",df.par[,"lag.from.center"],sep="")
  for (k in 1:nrow(df.par)){
    if(df.par$lag.from.center[k]==0){
      cov<-data[,paste(common_col_name,df.par[k,"center"],"days",sep="")]    
    }
    if(df.par$lag.from.center[k]>0){
      seq.days<-seq(df.par[k,"center"]-df.par[k,"lag.from.center"],df.par[k,"center"]+df.par[k,"lag.from.center"],by=1)
      h<-h.mlcv(seq.days)$h
      vector.lag<-seq.days-df.par[k,"center"]
      weights<-as.numeric(sapply(vector.lag,function(x) exp(-(x^2)/(2*(h^2)))))/sqrt(2*pi)
      cov<-rep(0,nrow(data))
      for(i in 1:length(vector.lag)){
        cov<-cov+weights[i]*data[,paste(common_col_name,seq.days[i],"days",sep="")]
      }
    }
    kernel.cov[,k]<-cov
  }
  kernel.cov
} 
