#############################################################################
# This script aims
#   - to obtain kernel covariates for observed data and evaluate 
#   - to create a ranking of kernel covariates per covariate and response
#############################################################################


# 1. Read the data
#    This step includes obtaining the covariates measured during a range of time according to some prior knowledge (Step 1. from the Supplementary Material S.4)
setwd("C:/Users/danye/OneDrive/Ambiente de Trabalho/PDMA/2º ano/Tese/Campanhas Sardinha 2000-2020/Github")
load("simulated_data2.rds")

# 2. Create the kernel covariate for each covariate measured during an past time interval
#    Apply the 'apply.kernel' for each covariate.
#    In the simulated dataset, there are 3 covariates 'X1', 'X2' and 'X3'
source("~/Script_1_Create_kernel_covariates.R")
cov.kernel.df1=apply.kernel(simulated.data,"X1_before_",7)
cov.kernel.df2=apply.kernel(simulated.data,"X2_before_",7)
cov.kernel.df3=apply.kernel(simulated.data,"X3_before_",7)

# 3. Attach the kernel covariates to the original dataset
simulated.data=cbind(simulated.data,cov.kernel.df1)
simulated.data=cbind(simulated.data,cov.kernel.df2)
simulated.data=cbind(simulated.data,cov.kernel.df3)

# 4. Create the responses for the two processes: presence-absense (named as 'pres') and strictly positive biomass ('pos.biom')
simulated.data$pres=ifelse(simulated.data$Y>0,1,0)
simulated.data$pos.biom=ifelse(simulated.data$Y>0,simulated.data$Y,NA)

