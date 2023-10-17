#############################################################################
# This script aims
#   - to create a ranking of kernel covariates per covariate and response
#############################################################################


# 5. Create a ranking of kernel covariate per covariate and per response (Step 2. from the Supplementary Material S.4)
#    This ranking:
#     - includes an analysis per type of smoothing basis function and considering different smoothing basis dimension
#     - is performed evaluating the AIC and the minimised generalized cross-validation (GCV) score after applying a simple
#     generalized additive model between each response and the kernel covariate
source("~/evaluate_knots.R")
rank.X1_4_PRES=apply.evaluate.knots(names(cov.kernel.df1),simulated.data,"pres",bs=c("bs","cr","tp"),nk=seq(5,20,by=2),family="bernoulli",link="log")
rank.X2_4_PRES=apply.evaluate.knots(names(cov.kernel.df2),simulated.data,"pres",bs=c("bs","cr","tp"),nk=seq(5,20,by=2),family="bernoulli",link="log")
rank.X3_4_PRES=apply.evaluate.knots(names(cov.kernel.df3),simulated.data,"pres",bs=c("bs","cr","tp"),nk=seq(5,20,by=2),family="bernoulli",link="log")
rank.X1_4_P.BIOM=apply.evaluate.knots(names(cov.kernel.df1),simulated.data,"pos.biom",bs=c("bs","cr","tp"),nk=seq(5,20,by=2),family="gamma",link="log")
rank.X2_4_P.BIOM=apply.evaluate.knots(names(cov.kernel.df2),simulated.data,"pos.biom",bs=c("bs","cr","tp"),nk=seq(5,20,by=2),family="gamma",link="log")
rank.X3_4_P.BIOM=apply.evaluate.knots(names(cov.kernel.df3),simulated.data,"pos.biom",bs=c("bs","cr","tp"),nk=seq(5,20,by=2),family="gamma",link="log")


# 6. Visualise the results from the previous step and select the kernel covariates to include in the model (Step 3 of Supplementary Material S.4)
head(rank.X1_4_PRES)
head(rank.X2_4_PRES)
head(rank.X3_4_PRES)
head(rank.X1_4_P.BIOM)
head(rank.X2_4_P.BIOM)
head(rank.X3_4_P.BIOM)

# 7. Save the rankings
save(rank.X1_4_PRES,rank.X2_4_PRES,rank.X3_4_PRES,rank.X1_4_P.BIOM,rank.X2_4_P.BIOM,rank.X3_4_P.BIOM,file="rankings.rds")

