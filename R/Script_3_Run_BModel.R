#############################################################################
# This script aims
#   - run the two-part spatio-temporal model considering the Barrier approach  
#############################################################################


# Libraries
library(INLA)
library(mgcv)

# 1. Read the data with kernel covariates and the evaluation of them
source("~/obtain_kernel_covariates_obser_data.R")

# 1. Load the polygon representing the study region
load("~/Data/area.rds")

# 2. Load the polygon representing the barrier
load("~/Data/barrier.rds")

# 3. Create year ID covariate. It must be from 1 to the dimension of period time. In our case, it should be from 1 to 6
simulated.data$year=as.numeric(substr(simulated.data$date,1,4))
simulated.data$year.id=simulated.data$year-2016

# 4. Create the mesh
segment=inla.sp2segment(barrier)
max.edge = 1
bound.outer = 1.5
mesh = inla.mesh.2d(boundary = area,
                    loc=cbind(simulated.data$lon, simulated.data$lat),
                    max.edge = c(1,5)*max.edge,
                    cutoff = max.edge/5,
                    offset = c(max.edge, bound.outer))


# 5. Construct observation weight matrix
#    group identifies the corresponding year id for each observation
#    n.group identifies the dimension of time period (in our case, 6 years)
Ast <- inla.spde.make.A(mesh = mesh, loc = as.matrix(simulated.data[,c("lon","lat")]),
                        group = simulated.data$year.id,
                        n.group = 6)


# 6. Specification of the spatio-temporal random effect which account for barriers

# - Number of triangles of the mesh
tl = length(mesh$graph$tv[,1])

# - Initialize a matrix containing the central coordinates of each triangle's 
posTri = matrix(0, tl, 2)
for (t in 1:tl){
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}
# - Convert to Spatial Points
posTri = SpatialPoints(posTri)

# - Compute the triangle positions
normal = over(barrier, posTri, returnList=T)

# - Checking which mesh triangles are inside the normal area
normal = unlist(normal)
barrier.triangles = setdiff(1:tl, normal)
poly.barrier = inla.barrier.polygon(mesh, barrier.triangles)

# - Define the priors for the hyperparameters
spde.barrier = inla.barrier.pcmatern(mesh, barrier.triangles = barrier.triangles, prior.range = c(50, 0.05),prior.sigma = c(1,0.05))
m <- spde.barrier$n.spde



# 7. Generates a list of named index vectors for an SPDE mode
field.z.idx <- inla.spde.make.index(name = 'x',n.spde = spde.barrier$f$n, n.group = 6)
field.y.idx <- inla.spde.make.index(name = 'u',n.spde = spde.barrier$f$n, n.group = 6)



# 8. Insert smoothing effects
# source("preliminary_eval_kernel_covariates.R")
# For presence and
# - For X1 covariate, it should be used the first covariate in the ranking:
# > X1_center_12_lag_0
# > nk = 7
# > bs = 'bs'
sm.X1<-smoothCon(s(X1_center_12_lag_0,bs="bs",k=7,fx=TRUE),data=simulated.data,absorb.cons=TRUE)[[1]]
X1.bs<-sm.X1$X
X1.df<-data.frame(X1.bs)
names(X1.df)<-paste("X1.bs",1:ncol(X1.bs),sep="")
lcs.X1<-inla.make.lincombs(X1.df)
names(lcs.X1)<-paste(names(lcs.X1),"X1",sep="")

# - For X2 covariate, it should be used:
# > X2_center_25_lag_3
# > nk = 9
# > bs = 'tp'
sm.X2<-smoothCon(s(X2_center_25_lag_3,bs="tp",k=9,fx=TRUE),data=simulated.data,absorb.cons=TRUE)[[1]]
X2.bs<-sm.X2$X
X2.df<-data.frame(X2.bs)
names(X2.df)<-paste("X2.bs",1:ncol(X2.bs),sep="")
lcs.X2<-inla.make.lincombs(X2.df)
names(lcs.X2)<-paste(names(lcs.X2),"X2",sep="")

# - For X3 covariate, it should be used:
# > X3_center_13_lag_1
# > nk = 9
# > bs = 'bs'
sm.X3<-smoothCon(s(X3_center_13_lag_1,bs="bs",k=9,fx=TRUE),data=simulated.data,absorb.cons=TRUE)[[1]]
X3.bs<-sm.X3$X
X3.df<-data.frame(X3.bs)
names(X3.df)<-paste("X3.bs",1:ncol(X3.bs),sep="")
lcs.X3<-inla.make.lincombs(X3.df)
names(lcs.X3)<-paste(names(lcs.X3),"X3",sep="")


# For biomass under presence and
# - For X1 covariate, it should be used:
# > X1_center_19_lag_0
# > nk = 5
# > bs = 'cr'
sm.X1.Y<-smoothCon(s(X1_center_19_lag_0,bs="cr",k=5,fx=TRUE),data=simulated.data,absorb.cons=TRUE)[[1]]
X1.bs.Y<-sm.X1.Y$X
X1.Y.df<-data.frame(X1.bs.Y)
names(X1.Y.df)<-paste("X1.bs.Y",1:ncol(X1.bs.Y),sep="")
lcs.X1.Y<-inla.make.lincombs(X1.Y.df)
names(lcs.X1.Y)<-paste(names(lcs.X1.Y),"X1.Y",sep="")

# - For X2 covariate, it should be used:
# > X2_center_28_lag_0
# > nk = 7
# > bs = 'tp'
sm.X2.Y<-smoothCon(s(X2_center_28_lag_0,bs="tp",k=7,fx=TRUE),data=simulated.data,absorb.cons=TRUE)[[1]]
X2.bs.Y<-sm.X2.Y$X
X2.Y.df<-data.frame(X2.bs.Y)
names(X2.Y.df)<-paste("X2.bs.Y",1:ncol(X2.bs.Y),sep="")
lcs.X2.Y<-inla.make.lincombs(X2.Y.df)
names(lcs.X2.Y)<-paste(names(lcs.X2.Y),"X2.Y",sep="")

# - For X3 covariate, it should be used:
# > X3_center_8_lag_0
# > nk = 5
# > bs = 'cr'
sm.X3.Y<-smoothCon(s(X3_center_8_lag_0,bs="cr",k=5,fx=TRUE),data=simulated.data,absorb.cons=TRUE)[[1]]
X3.bs.Y<-sm.X3.Y$X
X3.Y.df<-data.frame(X3.bs.Y)
names(X3.Y.df)<-paste("X3.bs.Y",1:ncol(X3.bs.Y),sep="")
lcs.X3.Y<-inla.make.lincombs(X3.Y.df)
names(lcs.X3.Y)<-paste(names(lcs.X3.Y),"X3.Y",sep="")



# 9. Combining data, effects and observation matrices

stk.z <- inla.stack(data=list(y=cbind(as.vector(simulated.data$pres), NA), link=1), 
                    A=list(Ast, 1),
                    effects=list(c(field.z.idx),list(a0=rep(1, length=length(simulated.data$pres)),
                                                     z.year=simulated.data$year.id,
                                                     z.X1=X1.df,
                                                     z.X2=X2.df,
                                                     z.X3=X3.df)), tag="zobs")

stk.y <- inla.stack(data=list(y=cbind( NA,as.vector(simulated.data$pos.biom)), link=2), 
                    A=list(Ast, 1),
                    effects=list(c(field.y.idx),list(b0=rep(1, length=length(simulated.data$pos.biom)),
                                                     y.year=simulated.data$year.id,
                                                     y.X1=X1.Y.df,
                                                     y.X2=X2.Y.df,
                                                     y.X3=X3.Y.df)), tag="yobs")

# 10. Read mesh data
load("~/Data/data_mesh.rds")
source("~/obtain_kernel_covariates_mesh.R")


# 11. Insert smoothing effects for mesh
XpX1 <- PredictMat(sm.X1, data.mesh)
XpX1.df<-data.frame(XpX1)
names(XpX1.df)<-paste("X1.bs",1:ncol(XpX1.df),sep="")

XpX1.y <- PredictMat(sm.X1.Y, data.mesh)
XpX1.y.df<-data.frame(XpX1.y)
names(XpX1.y.df)<-paste("X1.bs.Y",1:ncol(XpX1.y.df),sep="")

XpX2 <- PredictMat(sm.X2, data.mesh)
XpX2.df<-data.frame(XpX2)
names(XpX2.df)<-paste("X2.bs",1:ncol(XpX2.df),sep="")

XpX2.y <- PredictMat(sm.X2.Y, data.mesh)
XpX2.y.df<-data.frame(XpX2.y)
names(XpX2.y.df)<-paste("X2.bs.Y",1:ncol(XpX2.y.df),sep="")

XpX3 <- PredictMat(sm.X3, data.mesh)
XpX3.df<-data.frame(XpX3)
names(XpX3.df)<-paste("X3.bs",1:ncol(XpX3.df),sep="")

XpX3.y <- PredictMat(sm.X3.Y, data.mesh)
XpX3.y.df<-data.frame(XpX3.y)
names(XpX3.y.df)<-paste("X3.bs.Y",1:ncol(XpX3.y.df),sep="")


# 12. Create year ID  for mesh data
data.mesh$year=as.numeric(substr(data.mesh$date,1,4))
data.mesh$year.id=data.mesh$year-2016


# 13. Construct prediction weight matrix
A_pred <- inla.spde.make.A(mesh=mesh,loc=as.matrix(data.mesh[,1:2]),
                           group=data.mesh$year.id, 
                           n.group=6)


# 14. Combining data, effects and prediction matrices
stk.zp <- inla.stack(
  data = list(Y = matrix(NA, ncol(A_pred), 2), link = 1), 
  effects=list(field.z.idx,list(a0=rep(1, length=nrow(data.mesh)),
                                z.year=data.mesh$year.id,
                                z.X1=XpX1.df,
                                z.X2=XpX2.df,
                                z.X3=XpX3.df)),
  A = list(A_pred, 1),
  tag = 'zpred') 

stk.yp <- inla.stack(
  data = list(Y = matrix(NA, ncol(A_pred), 2), link = 2), 
  A = list(A_pred, 1),
  effects=list(field.y.idx,list(b0=rep(1, length=nrow(data.mesh)),
                    y.year=data.mesh$year.id,
                    y.X1=XpX1.y.df,
                    y.X2=XpX2.y.df,
                    y.X3=XpX3.y.df)), 
  tag = 'ypred')


# 15. Combining observation stack and prediction stack
stk.all <- inla.stack(stk.z, stk.y, stk.zp, stk.yp)


# 16. Concatenate the linear combinations 
Allcs<-c(lcs.X1,lcs.X2,lcs.X3,lcs.X1.Y,lcs.X2.Y,lcs.X3.Y)


# 17. Write the inla model formula
names(XpX1.df)
names(XpX2.df)
names(XpX3.df)
names(XpX1.y.df)
names(XpX2.y.df)
names(XpX3.y.df)
n=nrow(simulated.data)
nv=mesh$n

formula<- y ~ 0 + a0  + b0+ # intercepts
  X1.bs1+X1.bs2+X1.bs3+X1.bs4+X1.bs5+X1.bs6+ # X1 kernel covariate for Z process
  X2.bs1+X2.bs2+X2.bs3+X2.bs4+X2.bs5+X2.bs6+X2.bs7+X2.bs8+ # X2 kernel covariate for Z process
  X3.bs1+X3.bs2+X3.bs3+X3.bs4+X3.bs5+X3.bs6+X3.bs7+X3.bs8+ # X3 kernel covariate for Z process
  X1.bs.Y1+X1.bs.Y2+X1.bs.Y3+X1.bs.Y4+ # X1 kernel covariate for Y|(Z=1) process
  X2.bs.Y1+X2.bs.Y2+X2.bs.Y3+X2.bs.Y4+X2.bs.Y5+X2.bs.Y6+ # X2 kernel covariate for Y|(Z=1) process
  X3.bs.Y1+X3.bs.Y2+X3.bs.Y3+X3.bs.Y4+ # X3 kernel covariate for Y|(Z=1) process
  f(u, model = spde.barrier, group = u.group, control.group = list(model = 'ar1')) + # W_st for Y|(Z=1) process
  f(x, copy = "u", fixed = FALSE, group = x.group,control.group = list(model = 'ar1'))+ # W_st for Z process
  f(y.year,model="iid")+ # unstructured annual effect for Y|(Z=1) process
  f(z.year,model="iid") # unstructured annual effect for Z process


# 18. Run the inla model
model<-inla(formula, family=c('binomial', 'gamma'),
            control.family=list(list(link="logit"),list(link="log")),
            data = inla.stack.data(stk.all),
            lincomb=Allcs,
            control.predictor = list(A = inla.stack.A(stk.all),compute=TRUE,link=c(rep(1,n),rep(2,n),rep(1,nv*6),rep(2,nv*6))),
            control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
            control.inla = list(strategy = 'adaptive'), verbose=TRUE)


# 19. Save the model
#     More combinations of covariates should be considered.
#     For instance, removing some kernel covariates, adding other covariates and/or removing unstructured effects
save(model,file="model.rds")
save(stk.all,file="stack.rds")
save(mesh,file="mesh.rds")
