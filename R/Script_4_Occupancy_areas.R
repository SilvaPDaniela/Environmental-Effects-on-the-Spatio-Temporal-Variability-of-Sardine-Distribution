#############################################################################
# This script aims
#   - to classify areas according to species occupancy
#############################################################################

# 1. Load the model result, the stack, and the mesh
load("~/Data/model.rds")
load("~/Data/stack_all.rds")
load("~/Data/mesh.rds")


# 2. Save the index for each specific stack
idx.z <- inla.stack.index(stk.all,'zobs')$data
idx.y <- inla.stack.index(stk.all,'yobs')$data
idx.zp <- inla.stack.index(stk.all,'zpred')$data
idx.yp <- inla.stack.index(stk.all,'ypred')$data

# 3. Save the mesh coordinates

coords.utm=mesh$loc[,1:2]


# 4. Create the project prediction grid
(nxy <- round(c(diff(c(min(coords.utm[,1]), max(coords.utm[,1]))), diff(c(min(coords.utm[,2]), max(coords.utm[,2])))))*4)
projgrid <- inla.mesh.projector(mesh, xlim=c(min(coords.utm[,1]-1) , max(coords.utm[,1]+1)),ylim=c(min(coords.utm[,2])-1, max(coords.utm[,2])+1), dims=nxy)


# 5. Identify the grid points that are out of the prediction area
coords.pol<-area@polygons[[1]]@Polygons[[1]]@coords
xy.in1<-inout(projgrid$lattice$loc,as.matrix(coords.pol))


# 4. Save the marginal distribution of the estimated W_st on the mesh points for the prediction time
zone.index=model$marginals.random$x

# 5. Save the posterior mean of the marginal standard deviation (sigma)
std.s=exp(model$summary.hyperpar$mean[2])

# 6. Define the scale parameter
scale=4/5

# 5. Compute the probabilities P(|W_st|>scale*sigma)
prob.zone <- lapply(zone.index, function(x) {1-inla.pmarginal(std.s*scale, x)+inla.pmarginal(-std.s*scale,x)})


# 6. Save the posterior median of the predicted values of biomass
xmean_z_t1 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.zp[(1):(nv)],1]) # mean of Z for t1
xmean_y_t1 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.yp[(1):(nv)],4]) # median of Y|Z for t1
xmean_z_t1[!xy.in1] <- NA
xmean_y_t1[!xy.in1] <- NA
xmean_z_t2 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.zp[(1+nv):(2*nv)],1]) # mean of Z for t2
xmean_y_t2 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.yp[(1+nv):(2*nv)],4]) # median of Y|Z for t2
xmean_z_t2[!xy.in1] <- NA
xmean_y_t2[!xy.in1] <- NA
xmean_z_t3 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.zp[(1+2*nv):(3*nv)],1]) # mean of Z for t3
xmean_y_t3 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.yp[(1+2*nv):(3*nv)],4]) # median of Y|Z for t3
xmean_z_t3[!xy.in1] <- NA
xmean_y_t3[!xy.in1] <- NA
xmean_z_t4 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.zp[(1+3*nv):(4*nv)],1]) # mean of Z for t4
xmean_y_t4 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.yp[(1+3*nv):(4*nv)],4]) # median of Y|Z for t4
xmean_z_t4[!xy.in1] <- NA
xmean_y_t4[!xy.in1] <- NA
xmean_z_t5 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.zp[(1+4*nv):(5*nv)],1]) # mean of Z for t5
xmean_y_t5 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.yp[(1+4*nv):(5*nv)],4]) # median of Y|Z for t5
xmean_z_t5[!xy.in1] <- NA
xmean_y_t5[!xy.in1] <- NA
xmean_z_t6 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.zp[(1+5*nv):(6*nv)],1]) # mean of Z for t6
xmean_y_t6 <- inla.mesh.project(projgrid, model$summary.fitted.values[idx.yp[(1+5*nv):(6*nv)],4]) # median of Y|Z for t6
xmean_z_t6[!xy.in1] <- NA
xmean_y_t6[!xy.in1] <- NA

predicted_biom_data<-data.frame(lon=projgrid$lattice$loc[,1],
                           lat=projgrid$lattice$loc[,2],
                           occurrence_t1=as.numeric(xmean_z_t1),
                           occurrence_t2=as.numeric(xmean_z_t2),
                           occurrence_t3=as.numeric(xmean_z_t3),
                           occurrence_t4=as.numeric(xmean_z_t4),
                           occurrence_t5=as.numeric(xmean_z_t5),
                           occurrence_t6=as.numeric(xmean_z_t6),
                           biomassa_t1=as.numeric(xmean_z_t1*xmean_y_t1),
                           biomassa_t2=as.numeric(xmean_z_t2*xmean_y_t2),
                           biomassa_t3=as.numeric(xmean_z_t3*xmean_y_t3),
                           biomassa_t4=as.numeric(xmean_z_t4*xmean_y_t4),
                           biomassa_t5=as.numeric(xmean_z_t5*xmean_y_t5),
                           biomassa_t6=as.numeric(xmean_z_t6*xmean_y_t6))


# 6. Save the uncertainty values for the grid points
uncert_t1 <- inla.mesh.project(projgrid, as.numeric(prob.zone[1:nv]))
uncert_t1[!xy.in1] <- NA
uncert_t2 <- inla.mesh.project(projgrid, as.numeric(prob.zone[(1+nv):(2*nv)]))
uncert_t2[!xy.in1] <- NA
uncert_t3 <- inla.mesh.project(projgrid, as.numeric(prob.zone[(1+(2*nv)):(3*nv)]))
uncert_t3[!xy.in1] <- NA
uncert_t4 <- inla.mesh.project(projgrid, as.numeric(prob.zone[(1+(3*nv)):(4*nv)]))
uncert_t4[!xy.in1] <- NA
uncert_t5 <- inla.mesh.project(projgrid, as.numeric(prob.zone[(1+(4*nv)):(5*nv)]))
uncert_t5[!xy.in1] <- NA
uncert_t6 <- inla.mesh.project(projgrid, as.numeric(prob.zone[(1+(5*nv)):(6*nv)]))
uncert_t6[!xy.in1] <- NA

predicted_unc_data<-data.frame(lon=projgrid$lattice$loc[,1],
                               lat=projgrid$lattice$loc[,2],
                               uncert_t1=as.numeric(uncert_t1),
                               uncert_t2=as.numeric(uncert_t2),
                               uncert_t3=as.numeric(uncert_t3),
                               uncert_t4=as.numeric(uncert_t4),
                               uncert_t5=as.numeric(uncert_t5),
                               uncert_t6=as.numeric(uncert_t6))


# 7. Joint the two databases
predicted_all=cbind(predicted_biom_data,predicted_unc_data[,-c(1:2)])


# 8. Compute and save the median per unit time
med.per.year<-apply(predicted_all[,9:14],2,function(x) median(x,na.rm=T))


# 9. Define the threshold for uncertainty (in our paper, it was defined as 0.5)
p.cut=0.5


# 10. Classify the grid according to occupancy using predicted values and uncertainty
#     - 2 : Rare
#     - 3 : Low occasional
#     - 4 : High occasional
#     - 5 : Preferred

occupancy_data<-data.frame(lon=predicted_all$lon,lat=predicted_all$lat)
occupancy_data$occupancy_t1=NA
occupancy_data$occupancy_t1[which(predicted_all$biomassa_t1<=med.per.year[1] & predicted_unc_data$uncert_t1<=p.cut)]=rep(2,sum(predicted_all$biomassa_t1<=med.per.year[1] & predicted_unc_data$uncert_t1<=p.cut,na.rm=T))
occupancy_data$occupancy_t1[which(predicted_all$biomassa_t1<=med.per.year[1] & predicted_unc_data$uncert_t1>p.cut)]=rep(3,sum(predicted_all$biomassa_t1<=med.per.year[1] & predicted_unc_data$uncert_t1>p.cut,na.rm=T))
occupancy_data$occupancy_t1[which(predicted_all$biomassa_t1>med.per.year[1] & predicted_unc_data$uncert_t1<=p.cut)]=rep(5,sum(predicted_all$biomassa_t1>med.per.year[1] & predicted_unc_data$uncert_t1<=p.cut,na.rm=T))
occupancy_data$occupancy_t1[which(predicted_all$biomassa_t1>med.per.year[1] & predicted_unc_data$uncert_t1>p.cut)]=rep(4,sum(predicted_all$biomassa_t1>med.per.year[1] & predicted_unc_data$uncert_t1>p.cut,na.rm=T))
occupancy_data$occupancy_t2=NA
occupancy_data$occupancy_t2[which(predicted_all$biomassa_t2<=med.per.year[2] & predicted_unc_data$uncert_t2<=p.cut)]=rep(2,sum(predicted_all$biomassa_t2<=med.per.year[2] & predicted_unc_data$uncert_t2<=p.cut,na.rm=T))
occupancy_data$occupancy_t2[which(predicted_all$biomassa_t2<=med.per.year[2] & predicted_unc_data$uncert_t2>p.cut)]=rep(3,sum(predicted_all$biomassa_t2<=med.per.year[2] & predicted_unc_data$uncert_t2>p.cut,na.rm=T))
occupancy_data$occupancy_t2[which(predicted_all$biomassa_t2>med.per.year[2] & predicted_unc_data$uncert_t2<=p.cut)]=rep(5,sum(predicted_all$biomassa_t2>med.per.year[2] & predicted_unc_data$uncert_t2<=p.cut,na.rm=T))
occupancy_data$occupancy_t2[which(predicted_all$biomassa_t2>med.per.year[2] & predicted_unc_data$uncert_t2>p.cut)]=rep(4,sum(predicted_all$biomassa_t2>med.per.year[2] & predicted_unc_data$uncert_t2>p.cut,na.rm=T))
occupancy_data$occupancy_t3=NA
occupancy_data$occupancy_t3[which(predicted_all$biomassa_t3<=med.per.year[3] & predicted_unc_data$uncert_t3<=p.cut)]=rep(2,sum(predicted_all$biomassa_t3<=med.per.year[3] & predicted_unc_data$uncert_t3<=p.cut,na.rm=T))
occupancy_data$occupancy_t3[which(predicted_all$biomassa_t3<=med.per.year[3] & predicted_unc_data$uncert_t3>p.cut)]=rep(3,sum(predicted_all$biomassa_t3<=med.per.year[3] & predicted_unc_data$uncert_t3>p.cut,na.rm=T))
occupancy_data$occupancy_t3[which(predicted_all$biomassa_t3>med.per.year[3] & predicted_unc_data$uncert_t3<=p.cut)]=rep(5,sum(predicted_all$biomassa_t3>med.per.year[3] & predicted_unc_data$uncert_t3<=p.cut,na.rm=T))
occupancy_data$occupancy_t3[which(predicted_all$biomassa_t3>med.per.year[3] & predicted_unc_data$uncert_t3>p.cut)]=rep(4,sum(predicted_all$biomassa_t3>med.per.year[3] & predicted_unc_data$uncert_t3>p.cut,na.rm=T))
occupancy_data$occupancy_t4=NA
occupancy_data$occupancy_t4[which(predicted_all$biomassa_t4<=med.per.year[4] & predicted_unc_data$uncert_t4<=p.cut)]=rep(2,sum(predicted_all$biomassa_t4<=med.per.year[4] & predicted_unc_data$uncert_t4<=p.cut,na.rm=T))
occupancy_data$occupancy_t4[which(predicted_all$biomassa_t4<=med.per.year[4] & predicted_unc_data$uncert_t4>p.cut)]=rep(3,sum(predicted_all$biomassa_t4<=med.per.year[4] & predicted_unc_data$uncert_t4>p.cut,na.rm=T))
occupancy_data$occupancy_t4[which(predicted_all$biomassa_t4>med.per.year[4] & predicted_unc_data$uncert_t4<=p.cut)]=rep(5,sum(predicted_all$biomassa_t4>med.per.year[4] & predicted_unc_data$uncert_t4<=p.cut,na.rm=T))
occupancy_data$occupancy_t4[which(predicted_all$biomassa_t4>med.per.year[4] & predicted_unc_data$uncert_t4>p.cut)]=rep(4,sum(predicted_all$biomassa_t4>med.per.year[4] & predicted_unc_data$uncert_t4>p.cut,na.rm=T))
occupancy_data$occupancy_t5=NA
occupancy_data$occupancy_t5[which(predicted_all$biomassa_t5<=med.per.year[5] & predicted_unc_data$uncert_t5<=p.cut)]=rep(2,sum(predicted_all$biomassa_t5<=med.per.year[5] & predicted_unc_data$uncert_t5<=p.cut,na.rm=T))
occupancy_data$occupancy_t5[which(predicted_all$biomassa_t5<=med.per.year[5] & predicted_unc_data$uncert_t5>p.cut)]=rep(3,sum(predicted_all$biomassa_t5<=med.per.year[5] & predicted_unc_data$uncert_t5>p.cut,na.rm=T))
occupancy_data$occupancy_t5[which(predicted_all$biomassa_t5>med.per.year[5] & predicted_unc_data$uncert_t5<=p.cut)]=rep(5,sum(predicted_all$biomassa_t5>med.per.year[5] & predicted_unc_data$uncert_t5<=p.cut,na.rm=T))
occupancy_data$occupancy_t5[which(predicted_all$biomassa_t5>med.per.year[5] & predicted_unc_data$uncert_t5>p.cut)]=rep(4,sum(predicted_all$biomassa_t5>med.per.year[5] & predicted_unc_data$uncert_t5>p.cut,na.rm=T))
occupancy_data$occupancy_t6=NA
occupancy_data$occupancy_t6[which(predicted_all$biomassa_t6<=med.per.year[6] & predicted_unc_data$uncert_t6<=p.cut)]=rep(2,sum(predicted_all$biomassa_t6<=med.per.year[6] & predicted_unc_data$uncert_t6<=p.cut,na.rm=T))
occupancy_data$occupancy_t6[which(predicted_all$biomassa_t6<=med.per.year[6] & predicted_unc_data$uncert_t6>p.cut)]=rep(3,sum(predicted_all$biomassa_t6<=med.per.year[6] & predicted_unc_data$uncert_t6>p.cut,na.rm=T))
occupancy_data$occupancy_t6[which(predicted_all$biomassa_t6>med.per.year[6] & predicted_unc_data$uncert_t6<=p.cut)]=rep(5,sum(predicted_all$biomassa_t6>med.per.year[6] & predicted_unc_data$uncert_t6<=p.cut,na.rm=T))
occupancy_data$occupancy_t6[which(predicted_all$biomassa_t6>med.per.year[6] & predicted_unc_data$uncert_t6>p.cut)]=rep(4,sum(predicted_all$biomassa_t6>med.per.year[6] & predicted_unc_data$uncert_t6>p.cut,na.rm=T))



# 11. Map occupancy areas
layers <- c(1:6)
s2 <- list()

for (i in layers) {
  spg <- data.frame(x=occupancy_data$lon,y=occupancy_data$lat,Occupacy=occupancy_data[,2+i])
  coordinates(spg) <- ~ x + y
  gridded(spg) <- TRUE
  r <- raster(spg)
  r <- ratify(r)
  ## Levels
  rat <- levels(r)[[1]]
  lev<-c("Unfavourable", "Low occasional", "High occasional","Recurrent")
  if(length(setdiff(2:5,sort(unique(spg$Occupacy))))>0){
    rat$landcover<-lev[-setdiff(2:5,sort(unique(spg$Occupacy)))]
  }
  if(length(setdiff(2:5,sort(unique(spg$Occupacy))))==0){
    rat$landcover<-lev
  }
  rat.lev=data.frame(ID=2:5,landcover=c("Unfavourable", "Low occasional", "High occasional","Recurrent"))
  levels(r) <- rat.lev
  s2[[i]] <- r
  print(i)
}
s3<-stack(s2[[1]],s2[[2]],s2[[3]],s2[[4]],s2[[5]],s2[[6]])

levelplot(s3,xlab="", ylab="",col.regions=brewer.pal(n = 4, name = "BrBG"),
          names.attr=c("2017","2018","2019","2020","2021","2022"),
          par.strip.text=list(cex=1.5, lines=1))



# 12. Define according to persistent areas (considering the last of three years)
#     The idea was to maintain the preffered and rare areas, naming as favorable and unfavorable respectively
#      - 2 : Unfavorable
#      - 5 : Favorable
#      - 1 : Without persistence in rare and preferred areas

area_study<-occupancy_data
area_study$reg=1
area_study=area_study[!is.na(area_study$occupancy_t1),]

Fav.unf.zone<-occupancy_data
Fav.unf.zone$Fav.unf[!is.na(Fav.unf.zone$occupancy_t6)]<-rep(1,sum(!is.na(Fav.unf.zone$occupancy_t6)))

Fav.unf.zone[which(Fav.unf.zone$occupancy_t6==Fav.unf.zone$occupancy_t5 & Fav.unf.zone$occupancy_t6==2),ncol(Fav.unf.zone)]<-2
Fav.unf.zone[which(Fav.unf.zone$occupancy_t6==Fav.unf.zone$occupancy_t5 & Fav.unf.zone$occupancy_t6==5),ncol(Fav.unf.zone)]<-5
Fav.unf.zone[which(Fav.unf.zone[,ncol(Fav.unf.zone)]!=Fav.unf.zone$occupancy_t4),ncol(Fav.unf.zone)]<-1
Fav.unf.zone=Fav.unf.zone[,c(1:2,9)]




# 13. Map persistent areas

coordinates(Fav.unf.zone) <- ~ lon + lat
gridded(Fav.unf.zone) <- TRUE
r.fav_unf <- raster(Fav.unf.zone)
r.fav_unf <- ratify(r.fav_unf)
## Levels
rat <- levels(r.fav_unf)[[1]]
lev<-c("Non-persistence", "Favourable")
rat.lev=data.frame(ID=c(1,5),landcover=c("Non-persistence", "Favourable"))
levels(r.fav_unf) <- rat.lev

levelplot(r.fav_unf,xlab="", ylab="",col.regions=c("gray86",brewer.pal(n = 4, name = "BrBG")[4]),
          title="2020-2022",
          par.strip.text=list(cex=1.5, lines=1))

