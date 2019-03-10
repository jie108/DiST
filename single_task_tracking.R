################################################################################
#### Simulation for testing smoothing methods over new prelimilary fittings ####
################################################################################
rm(list=ls())
# need ncpu
ncpu=1

################
#### Source ####
################
source("dwi_fit.R")
source("dwi_basic.R")
source("dwi_track.R")

#######################
#### load data ####
#######################
load("pre_estimation.Rdata") ##this will load vovel-wise estimation from running single_task_estimation.R
load("smoothing.Rdata") ##this will load direction smoothing results from running single_task_smoothing.R

################################
#### Refit using updated v0 ####
################################
track.pre <- v.track(v.obj=pre, xgrid.sp=0.01875, ygrid.sp=0.01875,
                     zgrid=0.01875, braingrid=braingrid, elim=T, nproj=1) ##tracking based on voxel-wise estimation

track <- v.track(v.obj=v.obj$obj, xgrid.sp=0.01875, ygrid.sp=0.01875,
                     zgrid=0.01875, braingrid=braingrid, elim=T, nproj=1) ##tracking based on smoothed  estimation (DiST-mcv)


##save results
save(track.pre,track, file="tracking.Rdata")

############################
#### fiber realizations ####
############################
library(rgl) 
library(compositions)

tobj <- track
length(tobj$sorted.iinds[tobj$sorted.update.ind])
ndis <- length(tobj$sorted.iinds[tobj$sorted.update.ind])  #length(tobj$sorted.iinds[tobj$sorted.update.ind])  # number of fibers

open3d()
for (iind in (tobj$sorted.iinds[tobj$sorted.update.ind])[1:ndis]){
  cat(iind,"\n")
  # plot
  plot.fib(tobj$tracks1[[iind]]$inloc, tobj$tracks1[[iind]]$dir)
  plot.fib(tobj$tracks2[[iind]]$inloc, tobj$tracks2[[iind]]$dir)
}
