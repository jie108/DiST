##prepared by Jie Peng: 2019-03; jiepeng@ucdavis.edu
##ref: Wong, R.K.W., T.C.M. Lee, D. Paul, and J. Peng. Fiber direction estimation, smoothing and tracking in diffusion MRI (with Discussion)  (2016). The Annals of Applied Statistics, 10(3): 1137-1156.

rm(list=ls())
library(R.matlab)
library(xtable)
library(dwi.internals2)  ## package by Raymond Wong with ArrayIndex function

#### load library ####
library(rgl) 
library(compositions)
source("dwi_basic.R")
source("dwi_fit.R")
source("dwi_track_realdata.R")


x_r = 108:123  ## ROI used for tractography algorithm
y_r = 124:139
z_r = 37:42
voxel_r = paste0('x',toString(x_r[1]),'-',toString(x_r[length(x_r)]),'y',toString(y_r[1]),'-',toString(y_r[length(y_r)]),'z',toString(z_r[1]),'-',toString(z_r[length(z_r)]))

x_subr = 1:16 ## region to draw the tracking results
y_subr = 1:16
z_subr = 1:6
voxel_sub = paste0('x',toString(x_subr[1]),'-',toString(x_subr[length(x_subr)]),'y',toString(y_subr[1]),'-',toString(y_subr[length(y_subr)]),'z',toString(z_subr[1]),'-',toString(z_subr[length(z_subr)]))

temp = readMat("track.mat")
temp$n.fiber2 = c(temp$n.fiber2)

## run the tracking algorithm
nproj = 1  ## skip one voxles before termination
our.track <- v.track(v.obj=temp, xgrid.sp=temp$xgrid.sp, ygrid.sp=temp$ygrid.sp,
                     zgrid.sp=temp$zgrid.sp, braingrid=array(temp$braingrid,dim=c(3,length(x_subr),length(y_subr),length(z_subr))), elim=T, nproj=nproj,
                     vorient=c(1,1,1), elim.thres=5.5)



############################
#### fiber realizations ####
############################
tobj <- our.track

length(tobj$sorted.iinds[tobj$sorted.update.ind])
ndis <- length(tobj$sorted.iinds[tobj$sorted.update.ind])  #length(tobj$sorted.iinds[tobj$sorted.update.ind])  # number of fibers

open3d()
for (iind in (tobj$sorted.iinds[tobj$sorted.update.ind])[1:ndis]){
  cat(iind,"\n")
  # plot
  plot.fib(tobj$tracks1[[iind]]$inloc, tobj$tracks1[[iind]]$dir)
  plot.fib(tobj$tracks2[[iind]]$inloc, tobj$tracks2[[iind]]$dir)
}
decorate3d(xlim=range(temp$braingrid[1,,,]), ylim=range(temp$braingrid[2,,,]),
           zlim=range(temp$braingrid[3,,,]), box=F)

