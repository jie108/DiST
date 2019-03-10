################################################################################
### Part II: direction smoothing
################################################################################
rm(list=ls())
# need ncpu
ncpu=1

################
#### Source ####
################
source("dwi_fit.R")
source("dwi_basic.R")

#######################
#### load data ####
#######################
load("pre_estimation.Rdata") ##this will load vovel-wise estimation from running single_task_estimation.R

################################
#### Refit using updated v0 ####
################################
v.res <- v.smooth.bw(pre=pre, braingrid=braingrid, len=50,
                      range1=c(0.005,0.3), range2=c(0.005,0.3), cpus=ncpu,
                      xy.truncrange=100, z.truncrange=100, thres.v0.w=0.2, K=5,
                      cv.method="mad") ##mcv
v.obj <- update.v.obj(pre=pre, res=v.res)

##save results
save(v.obj, file="smoothing.Rdata") ##this will be used as inputs for fiber tracking 
