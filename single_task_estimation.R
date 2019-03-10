################################################################################
### Part I: direction estimation
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
#### Simulate data ####
#######################
source("sim-curve-new.R")

####################################
#### fit voxel level estimation ####
####################################
pre <- v.est(dwi.obs=dwi.obs, sigma=sigma, S0=S0const, b=1, grad.mat=grad.mat,
             braingrid=braingrid, cpus=ncpu, opt=list())

##save results
save(pre,braingrid, file="pre_estimation.Rdata") ##this will be used as inputs for the smoothing step.