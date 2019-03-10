# Example of two curve crossing
# SNR =20, FA=0.9

################
#### Source ####
################
source("dwi_fit.R")
source("dwi_basic.R")

##################
#### Settings ####
##################
# general
braindim <- c(16,16,4)

# for Rician noise model
sigma <- 0.5 # noise level

# for signal strength of DWI
S0const <- 10
reps <- 1

###################################
#### Generate gradient vectors ####
###################################
# generate gradient vectors (subj to transpose)
basis.verts <- tess.ico(2,initial="oct", random=F)$verts
basis.verts <- t(basis.verts) ## important

# generate the gradient matrix with replications
grad.mat <- matrix(rep(basis.verts,each=reps),nc=3)
######################
#### true tensors ####
######################
n.layer <- 2 # the maximum number of fiber in a voxel

n.voxel <- prod(braindim)
tensor.true <- array(0,dim=c(3,3,n.layer,n.voxel))
p.true <- array(0,dim=c(n.layer,n.voxel))
tensor.true.array <- array(0,dim=c(3,3,n.layer,braindim[1],braindim[2],braindim[3]))
p.true.array <- array(0,dim=c(n.layer,braindim[1],braindim[2],braindim[3]))
n.fiber <- array(0,dim=n.voxel)
v0.true <- array(dim=c(3, n.layer, n.voxel))

### laying down weights
p.true <- array(rep(c(1,rep(0,n.layer-1)),n.voxel),dim=dim(p.true))
p.true.array <- array(rep(c(1,rep(0,n.layer-1)),n.voxel),dim=dim(p.true.array))

### laying down the undirectional tensors
tensor.null <- diag(rep(1,3))
for (i in (1:braindim[1])){
  for (j in (1:braindim[2])){
    for (k in (1:braindim[3])){
      for (l in (1:n.layer)){
        voxelindex <- (k-1)*braindim[1]*braindim[2] + (j-1)*braindim[1] + i
        tensor.true[,,l,voxelindex] <- tensor.null
        tensor.true.array[,,l,i,j,k] <- tensor.null
      }
    }
  }
}

### laying down the first fiber
fib1.fun <- function(loc){
  width <- 0.15
  # using polar coordinate
  r <- sqrt(sum(loc[1:2]^2))
  ind <- (r>=(0.8-width))*(r<=(0.8+ width))*(loc[3]<=(0.5+width))*(loc[3]>=(0.5-width))
  # tangent
  if (ind){
    theta <- atan2(loc[2], loc[1]) + pi/2
    xy <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nr=2)%*%c(1,0)
    v <- c(as.vector(xy), 0); v <- v/sqrt(sum(v^2))
    tens <- get.sym.tens(lambdas=c(4, 0.3694233, 0.3694233), v1=v)
  } else {
    v <- NULL
    tens <- NULL
  }
  return(list(ind=as.logical(ind), v=v, tens=tens))
}

fib1.res <- fill.fiber(fun=fib1.fun, p.fib=0.7, braindim=braindim, tensor=tensor.true, tensor.array=tensor.true.array, p=p.true, p.array=p.true.array, v0=v0.true, n.fiber=n.fiber, final=F)

tensor.true <- fib1.res$tensor
tensor.true.array <- fib1.res$tensor.array
n.fiber <- fib1.res$n.fiber
p.true <- fib1.res$p
p.true.array <- fib1.res$p.array
v0.true <- fib1.res$v0

#### laying down the second fiber ####
fib2.fun <- function(loc){
  width <- 0.15
  # using polar coordinate
  r <- sqrt(sum((loc[1:2]-c(1,0))^2))
  ind <- (r>=(0.8-width))*(r<=(0.8+ width))*(loc[3]<=(0.5+width))*(loc[3]>=(0.5-width))
  # tangent
  if (ind){
    theta <- atan2(loc[2], loc[1]-1) + pi/2
    xy <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nr=2)%*%c(1,0)
    v <- c(as.vector(xy), 0); v <- v/sqrt(sum(v^2))
    tens <- get.sym.tens(lambdas=c(4, 0.3694233, 0.3694233), v1=v)
  } else {
    v <- NULL
    tens <- NULL
  }
  return(list(ind=as.logical(ind), v=v, tens=tens))
}

fib2.res <- fill.fiber(fun=fib2.fun, p.fib=0.3, braindim=braindim, tensor=tensor.true, tensor.array=tensor.true.array, p=p.true, p.array=p.true.array, v0=v0.true, n.fiber=n.fiber, final=T)

tensor.true <- fib2.res$tensor
tensor.true.array <- fib2.res$tensor.array
n.fiber <- fib2.res$n.fiber
p.true <- fib2.res$p
p.true.array <- fib2.res$p.array
v0.true <- fib2.res$v0

#plot.brain.ten(tensor.true, ps=p.true,braingrid=braingrid, n.layer=n.layer, n.fiber=n.fiber, minbraingrid=0.02)

###########################
#### tensor simulation ####
###########################

dwi.obs <- sim.mix.SB.whole(tensor.whole=tensor.true, p.whole=p.true, grad.mat=grad.mat, braindim=braindim, sigma=sigma, S0=S0const)


##############################################
#### creation of a three-dimensional grid ####
##############################################
# that reflects real distances between locations
xgrid.sp <- 0.01875 # voxel sidelength in x direction
ygrid.sp <- 0.01875 # voxel sidelength in x direction
zgrid.sp <- 0.01875 # voxel sidelength in x direction


braingrid <- array(0,dim=c(3,braindim))
for (i in (1:braindim[1])){
  for (j in (1:braindim[2])){
    for (k in (1:braindim[3])){
      braingrid[,i,j,k] <- c((i-braindim[1]/2-.5)*xgrid.sp, (j-braindim[2]/2-.5)*ygrid.sp, (k-braindim[3]/2-.5)*zgrid.sp)
    }
  }
}


