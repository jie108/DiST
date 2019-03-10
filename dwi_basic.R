##ref: Wong, R.K.W., T.C.M. Lee, D. Paul, and J. Peng. Fiber direction estimation, smoothing and tracking in diffusion MRI (with Discussion)  (2016). The Annals of Applied Statistics, 10(3): 1137-1156.

######################################################################
######## General Statistics/ functions ########
######################################################################
# load C rountines
library(dwi.internals2)
library(quadprog)

###########################################################
#### compute fractional anistropy of a 3*3 p.d. matrix ####
###########################################################
# input: a 3*3 positive definite matrix A
# output: fractional anisotropy of A and eigenvectors of A (as a matrix)

FA <- function(A, pd.ep=1e-5, eigval=NULL)
{
  if (is.null(eigval)){
    eigval <- eigen(A, symmetric=T, only.values=T)$values
  }
  if (eigval[1]>0){
    eigval <- pmax(eigval, eigval[1]*pd.ep)
    return(sqrt((sum(eigval^2)-3*mean(eigval)^2)/sum(eigval^2))*sqrt(3/2)) # FA index
  } else {
    return(0)
  }
}

##########################################################
#### Convert indices of an array to indices of a list ####
##########################################################
# input: dims = dimension of the bigger matrix
# xrange, yrange, zrange = indices in x, y and z directions
## inverse function of arrayInd (internal R function)
# has a C replacement in dwi.internals2

ArrayIndex.old <- function(dims, xrange, yrange, zrange)
{
  indexlist <- rep(0, length(xrange)*length(yrange)*length(zrange))
  count <- 0
  for (k in (1:length(zrange))){
    for (j in (1:length(yrange))){
      for (i in (1:length(xrange))){
        count <- count+1
        indexlist[count] <- (zrange[k]-1)*dims[1]*dims[2] + (yrange[j]-1)*dims[1] + xrange[i]
      }
    }
  }
  return(indexlist)
}

###################################################################
#### some rountine for interchanging between vector and matrix ####
###################################################################
# from vector to matrix

vec2mat <- function(vec)
{
  mat <- matrix(0,nr=3,nc=3); mat[1,1] <- vec[1]; mat[2,2] <- vec[2]; mat[3,3] <- vec[3]; mat[1,2] <- vec[4]; mat[1,3] <- vec[5]; mat[2,3] <- vec[6]; mat[2,1] <- vec[4]; mat[3,1] <- vec[5]; mat[3,2] <- vec[6]
  return(mat)
}

# from matrix to vector

mat2vec <- function(mat)
{
  vec <- c(mat[1,1],mat[2,2],mat[3,3],mat[1,2],mat[1,3],mat[2,3])
  return(vec)
}

#### from package Bessel
besselIasym <- function (x, nu, k.max = 10, expon.scaled = FALSE, log = FALSE) 
{
    stopifnot(k.max == round(k.max))
    d <- 0
    if (k.max >= 1) {
        m. <- 4 * nu^2
        x8 <- 8 * x
        for (k in k.max:1) {
            d <- (1 - d) * ((2 * (nu - k) + 1) * (2 * (nu + k) - 
                1))/(k * x8)
        }
    }
    if (log) {
        (if (expon.scaled) 
            log1p(-d)
        else x + log1p(-d)) - log(2 * pi * x)/2
    }
    else {
        (if (expon.scaled) 
            (1 - d)
        else exp(x) * (1 - d))/sqrt(2 * pi * x)
    }
}

############################
#### bessel function ####
############################
# has a regime to switch to asymptotic approximation when x is large

besselIm <- function(x, nu, expon.scaled=F, thres=100000){
  index <- (x<=thres)
  out <- array(dim=length(x))
  if (sum(index)){
    out[index] <- besselI(x[index], nu=nu, expon.scaled=expon.scaled)
  }
  if (sum(!index)){
    out[!index] <- besselIasym(x[!index], nu=nu, expon.scaled=expon.scaled)
  }
  return(out)
}


###############################
#### Quadratic programming ####
###############################
solve.QP2 <- function (Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE) 
{
    n <- nrow(Dmat)
    q <- ncol(Amat)
    if (missing(bvec)) 
        bvec <- rep(0, q)
    if (n != ncol(Dmat)) 
        stop("Dmat is not symmetric!")
    if (n != length(dvec)) 
        stop("Dmat and dvec are incompatible!")
    if (n != nrow(Amat)) 
        stop("Amat and dvec are incompatible!")
    if (q != length(bvec)) 
        stop("Amat and bvec are incompatible!")
    if ((meq > q) || (meq < 0)) 
        stop("Value of meq is invalid!")
    iact <- rep(0, q)
    nact <- 0
    r <- min(n, q)
    sol <- rep(0, n)
    lagr <- rep(0, q)
    crval <- 0
    work <- rep(0, 2 * n + r * (r + 5)/2 + 2 * q + 1)
    iter <- rep(0, 2)
    res1 <- .Fortran(quadprog:::.QP_qpgen2, as.double(Dmat), dvec = as.double(dvec), 
        as.integer(n), as.integer(n), sol = as.double(sol), lagr = as.double(lagr), 
        crval = as.double(crval), as.double(Amat), as.double(bvec), 
        as.integer(n), as.integer(q), as.integer(meq), iact = as.integer(iact), 
        nact = as.integer(nact), iter = as.integer(iter), work = as.double(work), 
        ierr = as.integer(factorized))
    list(solution = res1$sol, value = res1$crval, unconstrained.solution = res1$dvec, 
        iterations = res1$iter, Lagrangian = res1$lagr, iact = res1$iact[1:res1$nact], ierr=res1$ierr)
}

######################################################################
######################################################################
########  tensor estimation: single fiber ########
######################################################################
######################################################################

#####################################
#### Generation of design matrix ####
#####################################
# To generate the design matrix given the gradient matrix
# input: grad.mat is a n*3 gradient matrix
# output: 6*n matrix [b1^2,b2^2,b3^2,2b1*b2,2b1*b3,2b2*b3]

generate.design.noS0 <- function(grad.mat)
{
  design <- cbind(grad.mat[,1]^2, grad.mat[,2]^2, grad.mat[,3]^2,
                  2*grad.mat[,1]*grad.mat[,2], 2*grad.mat[,1]*grad.mat[,3],
                  2*grad.mat[,2]*grad.mat[,3])
  return(design)
}

###################################
#### Convert the tensor vector ####
###################################
# Convert the tensor vector (-[D11,D22,D33,D12,D13,D23]) to tensor matrix
vec2ten.noS0 <- function(vec)
{
  tensor <- matrix(0,nr=3,nc=3)
  vec <- -vec
  tensor[1,1] <- vec[1]
  tensor[2,2] <- vec[2]
  tensor[3,3] <- vec[3]
  tensor[1,2] <- vec[4]
  tensor[1,3] <- vec[5]
  tensor[2,3] <- vec[6]
  tensor[2,1] <- vec[4]
  tensor[3,1] <- vec[5]
  tensor[3,2] <- vec[6]
  return(tensor)
}

# from (-[D11,D22,D33,D12,D13,D23]) (6*N) to tensor matrix (3*3*N)
# input: vector tensor 6*N matrix, (-[D11,D22,D33,D12,D13,D23]): i.e., the coefficients for [b1^2,b2^2,b3^2,2b1*b2,2b1*b3,2b2*b3]
# output: 3*3*N matrix

vec2ten.onlyD <- function(vec)
{
  N <- dim(vec)[2]
  tensor.mat <- array(0,dim=c(3,3,N))
  for (i in (1:N)){
    tensor.mat[,,i] <- vec2ten.noS0(vec[,i])
  }
  return(tensor.mat)
}

#####################################################
#### make the tensor estimates positive definite ####
#####################################################
# Make a given matrix p.d. by adjusting its negative eigen-value(s)
make.pd <- function(mat, pd.ep=1e-5)
{
  eigval <- eigen(mat, symmetric=T, only.values=T)$values
  dd <- eigval[3]
  ep  <- eigval[1]*pd.ep
  if (dd>=ep){
    ad.mat <- mat
  } else {
    eigdecom <- eigen(mat, symmetric=T)
    eigval <- eigdecom$values
    eigvec <- eigdecom$vectors
    for (i in (1:3)){
      if (eigval[i]<ep){
        eigval[i] <- ep
      }
    }
    ad.mat <- eigval[1]*eigvec[,1]%*%t(eigvec[,1]) +
    eigval[2]*eigvec[,2]%*%t(eigvec[,2]) + eigval[3]*eigvec[,3]%*%t(eigvec[,3])
  }
  return(list(mat=ad.mat, eigval=eigval))
}

# Make the estimated tensors on whole pixels p.d.
# If not p.d., adjust the negative eigenvalues to be positive
# input: tensor.mat: the matrices to be adjusted
# output: the adjusted matrices
make.pd.whole.pixels <- function(tensor.mat,pd.ep=1e-5)
{
  N <- dim(tensor.mat)[3]
  ad.tensor.mat <- tensor.mat
  for (i in (1:N)){
    ad.tensor.mat[,,i] <- make.pd(tensor.mat[,,i],pd.ep)$mat
  }
  return(ad.tensor.mat)
}

########################################
#### Levenberg CMarquardt algorithm ####
########################################
# input: f= SB/S0.adv = exp(-t(b)%*% D %*% b)
# grad.mat is a n*3 gradient matrix, S0.adv: n*N S0
# pre.est 3*3*N estimates of tensors from the previous step: the first step uses results from linear regression
# output: 3*3*N from nonlinear regression

adv.reg.once <- function(f, grad.mat, braindim, pre.est, eps)
{
  n <- nrow(grad.mat)
  N <- dim(pre.est)[3]
  exp.pre <- matrix(0,nr=n,nc=N)
  for (i in (1:N)){
    exp.pre[,i] <- -diag(grad.mat%*%pre.est[,,i]%*%t(grad.mat)) # fitted dwi's from pre.est
  }
  exp.pre = exp(exp.pre)

  data <- f-exp.pre # residual between observed and fitted
  design.b <- generate.design.noS0(grad.mat)
  tt <- rep(1,6)
  diff <- matrix(0,nr=6,nc=N)
  new.est <- pre.est
  m.b <- ncol(design.b)
  for (k in (1:N)){
    design <- -sum(exp.pre[,k]*tt)*design.b
    diff[,k] <- solve(t(design)%*% design + eps*diag(rep(1,m.b)), t(design)%*%data[,k]) # get updated amount
  }
  diff.mat <- array(0,dim=c(3,3,N))
  for (i in (1:N)){
    diff.mat[1,1,i] <- diff[1,i]
    diff.mat[2,2,i] <- diff[2,i]
    diff.mat[3,3,i] <- diff[3,i]
    diff.mat[1,2,i] <- diff[4,i]
    diff.mat[1,3,i] <- diff[5,i]
    diff.mat[2,3,i] <- diff[6,i]
    diff.mat[2,1,i] <- diff[4,i]
    diff.mat[3,1,i] <- diff[5,i]
    diff.mat[3,2,i] <- diff[6,i]
  }
  new.est <- pre.est + diff.mat
  return(new.est)
}

###############################
#### non/Linear regression ####
###############################
# non/linear regression method for tensor estimation for each voxel
# S0 is known
# for linear regression: nstep=0
# for nonlinear regression: nstep>=1

dti.nlm <- function(grad.mat, SB, S0, braindim, nstep, eps=1e-5)
{
  n.voxel <- braindim[1]*braindim[2]*braindim[3]
  grad.mat <- diag(sqrt(1/diag(grad.mat%*%t(grad.mat)))) %*% grad.mat # standardization
  n <- nrow(grad.mat)
  S0.adv <- matrix(0,nr=n,nc=n.voxel) + S0
  design <- generate.design.noS0(grad.mat) # generate the design matrix
  dwi.mat <- log(SB)-log(S0.adv) # linear regression on log scale: n by # of voxels
  coefficient <- solve(t(design)%*%design,t(design)%*%dwi.mat)
  temp2 <- vec2ten.onlyD(coefficient) # convert the coefficients to the 3*3*N format
  pre.est <- make.pd.whole.pixels(temp2,eps)
  linear.est <- pre.est

  # for nonlinear
  f <- SB/S0.adv # exp(-t(b)%*%D%*%b)
  if (nstep>0){
    for (i in (1:nstep)){
      eps.c <- eps*i/nstep
      new.est <- adv.reg.once(f, grad.mat, braindim, pre.est, eps.c)
      new.est <- make.pd.whole.pixels(new.est,eps)
      pre.est <- new.est
    }
  } else {
    new.est <- NA
  }
  return(list(linear.est=linear.est, new.est=new.est))
}

######################################################################
######################################################################
######## Data Simulation ########
######################################################################
######################################################################

#######################################
#### data simulation: Rician noise ####
#######################################
# To simulate DWI data of a given pixel, given the underlying tensor and the gradient direction (standardized to norm 1).
# Gaussian tensor model and rician noise model are used.
# input: the true tensor, the gradient and the simulated error. error is a 2*1 vector. grad is a 1*3 vector
# output: the observed DWI
# model for simulation: S(grad) = S0*exp(-grad %*% tensor %*% t(grad))
sim.dwi.per.pixel <- function(tensor, grad, error, S0)
{
  true.dwi <- S0*exp(-as.vector(t(grad) %*% tensor %*% grad)) # note that grad is interpreted as the column vector object in R
  dwi <- sqrt((true.dwi + error[1])^2 + error[2]^2) # for gaussian model, the complex part is zero and thus it is white noise
  return(dwi)
}

# lso note: tensor.true and norm of gradient should be in a comparable scale sucht that the resulting dwi is not too small/large
# e.g., if S0=1000, then Si should be in hundreds
# input: dt6.smooth is the underlying tensor, sigma is the std of the error. grad is a 1*3 vector. S0 is a 311296 vector.

sim.whole.pixels <- function(dt6.smooth, grad, braindim, sigma, S0)
{
  N <- braindim[1]*braindim[2]*braindim[3]
  dwi <- array(dim=N)
  errors <- sigma*matrix(rnorm(2*N),nr=2,nc=N)
  for (i in (1:N)){
    dwi[i] <- sim.dwi.per.pixel(dt6.smooth[,,i],grad,errors[,i],S0[i])
  }
  return(dwi)
}


##################################################
#### To simulate DWI data of the whole pixels ####
##################################################
# given a nx3 gradient matrix , the std sigma, and the S0.
# input: dt6.smooth is the underlying tensor, sigma is the std of the error.
# grad.mat is a n*3 gradient matrix
# output: n*N matrix. the row is the simulated DWI data for a particular
# gradient. And the column is the simulated DWI data for a particular pixel.
# Note: grad_matrix should be standardized to have norm 1 for each row (gradient direction)
# also, sigma should be in a  comparable scale as to S0

adv.multi.grad.organize <- function(dt6.smooth, grad.mat, braindim, sigma, S0)
{
  n <- nrow(grad.mat)
  N <- braindim[1]*braindim[2]*braindim[3]
  dwi.mat <- matrix(0,nr=n,nc=N)
  for (i in (1:n)){
    dwi.mat[i,] <- sim.whole.pixels(dt6.smooth, grad.mat[i,],braindim,sigma,S0)
  }
  return(dwi.mat)
}

###############
#### getSB ####
###############
getSB <- function(dt6.smooth,grad.mat, braindim, sigma, S0)
{
  n.voxel <- braindim[1]*braindim[2]*braindim[3]
  grad.mat <- diag(sqrt(1/diag(grad.mat%*%t(grad.mat)))) %*% grad.mat # standardization to make gradients norm 1

  n <- dim(grad.mat)
  S0.adv <- array(0,dim=c(n[1],n.voxel)) + S0 # S0 for each voxel
  SB <- adv.multi.grad.organize(dt6.smooth, grad.mat, braindim, sigma, S0.adv[1,])
  return(SB)
}

#######################################
#### simulation of tensors (whole) ####
#######################################
# For one voxel
sim.mix.SB <- function(tens, p, design, error, S0, b=1){
  #grad.mat <- diag(sqrt(1/diag(grad.mat%*%t(grad.mat)))) %*% grad.mat # standardization
  #design <- generate.design.noS0(grad.mat)
  tdwis <- (S0*exp(-b*design%*% tens))%*%p
  dwis <- sqrt((tdwis + error[,1])^2 + error[,2]^2) # for gaussian model, the complex part is zero and thus it is white noise
  return(dwis)
}

sim.mix.SB.whole <- function(tensor.whole, p.whole, grad.mat, braindim, sigma, S0, b=1){
  grad.mat <- diag(sqrt(1/diag(grad.mat%*%t(grad.mat)))) %*% grad.mat # standardization
  design <- generate.design.noS0(grad.mat)
  n.voxel <- prod(braindim)
  n <- nrow(grad.mat)
  dwi.mat <- matrix(0,nr=n,nc=n.voxel)
  errors <- sigma*array(rnorm(n*2*n.voxel),dim=c(n,2,n.voxel))
  n.layer <- dim(tensor.whole)[3]
  for (i in (1:n.voxel)){
    temp.tens <- matrix(0,nr=6,nc=n.layer)
    for (j in (1:n.layer)){
      temp.tens[,j] <- mat2vec(tensor.whole[,,j,i])
    }
    dwi.mat[,i] <- sim.mix.SB(temp.tens,p.whole[,i],design,errors[,,i],S0,b)
  }
  return(dwi.mat)
}

###########################
#### For fiber example ####
###########################
gram.schmidt.3d <- function(mat){
  mat[,1] <- mat[,1]/sqrt(sum(mat[,1]^2))
  mat[,2] <- mat[,2] - sum(mat[,1]*mat[,2])*mat[,1]
  mat[,2] <- mat[,2]/sqrt(sum(mat[,2]^2))
  mat[,3] <- mat[,3] - sum(mat[,1]*mat[,3])*mat[,1] - sum(mat[,2]*mat[,3])*mat[,2]
  mat[,3] <- mat[,3]/sqrt(sum(mat[,3]^2))
  return(mat)
}

get.sym.tens <- function(lambdas,v1){
  if (qr(cbind(v1,c(1,0,0),c(0,1,0)))$rank==3){
    mat <- cbind(v1,c(1,0,0),c(0,1,0))
  } else if (qr(cbind(v1,c(1,0,0),c(0,0,1)))$rank==3){
    mat <- cbind(v1,c(1,0,0),c(0,0,1))
  } else if (qr(cbind(v1,c(0,1,0),c(0,0,1)))$rank==3){
    mat <- cbind(v1,c(0,1,0),c(0,0,1))
  }

  eig.vecs <- gram.schmidt.3d(mat)
  res <- eig.vecs%*%diag(lambdas)%*%t(eig.vecs)
  return(res)
}

sym.tens.from.fit <- function(j,res){
  mat2vec(get.sym.tens(res$ten.lambdas,res$basis.verts[j,]))
}

# fun is a function taking a three dimensional vectors and return either in or out of this fiber
# n.fiber is the current n.fiber
fill.fiber <- function(fun, p.fib, braindim, tensor, tensor.array, p, p.array,
                       v0, n.fiber, final, xrange=c(0,1), yrange=c(0,1),
                       zrange=c(0,1)){
  xgrid <- seq(xrange[1], xrange[2], len=braindim[1]+2)[-c(1,braindim[1]+2)]
  ygrid <- seq(yrange[1], yrange[2], len=braindim[2]+2)[-c(1,braindim[2]+2)]
  zgrid <- seq(zrange[1], zrange[2], len=braindim[3]+2)[-c(1,braindim[3]+2)]

  for (i in (1:braindim[1])){
    for (j in (1:braindim[2])){
      for (k in (1:braindim[3])){
        out <- fun(c(xgrid[i], ygrid[j], zgrid[k]))
        voxindex <- ArrayIndex(braindim, i, j, k)
        if (out$ind){
          tensor[,,n.fiber[voxindex]+1, voxindex] <- out$tens
          tensor.array[,,n.fiber[voxindex]+1, i, j, k] <- out$tens
          p[n.fiber[voxindex]+1, voxindex] <- p.fib
          p.array[n.fiber[voxindex]+1, i, j, k] <- p.fib
          v0[,n.fiber[voxindex]+1,voxindex] <- out$v
          n.fiber[voxindex] <- n.fiber[voxindex]+1
        }
        if (final){
          p[,voxindex] <- p[,voxindex]/sum(p[,voxindex])
          p.array[,i,j,k] <- p.array[,i,j,k]/sum(p.array[,i,j,k])
        }
      }
    }
  }
  return(list(tensor=tensor, tensor.array=tensor.array, p=p, p.array=p.array, v0=v0, n.fiber=n.fiber))
}

##########################################################
##########################################################
######## plotting ########
##########################################################
##########################################################
# for ploting examples generated from sim-curve.R
plot.brain.ten <- function(tens, ps, braingrid, n.layer, n.fiber, minbraingrid){
  braindim <- dim(braingrid)[-1]
  n.voxel <- prod(braindim)

  # find a const
  const <- (minbraingrid/(max(abs(tens))*4))^2

  #   plot3d.new(xlim=c(0,max(braingrid[1,,,])), ylim=c(0,max(braingrid[2,,,])), zlim=c(0,max(braingrid[3,,,])), aspect=T, xlab="", ylab="", zlab="", box=F)
  open3d()
  for (i in (1:n.voxel)){
    voxindex <- as.vector(arrayInd(i, braindim))
    mu <- braingrid[, voxindex[1], voxindex[2], voxindex[3]]
    if (n.fiber[i]>0){
      for (j in (1:n.fiber[i])){
        plot3d(ellipse3d(tens[,,j, i]*const, centre=mu), add=T, alpha=0.8*ps[j,i], col=get.color(tens[,,j, i]), aspect=T, box=T)
      }
    }
  }
}

# for ploting the results (not useful since we are fitting tensor anymore)
# n.fiber is n.fiber2
plot.brain.ten.fit <- function(tens, ps, locs, n.fiber, minbraingrid){
  # find a const
  const <- (minbraingrid/(max(abs(tens))*4))^2
  open3d()
  for (iind in (1:length(n.fiber))){
    if (n.fiber[iind]>0){
      plot3d(ellipse3d(tens[iind,,]*const, centre=locs[iind,]), add=T, alpha=0.8*ps[iind], col=get.color(tens[iind,,]), aspect=T)
    }
  }
}

get.color <- function(mat){
  temp <- diag(mat)/sum(diag(mat))
  return(rgb(red=temp[1], green=temp[2], blue=temp[3], alpha=1))
}

