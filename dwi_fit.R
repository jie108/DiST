## ref: Wong, R.K.W., T.C.M. Lee, D. Paul, and J. Peng. Fiber direction estimation, smoothing and tracking in diffusion MRI (with Discussion)  (2016). The Annals of Applied Statistics, 10(3): 1137-1156.
#########################################
### functions for direction estimation and smoothing
#########################################

###########################
#### loading libraries ####
###########################
library(doParallel)
library(rgenoud)
library(cluster)


################################
#### tessellation of sphere ####
################################
tess.ico <- function(depth, no.rev=T, initial="ico", random=F)
{
  if (initial=="ico"){
    # construct icosahedron
    t <- (1+sqrt(5))/2; tau <- t/sqrt(1+t^2) # tau <- 0.8506508084
    one <- 1/sqrt(1+t^2) # one <- 0.5257311121

    za <- c(tau, one, 0)
    zb <- c(-tau, one, 0)
    zc <- c(-tau, -one, 0)
    zd <- c(tau, -one, 0)
    ya <- c(one, 0, tau)
    yb <- c(one, 0, -tau)
    yc <- c(-one, 0, -tau)
    yd <- c(-one, 0, tau)
    xa <- c(0, tau, one)
    xb <- c(0, -tau, one)
    xc <- c(0, -tau, -one)
    xd <- c(0, tau, -one)

    vert <- cbind(za,zb,zc,zd,ya,yb,yc,yd,xa,xb,xc,xd)
    vert <- apply(vert,2,function(x){x/sqrt(sum(x^2))})
    #face <- cbind(c(5 9 8), c(5 8 10), c(6 7 12), c(6 11 7), c(1 5 4), c(1 4 6), c(3 8 2), c(3 2 7), c(9 1 12), c(9 12 2), c(10 11 4), c(10 3 11), c(9 5 1), c(12 1 6), c(5 10 4), c(6 4 11), c(8 9 2), c(7 2 12), c(8 3 10), c(7 11 3))
    tris <- array(dim=c(3,3,20))
    tris[,,1] <- cbind(ya, xa, yd)
    tris[,,2] <- cbind(ya, yd, xb)
    tris[,,3] <- cbind(yb, yc, xd)
    tris[,,4] <- cbind(yb, xc, yc)
    tris[,,5] <- cbind(za, ya, zd)
    tris[,,6] <- cbind(za, zd, yb)
    tris[,,7] <- cbind(zc, yd, zb)
    tris[,,8] <- cbind(zc, zb, yc)
    tris[,,9] <- cbind(xa, za, xd)
    tris[,,10] <- cbind(xa, xd, zb)
    tris[,,11] <- cbind(xb, xc, zd)
    tris[,,12] <- cbind(xb, zc, xc)
    tris[,,13] <- cbind(xa, ya, za)
    tris[,,14] <- cbind(xd, za, yb)
    tris[,,15] <- cbind(ya, xb, zd)
    tris[,,16] <- cbind(yb, zd, xc)
    tris[,,17] <- cbind(yd, xa, zb)
    tris[,,18] <- cbind(yc, zb, xd)
    tris[,,19] <- cbind(yd, zc, xb)
    tris[,,20] <- cbind(yc, xc, zc)
    ntris <- 20
  } else if (initial=="tetra"){
    cat("extra points... from tetra options\n")
    ppp <- c(1,1,1)
    mmp <- c(-1,-1,1)
    mpm <- c(-1,1,-1)
    pmm <- c(1,-1,-1)

    vert <- cbind(ppp,mmp,mpm,pmm)
    vert <- apply(vert,2,function(x){x/sqrt(sum(x^2))})

    tris <- array(dim=c(3,3,4))
    tris[,,1] <- cbind(ppp,mmp,mpm)
    tris[,,2] <- cbind(ppp,pmm,mmp)
    tris[,,3] <- cbind(mpm,mmp,pmm)
    tris[,,4] <- cbind(pmm,ppp,mpm)
    ntris <- 4
  } else if (initial=="oct"){
    xplus <- c(1,0,0)
    xmin <- c(-1,0,0)
    yplus <- c(0,1,0)
    ymin <- c(0,-1,0)
    zplus <- c(0,0,1)
    zmin <- c(0,0,-1)

    vert <- cbind(xplus,xmin,yplus,ymin,zplus,zmin)
    vert <- apply(vert,2,function(x){x/sqrt(sum(x^2))})

    tris <- array(dim=c(3,3,8))
    tris[,,1] <- c(xplus,zplus,yplus)
    tris[,,2] <- c(yplus,zplus,xmin)
    tris[,,3] <- c(xmin,zplus,ymin)
    tris[,,4] <- c(ymin,zplus,xplus)
    tris[,,5] <- c(xplus,yplus,zmin)
    tris[,,6] <- c(yplus,xmin,zmin)
    tris[,,7] <- c(xmin,ymin,zmin)
    tris[,,8] <- c(ymin,xplus,zmin)
    ntris <- 8
  }

  if (depth>0){
    for (d in (1:depth)){
      # construct normalized midpoints
      midA <- apply(tris[,c(1,3),],c(1,3),mean)
      midA <- apply(midA,2,function(x){x/sqrt(sum(x^2))})
      midB <- apply(tris[,c(1,2),],c(1,3),mean)
      midB <- apply(midB,2,function(x){x/sqrt(sum(x^2))})
      midC <- apply(tris[,c(2,3),],c(1,3),mean)
      midC <- apply(midC,2,function(x){x/sqrt(sum(x^2))})

      # construct 4 new triangles for every existing triangle
      tris1 <- tris
      tris <- array(dim=c(3,3,4*ntris))
      for (i in (1:ntris)){
        tris[,,i] <- cbind(tris1[,1,i], midB[,i], midA[,i])
        tris[,,ntris+i] <- cbind(midB[,i],tris1[,2,i],midC[,i])
        tris[,,2*ntris+i] <- cbind(midA[,i], midB[,i],midC[,i])
        tris[,,3*ntris+i] <- cbind(midA[,i],midC[,i],tris1[,3,i])
      }
      ntris <- ntris*4
    }
  }
  verts <- (unique(t(matrix(tris,nr=3))))
  if (no.rev){
    verts <- verts[verts[,3]>=0,]
    verts <- verts[((verts[,3]==0)*(verts[,2]<0))==0,]
    verts <- verts[((verts[,3]==0)*(verts[,2]==0)*(verts[,1]<0))==0,]
  }
  verts <- t(verts)

  if (random){
    angx <- runif(1,0,2*pi)
    angy <- runif(1,0,2*pi)
    angz <- runif(1,0,2*pi)
    rot.mat <- sphere.rotate(angx,angy,angz)
    verts <- rot.mat%*%verts
  }
  return(list(verts=verts,tris=tris))
}

subtess <- function(tess.obj,n,include=NULL,nthres=10,seed=1234,Nsim=5000)
{
  verts <- tess.obj$verts
  nverts <- ncol(verts)
  if (n >= nverts){
    stop("n>=nverts!\n")
  }

  set.seed(seed)
  cand <- matrix(nr=Nsim,nc=n)
  ninclude <- length(include)
  if (ninclude>=1){
    cand[,1:ninclude] <- rep(include,each=Nsim)
    ind1 <- (ninclude+1):n
  } else {
    ind1 <- 1:n
  }
  ind0 <- setdiff(1:nverts,include)
  crit <- array(dim=Nsim)
  for (i in (1:Nsim)){
    cand[i,ind1] <- sample(ind0,size=n-ninclude,replace=F)
    temp <- verts[,cand[i,]]
    temptemp <- abs(t(temp)%*%temp)
    diag(temptemp) <- 1
    #crit[i] <- sum((acos((temptemp))))
    temptemptemp <- acos(temptemp)
    crit[i] <- sum(apply(temptemptemp,2,function(x){sum(sort(x)[1:nthres])}))
  }
  out.verts <- verts[,cand[which.max(crit),]]
  return(out.verts)
}

#########################
#### rotation matrix ####
#########################
sphere.rotate <- function(angx, angy, angz)
{
  return(matrix(
                c( cos(angy)*cos(angz), -cos(angx)*sin(angz)+sin(angx)*sin(angy)*cos(angz),
                  sin(angx)*sin(angz) + cos(angx)*sin(angy)*cos(angz),
                  cos(angy)*sin(angz), cos(angx)*cos(angz)+sin(angx)*sin(angy)*sin(angz),
                  -sin(angx)*cos(angz)+ cos(angx)*sin(angy)*sin(angz),
                  -sin(angy), sin(angx)*cos(angy), cos(angx)*cos(angy)
                  ),nr=3,byrow=T))
}

###############################
#### L1 fit (racian noise) ####
###############################

#### an internal function for solving lambda when the equality constraint is implement
dwi.fit.l1.constraint.internal <- function(lasso.lam, betas, y, S, gam,
                                           betas.w, sigma, Nmiddle.max=100,
                                           Ninner=25)
{
  internal.fit <- dwi.fit.l1.internal(lasso.lam, betas, y, S, gam, betas.w,
                                      sigma, Nmiddle.max=Nmiddle.max,
                                      Ninner=Ninner, tol=1e-8, beta.max=100)
  return(sum(internal.fit$betas)-1)
}

#####
# if lasso.lam is given, the equality constraint is used to select lasso.lam
dwi.fit.l1.constraint <- function(dwis, sigma, S, gam=0,
                                  betas.w=rep(1,ncol(S)), betas0=NULL,
                                  Ninner=25, Nmiddle.max=100,
                                  lasso.lam.max=NULL,lasso.lam.min=NULL,
                                  lasso.lam=NULL, betas0.equal=T)
{
  # initial constructions
  y <- dwis
  p <- ncol(S)
  n <- nrow(S)
  lasso.lam.given <- T

  if (is.null(betas0)){
    ridge.lam <- 1e-7
    Dmat <- t(S)%*%S+ridge.lam*diag(ncol(S))
    dvec <- t(S)%*%y -1/2*gam*betas.w
    if (betas0.equal){
      Amat <- cbind(1,diag(ncol(S)))
      bvec <- c(1,rep(0,ncol(S)))
      meq <- 1
    } else {
      Amat <- cbind(diag(ncol(S)))
      bvec <- rep(0,ncol(S))
      meq <- 0
    }
    qpres <- solve.QP2(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=meq)
    icount <- 0
    while ((qpres$ierr>0)&&(icount<=5)){
      ridge.lam <- ridge.lam*10
      Dmat <- t(S)%*%S+ridge.lam*diag(ncol(S))
      qpres <- solve.QP2(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=meq)
      icount <- icount + 1
    }
    betas0 <- qpres$solution
  }
  betas <- betas0

  if (is.null(lasso.lam)){
    # getting the maximum lambda
    if (is.null(lasso.lam.max)||is.null(lasso.lam.min)){
      temp <- as.vector(S%*%betas)/sigma^2
      temp2 <- temp*y
      temp3 <- -temp + besselIm(temp2,1,T)/besselIm(temp2,0,T)*y/sigma^2
      A1 <- matrix(temp3,nr=n,nc=p)*S #n*p matrix

      Aw.sumi <- as.vector(A1)
      Bw.sumi <- (t(S)%*%S)/sigma^2
      denom <- diag(Bw.sumi)/n

      fitrss <- Aw.sumi/n
      if (is.null(lasso.lam.max)){
        lasso.lam.max <- max(fitrss+denom*betas-gam*betas.w-1/(2*p)*denom)
      }
      if (is.null(lasso.lam.min)){
        lasso.lam.min <- min(fitrss+denom*betas-gam*betas.w-2/p*denom)
      }
    }

    uniroot.fit <- myuniroot(f=dwi.fit.l1.constraint.internal,
                             interval=c(lasso.lam.min, lasso.lam.max),
                             betas=betas0, y=y, S=S, gam=gam,
                             betas.w=betas.w, sigma=sigma,
                             Nmiddle.max=Nmiddle.max, Ninner=Ninner, tol=1e-20)
    if (abs(uniroot.fit$f.root)>1e-4){
      cat("problem in uniroot:", uniroot.fit$f.root,"\n")
    }

    lasso.lam <- uniroot.fit$root
    lasso.lam.given <- F
  }
  internal.fit <- dwi.fit.l1.internal(lasso.lam, betas0, y, S, gam,
                                      betas.w, sigma, Nmiddle.max=Nmiddle.max,
                                      Ninner=Ninner, tol=1e-8)
  betas <- internal.fit$betas

  return(list(betas=betas, lasso.lam=lasso.lam, S=S, dwis=dwis, sigma=sigma,
              gam=gam, betas.w=betas.w, betas0=betas0,
              lasso.lam.given=lasso.lam.given))
}

#  modification of uniroot
myuniroot <- function(f, interval, ..., lower = min(interval), upper =
                      max(interval),f.lower = f(lower, ...), f.upper = f(upper,
                                                                         ...),
                      tol = .Machine$double.eps^0.25, maxiter = 1000,
                      max.count=10)
{
  len <- upper-lower
  count <- 0
  while ((f.lower * f.upper > 0)&&(count<=max.count)){
    lower <- lower-len
    upper <- upper+len
    f.lower <- f(lower,...)
    f.upper <- f(upper,...)
    count <- count + 1
  }
  if (count>max.count){
    stop('Cannot get the right limits for uniroot!\n')
  }
  return(uniroot(f=f, interval=interval, ..., lower = lower, upper = upper, f.lower = f.lower, f.upper = f.upper, tol = tol, maxiter = maxiter))
}

##############################
#### Fit isotropic tensor ####
##############################
isotropic.mle.fn <- function(tdwi, dwis, sigma, S0){
  temp <- dwis*tdwi/sigma^2
  nll <- sum(tdwi^2/(2*sigma^2) - log(besselIm(temp,0,T)) - temp)
  return(nll)
}

isotropic.mle <- function(tdwi0, dwis, b, sigma, S0)
{
  res <- optim(par=tdwi0, fn=isotropic.mle.fn, dwis=dwis, sigma=sigma, S0=S0,
               method="L-BFGS-B", lower=0, upper=Inf)
  return(list(eig=log(res$par/S0)/(-b),optim=res,ll=-res$val))
}


#############################################################################################
#### new parametrization ####
#############################################################################################
# alpha = lam1-lam2
# tau = wj exp(-b * lam2) # assuming q as unit vector

# likelihood function
mle.cs2.fn <- function(pars, J, dwis, sigma, S0, design)
{
  angs <- matrix(pars[(1:(J*2))], nrow=J)
  alphas <- pars[(J*2+1):(J*3)]
  taus <- pars[(J*3+1):(J*4)]
  coefficient <- matrix(nr=6,nc=J)
  for (j in (1:J)){
    direction <- c(cos(angs[j,1])*sin(angs[j,2]), sin(angs[j,1])*sin(angs[j,2]), cos(angs[j,2]))
    coefficient[,j] <- mat2vec(alphas[j]*direction%o%direction)
  }
  tdwi <- as.vector((S0*exp(-design %*% coefficient))%*%taus)
  temp <- dwis*tdwi/sigma^2
  nll <- sum(tdwi^2/(2*sigma^2) - log(besselIm(temp,0,T)) - temp)
  return(nll)
}


mle.cs2 <- function(J, dwis, sigma, S0, design, alphas, taus, directions,
                    control=list(), method="optim"){

  # expect: the scale of alpha*coefficient ~ 1
  if (method=="optim"){
    directions <- matrix(directions,nrow=J) # to safeguard the degeneracy when J=1
    angs <- matrix(nr=J, nc=2)
    for (j in (1:J)){
      angs[j,1] <- atan2(directions[j,2], directions[j,1])
      angs[j,2] <- acos(directions[j,3])
    }

    par0 <- c(as.vector(angs), alphas, taus)
    lower <- c(rep(-Inf,2*J),rep(1e-8,J),rep(1e-8,J))
    upper <- c(rep(Inf, 2*J), rep(Inf,J),rep(1-1e-8,J))
    res <- optim(par0, mle.cs2.fn, J=J, dwis=dwis, sigma=sigma, S0=S0,
                 design=design, method="L-BFGS-B", lower=lower,
                 upper=upper, control=control, hessian=F)
  } else if (method=="genoud") {
    lower <- c(rep(-2*pi,2*J),rep(1e-8,J),rep(1e-8,J))
    upper <- c(rep(4*pi, 2*J), rep(20,J),rep(1-1e-8,J))
    res <- genoud(mle.cs2.fn, nvars=length(lower), J=J, dwis=dwis,
                  sigma=sigma, S0=S0, design=design,
                  Domains=cbind(lower, upper), boundary.enforcement = 2,
                  print.level=1)
  }

  angs <- matrix(res$par[(1:(J*2))], nrow=J)
  alphas <- res$par[(J*2+1):(J*3)]
  taus <- res$par[(J*3+1):(J*4)]
  coefficient <- matrix(nr=6,nc=J)
  directions <- matrix(nr=J, nc=3)
  for (j in (1:J)){
    directions[j,] <- c(cos(angs[j,1])*sin(angs[j,2]),
                        sin(angs[j,1])*sin(angs[j,2]), cos(angs[j,2]))
    coefficient[,j] <- mat2vec(alphas[j]*directions[j,]%o%directions[j,])
  }

  return(list(coefficient=coefficient, alphas=alphas, taus=taus,
              directions=directions, optim=res))
}

mle.cs2.model <- function(dwis, sigma, S0, design, maxJ=4){
  tempfit <- list()
  ll <- array(dim=maxJ)
  aic <- array(dim=maxJ)
  bic <- array(dim=maxJ)

  for (j in (1:maxJ)){
    cat(j,"...\n")
    tempfit[[j]] <- mle.cs2(J=j, dwis=dwis, sigma=sigma, S0=S0,
                            design=design, alphas=NULL, taus=NULL,
                            directions=NULL, method="genoud")
    ll[j] <- tempfit[[j]]$optim$value
    aic[j] <- 2*ll[j] + 2*length(tempfit[[j]]$optim$par)
    bic[j] <- 2*ll[j] + log(length(dwis))*length(tempfit[[j]]$optim$par)
  }
  list(ll=ll, aic.mod=tempfit[[which.min(aic)]],
       bic.mod=tempfit[[which.min(bic)]], fits=tempfit, aicJ=which.min(aic),
       bicJ=which.min(bic))
}

######################################################
#### Grid methods for initial point of directions ####
######################################################

## get design matrix
get.design.for.basis3 <- function(design, S0, b, basis.verts, alphas){
  # alphas is a vector of alpha for each direction # length(alphas)=nverts
  nverts <- nrow(basis.verts)
  n <- nrow(design)
  S <- matrix(nr=n, nc=nverts)
  for (k in (1:nverts)){
    ten <- mat2vec(alphas[k]*basis.verts[k,]%o%basis.verts[k,])
    S[,k] <- S0 * exp(-b*design%*%ten)
  }
  return(list(S=S, alphas=alphas, basis.verts=basis.verts, nverts=nverts))
}

## linear models for direction
dir.lm.l1 <- function(dwis, sigma, S0, design, b=1, def.alphas=2,
                   betas.w=NULL, gam=0, basis.depth=3,
                   basis.random=F, basis.init="ico", Ninner=25,
                   Nmiddle.max=100, basis.verts=NULL, like=F){

  # generate basis
  #grad.mat <- diag(sqrt(1/diag(grad.mat%*%t(grad.mat)))) %*% grad.mat # standardization
  #design <- generate.design.noS0(grad.mat)
  
  if (is.null(basis.verts)){
    basis.verts <- t(tess.ico(basis.depth, no.rev=T, random=basis.random,
                              initial=basis.init)$verts)
  }
  nverts <- nrow(basis.verts)

  # initialize the lambda1 and fa
  alphas <- rep(def.alphas, nverts)

  # generate the design matrix
  design.obj <- get.design.for.basis3(design=design, S0=S0, b=b,
                                      basis.verts=basis.verts, alphas=alphas)
  S <- design.obj$S
  if (is.null(betas.w)){
    betas.w <- rep(1, ncol(S))
  }

  # fitting linear models: lasso.lam=0, betas.equal=F
  # lasso.lam.max, lasso.lam.min are useless since lasso.lam=0
  fit <- dwi.fit.l1.constraint(dwis=dwis, sigma=sigma, S=S, gam=gam,
                               betas.w=betas.w, betas0=NULL,
                               Ninner=Ninner, Nmiddle.max=Nmiddle.max,
                               lasso.lam.max=NULL, lasso.lam.min=NULL,
                               lasso.lam=0, betas0.equal=F)
  if (like){
    tdwi <- as.vector(S%*%fit$betas)
    temp <- dwis*tdwi/sigma^2
    nll <- sum(tdwi^2/(2*sigma^2) - log(besselIm(temp,0,T)) - temp)
    bic <- 2*nll + sum(fit$betas>1e-15)*log(length(dwis))
  } else {
    nll <- NULL
    bic <- NULL
  }
  return(list(fit=fit, basis.verts=basis.verts, def.alphas=def.alphas, nll=nll, bic=bic))
}

################################
#### Biased (direction) MLE ####
################################
# with clustering embedded for getting initial for the MLE
biased.mle.cs2 <- function(dwis, design, sigma, S0, b=1, maxJ=4,
                           basis.depth=3, basis.random=F, def.alphas=2,
                           betas0=NULL, Ninner=10, Nmiddle.max=40, alphas.const=NULL)
{
  if (is.null(alphas.const)){
    alphas.const <- b
  }
  def.alphas <- def.alphas/alphas.const

  tol <- 1e-15 # tolerance for claiming the betas as 0

  # linear model for directions
  res <- dir.lm.l1(dwis=dwis, sigma=sigma, S0=S0, design=design, b=b,
                   def.alphas=def.alphas, betas.w=NULL, gam=0,
                   basis.depth=basis.depth, basis.random=basis.random,
                   basis.init="ico", Ninner=Ninner, Nmiddle.max=Nmiddle.max,
                   basis.verts=NULL, like=F)

  dfit <- res$fit
  basis.verts <- res$basis.verts

  if (max(dfit$betas)>0){
    J1 <- sum(dfit$betas>0)
    #cat("Got", J1, "directions...\n")
    directions <- basis.verts[dfit$betas>tol,]
    betas <- dfit$betas[dfit$betas>tol]
    dis <- get.similar(directions)
    J <- min(maxJ, J1)
    use0 <- F
  } else if (max(dfit$betas0)>0){
    J1 <- sum(dfit$betas0>0)
    #cat("Got", J1, "directions...\n")
    directions <- basis.verts[dfit$betas0>tol,]
    betas <- dfit$betas[dfit$betas0>tol]
    dis <- get.similar(directions)
    J <- min(c(3, maxJ, J1)) # if betas0 is used, it usually returns a lot, so J is cap at 3 for computational purpose
    use0 <- T
  }

  # getting initial directions
  pobj.list <- list()
  dir.list <- list()
  if (J1==0){
    # try single fiber model
    # use betas0
    J <- 1
    pobj.list[[1]] <- list(clustering=1)
    dir.list[[1]] <- matrix(c(1,0,0),nc=3) # just choose a stupid start
  } else if (J1==1){
    pobj.list[[1]] <- list(clustering=1)
    dir.list[[1]] <- matrix(directions,nc=3)
  } else {
    for (jj in (1:J)){
      if (jj==J1){
        pobj.list[[jj]] <- list(clustering=1:J1)
      } else {
        pobj.list[[jj]] <- pam(x=dis, k=jj, diss=T)
      }
      dir.list[[jj]] <- matrix(nr=jj, nc=3)
      for (kk in (1:jj)){
        ind <- pobj.list[[jj]]$clustering==kk
        if (sum(ind)>1){
          tempdirect <- trans(directions[ind,][which.max(betas[ind]),], directions[ind,])
          dir.list[[jj]][kk,] <- get.w.i.mean2(tempdirect, w=rep(1, sum(ind)))$v0
        } else {
          dir.list[[jj]][kk,] <- directions[ind,]
        }
      }
    }
  }

  # mle fits
  nll <- array(dim=J)
  aic <- array(dim=J)
  bic <- array(dim=J)
  mle.fits <- list()
  for (j in (1:J)){
    # note that mle.cs2 does not take care of b naturally, but it is equivalence to b*design
    mle.fits[[j]] <- mle.cs2(J=j, dwis=dwis, sigma=sigma, S0=S0,
                              design=design*(b/alphas.const), alphas=rep(def.alphas,j),
                              taus=rep(1/j, j), directions=dir.list[[j]],
                              control=list(), method="optim")
    mle.fits[[j]]$alphas <- mle.fits[[j]]$alphas/alphas.const
    nll[j] <- mle.fits[[j]]$optim$value
    aic[j] <- 2*nll[j] + 2*length(mle.fits[[j]]$optim$par)
    bic[j] <- 2*nll[j] + log(length(dwis))*length(mle.fits[[j]]$optim$par)
  
  }

  return(list(lm.fit=dfit, lm.dir=directions, J1=J1, J=J, pobj=pobj.list,
              dirs=dir.list, mle.fit=mle.fits, aic=aic, bic=bic, use0=use0,
              alphas.const=alphas.const))
}

############################################################################################
#### Profile likelihood for demonstration of multiple maximizers of likelihood function ####
############################################################################################

# likelihood function
mle.cs2.profile.fn <- function(pars,J,angs,dwis,sigma,S0,design){
  alphas <- pars[(1):(J)]
  taus <- pars[(J+1):(J*2)]
  nnbhd <- ncol(dwis)
  dwis <- as.vector(dwis)
  coefficient <- matrix(nr=6,nc=J)
  for (j in (1:J)){
    direction <- c(cos(angs[j,1])*sin(angs[j,2]), sin(angs[j,1])*sin(angs[j,2]), cos(angs[j,2]))
    coefficient[,j] <- mat2vec(alphas[j]*direction%o%direction)
  }
  tdwi <- (S0*exp(-design %*% coefficient))%*%taus
  tdwi <- rep(tdwi,nnbhd)
  temp <- dwis*tdwi/sigma^2
  nll <- as.vector(rep(-1,nrow(design))%*%matrix(-tdwi^2/(2*sigma^2) + log(besselIm(temp,0,T)) + temp,nc=nnbhd))
  return(sum(nll))
}


mle.cs2.profile <- function(J, dwis, sigma, S0, design, alphas, taus, directions,
                    control=list(), method="optim"){

  directions <- matrix(directions,nrow=J) # to safeguard the degeneracy when J=1
  angs <- matrix(nr=J, nc=2)
  for (j in (1:J)){
    angs[j,1] <- atan2(directions[j,2], directions[j,1])
    angs[j,2] <- acos(directions[j,3])
  }
  if (method=="optim"){

    #par0 <- c(as.vector(angs), alphas, taus)
    par0 <- c(alphas, taus)
    lower <- c(rep(1e-8,J),rep(1e-8,J))
    upper <- c(rep(Inf,J),rep(1-1e-8,J))
    res <- optim(par0, mle.cs2.profile.fn, J=J, angs=angs, dwis=dwis, sigma=sigma, S0=S0,
                 design=design, method="L-BFGS-B", lower=lower,
                 upper=upper, control=control, hessian=F)
  } else if (method=="genoud") {
    lower <- c(rep(1e-8,J),rep(1e-8,J))
    upper <- c(rep(20,J),rep(1-1e-8,J))
    res <- genoud(mle.cs2.profile.fn, nvars=length(lower), J=J, angs=angs, dwis=dwis,
                  sigma=sigma, S0=S0, design=design,
                  Domains=cbind(lower, upper), boundary.enforcement = 2,
                  print.level=1)
  }

  alphas <- res$par[(1):(J)]
  taus <- res$par[(J+1):(J*2)]
  dwis <- as.vector(dwis)
  coefficient <- matrix(nr=6,nc=J)
  directions <- matrix(nr=J, nc=3)
  for (j in (1:J)){
    directions[j,] <- c(cos(angs[j,1])*sin(angs[j,2]), sin(angs[j,1])*sin(angs[j,2]), cos(angs[j,2]))
    coefficient[,j] <- mat2vec(alphas[j]*directions[j,]%o%directions[j,])
  }

  return(list(coefficient=coefficient, alphas=alphas, taus=taus,
              directions=directions, optim=res))
}


######################################################################################
#### two-stage procedure for switching between isotropic and anisotropic tensor ####
######################################################################################
####
lm.log <- function(dwis, design, S0, b=1, pd.correct=F, pd.ep=1e-5){
  #design <- generate.design.noS0(grad.mat)
  y <- log(dwis)-log(S0)
  res <- lm(y~design-1)
  mat  <- vec2mat(-res$coefficients)/b
  return(make.pd(mat, pd.ep=1e-5))
}

# note def.alphas will be adjusted to be def.alphas/alphas.const
biased.mle.cs2.iso <- function(dwis, design, sigma, S0, b=1, maxJ=4,
                               basis.depth=3, basis.random=F, def.alphas=2,
                               betas0=NULL, Ninner=10, Nmiddle.max=40,
                               display=F, pre.thres=0.15, alphas.const=NULL){

  # prelimilary fit: linear regression on log scale
  lmres <- lm.log(dwis=dwis, design=design, S0=S0, b=b, pd.correct=T)
  pre.fa <- FA(A=lmres$mat, eigval=lmres$eigval)

  # get alphas.const
  if (is.null(alphas.const)){
    if (lmres$eigval[1]>0){
      alphas.const <- 10^(round(-log10(lmres$eigval[1])))
      alphas.const <- 1/(lmres$eigval[1])*def.alphas
    } else {
      alphas.const <- b
    }
  }
  if (display){
    cat("Prelimilary fit with FA:", pre.fa, ".\n")
  }

  # isotropic estmation
  if (display){
    cat("Isotropic estimation...\n")
  }
  fit1 <- isotropic.mle(tdwi0=mean(dwis), dwis=dwis, b=b, sigma=sigma, S0=S0)

  # anisotropic estimation
  if (pre.fa<pre.thres){
    if (display){
      cat("Anisotropic estimation is not implemented because of small FA (prelimilary fit)!\n")
    }
    n.fiber <- 0
    fit2 <- NULL
  } else {
    if (display){
      cat("Anisotropic estimation...\n")
    }
    fit2 <- biased.mle.cs2(dwis=dwis, design=design, sigma=sigma, S0=S0,
                           b=b, maxJ=maxJ, basis.depth=basis.depth,
                           basis.random=basis.random, def.alphas=def.alphas,
                           betas0=betas0, Ninner=Ninner,
                           Nmiddle.max=Nmiddle.max, alphas.const=alphas.const)
    bic <- c(-2*fit1$ll + log(length(dwis)), fit2$bic)
    n.fiber <- which.min(bic)-1
    if (display){
      cat("BIC:", bic, "\n")
      cat("Number of fibers:", n.fiber, "\n")
    }
  }

  if (n.fiber == 0){
    vs <- NULL
  } else {
    vs <- fit2$mle.fit[[which.min(fit2$bic)]]$directions
  }

  return(list(vs=vs, iso.fit=fit1, aniso.fit=fit2, n.fiber=n.fiber))

}

###################################
#### similarity matrix for PAM ####
###################################
get.similar <- function(vecs){
  vecs <- matrix(vecs,nc=3)
  n <- nrow(vecs)
  angs <- abs(acos(pmin(abs(vecs %*% t(vecs)), 1)))
  return(angs)
}

##########################################
#### trimmed weighted instrinsic mean ####
##########################################
get.ang.ss <- function(angs, vecs, w, trim){
  v0 <- c(cos(angs[1])*sin(angs[2]), sin(angs[1])*sin(angs[2]), cos(angs[2]))
  ss <- (acos(pmin(abs(vecs%*%v0), 1)))^2
  return(sum(ss*w*(rank(-ss)>(trim*length(w)))))
}

# note that the vecs are assumed to be aligned in the similar directions
get.w.i.mean <- function(vecs, w, trim=0.05){
  vecs <- matrix(vecs, nc=3)
  # min arc distance estimate
  est0 <- apply(vecs, 2, function(x,w){mean(w*x)}, w=w)
  est0 <- est0/sqrt(sum(est0^2))
  angs0 <- spher.coord1(est0)[-1]

  res <- optim(par=angs0, fn=get.ang.ss, vecs=vecs, w=w, trim=trim, method="L-BFGS-B") # no boundary is needed
  angs <- res$par
  v0 <- c(cos(angs[1])*sin(angs[2]), sin(angs[1])*sin(angs[2]), cos(angs[2]))
  return(list(v0=v0, angs=angs, optim=res, est0=est0))
}


#### without trimming to fasten the speed
get.ang.ss2 <- function(angs, vecs, w){
  v0 <- c(cos(angs[1])*sin(angs[2]), sin(angs[1])*sin(angs[2]), cos(angs[2]))
  ss <- (acos(pmin(abs(vecs%*%v0), 1)))^2
  return(sum(ss*w))
}

# note that the vecs are assumed to be aligned in the similar directions
get.w.i.mean2 <- function(vecs, w){
  vecs <- matrix(vecs, nc=3)
  # min arc distance estimate
  est0 <- apply(vecs, 2, function(x,w){mean(w*x)}, w=w)
  est0 <- est0/sqrt(sum(est0^2))
  angs0 <- spher.coord1(est0)[-1]

  res <- optim(par=angs0, fn=get.ang.ss2, vecs=vecs, w=w, method="L-BFGS-B") # no boundary is needed
  angs <- res$par
  v0 <- c(cos(angs[1])*sin(angs[2]), sin(angs[1])*sin(angs[2]), cos(angs[2]))
  return(list(v0=v0, angs=angs, optim=res, est0=est0))
}


####################################
#### aligning direction vectors ####
####################################

#### for transforming the directions to the closet of the true ####
trans <- function(v0s, vecs){
  n <- nrow(vecs)
  temp <- vecs%*%v0s
  ind <- apply(temp, 1, function(x){x[which.max(abs(x))]<0})
  vecs[ind,] <- -vecs[ind,]
  return(vecs)
}

#### for only one v0
trans1 <- function(v0, vecs){
  n <- nrow(vecs)
  temp <- vecs%*%v0
  ind <- as.vector(temp<0)
  vecs[ind,] <- -vecs[ind,]
  return(vecs)
}

###################################
#### get spherical cooridnates ####
###################################
# v: a three dimensional unit vector

spher.coord <- function(v){
  r <- sqrt(sum(v^2))
  theta <- acos(v[3]/r)
  psi <- atan2(v[2], v[1])
  #psi <- psi+2*pi*(psi<0)
  return(c(r, theta, psi))
}

spher.coord1 <- function(v){
  r <- sqrt(sum(v^2))
  theta <- atan2(v[2],v[1])
  psi <- acos(v[3]/r)
  #psi <- psi+2*pi*(psi<0)
  return(c(r, theta, psi))
}

########################################################################################
########################################################################################
################################# direction smoothing ##################################
########################################################################################
########################################################################################

###########################################
#### voxel estimation for every voxels ####
###########################################
v.est <- function(dwi.obs, sigma, S0, b, grad.mat, braingrid, cpus=4, opt=list())
{
  if (length(S0)==1){
    S0 <- rep(S0, ncol(dwi.obs))
  }

  if (length(sigma)==1){
    sigma <- rep(sigma, ncol(dwi.obs))
  }

  # Fit all directions individually (voxel-wise)

  if (is.null(opt$maxJ)){
    opt$maxJ <- 4
  }
  if (is.null(opt$basis.depth)){
    opt$basis.depth <- 3
  }
  if (is.null(opt$basis.random)){
    opt$basis.random <- T
  }
  if (is.null(opt$def.alphas)){
    opt$def.alphas <- 2
  }
  if (is.null(opt$Ninner)){
    opt$Ninner <- 10
  }
  if (is.null(opt$Nmiddle.max)){
    opt$Nmiddle.max <- 50
  }
  if (is.null(opt$pre.thres)){
    opt$pre.thres <- 0.2
  }

  grad.mat <- diag(sqrt(1/diag(grad.mat%*%t(grad.mat)))) %*% grad.mat # standardization
  design <- generate.design.noS0(grad.mat)

  # container
  braindim <- dim(braingrid)[-1]
  n.voxel <- prod(braindim)
  templen <- n.voxel*opt$maxJ
  map <- array(dim=templen)
  rmap <- array(dim=n.voxel)
  n.fiber <- array(dim=n.voxel)
  n.fiber2 <- array(dim=templen) # replicated version of n.fiber
  vec <- array(dim=c(templen, 3))
  loc <- array(dim=c(templen, 3))


  # fit (for each voxel)
  registerDoParallel(cores=cpus)
  fitlist <- foreach(i=1:n.voxel, .verbose=T, .packages=c("quadprog",
                                                          "dwi.internals2",
                                                          "cluster")) %dopar%{
    icount <- 0
    out <- 1; class(out) <- "try-error"
    while ((inherits(out,"try-error"))&&(icount<=2)){
      out <- try(biased.mle.cs2.iso(dwis=dwi.obs[,i], design=design,
                                    sigma=sigma[i], S0=S0[i], b=b, maxJ=opt$maxJ,
                                    basis.depth=opt$basis.depth,
                                    basis.random=opt$basis.random,
                                    def.alphas=opt$def.alphas, betas0=NULL,
                                    Ninner=opt$Ninner, Nmiddle.max=opt$Nmiddle.max,
                                    display=F, pre.thres=opt$pre.thres,
                                    alphas.const=opt$alphas.const))
      icount <- icount + 1
    }
    if(inherits(out,"try-error")){
      cat(i,"\n")
    }
    out
  }

  # organize the result
  iind <- 1
  for (i in (1:n.voxel)){
    voxindex <- as.vector(arrayInd(i, braindim))
    n.fiber[i] <- fitlist[[i]]$n.fiber
    rmap[i] <- iind
    if (n.fiber[i]>=1){
      for (j in (1:n.fiber[i])){
        map[iind] <- i
        n.fiber2[iind] <- n.fiber[i]
        loc[iind,] <- braingrid[,voxindex[1], voxindex[2], voxindex[3]]
        vec[iind,] <- fitlist[[i]]$vs[j,]
        iind <- iind + 1
      }
    } else {
      map[iind] <- i
      n.fiber2[iind] <- n.fiber[i]
      loc[iind,] <- braingrid[,voxindex[1], voxindex[2], voxindex[3]]
      iind <- iind + 1
    }
  }

  # reduce the size to temp.n.layer
  m <- iind - 1
  map <- map[1:m]
  n.fiber2 <- n.fiber2[1:m]
  vec <- vec[1:m,]
  loc <- loc[1:m,]

  return(list(map=map, rmap=rmap, n.fiber2=n.fiber2, n.fiber=n.fiber, vec=vec, loc=loc, fitslist=fitlist, opt=opt))  
}

#################
#### weights ####
#################

# s0 = matrix of locations, each row is a vector of 3 dimension
# v0 = matrix of directions, each row is a vector of 3 dimension
# tens.sp = ?x3x3, array of anisotropic matrix used in anisotropic weight
#       e.g. s0=pre$loc, v0=pre$eig, ten.sp=pre$ten
# ind.not = matrix with number of row = length(mind), containing index that are not used
# pre = obj from tensor.est

get.raw.w <- function(s0, v0, ten.sp, pre, braingrid, xy.truncrange=1, z.truncrange=1, ind.not=NULL, dir.weight=T){

  s0 <- matrix(s0,nc=3)
  v0 <- matrix(v0,nc=3)
  if (nrow(s0)!=nrow(v0)){
    stop("s0 does not match with v0!\n")
  }

  ten.sp <- array(ten.sp, dim=c(length(ten.sp)/9,3,3))

  m <- length(pre$map)
  n <- nrow(s0)
  braindim <- dim(braingrid)[-1]
  n.voxel <- prod(braindim)

  if (is.vector(ind.not)){
    ind.not <- matrix(ind.not, nr=n, nc=length(ind.not))
  }
  tempxy <- -xy.truncrange: xy.truncrange
  tempz <- -z.truncrange: z.truncrange

  weight <- list()
  for (i in (1:n)){
    # create an index vector for considered tensors
    consider <- 1:m
    consider <- setdiff(consider, ind.not[i,])

    weight[[i]] <- list()
    
    # get neighborhood
    dis0 <- t(matrix(pre$loc, nc=3)) - s0[i,]
    ii <- which.min(apply(dis0^2, 2, sum))
    voxindex <- as.vector(arrayInd(pre$map[ii], braindim))

    xrange <- voxindex[1] + tempxy
    yrange <- voxindex[2] + tempxy
    zrange <- voxindex[3] + tempz

    xrange <- xrange[((xrange>=1)*(xrange<=braindim[1]))==1]
    yrange <- yrange[((yrange>=1)*(yrange<=braindim[2]))==1]
    zrange <- zrange[((zrange>=1)*(zrange<=braindim[3]))==1]

    index.nbhd <- ArrayIndex(braindim, xrange, yrange, zrange)
    index.nbhd2 <- get.map(index.nbhd, pre$rmap, pre$n.fiber)
    consider <- intersect(consider, index.nbhd2)
    consider1 <- consider[pre$n.fiber2[consider]>0]
    consider2 <- consider[pre$n.fiber2[consider]==0]
    aniso <- prod(!is.na(v0[i,]))
    weight[[i]]$aniso <- aniso
    if (aniso){
      weight[[i]]$consider <- consider1
      consider <- consider1
    } else {
      weight[[i]]$consider <- consider2
      consider <- consider2
    }

    # get spatial weight
    ten <- ten.sp[i,,]
    if (length(consider)){
      #dis <- t(matrix(pre$loc[consider,], nc=3)) - s0[i,]
      dis <- dis0[,consider]
      weight[[i]]$sqdists <- as.vector(rep(1,3)%*%((dis)*solve(ten,dis))*sum(diag(ten)))
    } else {
      weight[[i]]$sqdists <- numeric(0)
    }

    # get directional weight
    if (dir.weight){
      if (aniso){
        if (length(consider)){
          ang <- apply(t(matrix(pre$vec[consider,], nc=3))*v0[i,], 2, function(x){acos(min(abs(sum(x)),1))})
          weight[[i]]$ang <- ang
        } else {
          weight[[i]]$ang <- numeric(0)
        }
      }
    }
  }
  return(weight)
}

#################
#### get.map ####
#################

# get the internal index from the location index
get.map <- function(ind, rmap, n.fiber){
  n.ten <- n.fiber + (n.fiber==0)
  rep(rmap[ind], n.ten[ind]) + as.vector(unlist(sapply(n.ten[ind], function(x){1:x-1})))
}

#######################
#### Compute error ####
#######################
get.ten.err <- function(tens1, tens2, dis="logE"){
  p1 <- ncol(tens1)
  p2 <- ncol(tens2)
  temp <- matrix(nr=p1, nc=p2)
  for (i in (1:p1)){
    for (j in (1:p2)){
      if (dis=="logE"){
        temp[i,j] <- log.dist(vec2mat(tens1[,i]), vec2mat(tens2[,j]))
      } else if (dis=="affine"){
        temp[i,j] <- affine.dist(vec2mat(tens1[,i]), vec2mat(tens2[,j]))
      } else if (dis=="Euclid"){
        temp[i,j] <- norm(vec2mat(tens1[,i])-vec2mat(tens2[,j]), type="F")
      }
    }
  }
  ses <- array(dim=p1)
  min.inds <- matrix(nr=p1, nc=2)
  index1 <- 1:p1
  index2 <- 1:p2
  for (i in (1:p1)){
    min.ind <- arrayInd(which.min(temp), c(p1-i+1, p2-i+1))
    min.inds[i,] <- c(index1[min.ind[1]], index2[min.ind[2]])
    index1 <- index1[-min.ind[1]]
    index2 <- index2[-min.ind[2]]
    ses[i] <- temp[min.ind[1], min.ind[2]]^2
    if (p2-i>0){
      temp <- matrix(temp[-min.ind[1],-min.ind[2]], nc=p2-i)
    }
  }
  return(list(err=c(sum(ses), mean(ses), max(ses)), min.inds=min.inds, ses=ses))
}

logE.eval <- function(tensor.true, n.fiber, pre, out){
  mnf <- max(n.fiber)
  correct <- array(dim=(mnf+1))
  mse <- array(0,dim=(mnf+1))
  err.detail <- array(dim=c(3,length(n.fiber)))
  for (i in (1:(mnf+1))){
    nf <- i-1
    tempind <- which(n.fiber==nf)
    tempind2 <- which(pre$n.fiber==nf)
    tempind3 <- intersect(tempind, tempind2)
    correct[i] <- length(tempind3)/length(tempind)
    nf2 <- nf+(nf==0)
    iind <- get.map(tempind3, pre$rmap, pre$n.fiber)
    nse <- length(tempind3)*nf2
    ii <- 1
    temptens1 <- matrix(nr=6,nc=nf2)
    temptens2 <- matrix(nr=6,nc=nf2)
    for (j in (tempind3)){
      for (k in (1:nf2)){
        temptens1[,k] <- mat2vec(tensor.true[,,k,j])
        temptens2[,k] <- mat2vec(out$ten[iind[ii],,])
        ii <- ii + 1
      }
      err.detail[,j] <- get.ten.err(tens1=temptens1, tens2=temptens2, dis="logE")$err
      mse[i] <- mse[i] + err.detail[1,j]
    }
    mse[i] <- mse[i]/nse
  }
  return(list(mse=mse, correct=correct, err.detail=err.detail))
}

get.v0.err <- function(v1, v2){
  p1 <- ncol(v1)
  p2 <- ncol(v2)
  temp <- matrix(nr=p1, nc=p2)
  for (i in (1:p1)){
    for (j in (1:p2)){
      temp[i,j] <- (acos(min(abs(sum(v1[,i]*v2[,j])),1)))^2 # squared error
    }
  }
  ses <- array(dim=p1)
  min.inds <- matrix(nr=p1, nc=2)
  index1 <- 1:p1
  index2 <- 1:p2
  for (i in (1:p1)){
    min.ind <- arrayInd(which.min(temp), c(p1-i+1, p2-i+1))
    min.inds[i,] <- c(index1[min.ind[1]], index2[min.ind[2]])
    index1 <- index1[-min.ind[1]]
    index2 <- index2[-min.ind[2]]
    ses[i] <- temp[min.ind[1], min.ind[2]]
    if (p2-i>0){
      temp <- matrix(temp[-min.ind[1],-min.ind[2]], nc=p2-i)
    }
  }
  return(list(err=c(sum(ses), mean(ses), max(ses)), min.inds=min.inds, ses=ses))
}

v0.eval <- function(v0.true, n.fiber, pre, out){
  mnf <- max(n.fiber)
  correct <- array(dim=(mnf+1))
  mse <- array(0,dim=(mnf+1))
  err.detail <- array(dim=c(3,length(n.fiber)))
  for (i in (2:(mnf+1))){
    nf <- i-1
    tempind <- which(n.fiber==nf)
    tempind2 <- which(pre$n.fiber==nf)
    tempind3 <- intersect(tempind, tempind2)
    correct[i] <- length(tempind3)/length(tempind)
    nf2 <- nf+(nf==0)
    iind <- get.map(tempind3, pre$rmap, pre$n.fiber)
    nse <- length(tempind3)*nf2
    ii <- 1
    tempv1 <- matrix(nr=3,nc=nf2)
    tempv2 <- matrix(nr=3,nc=nf2)
    for (j in (tempind3)){
      for (k in (1:nf2)){
        tempv1[,k] <- v0.true[,k,j]
        tempv2[,k] <- eigen(out$ten[iind[ii],,])$vectors[,1]
        ii <- ii + 1
      }
      err.detail[,j] <- get.v0.err(v1=tempv1, v2=tempv2)$err
      mse[i] <- mse[i] + err.detail[1,j]
    }
    mse[i] <- mse[i]/nse
  }
  return(list(mse=mse, correct=correct, err.detail=err.detail))
}

v0.eval2 <- function(v0.true, n.fiber, pre, v0s){
  mnf <- max(n.fiber)
  correct <- array(dim=(mnf+1))
  mse <- array(0,dim=(mnf+1))
  err.detail <- array(dim=c(3,length(n.fiber)))
  for (i in (2:(mnf+1))){
    nf <- i-1
    tempind <- which(n.fiber==nf)
    tempind2 <- which(pre$n.fiber==nf)
    tempind3 <- intersect(tempind, tempind2)
    correct[i] <- length(tempind3)/length(tempind)
    nf2 <- nf+(nf==0)
    iind <- get.map(tempind3, pre$rmap, pre$n.fiber)
    nse <- length(tempind3)*nf2
    ii <- 1
    tempv1 <- matrix(nr=3,nc=nf2)
    tempv2 <- matrix(nr=3,nc=nf2)
    for (j in (tempind3)){
      for (k in (1:nf2)){
        tempv1[,k] <- v0.true[,k,j]
        tempv2[,k] <- v0s[iind[ii],]
        ii <- ii + 1
      }
      err.detail[,j] <- get.v0.err(v1=tempv1, v2=tempv2)$err
      mse[i] <- mse[i] + err.detail[1,j]
    }
    mse[i] <- mse[i]/nse
  }
  return(list(mse=mse, correct=correct, err.detail=err.detail))
}

##################
#### graphics ####
##################

plotv <- function(vecs, nonrandom=F, ind=rep(2,nrow(vecs)), xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1), const=0.1,...){
  if (nonrandom){
    dat <- data.frame(vecs)
    temp <- aggregate(rep(1, nrow(dat)), by=as.list(dat), FUN=sum)
    count <- (temp[,4]/sum(temp[,4]))^(1/3)*const
    color <- "black"
  } else {
    temp <- cbind(vecs,1)
    count <- (temp[,4]/sum(temp[,4]))^(1/3)*const
    color <- c("yellow", "black", "green")[ind]
  }

  plot3d(temp[,1:3], xlim=xlim, ylim=ylim, zlim=zlim, ...)
  spheres3d(temp[,1:3], radius=count, color=color)
}

#############################
#### Direction smoothing ####
#############################

#### for Silhouette: one cluster
get.connect2 <- function(dis, ang=0.1){
  if (nrow(dis)<=4){
    return(F)
  } else {
    cres <- pam(x=dis, k=2, diss=T)
    mdis <- min(dis[cres$clustering==1,cres$clustering==2])
    return(mdis>ang)
  }
}

get.connect3 <- function(dis, ang=0.1){
  n <- nrow(dis)
  if (n==1){
    return(list(ind=F, k=1, clustering=1, id.med=1))
  } else if (n==2){
    if (dis[1,2]<=ang){
      return(list(ind=F, k=1, clustering=c(1,1), id.med=1))
    } else {
      return(list(ind=F, k=2, clustering=c(1,2), id.med=c(1,2)))
    }
  } else if (n==3){
    cres <- pam(x=dis, k=2, diss=T)
    mdis <- min(dis[cres$clustering==1,cres$clustering==2])
    if (mdis<=ang){
      ccres <- pam(x=dis, k=1, diss=T)
      return(list(ind=F, k=1, clustering=c(1,1,1), id.med=ccres$id.med))
    } else {
      if (min(c(dis[1,2], dis[1,3], dis[2,3]))<=ang){
        return(list(ind=F, k=2, clustering=cres$clustering, id.med=cres$id.med))
      } else {
        return(list(ind=F, k=3, clustering=1:3, id.med=1:3))
      }
    }
  } else {
    cres <- pam(x=dis, k=2, diss=T)
    mdis <- min(dis[cres$clustering==1,cres$clustering==2])
    if (mdis<=ang){
      ccres <- pam(x=dis, k=1, diss=T)
      return(list(ind=F, k=1, clustering=rep(1,n), id.med=ccres$id.med))
    } else {
      return(list(ind=T))
    }
  }
}

#################################
#### Assign group to centers ####
#################################

#### for transforming the directions to the closet of id ####
labeling <- function(v0s, vecs){
  n <- nrow(vecs)
  temp <- vecs%*%t(v0s)
  label <- apply(temp, 1, function(x){which.max(abs(x))})
  return(label)
}

#############################
#### direction smoothing ####
#############################

v.smooth <- function(pre, h1, braingrid, xy.truncrange, z.truncrange,
                     thres.v0.w=0.2, K=5, leave=F, wobjs=NULL,
                     update.vox=NULL){
  n.voxel <- prod(dim(braingrid)[-1])
  if (is.null(update.vox)){
    update.vox <- 1:n.voxel
  }

  # ten.sp is used for anisotropic weight; Here, we set it to identity -> isotropic weight
  ten.sp <- array(0,dim=c(1,3,3))
  ten.sp[,1,1] <- 1
  ten.sp[,2,2] <- 1
  ten.sp[,3,3] <- 1

  v0s <- pre$vec

  update.ind <- array(F, dim=length(pre$n.fiber2))
  
  for (vox in (update.vox)){
    if (pre$n.fiber[vox]>0){
      #cat(vox,"\n")
      iind <- pre$rmap[vox] # if more than one tensor, the first one is used; Since, they share the same spatial weight, so it won't change the result. 

      ##
      if (!is.null(wobjs)){
        wobj <- wobjs[[vox]]
      } else {
        if (leave){
          wobj <- get.raw.w(s0=pre$loc[iind,], v0=pre$vec[iind,], ten.sp=ten.sp,
                            pre=pre, braingrid=braingrid,
                            xy.truncrange=xy.truncrange, z.truncrange=z.truncrange,
                            ind.not=iind+(0:(pre$n.fiber[vox]-1)), dir.weight=F)[[1]]
        } else {
          wobj <- get.raw.w(s0=pre$loc[iind,], v0=pre$vec[iind,], ten.sp=ten.sp,
                            pre=pre, braingrid=braingrid,
                            xy.truncrange=xy.truncrange, z.truncrange=z.truncrange,
                            ind.not=NULL, dir.weight=F)[[1]]
        }
      }
      ## only the spatial weight is needed
      vecs <- matrix(nr=length(wobj$consider), nc=3)
      locs <- matrix(nr=length(wobj$consider), nc=3)
      for (i in (1:length(wobj$consider))){
        vecs[i,] <- pre$vec[wobj$consider[i],]
        locs[i,] <- pre$loc[wobj$consider[i],]
      }

      ## get the weights
      w1 <- exp(-wobj$sqdists/(2*h1[iind]^2)) 
      # may want to adjust w3
      w <- w1
      w <- w/sum(w)

      ## thresholding the weights
      thres.rank <- sum(cumsum(sort(w)) <= thres.v0.w)
      thres.ind <- (rank(w)>thres.rank)
      vecs1 <- vecs[thres.ind,]
      vecs1 <- matrix(vecs1, nc=3)
      w1 <- w[thres.ind]
      consider1 <- wobj$consider[thres.ind]

      ## decide if more than 1 cluster is needed
      ### using connectivity
      dis <- get.similar(vecs1)
      connect.obj <- get.connect3(dis)
      ind1 <- connect.obj$ind

      if (ind1){
        K1 <- min(K, length(w1)-1)
        #dis <- get.similar(vecs1)
        cres.list <- list()
        sil.widths <- array(dim=K1)
        for (k in (2:K1)){
          cres.list[[k]] <- pam(x=dis, k=k, diss=T)
          sil.widths[k] <- cres.list[[k]]$silinfo$avg.width
        }
        sel.k <- which.max(sil.widths)
        ocres <- cres.list[[sel.k]]
        clustering <- ocres$clustering
        id.med <- ocres$id.med
      } else {
        sel.k <- connect.obj$k
        clustering <- connect.obj$clustering
        id.med <- connect.obj$id.med
      }

      ## labeling all directions
      #clustering0 <- labeling(vecs1[id.med,], vecs)

      ## get weighted trimmed geodesic mean
      nv0s <- matrix(nr=sel.k, nc=3)
      for (k in (1:sel.k)){
        lab <- (clustering==k)
        vecs2 <- trans1(vecs1[id.med[k],], matrix(vecs1[lab,],nc=3)) #####
        nv0s[k,] <- get.w.i.mean2(vecs2, w1[lab])$v0 # require "trans" as should be done beforehands
        #lab <- (clustering0==k)
        #vecs2 <- trans1(vecs1[id.med[k],], matrix(vecs[lab,],nc=3)) #####
        #nv0s[k,] <- get.w.i.mean2(vecs2, w[lab])$v0 # require "trans" as should be done beforehands
      }

      ## assign group
      tempdot <- abs(nv0s%*%t(matrix(pre$vec[iind+(0:(pre$n.fiber[vox]-1)),], nc=3)))
      index1 <- 1:nrow(tempdot)
      index2 <- 1:ncol(tempdot)
      while(length(tempdot)>0){
        tind <- arrayInd(which.max(tempdot),dim(tempdot))
        tind2 <- c(index1[tind[1]], index2[tind[2]])
        v0s[iind+tind2[2]-1,] <- nv0s[tind2[1],]
        update.ind[iind+tind2[2]-1] <- T
        tempdot <- tempdot[-tind[1], -tind[2]]
        index1 <- index1[-tind[1]]
        index2 <- index2[-tind[2]]
        if (length(tempdot)>0){
          tempdot <- matrix(tempdot,nr=length(index1))
        }
      }
    }
  }
  return(list(vec=v0s, update.ind=update.ind))
}


# with automatic bandwidth selection
# update.v0.clust2.bw
v.smooth.bw <- function(pre, braingrid, len=30, range1=c(0.01,0.1),
                        range2=c(0.1,1), cpus=1, xy.truncrange=100,
                        z.truncrange=100, thres.v0.w=0.2, K=5,
                        cv.method="trim-cv", exp.range=F){
  n.voxel <- prod(dim(braingrid)[-1])
  find1 <- (pre$n.fiber2==1)
  find2 <- (pre$n.fiber2>=2)

  if (!exp.range){
    h11s <- seq(range1[1], range1[2], len=len)
    h12s <- seq(range2[1], range2[2], len=len)
  } else {
    h11s <- exp(seq(log(range1[1]), log(range1[2]), len=len))
    h12s <- exp(seq(log(range2[1]), log(range2[2]), len=len))
  }
  err1 <- array(dim=len)
  err2 <- array(dim=len)

  # get wobjs
  wobjs <- list()
  ten.sp <- array(0,dim=c(1,3,3))
  ten.sp[,1,1] <- 1
  ten.sp[,2,2] <- 1
  ten.sp[,3,3] <- 1
  for (vox in (1:n.voxel)){
    if (pre$n.fiber[vox]>0){
      iind <- pre$rmap[vox] # if more than one tensors, the first one is used; Since, they share the same spatial weight, so it won't change the result. 
      wobjs[[vox]] <- get.raw.w(s0=pre$loc[iind,], v0=pre$vec[iind,], ten.sp=ten.sp,
                                pre=pre, braingrid=braingrid,
                                xy.truncrange=xy.truncrange, z.truncrange=z.truncrange,
                                ind.not=iind+(0:(pre$n.fiber[vox]-1)), dir.weight=F)[[1]]
    } 
  }

  registerDoParallel(cores=cpus)
  cv.res <- foreach(k=1:len, .verbose=T, .combine=cbind)%dopar%{
    h1 <- array(dim=length(pre$n.fiber2))
    h1[find1] <- h11s[k]
    h1[find2] <- h12s[k]

    clust2.res <- v.smooth(pre=pre, h1=h1, braingrid=braingrid, xy.truncrange=xy.truncrange, z.truncrange=z.truncrange, thres.v0.w=thres.v0.w, K=K, leave=T, wobjs=wobjs)
    ang.err <- apply(clust2.res$vec*pre$vec, 1, function(x){acos(min(abs(sum(x)),1))})
    find11 <- (find1*clust2.res$update.ind)==1
    find22 <- (find2*clust2.res$update.ind)==1
    if (cv.method=="cv"){
      eachout <- c(mean(ang.err[find11]^2),mean(ang.err[find22]^2), sum(find11), sum(find22))
    } else if (cv.method=="trim-cv") {
      eachout <- c(mean(ang.err[find11]^2, trim=0.05), mean(ang.err[find22]^2, trim=0.05), sum(find11), sum(find22))
    } else if (cv.method=="mad"){
      eachout <- c(median(ang.err[find11]), median(ang.err[find22]), sum(find11), sum(find22))
    }
    eachout
  }
  min.ind <- apply(cv.res[1:2,],1,which.min)
  if (!length(min.ind[[1]])){
    min.ind[[1]] <- 1
  }
  if (!length(min.ind[[2]])){
    min.ind[[2]] <- 1
  }
  min.ind <- unlist(min.ind)
  h1 <- array(dim=length(pre$n.fiber2))
  h1[find1] <- h11s[min.ind[1]]
  h1[find2] <- h12s[min.ind[2]]
  v0s.res <- v.smooth(pre=pre, h1=h1, braingrid=braingrid, xy.truncrange=xy.truncrange, z.truncrange=z.truncrange, thres.v0.w=thres.v0.w, K=K, leave=F, wobjs=NULL)
  return(list(vec=v0s.res$vec, update.ind=v0s.res$update.ind, h1=h1,
              min.ind=min.ind, h11s=h11s, h12s=h12s, cv.res=cv.res,
              warn1=(min.ind[1]==1)+(min.ind[1]==len),
              warn2=(min.ind[2]==1)+(min.ind[2]==len)))
}

##############################################################
##### Update number of fiber directions after v.smooth.bw ####
##############################################################
# pre from v.est
# res from v.smooth.bw or v.smooth or any list(vec, update.ind)
update.v.obj <- function(pre, res){
  iinds <- which(((!res$update.ind)*(pre$n.fiber2>0))==1)
  if (length(iinds)>0){
    nvox <- length(pre$rmap)
    tt <- tabulate(pre$map[iinds],nvox)
    voxs <- which(tt>0)
    nvoxs <- tt[voxs]


    n.fiber <- pre$n.fiber
    n.fiber[voxs] <- n.fiber[voxs] - nvoxs
    n.fiber2 <- pre$n.fiber2

    rmap <- pre$rmap
    for (i in (1:length(voxs))){
      vox <- voxs[i]
      tiind <- pre$rmap[vox]
      ii <- tiind:(tiind + max(0,pre$n.fiber[vox]-1))
      n.fiber2[ii] <- n.fiber2[ii] - nvoxs[i]
      if ((vox+1) <= nvox){
        if (n.fiber[vox]!=0){
          rmap[(vox+1):nvox] <- rmap[(vox+1):nvox] - nvoxs[i]
        } else {
          rmap[(vox+1):nvox] <- rmap[(vox+1):nvox] - nvoxs[i] + 1
        }
      }
    }

    iinds1 <- pre$rmap[voxs[n.fiber[voxs]==0]]
    iinds2 <- setdiff(iinds, iinds1)

    n.fiber2 <- n.fiber2[-iinds2]

    map <- pre$map[-iinds2]
    vec <- res$vec
    vec[iinds1,] <- NA
    vec <- vec[-iinds2,] ###
    loc <- pre$loc[-iinds2,]

    obj <- list(map=map, rmap=rmap, n.fiber2=n.fiber2, n.fiber=n.fiber,
                vec=vec, loc=loc)
    return(list(obj=obj, up=T, remove.iinds=iinds))
  } else {
    obj <- list(map=pre$map, rmap=pre$rmap, n.fiber2=pre$n.fiber2,
                n.fiber=pre$n.fiber, vec=res$vec, loc=pre$loc)
    return(list(obj=obj, up=F, remove.iinds=NULL))
  }
}

##########################
#### get S0 and sigma ####
##########################
# suppose to be voxel-wise
# use only u=(0,0,0)
# dwis0 is a vector
get.Ss <- function(dwis0){
  if (sum(dwis0)==0){
    return(list(S0=0, sigma=0))
  } else {
    diws0 <- as.vector(dwis0)
    n <- length(dwis0)
    par0 <- c(mean(dwis0), sd(dwis0))
    res <- optim(par0, get.Ss.fn, dwis0=dwis0, method="L-BFGS-B", lower=c(1e-10,1e-10))
    return(list(S0=res$par[1], sigma=res$par[2], optim=res))
  }
}

get.Ss.fn <- function(pars, dwis0){
  S0 <- pars[1]
  sigma2 <- pars[2]^2
  temp <- S0*dwis0/sigma2
  return(-sum(log(dwis0/sigma2) - (dwis0^2 + S0^2)/(2*sigma2) + log(besselIm(temp, 0, T)) + temp))
}

#################################
#### get S0 with given sigma ####
#################################
# suppose to be voxel-wise
# use only u=(0,0,0)
# dwis0 is a vector
get.Ss2 <- function(dwis0, sigma){
  if ((sum(dwis0)==0)||(sigma==0)){
    return(0)
  } else {
    diws0 <- as.vector(dwis0)
    n <- length(dwis0)
    par0 <- mean(dwis0)
    res <- optim(par0, get.Ss2.fn, dwis0=dwis0, sigma2=sigma^2, method="L-BFGS-B", lower=c(1e-10,1e-10))
    return(res$par)
  }
}

get.Ss2.fn <- function(S0, dwis0, sigma2){
  temp <- S0*dwis0/sigma2
  return(-sum(log(dwis0/sigma2) - (dwis0^2 + S0^2)/(2*sigma2) + log(besselIm(temp, 0, T)) + temp))
}


#########################
##### subset pre obj ####
#########################
# pre from v.est
# retain: an array containning the index of voxels that will be retained, mostly from ArrayIndex
subset.v.est <- function(pre, retain)
{
  n.fiber <- pre$n.fiber[retain]
  n <- length(retain)
  n1 <- sum(n.fiber+ (n.fiber==0))
  rmap <- cumsum(c(1,n.fiber+ (n.fiber==0)))[-(1+n)]
  map <- array(dim=n1)
  retain.iind <- array(dim=n1)
  retain.iind[rmap] <- pre$rmap[retain]
  map[rmap] <- 1:n
  for (iind in (1:n1)){
    if (is.na(map[iind])){
      map[iind] <- map[iind-1]
      retain.iind[iind] <- retain.iind[iind-1] + 1
    }
  }

  vec <- pre$vec[retain.iind,]
  loc <- pre$loc[retain.iind,]
  n.fiber2 <- pre$n.fiber2[retain.iind]
  fitslist <- lapply(retain, function(x){pre$fitslist[[x]]})

  return(list(map=map, rmap=rmap, n.fiber2=n.fiber2, n.fiber=n.fiber,
              vec=vec, loc=loc, fitslist=fitslist))
}
