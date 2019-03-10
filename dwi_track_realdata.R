##ref: Wong, R.K.W., T.C.M. Lee, D. Paul, and J. Peng. Fiber direction estimation, smoothing and tracking in diffusion MRI (with Discussion)  (2016). The Annals of Applied Statistics, 10(3): 1137-1156.

########################################
### functions for fiber tracking (slightly modified version to fit the format of the real data outputs)
#######################################
fiber.track <- function(iind, eig, loc, map, rmap, n.fiber, xgrid.sp, ygrid.sp,
                        zgrid.sp, braingrid, max.line=1000, nproj=1, thres.ang=0.5235988,
                        vorient=c(1,1,1))
{

  braindim <- dim(braingrid)[-1]
  nvox <- prod(braindim)
  dimens <- c(xgrid.sp, ygrid.sp, zgrid.sp)

  path.voxel <- array(dim=max.line)
  path.dir <- array(dim=c(max.line, 3))
  path.in <- array(dim=c(max.line, 3))
  path.change <- array(dim=max.line)
  path.iind <- array(dim=max.line)
  pass.vox <- NULL
  pass.dir <- NULL
  pass.dis <- NULL
  pass.pdis <- NULL # perpendicular distance
  

  # initialization of the tract
  path.voxel[1] <- map[iind]  ## tract starting voxel
  path.dir[1,] <- eig[iind,]  ## tract starting direction
  path.in[1,] <- loc[iind,]   ## tract starting location
  path.change[1] <- T     
  path.iind[1] <- iind  ## tract starting index (index of all fiber directions in the ROI)

  ii <- 1
  while ((ii<max.line)){
    
    fio <- fiber.in.out(inc=path.in[ii,]-loc[path.iind[ii],], direct=path.dir[ii,], dimens=dimens)
    path.in[ii+1,] <- fio$outc + loc[path.iind[ii],]

    # for previous pass.dis and pass.pdis, using the previous "change"
    if ((!path.change[ii])&&(n.fiber[path.voxel[ii]]>0)){
      pass.pdis <- c(pass.pdis, dist.line(loc[path.iind[ii],], path.in[ii,], path.in[ii+1,]))
      pass.dis <- c(pass.dis, sqrt(sum((path.in[ii,]-path.in[ii+1,])^2)))
    }

    # determine which voxel it is going to
    next.vox <- get.out.vox(fio$index, path.voxel[ii], braindim=braindim, vorient=vorient)

    if (is.na(next.vox)){
      break
    }

    # determine if we should stop
    pro.res <- project.proceed(inc0=path.in[ii+1,], vox0=next.vox,
                               dir0=path.dir[ii,], loc, eig, rmap, n.fiber,
                               braindim, dimens, nproj=nproj,
                               thres.ang=thres.ang, vorient=vorient)
    change <- pro.res$first
    good <- pro.res$last

    if (!good){
      break
    }

    # update voxel
    path.voxel[ii+1] <- next.vox

    # update dir, iind and change
    if (n.fiber[next.vox]<=1){
      path.iind[ii+1] <- rmap[next.vox]
      path.change[ii+1] <- change
      if (change){
        path.dir[ii+1,] <- eig[path.iind[ii+1],]
      } else {
        path.dir[ii+1,] <- path.dir[ii,]
        if (n.fiber[next.vox]==1){
          pass.vox <- c(pass.vox,next.vox)
          pass.dir <- rbind(pass.dir, path.dir[ii,])
        }
      }
    } else {
      # thresholding rule -> determine stop or not, and within the thresholding rule, choose the closest

      if (change){
        # decide which directions
        tiind <- rmap[next.vox]
        chosen <- which.max(abs(eig[tiind+(0:(n.fiber[next.vox]-1)),]%*%path.dir[ii,]))
        path.iind[ii+1] <- tiind+chosen-1
        path.dir[ii+1,] <- eig[path.iind[ii+1],]
        path.change[ii+1] <- T
      } else {
        path.iind[ii+1] <- rmap[next.vox]
        path.change[ii+1] <- F
        path.dir[ii+1,] <- path.dir[ii,]
        pass.vox <- c(pass.vox,next.vox)
        pass.dir <- rbind(pass.dir, path.dir[ii,])
      }
    }

    # align directions
    path.dir[ii+1,] <- sign(sum(path.dir[ii+1,]*path.dir[ii,]))*path.dir[ii+1,]

    ii <- ii+1
  }

  if (ii<max.line){
    path.in <- path.in[1:(ii+1),]
    path.iind <- path.iind[1:ii]
    path.dir <- path.dir[1:ii,]
    path.change <- path.change[1:ii]
  }
  return(list(inloc=path.in, dir=path.dir, iinds=path.iind, change=path.change, pvox=pass.vox, pdir=pass.dir, pdis=pass.dis, ppdis=pass.pdis))
}

fiber.in.out <- function(inc, direct, dimens){
  # assuming inc, outc are coordinates with the center of the voxel being (0,0,0)
  if (sum(dimens==0)){
    stop("directions has zero component, not yet supported! Please modify fiber.in.out\n")
  }
  # compute the distance of the current fiber directon to each face of the current voxel
  # create a box with x,y,z dimens and divide by travel speed(direct)
  tempdiff <- (round(cbind(dimens/2-inc,-inc-dimens/2),5)/direct)  ## Hao: add round5
  
  # 
  tempdiff[tempdiff==-Inf]=Inf
  
  
  # get which axis is the current fiber direction hitting face of the current voxel first
  # 1:x  2:y  3:z
  index1 <- which.min(diag(tempdiff[,2-(direct>=0)]))  
  # which direction it is hitting 1:positive  2:negative
  index <- c(index1, (2-(direct>0))[index1])
  const <- tempdiff[index[1],index[2]]
  outc <- round(inc + const*direct,5)     
  
  return(list(outc=outc, index=as.vector(index)))
}

get.out.vox <- function(index, cvox, braindim, vorient){
  # convert current voxel vector index to the 3d index in ROI
  cvoxindex <- as.vector(arrayInd(cvox, braindim))
  if (index[2]==1){
    # positive sides
    cvoxindex[index[1]] <- cvoxindex[index[1]] + vorient[index[1]]
  } else {
    # negative sides
    cvoxindex[index[1]] <- cvoxindex[index[1]] - vorient[index[1]]
  }
  if ((cvoxindex[index[1]]<1)||(cvoxindex[index[1]]>braindim[index[1]])){
    return(NA)
  } else {
    return(ArrayIndex(braindim, cvoxindex[1], cvoxindex[2], cvoxindex[3]))
  }
}

proceed <- function(vox0, dir0, eig, rmap, n.fiber, thres.ang=0.5235988){
  good <- T
  if (n.fiber[vox0]==0){
    good <- F
  } else if (n.fiber[vox0]==1) {
    good <- acos(min(abs(eig[rmap[vox0],]%*%dir0),1))<thres.ang
  } else {
    good <- as.logical(sum(as.vector(acos(pmin(abs(eig[rmap[vox0]+(0:(n.fiber[vox0]-1)),]%*%dir0),1)))<thres.ang))
  }
  return(good)
}


project.proceed <- function(inc0, vox0, dir0, loc, eig, rmap, n.fiber,
                            braindim, dimens, nproj=2, thres.ang=0.5235988,
                            vorient=c(1,1,1)){
  # first
  first <- proceed(vox0, dir0, eig, rmap, n.fiber, thres.ang)

  cont <- !first
  num <- 1
  tinc <- inc0
  vox <- vox0
  last <- first
  iind <- rmap[vox]
  while (cont && (num <= nproj)){
    fio <- fiber.in.out(inc=tinc-loc[iind,], direct=dir0, dimens=dimens)
    tinc <- fio$outc + loc[iind,]

    # determine which voxel it is going to
    vox <- get.out.vox(fio$index, vox, braindim=braindim, vorient=vorient)
    if (is.na(vox)){
      last <- F
      break
    }
    iind <- rmap[vox]
    last <- proceed(vox, dir0, eig, rmap, n.fiber, thres.ang)
    if (last){
      break
    }
    num <- num + 1
  }
  return(list(first=first, last=last))
}

dist.line <- function(x, y, z){
  # the distance between x and line(y,z)
  v <- y-z; v <- v/sqrt(sum(v^2))
  w <- y-x
  return(sqrt(sum((w - sum(w*v)*v)^2)))
}

############################################
#### tracking and updating of the brain ####
############################################
# vorient: if the voxel orientation does not match the physical orientation, make it -1
v.track <- function(v.obj, xgrid.sp, ygrid.sp, zgrid.sp, braingrid,
                    max.line=100, nproj=1, elim=T, elim.thres=1,
                    thres.ang=0.5235988, vorient=rep(1,1,1))
{
  tracks1 <- list()
  tracks2 <- list()
  all.pvox <- NULL
  all.pdir <- NULL
  all.pdis <- NULL
  all.ppdis <- NULL
  n.use.iind <- array(0, dim=length(v.obj$n.fiber2))
  n.iinds <- array(0,dim=length(v.obj$n.fiber2))
  lens <- array(0, dim=length(v.obj$n.fiber2))
  
  for (iind in which(v.obj$n.fiber2>0)){
    cat(iind,"\n")
    tracks1[[iind]] <- fiber.track(iind=iind, eig=v.obj$vec, loc=v.obj$loc,
                                   map=v.obj$map, rmap=v.obj$rmap,
                                   n.fiber=v.obj$n.fiber, xgrid.sp=xgrid.sp,
                                   ygrid.sp=ygrid.sp, zgrid.sp=zgrid.sp, braingrid=braingrid,
                                   max.line=max.line, nproj=nproj, thres.ang=thres.ang,
                                   vorient=vorient)

    tracks2[[iind]] <- fiber.track(iind=iind, eig=-v.obj$vec, loc=v.obj$loc,
                                   map=v.obj$map, rmap=v.obj$rmap,
                                   n.fiber=v.obj$n.fiber, xgrid.sp=xgrid.sp,
                                   braingrid=braingrid,
                                   ygrid.sp=ygrid.sp, zgrid.sp=zgrid.sp,
                                   max.line=max.line, nproj=nproj, thres.ang=thres.ang,
                                   vorient=vorient)
    
    n.use.iind[tracks1[[iind]]$iinds] <- n.use.iind[tracks1[[iind]]$iinds] + 1
    n.use.iind[tracks2[[iind]]$iinds] <- n.use.iind[tracks2[[iind]]$iinds] + 1
    n.use.iind[iind] <- n.use.iind[iind] - 1
    n.iinds[iind] <- length(union(tracks1[[iind]]$iinds, tracks2[[iind]]$iinds))
    lens[iind] <- get.fdis(tracks1[[iind]]$inloc) + get.fdis(tracks2[[iind]]$inloc)

    if (length(all.pdis)!=length(all.pvox)){
      break
    }
  }
  #len.ord <- order(n.iinds, decreasing=T)
  len.ord <- order(lens, decreasing=T)
  if (max(lens[n.iinds<=1])> elim.thres){
    cat("elim.thres is too small: it should be set at least", max(lens[n.iinds<=1]),"\n") 
  }

  if (elim){
    update.ind <- rep(T, length(v.obj$n.fiber2))
    update.ind[as.logical((v.obj$n.fiber2==0)+(lens<=elim.thres))] <- F
    nv.obj <- update.v.obj(v.obj, list(vec=v.obj$vec, update.ind=update.ind))$obj
  } else {
    nv.obj <- v.obj
    update.ind <- rep(T, length(v.obj$n.fiber2))
    update.ind[as.logical((v.obj$n.fiber2==0)+(lens<=elim.thres))] <- F
  }
  sorted.iinds <- (1:length(v.obj$n.fiber2))[len.ord]
  sorted.update.ind <- update.ind[len.ord]
  return(list(v.obj=nv.obj, tracks1=tracks1, tracks2=tracks2, n.iinds=n.iinds,
              n.use.iind=n.use.iind, update.ind=update.ind, sorted.iinds=sorted.iinds,
              sorted.update.ind=sorted.update.ind, lens=lens))
}

#### get fiber distance ####
get.fdis <- function(inloc){
  if (nrow(inloc)==2){
    return(sum((apply(inloc,2,diff))^2))
  } else {
    return(sum(apply(apply(inloc,2,diff),1,function(x){sum(x^2)})))
  }
}

#############################
#### additional function ####
#############################

#### color-coded drawing of the function ####
get.fib.col <- function(vecs)
{
  vecs <- abs(vecs)
  if (length(vecs)==3){
    vecs <- matrix(vecs, nc=3)
  }
  return(rgb(red=vecs[,1], green=vecs[,2], blue=vecs[,3], alpha=1))
}

plot.fib <- function(loc, vec, x_lim=NULL, y_lim=NULL, z_lim=NULL)
{
  if(is.null(x_lim)&&is.null(y_lim)&&is.null(z_lim)){
    # repeat location for get color for each line segment
    loc1 <- matrix(rep(loc, each=2), nc=3)
    loc1 <- loc1[-c(1, nrow(loc1)), ]
    
    # get color
    col1 <- get.fib.col(vec)
    col1 <- rep(col1, each=2)
    
    lines3d(loc1, lwd=1, col=col1)
  }else{
    epsi = 1e-5
    index_t = 1:dim(loc)[1]
    dele_t = which((((loc[,2]>y_lim[1]-epsi)&(loc[,2]<y_lim[2]+epsi))==0)|(((loc[,1]>x_lim[1]-epsi)&(loc[,1]<x_lim[2]+epsi))==0)|(((loc[,3]>z_lim[1]-epsi)&(loc[,3]<z_lim[2]+epsi))==0))
    if(length(dele_t)==length(index_t)){

    }else if((length(dele_t)>0)&&(length(index_t)-length(dele_t)>1)){
      index_u = index_t[-dele_t]
      loc = loc[index_u,]
      
      index_v = integer(0)
      for(i in 1:(length(index_u)-1)){
        if((index_u[i+1]-index_u[i])==1){
          index_v = c(index_v,index_u[i])
        }
      }
      vec = vec[index_v,]
      # repeat location for get color for each line segment
      loc1 <- matrix(rep(loc, each=2), nc=3)
      loc1 <- loc1[-c(1, nrow(loc1)), ]
      
      # get color
      col1 <- get.fib.col(vec)
      col1 <- rep(col1, each=2)
      
      lines3d(loc1, lwd=1, col=col1)
    } else if(length(dele_t)==0) {
      loc1 <- matrix(rep(loc, each=2), nc=3)
      loc1 <- loc1[-c(1, nrow(loc1)), ]
      
      # get color
      col1 <- get.fib.col(vec)
      col1 <- rep(col1, each=2)
      
      lines3d(loc1, lwd=1, col=col1)
    }
  }
}

#### generate pre object from the eigen-vector of single fiber model ####
gen.pre.single <- function(vec, loc, fa, fa.thres)
{
  ind  <- (fa>fa.thres)
  n <- nrow(vec)
  vec[!ind,] <- NA
  nfib <- rep(1,n)
  nfib[!ind] <- 0
  return(list(map=1:n, rmap=1:n, n.fiber2=nfib, n.fiber=nfib, vec=vec, loc=loc))
}
