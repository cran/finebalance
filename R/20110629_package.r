finebalance<-function(distance.matrix, z, f, type=c("total","order","chisq"),
    abs.diff=NULL, perc.diff=NULL, force=NULL){
    # distance.matrix       nt*nc matrix
    # z                     treatment (0,1 valued)
    # f                     vector of length nt+nc, factor to be balanced
    # type                  which criterion to minimize
    # abs.diff              scaler, the maximum absolute difference for each level of f
    # perc.diff             scaler, the maximum percentage difference for each level of f
    # force                 vector of length nt+nc, 0, 1 valued, 
    #                       1 means the individual must be matched,
    #                       0 means the individual can be either matched or not.
    # if force is provided, all the other will be ignored
    
    call <- match.call()
    type <- match.arg(type)

    # to make the entries in the distance matrix integer-valued
    reso<-(.Machine$integer.max/64 -2)/max(distance.matrix[is.finite(distance.matrix)])
    distance.matrix.original<-distance.matrix
    distance.matrix<-round(distance.matrix*reso)+1

    # number of treated and control
    nt <- dim(distance.matrix)[1]
    nc <- dim(distance.matrix)[2]
    if(nt!=sum(z==1)|nc!=sum(z==0)){
      stop("Distance matrix and treatment variable does not agree!")
    }
    
    # analyze each factor
    f <- as.factor(as.integer(f))
    n.level <- nlevels(f)
    n.treated<-table(f[z==1])
    n.control<-table(f[z==0])
    if(!is.null(force)){
      n.control.force<-table(f[z==0&force==1])
      n.control.not.force<-table(f[z==0&force!=1])
    }

    # keep record of the names of treated and control in the original matrix
    code<-c(rownames(distance.matrix),colnames(distance.matrix))
    rownames(distance.matrix)<-1:nt
    colnames(distance.matrix)<-nt+1:nc

    # calculate the upper and lower bound of number of matched controls
    if(is.null(force)){
      if(type=="total"){
        n.matched<-mintotal(n.treated,n.control, abs.diff, perc.diff)
      }
      if(type=="order"){
        n.matched<-minorder(n.treated,n.control, abs.diff, perc.diff)
      }
      if(type=="chisq"){
        n.matched<-minchisq(n.treated,n.control, abs.diff, perc.diff)
      }
      n.matched.up<-n.matched$n.matched.up
      n.matched.low<-n.matched$n.matched.low
    } else {
      if(type=="total"){
        n.matched<-mintotal.force(n.treated,n.control, abs.diff, perc.diff, n.control.force)
      }
      if(type=="order"){
        n.matched<-minorder.force(n.treated,n.control, abs.diff, perc.diff, n.control.force)
      }
      if(type=="chisq"){
        n.matched<-minchisq.force(n.treated,n.control, abs.diff, perc.diff, n.control.force)
      }
      n.matched.force.low<-n.matched$low
      excess<-n.matched$excess
    }

    # set up the starting node, end node, distance, and upper capacity for each edge
    # note that lower capacity must be 0 for each edge (have to adjust b accordingly)
    # b=out-in for each node
    # first layer: from treated to control
    startn <- rep(c(1:nt), nc)
    endn <- nt + rep(c(1:nc), rep(nt ,nc))
    dists <- as.vector(distance.matrix)
    ucap <- rep(1, (nt*nc))
    b <- c(rep(1, nt), rep(0, nc))
    if(!is.null(force)){
      b[nt+1:nc][force[z==0]==1] <- -1
    }
    # second layer: from control to stratum
    if(is.null(force)){
      startn <- c(startn, nt + 1:nc)
      endn <- c(endn, as.integer(f[z==0])+nt+nc)
      dists <- c(dists, rep(0, nc))
      ucap <- c(ucap, rep(1, nc))
      b <- c(b, -n.matched.low)
    } else {
      startn <- c(startn, (nt + 1:nc)[force[z==0]!=1])
      endn <- c(endn, as.integer(f[force!=1&z==0])+nt+nc)
      dists <- c(dists, rep(0, sum(force!=1&z==0)))
      ucap <- c(ucap, rep(1, sum(force!=1&z==0)))
      b <- c(b, -n.matched.force.low)
    }
    # third layer: from stratum to sink
    dists <- c(dists, rep(0, n.level))
    startn <- c(startn, nt+nc+1:n.level)
    endn <- c(endn, rep(nt+nc+n.level+1,n.level))
    if(is.null(force)){
      ucap <- c(ucap, n.matched.up-n.matched.low)
      b <- c(b, -nt+sum(n.matched.low))
    } else {
      ucap <- c(ucap, excess)
      b <- c(b, -nt+sum(n.matched.force.low)+sum(z==0&force==1))
    }

    fop <- .Fortran("relaxalg",
            as.integer(nc+nt+n.level+1),  #number of total nodes
            as.integer(length(startn)),
            as.integer(startn),
            as.integer(endn),
            as.integer(dists),
            as.integer(ucap),
            as.integer(b),
            x1=integer(length(startn)),
            crash1=as.integer(0),
            large1=as.integer(.Machine$integer.max/4),
            feasible1=integer(1),
            NAOK = FALSE,
            DUP = TRUE,
            PACKAGE = "optmatch")
    x <-fop$x1
    st=startn[x==1][1:nt]
    en=endn[x==1][1:nt][order(st)]
    matchpairnames=code[en]
    allnames<-rep(NA,nc+nt)
    allnames[z==1]<-code[1:nt]
    allnames[z==0]<-code[nt+1:nc]
    matchpair=match(matchpairnames,allnames)
    total.distance=sum(diag(distance.matrix.original[,matchpairnames]))
    # the matchpire[i]th person is matched to the ith treated
    # matchpairnames[i] is matched to the ith treated
    # matchtable gives the number of treated, number of control, number of matched control
    list("matchpair"=matchpair, matchpairnames=matchpairnames, 
      total.distance=total.distance, 
      "matchtable"=rbind(n.treated,n.control,table(f[matchpair])))
}

smahal<-function(z,X){
# Rank based Mahalanobis distance.  Prevents an outlier from
# inflating the variance for a variable, thereby decreasing its importance.
# Also, the variances are not permitted to decrease as ties
# become more common, so that, for example, it is not more important
# to match on a rare binary variable than on a common binary variable
# z is a vector, length(z)=n, with z=1 for treated, z=0 for control
# X is a matrix with n rows containing variables in the distance
  X<-as.matrix(X)
  n<-dim(X)[1]
  rownames(X)<-1:n
  k<-dim(X)[2]
  m<-sum(z)
  for (j in 1:k) X[,j]<-rank(X[,j])
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  out<-matrix(NA,m,n-m)
  Xc<-X[z==0,]
  Xt<-X[z==1,]
  rownames(out)<-rownames(X)[z==1]
  colnames(out)<-rownames(X)[z==0]
  library(MASS)
  icov<-ginv(cv)
  for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=TRUE)
  out
}

bound <- function(n.treated, n.control, abs.diff=NULL, perc.diff=NULL){
  differ <- n.control-n.treated
  sel <- differ>0
  less <- -sum(differ[!sel])       # how many we are in short of in total
  differ.sel <- differ[sel]        # how many more we have in each stratum
  n.treated.sel <- n.treated[sel]
  if(!is.null(abs.diff)&is.null(perc.diff)){
    if(sum(pmin(differ.sel,abs.diff))<less){
      warning("abs.diff invalid, abs.diff ignored")
    } else {
      differ.sel<-pmin(differ.sel,abs.diff)
    }
  }
  if(!is.null(perc.diff)&is.null(abs.diff)){
    if(sum(pmin(differ.sel,floor(n.treated.sel*perc.diff)))<less){
      warning("perc.diff invalid, perc.diff ignored")
    } else {
      differ.sel<-pmin(differ.sel, floor(n.treated.sel*perc.diff))
    }
  }
  if(!is.null(perc.diff)&!is.null(abs.diff)){
    if(sum(pmin(differ.sel,floor(n.treated.sel*perc.diff),abs.diff))<less){
      warning("the combination of abs.diff and perc.diff is invalid, both are ignored")
    } else {
      differ.sel<-pmin(differ.sel, floor(n.treated.sel*perc.diff), abs.diff)
    }
  }
  differ.sel
}

bound.force <- function(n.treated, n.control, abs.diff=NULL, perc.diff=NULL, n.force){
  differ <- n.control-n.treated
  sel <- differ>0
  if(sum(n.force)>sum(n.treated)){
    stop("two many forced")
  }
  if(!is.null(abs.diff)|!is.null(perc.diff)){
    warning("when force is present, abs.diff and perc.diff ignored")
  }
  if(all(n.force[sel]<=n.treated[sel])){
    less <- -sum(differ[!sel])
    low <- n.control
    low[!sel] <- n.control[!sel] - n.force[!sel]
    low[sel] <- n.treated[sel]-n.force[sel]
    excess <- differ
    excess[!sel] <- 0
  } else {
    less <- -sum(differ[!sel]) - sum((n.force-n.treated)[n.force>n.treated])
    if(less<0){
      stop("two many forces")
    }
    low <- n.control
    low[sel] <- pmax(n.treated[sel]-n.force[sel],0)
    low[!sel] <- n.control[!sel] - n.force[!sel]
    excess <- differ
    excess[!sel] <- 0
    excess[sel] <- n.control[sel] - pmax(n.treated[sel], n.force[sel])
  }
  list(low=low, excess=excess, less=less)
}

mintotal <- function(n.treated, n.control, abs.diff=NULL, perc.diff=NULL){
  differ <- n.control-n.treated
  sel <- differ>0
  less <- -sum(differ[!sel])       # how many we are in short of in total
  n.treated.sel <- n.treated[sel]
  n.control.sel <- n.control[sel]
  differ.sel <- bound(n.treated, n.control, abs.diff=abs.diff, perc.diff=perc.diff)
  n.matched.up <- n.control
  n.matched.low <- n.control
  n.matched.up[sel] <- n.treated.sel+differ.sel
  n.matched.low[sel] <- n.treated.sel
  list(n.matched.up=n.matched.up, n.matched.low=n.matched.low)
}

mintotal.force <- function(n.treated, n.control, abs.diff=NULL, perc.diff=NULL, n.force){
  out <- bound.force(n.treated, n.control, abs.diff=abs.diff, perc.diff=perc.diff, n.force)
  list(low=out$low, excess=out$excess)
}

minorder.force <- function(n.treated, n.control, abs.diff=NULL, perc.diff=NULL, n.force){
  out <- bound.force(n.treated, n.control, abs.diff=abs.diff, perc.diff=perc.diff, n.force)
  less <- out$less
  excess <- out$excess
  low <- out$low
  ans <- rep(0,length(excess))
  i=0
  j=0
  ind <- 0 #indicator
  while(i<less){
    loc <- (excess>0) & (n.force+low+ans==n.treated+j)
    if((i+sum(loc))>less){
      ind <- 1
      break
    } else{
      ans[loc] <- ans[loc]+1
      excess[loc] <- excess[loc]-1
      i <- i+sum(loc)
    }
    j=j+1
  }
  if(ind==1){
    ans.low <- ans
    ans.up <- ans
    ans.up[loc] <- ans.up[loc]+1
  } else{
    ans.low <- ans
    ans.up <- ans
  }
  list(low=low+ans.low, excess=ans.up-ans.low)
}

minchisq.force <- function(n.treated, n.control, abs.diff=NULL, perc.diff=NULL, n.force){
  out <- bound.force(n.treated, n.control, abs.diff=abs.diff, perc.diff=perc.diff, n.force)
  less <- out$less
  excess <- out$excess
  low <- out$low
  ans <- rep(0,length(excess))
  i=0
  ind <- 0 #indicator
  grad <- (2*low+2*n.force-2*n.treated+1)/n.treated
  grad[excess==0] <- Inf
  while(i<less){
    loc <- which(grad==min(grad))
    if((i+length(loc))>less){
      ind <- 1
      break
    } else{
      ans[loc] <- ans[loc]+1
      grad[loc] <- (2*ans[loc]+2*n.force[loc]+2*low[loc]-2*n.treated[loc]+1)/n.treated[loc]
      if(any(ans[loc]>=excess[loc])){
        grad[loc][ans[loc]>=excess[loc]] <- Inf
      }
      i <- i+length(loc)
    }
  }
  if(ind==1){
    ans.low <- ans
    ans.up <- ans
    ans.up[loc] <- ans.up[loc]+1
  } else{
    ans.low <- ans
    ans.up <- ans
  }
  list(low=low+ans.low, excess=ans.up-ans.low)
}

minorder <- function(n.treated, n.control, abs.diff=NULL, perc.diff=NULL){
  differ <- n.control-n.treated
  sel <- differ>0
  less <- -sum(differ[!sel])       # how many we are in short of in total
  n.treated.sel <- n.treated[sel]
  n.control.sel <- n.control[sel]
  differ.sel <- bound(n.treated, n.control, abs.diff=abs.diff, perc.diff=perc.diff)
  ans <- rep(0,length(differ.sel))
  i=0
  ind <- 0 #indicator
  while(i<less){
    loc <- differ.sel>0
    if((i+sum(loc))>less){
      ind <- 1
      break
    } else{
      ans[loc] <- ans[loc]+1
      differ.sel[loc] <- differ.sel[loc]-1
      i <- i+sum(loc)
    }
  }
  if(ind==1){
    ans.low <- ans
    ans.up <- ans
    ans.up[loc] <- ans.up[loc]+1
  } else{
    ans.low <- ans
    ans.up <- ans
  }
  n.matched.up <- n.control
  n.matched.low <- n.control
  n.matched.up[sel] <- n.treated.sel+ans.up
  n.matched.low[sel] <- n.treated.sel+ans.low
  list(n.matched.up=n.matched.up, n.matched.low=n.matched.low,max.diff=max(ans.up))
}

minchisq <- function(n.treated, n.control, abs.diff=NULL, perc.diff=NULL){
  differ <- n.control-n.treated
  sel <- differ>0
  less <- -sum(differ[!sel])       # how many we are in short of in total
  n.treated.sel <- n.treated[sel]
  n.control.sel <- n.control[sel]
  differ.sel <- bound(n.treated, n.control, abs.diff=abs.diff, perc.diff=perc.diff)
  ans <- rep(0,length(differ.sel))
  grad <- 1/n.treated.sel
  i=0
  ind <- 0
  while(i<less){
    loc <- which(grad==min(grad))
    if((i+length(loc))>less){
      ind <- 1
      break
    } else{
      ans[loc] <- ans[loc]+1
      grad[loc] <- (2*ans[loc]+1)/n.treated.sel[loc]
      if(any(ans[loc]>=differ.sel[loc])){
        grad[loc][ans[loc]>=differ.sel[loc]] <- Inf
      }
      i <- i+length(loc)
    }
  }
  if(ind==1){
    ans.low <- ans
    ans.up <- ans
    ans.up[loc] <- ans.up[loc]+1
  } else{
    ans.low <- ans
    ans.up <- ans
  }
  n.matched.up <- n.control
  n.matched.low <- n.control
  n.matched.up[sel] <- n.treated.sel+ans.up
  n.matched.low[sel] <- n.treated.sel+ans.low
  list(n.matched.up=n.matched.up, n.matched.low=n.matched.low)
}


