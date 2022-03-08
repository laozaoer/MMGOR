Data_object=function(t,X,Z,b,r,beta,gamma,theta,U,H){
  Ht=H(t)
  if(r>0){
    S=(1+r*Ht*exp(sum(beta*X)+sum(gamma*Z)+theta*b))^(-1/r)
  }
  else{
    S=exp(-Ht*exp(sum(beta*X)+sum(gamma*Z)+theta*b))
  }
  
  return(S-U)
}
Generate_singleT=function(X,Z,b,r,beta,gamma,theta,H){
  U=runif(1,0,1)
  tresult=nleqslv::nleqslv(x=0.1,fn=Data_object,X=X,Z=Z,b=b,r=r,beta=beta,gamma=gamma,theta=theta,U=U,H=H)
  return(tresult$x)
}
Generate_T=function(X,Z,b,r,beta,gamma,theta,n,ni,H){
  result=list()
  length(result)=n
  for (i in 1:n) {
    result[[i]]=rep(0,ni[i])
    for (j in 1:ni[i]) {
      result[[i]][j]=Generate_singleT(X[[i]][j,],Z[i,],b[i],r,beta,gamma,theta,H)
    }
  }
  return(result)
}
data_for_est=function(r,beta,gamma,theta,n,H){
  betadim=length(beta)
  gammadim=length(gamma)
  Z=matrix(runif(n*gammadim,-1,1),nrow = n,ncol = gammadim)
  mi=rep(0,n)
  b=rnorm(n,0,1)
  for(i in 1:n){
    mi[i]=extraDistr::rtpois(1,exp(1.7),a=1,b=8)
  }
  C=list()
  length(C)=n+1
  for(i in 1:n){
    C[[i]]=runif(mi[i],0,1)
  }
  X=list()
  length(X)=n
  for (i in 1:n) {
    X[[i]]=matrix(runif(mi[i]*betadim,-1,1),nrow = mi[i],ncol=betadim)
  }
  
  
  
  rawC=suppressWarnings(Generate_T(X,Z,b,r,beta,gamma,theta,n,mi,H))
  lowC=0
  upC=quantile(as.numeric(as.character(unlist(rawC))),probs = 0.85)
  
  for(i in 1:n){
    C[[i]]=runif(mi[i],lowC,upC)
  }
  C[[n+1]]=upC
  Delta=list()
  length(Delta)=n
  for (i in 1:n) {
    Delta[[i]]=rep(0,mi[i])
    for (j in 1:mi[i]) {
      if(rawC[[i]][j]<=C[[i]][j]){
        Delta[[i]][j]=1
      }
    }
  }
  
  
  return(list(X=X,Z=Z,n=n,ni=mi,r=r,Delta=Delta,C=C))
}


GOR_MM=function(Delta,X,Z,n,ni,r,C,knotsnum,order,cluster.ind=TRUE,pen.ind=FALSE,lambda = 0,itermax=500,tol=1e-7,quadnum=30){
  
  if(cluster.ind){
    if(is.null(X)){
      betadim=0
      gammadim=ncol(Z)
      lowC=0
      upC=C[[n+1]]
      blC <- list()
      length(blC) <- n
      knots <- seq(0,1  , length.out = (knotsnum + 2))
      knots=knots[3:length(knots)-1]
      for (i in 1:n) {
        blC[[i]]=t(splines2::ibs((C[[i]]-lowC)/(upC-lowC),knots = knots,degree=order,Boundary.knots = c(0,1),intercept = TRUE))
      }
      bspl1.1 <- fda::create.bspline.basis(rangeval = c(0,1),breaks = c(0,knots,1),norder = (order+1))
      R <- fda::bsplinepen(bspl1.1, Lfdobj = 1)
      
      myrules=gaussquad::hermite.h.quadrature.rules(quadnum,normalized=FALSE)
      myrules=as.matrix(myrules[[quadnum]])
      initial_value=rep(0,(betadim+gammadim+knotsnum+order+2))
      X=list()
      length(X)=n
      for (i in 1:n) {
        X[[i]]=as.matrix(runif(ni[i],-1,1))
      }
      parest=MainFunc(initial_value,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,itermax,tol,pen.ind,lambda,R)
      hessian1=numDeriv::hessian(func=testquadrature1current,x=parest[,1],rules=myrules,Delta=Delta,
                                 X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=FALSE,lambda=0,R=R)
      
      hessian=numDeriv::hessian(func=testquadrature1current,x=parest[,1],rules=myrules,Delta=Delta,
                                X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=pen.ind,lambda=lambda,R=R)
      df=sum(diag(hessian1%*%(solve(hessian,tol=1e-40))))
      var=diag(solve(hessian,tol=1e-40))
      result=cbind(parest[(1:(betadim+gammadim+1)),1],var[1:(betadim+gammadim+1)])
      result[(betadim+gammadim+1),1]=exp(result[(betadim+gammadim+1),1])
      result[(betadim+gammadim+1),2]=(result[(betadim+gammadim+1),1])^2*result[(betadim+gammadim+1),2]
      result=as.data.frame(result)
      colnames(result)=c("Est","SE")
      log_likelihood=testquadrature1current(parest[,1],rules=myrules,Delta=Delta,
                                            X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=FALSE,lambda=0,R=R)
      
      
      for (i in ((betadim+1):(gammadim+betadim))) {
        rownames(result)[i]=paste("gamma",i-betadim,sep="_")
      }
      rownames(result)[betadim+gammadim+1]="theta"
      result[,2]=sqrt(result[,2])
      result=round(result,digits = 3)
      AIC=-log_likelihood+(1+df/n)/(1-(df+2)/n)
      
    }else if(is.null(Z)){
      betadim=ncol(X[[1]])
      gammadim=0
      lowC=0
      upC=C[[n+1]]
      blC <- list()
      length(blC) <- n
      knots <- seq(0,1  , length.out = (knotsnum + 2))
      knots=knots[3:length(knots)-1]
      for (i in 1:n) {
        blC[[i]]=t(splines2::ibs((C[[i]]-lowC)/(upC-lowC),knots = knots,degree=order,Boundary.knots = c(0,1),intercept = TRUE))
      }
      bspl1.1 <- fda::create.bspline.basis(rangeval = c(0,1),breaks = c(0,knots,1),norder = (order+1))
      R <- fda::bsplinepen(bspl1.1, Lfdobj = 1)
      
      myrules=gaussquad::hermite.h.quadrature.rules(quadnum,normalized=FALSE)
      myrules=as.matrix(myrules[[quadnum]])
      initial_value=rep(0,(betadim+gammadim+knotsnum+order+2))
      Z=as.matrix(runif(n,-1,1))
      parest=MainFunc(initial_value,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,itermax,tol,pen.ind,lambda,R)
      hessian1=numDeriv::hessian(func=testquadrature1current,x=parest[,1],rules=myrules,Delta=Delta,
                                 X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=FALSE,lambda=0,R=R)
      
      hessian=numDeriv::hessian(func=testquadrature1current,x=parest[,1],rules=myrules,Delta=Delta,
                                X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=pen.ind,lambda=lambda,R=R)
      df=sum(diag(hessian1%*%(solve(hessian,tol=1e-40))))
      var=diag(solve(hessian,tol=1e-40))
      result=cbind(parest[(1:(betadim+gammadim+1)),1],var[1:(betadim+gammadim+1)])
      result[(betadim+gammadim+1),1]=exp(result[(betadim+gammadim+1),1])
      result[(betadim+gammadim+1),2]=(result[(betadim+gammadim+1),1])^2*result[(betadim+gammadim+1),2]
      result=as.data.frame(result)
      colnames(result)=c("Est","SE")
      log_likelihood=testquadrature1current(parest[,1],rules=myrules,Delta=Delta,
                                            X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=FALSE,lambda=0,R=R)
      for(i in 1:betadim){
        rownames(result)[i]=paste("beta",i,sep = "_")
      }
      
      rownames(result)[betadim+gammadim+1]="theta"
      result[,2]=sqrt(result[,2])
      result=round(result,digits = 3)
      AIC=-log_likelihood+(1+df/n)/(1-(df+2)/n)
      
      
    }else{
      betadim=ncol(X[[1]])
      gammadim=ncol(Z)
      
      lowC=0
      upC=C[[n+1]]
      blC <- list()
      length(blC) <- n
      knots <- seq(0,1  , length.out = (knotsnum + 2))
      knots=knots[3:length(knots)-1]
      for (i in 1:n) {
        blC[[i]]=t(splines2::ibs((C[[i]]-lowC)/(upC-lowC),knots = knots,degree=order,Boundary.knots = c(0,1),intercept = TRUE))
      }
      bspl1.1 <- fda::create.bspline.basis(rangeval = c(0,1),breaks = c(0,knots,1),norder = (order+1))
      R <- fda::bsplinepen(bspl1.1, Lfdobj = 1)
      
      myrules=gaussquad::hermite.h.quadrature.rules(quadnum,normalized=FALSE)
      myrules=as.matrix(myrules[[quadnum]])
      initial_value=rep(0,(betadim+gammadim+knotsnum+order+2))
      parest=MainFunc(initial_value,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,itermax,tol,pen.ind,lambda,R)
      
      hessian1=numDeriv::hessian(func=testquadrature1current,x=parest[,1],rules=myrules,Delta=Delta,
                                 X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=FALSE,lambda=0,R=R)
      
      hessian=numDeriv::hessian(func=testquadrature1current,x=parest[,1],rules=myrules,Delta=Delta,
                                X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=pen.ind,lambda=lambda,R=R)
      df=sum(diag(hessian1%*%(solve(hessian,tol=1e-40))))
      var=diag(solve(hessian,tol=1e-40))
      result=cbind(parest[(1:(betadim+gammadim+1)),1],var[1:(betadim+gammadim+1)])
      result[(betadim+gammadim+1),1]=exp(result[(betadim+gammadim+1),1])
      result[(betadim+gammadim+1),2]=(result[(betadim+gammadim+1),1])^2*result[(betadim+gammadim+1),2]
      result=as.data.frame(result)
      colnames(result)=c("Est","SE")
      log_likelihood=testquadrature1current(parest[,1],rules=myrules,Delta=Delta,
                                            X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=FALSE,lambda=0,R=R)
      for(i in 1:betadim){
        rownames(result)[i]=paste("beta",i,sep = "_")
      }
      for (i in ((betadim+1):(gammadim+betadim))) {
        rownames(result)[i]=paste("gamma",i-betadim,sep="_")
      }
      rownames(result)[betadim+gammadim+1]="theta"
      result[,2]=sqrt(result[,2])
      result=round(result,digits = 3)
      AIC=-log_likelihood+(1+df/n)/(1-(df+2)/n)
      
    }
  }else{
    
    
    betadim=ncol(X[[1]])
    gammadim=ncol(Z)
    
    lowC=0
    upC=C[[n+1]]
    blC <- list()
    length(blC) <- n
    knots <- seq(0,1  , length.out = (knotsnum + 2))
    knots=knots[3:length(knots)-1]
    for (i in 1:n) {
      blC[[i]]=t(splines2::ibs((C[[i]]-lowC)/(upC-lowC),knots = knots,degree=order,Boundary.knots = c(0,1),intercept = TRUE))
    }
    bspl1.1 <- fda::create.bspline.basis(rangeval = c(0,1),breaks = c(0,knots,1),norder = (order+1))
    R <- fda::bsplinepen(bspl1.1, Lfdobj = 1)
    

    myrules=cbind(c(0,0),c(1/2,1/2))
    initial_value=rep(0,(betadim+gammadim+knotsnum+order+1))
    parest=MainFuncnocluster(initial_value,myrules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,itermax,tol,pen.ind,lambda,R)
    
    hessian1=numDeriv::hessian(func=testquadrature1currentnocluster,x=parest[,1],rules=myrules,Delta=Delta,
                               X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=FALSE,lambda=0,R=R)
    
    hessian=numDeriv::hessian(func=testquadrature1currentnocluster,x=parest[,1],rules=myrules,Delta=Delta,
                              X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=pen.ind,lambda=lambda,R=R)
    df=sum(diag(hessian1%*%(solve(hessian,tol=1e-40))))
    var=diag(solve(hessian,tol=1e-40))
    result=cbind(parest[(1:(betadim+gammadim)),1],var[1:(betadim+gammadim)])
    result=as.data.frame(result)
    colnames(result)=c("Est","SE")
    log_likelihood=testquadrature1currentnocluster(parest[,1],rules=myrules,Delta=Delta,
                                                   X=X,Z=(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,penind=FALSE,lambda=0,R=R)
    for(i in 1:betadim){
      rownames(result)[i]=paste("beta",i,sep = "_")
    }
    for (i in ((betadim+1):(gammadim+betadim))) {
      rownames(result)[i]=paste("gamma",i-betadim,sep="_")
    }
    result[,2]=sqrt(result[,2])
    result=round(result,digits = 3)
    AIC=-log_likelihood+(1+df/n)/(1-(df+2)/n)
    
  }
  
  
  
  
  return(list(result=result,AIC=AIC))
}
