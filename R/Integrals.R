getliu <- function(covs,R){
  ttt  <- R%*%covs
  ttt2 <- ttt%*%ttt#제곱
  ttt3 <- ttt2%*%ttt#세제곱
  ttt4 <- ttt2%*%ttt2#
  c1   <- sum(diag(ttt))
  c2   <- sum(diag(ttt2))
  c3   <- sum(diag(ttt3))
  
  
  c4   <- sum(diag(ttt4))
  
  muQ  <- c1
  sigQ <- sqrt(2*c2)
  s1   <- c3/(c2*sqrt(c2))
  s2   <- c4/c2^2
  s1_2 <- s1*s1
  
  a    <- ifelse(s1_2>s2, 1/(s1-sqrt(s1_2 - s2)), 1/sqrt(s2))
  delta<- ifelse(s1_2>s2, s1*a^3 - a^2, 0)
  l    <- a^2 - 2*delta
  muX  <- l + delta
  sigX <- sqrt(2*(l+2*delta))
  
  list(muX=muX, sigX=sigX, muQ=muQ, sigQ=sigQ, l=l, delta=delta)
}

qliu   <- function(pval,lius,Lower=FALSE){
  (qchisq(pval,df=lius$l,ncp=lius$delta,lower=Lower)-lius$muX)*lius$sigQ/lius$sigX + lius$muQ
}

pliu   <- function(quan,lius,Lower=FALSE){
  tt <- (quan-lius$muQ)/lius$sigQ*lius$sigX + lius$muX
  pchisq(tt,df=lius$l,ncp=lius$delta,lower=Lower)
}

get_K<-function(distance){
  Ds<-distance^2
  p<-dim(Ds)[1]
  agr1_K<-(diag(p)-matrix(1/p,ncol=p,nrow=p))
  K<- (-1/2)*agr1_K%*%Ds%*%agr1_K
}


res_hom <- try(unlist(getHomFARVAT(nphe,nOTU,rhos,S,W,SIGMA,yAy,fileout="log/")),TRUE)
# W<-1
# i<-1
getHomFARVAT <- function(nphe,nOTU,rhos,S,W,SIGMA,yAy,fileout){  
  getR <-function(rho,I,II) W%*%((1-rho)*I+rho*II%*%t(II))%*%W
  IIq <- matrix(rep(1,nphe),ncol=1)
  Im <- diag(nOTU)
  IIm <- matrix(rep(1,nOTU),ncol=1)
  
  rRho <- qmins <- pvals <- stats <- seq(along=rhos)
  lius <- as.list(stats)
  for(i in seq(along=rhos)){
    R <- getR(rhos[i],Im,IIm)
    lius[[i]] <- getliu(SIGMA,R)
    tmp_llqS<-t(IIq)%*%S
    stats[i] <- tmp_llqS%*%R%*%t(tmp_llqS)/sum(yAy) 
    pvals[i] <- pliu(stats[i],lius[[i]])
  }  					
  
  e <- eigen(SIGMA,T)
  evectors <- Re(e$vectors)
  evalues <-Re(e$values)
  o <- evalues<0
  evalues[o] <- 0
  Z <- evectors%*%diag(sqrt(evalues),length(evalues),length(evalues))%*%t(evectors)%*%W
  zbar <- matrix(rowMeans(Z), ncol=1)
  zbar2 <- t(zbar)%*%zbar
  rRho <- nOTU^2*rhos*zbar2 + (1-rhos)/zbar2*sum((t(zbar)%*%Z)^2)
  M <- zbar%*%solve(t(zbar)%*%zbar)%*%t(zbar)
  t_lmZ<-(Im-M)%*%Z
  IMV <-  t_lmZ%*%t(t_lmZ)
  IMV2 <- IMV%*%IMV
  muQ <- sum(diag(IMV))
  muQ2 <- sum(diag(IMV2))
  sigXi <- 2*sqrt( sum( diag(t(Z)%*%M%*%Z%*%t(Z)%*%(Im-M)%*%Z) ) )
  sigQ <- sqrt( 2*muQ2 + sigXi^2 )
  id <- which.min(pvals)
  rhoMin <- rhos[id]
  pMin <- pvals[id]
#   write(paste0("zbar2 = ",zbar2,", sigXi = ",sigXi,", sigQ = ",sigQ),paste("./",fileout,'res.log',sep=''),append=TRUE)
#   write(paste0("rhoMin = ",rhoMin,", pMin = ",pMin),paste("./",fileout,'res.log',sep=''),append=TRUE)
  
  for(i in seq(along=rhos))  qmins[i] <- qliu(pMin,lius[[i]])
  #browser() 
  ff <- function(x,qmins,rRho,rhos,muQ,sigQ,sigXi,lius){
    fff <- function(xxx,qmins,rRho,rhos) min((qmins-rRho*xxx)/(1-rhos))
    dx <- unlist(lapply(x,fff,qmins=qmins,rRho=rRho,rhos=rhos))
    deltax <- (dx-muQ)*sqrt(sigQ^2-sigXi^2)/sigQ+muQ
    pliu(deltax,lius,Lower=TRUE)*dchisq(x,df=1)
  }
  
  liu <- getliu(Im,IMV)
  pvalSKATO <- 1-unlist(integrate(ff,0,Inf,qmins=qmins,rRho=rRho,rhos=rhos,muQ=muQ,sigQ=sigQ,sigXi=sigXi,lius=liu))$value
  res <- c(stats[c(1,7)],pvals[c(1,7)],rhoMin,pvalSKATO)
  names(res)<-c(paste0("stat",c(1,7)),paste0("pval",c(1,7)),"rhoMin","pvalSKATO")
  return(res)
}

