MiRKAT2<-function (y, X = NULL, Ks, out_type = "C", nperm = 999, method = "davies",distH0=NULL,returnH0=F) {

    n = length(y)
    if (any(is.na(y))) {
        ids = which(is.na(y))
        stop(paste("subjects", ids, "has missing response, please remove before proceed \n"))
    }
    if (is.null(X) == FALSE) {
        if (NROW(X) != length(y)) 
            stop("Dimensions of X and y don't match.")
    }
    if (class(Ks) == "matrix") {
        Ks = list(Ks)
    }
    if (class(Ks) == "list") {
        if ((any(lapply(Ks, "nrow") != n)) | (any(lapply(Ks, 
            "ncol") != n))) {
            stop("distance matrix need to be n x n, where n is the sample size \n ")
        }
        if (class(Ks) != "list") {
            stop("Distance needs to be a list of n x n matrices or a single n x n matrix \n")
        }
    }
    if (!is.null(X)) {
        if (any(is.na(X))) {
            stop("NAs in  covariates X, please impute or remove subjects which has missing covariates values")
        }
    }
    if (method == "moment" & n < 100 & out_type == "C") {
        warning("Continuous outcome: sample size < 100, p-value using moment matching can be inaccurate at tails, davies or permutation is recommended")
    }
    if (method == "moment" & n < 200 & out_type == "D") {
        warning("Continuous outcome: sample size < 200, p-value using moment matching can be inaccurate at tails, davies or permutation is recommended")
    }
    if (!(out_type %in% c("C", "D"))) {
        stop("Currently only continuous and Binary outcome are supported. Please choose out_type = \"C\" or \"D\" ")
    }
    if (out_type == "C") {
        re = MiRKAT_continuous(y, X = X, Ks = Ks, method = method, 
            nperm = nperm)
    }
    if (out_type == "D") {
	if(is.null(distH0)){
		if(returnH0){
			print("get H0 for final MiRKAT")
		}

		re = MiRKAT_binary2(y, X = X, Ks = Ks, method = method, nperm = nperm,returnH0=returnH0)
	}else{
		re = MiRKAT_binary2_final(y, X = X, Ks = Ks, method = method, nperm = nperm,distH0=distH0)
	}
    }
    return(re)
}


# For continuous outcome
permuted.index = function (n){
  out = sample.int(n, n)
  return(out)
}

getQ = function(K, res, s2){    
  Q = 1/s2*res %*% K %*% res
}

getLambda_davies = function(K, P0){
  PKP = P0 %*% K %*% P0
  ee = eigen(PKP, symmetric = T)         
  lambda0 = ee$values[ee$values > 1e-10]
  return(lambda0)    
}

getindivP_davies = function(Q, lambda0, n, px){
  
  if (length(lambda0) >= n-px){ 
    # In fact, it can not go bigger than n-p because P0 has rank n-p
    lambda = c(lambda0 - Q/(n-px))
    k = length(lambda)
    h = rep(1,k)
  }else{
    lambda = c(lambda0 - Q/(n-px), -Q/(n-px))
    k = length(lambda0)
    h = c(rep(1, k), n-px - k)
  }
  
  p_davies = davies(q = 0, lambda = lambda, h = h, acc = 1e-6)$Qq  
  p_davies = ifelse(p_davies < 0, 0, p_davies) 
  
  return(p_davies)  
}

getSat = function(Q, keppa_tlt, niu_tlt){
  p_sat= 1 - pchisq(Q/keppa_tlt, df = niu_tlt)
  p_sat = ifelse(p_sat<0, 0, p_sat)
  return(p_sat)
}

getParamSat = function(K, P0,px){
  n        = nrow(K)
  POK      = P0%*%K 
  e_tlt    = sum(diag(POK))/2
  Iss      = 0.5*(n-px)#.5*sum(diag(P0 %*% P0))
  W        = POK %*% P0
  Its      = .5*sum(diag(W))
  Itt      = .5*sum(diag(W %*% K))
  Itt_tlt  = Itt - Its^2/Iss
  niu_tlt  = 2*e_tlt^2/Itt_tlt
  keppa_tlt= Itt_tlt/e_tlt
  return(list(niu_tlt = niu_tlt, keppa_tlt = keppa_tlt))
}
# This is to permute the residuals 
getQsim_continuous = function(mod,  nperm, X1, Ks){
  res = resid(mod)
  n = length(res)
  perm    = sapply(1:nperm, function(x) permuted.index(n))
  y_star  = mod$fitted + matrix(res[perm],n,nperm)
  
  res_sim = qr.resid(qr(X1), y_star) 
  px = NCOL(X1)
  modelVar= function(x, px){sum((x - mean(x))^2)/(n - px)}  
  
  sigma2_sim= apply(res_sim, 2, modelVar, px) # Already the variance
  
  Q_sim   = sapply(1:length(Ks), function(j){
    sapply(1:nperm, function(i){
      res_sim[,i] %*% Ks[[j]] %*%res_sim[,i]/sigma2_sim[i]})
  })
  return(Q_sim)    
}

# For binary outcome 
Get_Var_Elements =function(m4,u1,u2){
  temp1 = u1^2 * u2^2
  a1    = sum(m4 * temp1)
  a2    = sum(u1^2) * sum(u2^2) - sum(temp1)
  a3    = sum(u1*u2)^2 - sum(temp1)
  return(a1+a2+2*a3)
}
getHm = function(Q,muQ, varQ, df){
  Q_corrected= (Q - muQ)*sqrt(2*df)/sqrt(varQ) + df
  p = 1 - pchisq(Q_corrected ,df = df)
  p = ifelse(p<0, 0, p)
  return(p)
}
getIndivP_hm = function(K, res, mu, D0, P0){
  Q   = t(res)%*% K %*% res
  # K1  = 1/D0 * P01 %*% K %*% t(1/D0 * P01 )
  K1 = P0 %*% (D0*t(D0*K)) %*% P0
  eK  = eigen(K1, symmetric = T)
  # Instead of matching the first two moments, match to the fourth moment
  # Code adapted from SKAT package
  lambda = eK$values[eK$values > 1e-10]
  U   = as.matrix(eK$vectors[,eK$values > 1e-10])
  p.m = length(lambda)
  m4  = (3*mu^2-3*mu +1)/(mu*(1-mu))
  
  zeta =rep(0,p.m)
  var_i=rep(0,p.m)
  varQ = 0
  
  for(i in 1:p.m){   # The diagonals
    temp.M1 = sum(U[,i]^2)^2 - sum(U[,i]^4)
    zeta[i] = sum(m4 * U[,i]^4) + 3* temp.M1 # because ( \sum .)^4, not ^2
    var_i[i]= zeta[i] - 1
  }
  
  if(p.m == 1){
    Cov_Mat = matrix(zeta* lambda^2, ncol=1,nrow=1)
  } else if(p.m > 1){
    Cov_Mat = diag(zeta* lambda^2)
    for(i in 1:(p.m-1)){
      for(j in (i+1):p.m){
        Cov_Mat[i,j] = Get_Var_Elements(m4,U[,i],U[,j])
        Cov_Mat[i,j] = Cov_Mat[i,j]* lambda[i]* lambda[j]
      }
    }
  }
  Cov_Mat       = Cov_Mat + t(Cov_Mat)
  diag(Cov_Mat) = diag(Cov_Mat)/2
  
  varQ = sum(Cov_Mat) - sum(lambda)^2
  muQ  = sum(lambda)
  lambda.new = lambda * sqrt(var_i)/sqrt(2)
  df         =  sum(lambda.new^2)^2/sum(lambda.new^4)
  Q_corrected= (Q - muQ)*sqrt(2*df)/sqrt(varQ) + df
  p_corrected= 1 - pchisq(Q_corrected ,df = df)
  
  p_corrected = ifelse(p_corrected <0, 0, p_corrected)
  return(list(p_hm= p_corrected, Q = Q, muQ = muQ, varQ = varQ, df = df))
}

getIndivP_binary = function(K, res, D0, px, P0){
  n = length(res)
  Q <- as.numeric(res %*% K %*% res) 
  PKP = P0 %*% (D0*t(D0 * K)) %*% P0 
  eP0 = c(rep(1, n-px), rep(0, px))
  ePKP = eigen(PKP, symmetric = T)$values
  lambda0 = ePKP - Q*eP0/n  # the MLE of s2 is 1, therefore, we should divide by .        
  lambda0 = lambda0[abs(lambda0) >= 1e-10]
  p = davies(0, lambda=lambda0, acc = 1e-6)$Qq
  p = ifelse(p < 0, 0, p)
  return(list(Q = Q, ePKP = ePKP, p = p ))
}

davies <- function(q,lambda,h = rep(1,length(lambda)),delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {
  
  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  
  out <- .C("qfc", lambdas=as.double(lambda), noncentral=as.double(delta), df=as.integer(h), r=as.integer(r), sigma=as.double(sigma), q=as.double(q), lim=as.integer(lim), acc=as.double(acc), trace=as.double(rep(0,7)), ifault=as.integer(0), res=as.double(0), PACKAGE="MiRKAT")
  
  out$res <- 1 - out$res
  
  return(list(trace=out$trace,ifault=out$ifault,Qq=out$res))
  
}


getOmnibous_p<-function(Ks,nperm,n,res,Qs){
  q_sim = sapply(1:nperm, function(i){
    ind <- sample(n)
    p1 = sapply(1:length(Ks), function(j) {
      Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res) # the adjusted is zero in this case 
      return(Q1)
    })  
  }) 
  
  q_sim = t(q_sim)
  Q_all = rbind(unlist(Qs), q_sim)
  p_all = 1 - (apply(Q_all, 2, base::rank)-1)/(nperm + 1)  # The smallest Q gets pvalue = 0 and the biggest one gets p value = 1
  p_perm = p_all[1,]
  minP_all= apply(p_all,1, min)
  p_final = base::rank(minP_all)[1]/(nperm  + 1)  
}

MiRKAT_binary2<-function (y, X = NULL, Ks, family = "binomial", nperm = 999, 
    method = "davies",returnH0=F) 
{
        #browser()
    n <- length(y)
    if (is.null(X)) {
        X1 <- matrix(rep(1, length(y)), ncol = 1)
    }
    else {
        X1 <- model.matrix(~., as.data.frame(X))
    }
    qX1 <- qr(X1)
    X1 <- X1[, qX1$pivot, drop = FALSE]
    X1 <- X1[, 1:qX1$rank, drop = FALSE]
    options(warn = 2)
    mod <- glm(y ~ X1 - 1, family = binomial)
    options(warn = 1)
    px = NCOL(X1)
    mu = mod$fitted.values
    res = y - mu
    w = mu * (1 - mu)
    D0 = sqrt(w)
    DX12 = D0 * X1
    P0 = diag(n) - DX12 %*% solve(t(DX12) %*% (DX12)) %*% t(DX12)
    if (method == "davies") {
        if (n < 50) {
            warning("For binary outcome and n < 50, p-value using davies method can be inaccurate at tails, permutation is recommended.")
        }
        S = sapply(Ks, getIndivP_binary, res, D0, px, P0)
        ps = as.numeric(unlist(S[3, ]))
        if (length(Ks) == 1) {
            return(indivP = ps)
        }
        eP0 = c(rep(1, n - px), rep(0, px))
        Qs = unlist(S[1, ])
        q_sim = sapply(1:nperm, function(i) {
            ind <- sample(n)
            p1 = sapply(1:length(Ks), function(j) {
                Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% 
                  res)
                return(Q1)
            })
        })
        q_sim = t(q_sim)
        Q_all = rbind(unlist(Qs), q_sim)
        p_all = 1 - (apply(Q_all, 2, base::rank) - 1)/(nperm + 1)
	#p_all<-apply(Q_all, 2, function(data){
	#	permu_p_MiRKATmethod(data[1],data[-1],lower=F)
	#})
        p_perm = p_all[1, ]
        minP_all = apply(p_all, 1, min)

        #p_final = base::rank(minP_all)[1]/(nperm + 1)
	p_final<-permu_p_MiRKATmethod(minP_all[1],minP_all[-1],lower=T)
        return(list(indivP = ps, omnibus_p = p_final))
    }
    if (method == "moment") {
        S = sapply(Ks, getIndivP_hm, res, mu, D0, P0)
        ps = as.numeric(unlist(S[1, ]))
        Qs = unlist(S[2, ])
        if (length(Ks) == 1) {
            return(indivP = ps)
        }
        q_sim = sapply(1:nperm, function(i) {
            ind <- sample(n)
            p1 = sapply(1:length(Ks), function(j) {
                Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% 
                  res)
                return(Q1)
            })
        })
        q_sim = t(q_sim)
        Q_all = rbind(unlist(Qs), q_sim)
        p_all = 1 - (apply(Q_all, 2, base::rank) - 1)/(nperm + 1)
	#p_all<-apply(Q_all, 2, function(data){
	#	permu_p_MiRKATmethod(data[1],data[-1],lower=F)
	#})
        p_perm = p_all[1, ]
        minP_all = apply(p_all, 1, min)
        #p_final = base::rank(minP_all)[1]/(nperm + 1)
	p_final<-permu_p_MiRKATmethod(minP_all[1],minP_all[-1],lower=T)
        return(list(indivP = ps, omnibus_p = p_final))
    }
    if (method == "permutation") {
        Qs = lapply(Ks, getQ, res, s2 = 1)
        q_sim = sapply(1:nperm, function(i) {
            ind <- sample(n)
            p1 = sapply(1:length(Ks), function(j) {
                Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% 
                  res)
                return(Q1)
            })
        })
        if (length(Ks) == 1) {
            p_perm = (sum(q_sim > Qs) + 1)/(nperm + 1)
            return(indivP = p_perm)
        }
        q_sim = t(q_sim)
        Q_all = rbind(unlist(Qs), q_sim)
        p_all = 1 - (apply(Q_all, 2, base::rank) - 1)/(nperm + 1)
	#p_all<-apply(Q_all, 2, function(data){
	#	permu_p_MiRKATmethod(data[1],data[-1],lower=F)
	#})
        p_perm = p_all[1, ]
        minP_all = apply(p_all, 1, min)
        #p_final = base::rank(minP_all)[1]/(nperm + 1)
	p_final<-permu_p_MiRKATmethod(minP_all[1],minP_all[-1],lower=T)

	  if(returnH0){
		distH0<-list(q_sim=q_sim)
	  }else{
		distH0<-NULL
	  }

        return(list(indivP = p_perm, omnibus_p = p_final,distH0=distH0))
    }
}


MiRKAT_binary2_final<-function (y, X = NULL, Ks, family = "binomial", nperm = 999, 
    method = "davies",distH0) 
{

  #browser()
  q_sim<-distH0$distH0$q_sim
    n <- length(y)
    if (is.null(X)) {
        X1 <- matrix(rep(1, length(y)), ncol = 1)
    }
    else {
        X1 <- model.matrix(~., as.data.frame(X))
    }
    qX1 <- qr(X1)
    X1 <- X1[, qX1$pivot, drop = FALSE]
    X1 <- X1[, 1:qX1$rank, drop = FALSE]
    options(warn = 2)
    mod <- glm(y ~ X1 - 1, family = binomial)
    options(warn = 1)
    px = NCOL(X1)
    mu = mod$fitted.values
    res = y - mu
    w = mu * (1 - mu)
    D0 = sqrt(w)
    DX12 = D0 * X1
    P0 = diag(n) - DX12 %*% solve(t(DX12) %*% (DX12)) %*% t(DX12)
    if (method == "davies") {
        if (n < 50) {
            warning("For binary outcome and n < 50, p-value using davies method can be inaccurate at tails, permutation is recommended.")
        }
        S = sapply(Ks, getIndivP_binary, res, D0, px, P0)
        ps = as.numeric(unlist(S[3, ]))
        if (length(Ks) == 1) {
            return(indivP = ps)
        }
        eP0 = c(rep(1, n - px), rep(0, px))
        Qs = unlist(S[1, ])
        Q_all = rbind(unlist(Qs), q_sim)
        p_all = 1 - (apply(Q_all, 2, base::rank) - 1)/(nperm + 1)
        p_perm = p_all[1, ]
        minP_all = apply(p_all, 1, min)
        #p_final = base::rank(minP_all)[1]/(nperm + 1)
	p_final<-permu_p_MiRKATmethod(minP_all[1],minP_all[-1],lower=T)
        return(list(indivP = ps, omnibus_p = p_final))
    }
    if (method == "moment") {
        S = sapply(Ks, getIndivP_hm, res, mu, D0, P0)
        ps = as.numeric(unlist(S[1, ]))
        Qs = unlist(S[2, ])
        if (length(Ks) == 1) {
            return(indivP = ps)
        }
        Q_all = rbind(unlist(Qs), q_sim)
        p_all = 1 - (apply(Q_all, 2, base::rank) - 1)/(nperm + 1)
        p_perm = p_all[1, ]
        minP_all = apply(p_all, 1, min)
        #p_final = base::rank(minP_all)[1]/(nperm + 1)
	p_final<-permu_p_MiRKATmethod(minP_all[1],minP_all[-1],lower=T)
        return(list(indivP = ps, omnibus_p = p_final))
    }
    if (method == "permutation") {
        Qs = lapply(Ks, getQ, res, s2 = 1)
        if (length(Ks) == 1) {
            p_perm = (sum(q_sim > Qs) + 1)/(nperm + 1)
            return(indivP = p_perm)
        }
        Q_all = rbind(unlist(Qs), q_sim)
        p_all = 1 - (apply(Q_all, 2, base::rank) - 1)/(nperm + 1)
        p_perm = p_all[1, ]
        minP_all = apply(p_all, 1, min)
        #p_final = base::rank(minP_all)[1]/(nperm + 1)
	p_final<-permu_p_MiRKATmethod(minP_all[1],minP_all[-1],lower=T)
        return(list(indivP = p_perm, omnibus_p = p_final))
    }
}