## MCMC FUNCTION

CensArealCount <- function(Yna, X, C, Adj,cl, 
                           beta_init=NULL, sigma2_init = NULL, rhos_init = NULL, rhot_init = NULL,
                           jitter=0,
                           beta0mn = 0, beta0sd = 10,
                           beta1mn = 0, beta1sd = 10,
                           a_s = 9, b_s = 1,
                           a_t = 5, b_t = 1,
                           a_sig = 0.1, b_sig = 0.1,
                           niter=1000, nburn=500, nthin=2){
  
  ns <- ncol(Yna)
  nt <- nrow(Yna)
  
  M <- diag(rowSums(Adj))
  
  Scindsl <- lapply(1:nt,function(i){which(C[i,]==1)})
  Sncindsl <- lapply(1:nt,function(i){which(C[i,]==0)})
  NAindsl <- lapply(1:nt,function(i){which(is.na(Yna[i,]))})
  
  cnvec <- as.numeric(summary(Scindsl)[,1])
  ncnvec <- as.numeric(summary(Sncindsl)[,1])
  nanvec <- as.numeric(summary(NAindsl)[,1])
  
  Y_c <- lapply(1:nt,function(i){Yna[i,Scindsl[[i]]]})
  Y_nc <- lapply(1:nt,function(i){Yna[i,Sncindsl[[i]]]})
  
  X_c <- lapply(1:nt,function(i){X[i,Scindsl[[i]]]})
  X_nc <- lapply(1:nt,function(i){X[i,Sncindsl[[i]]]})
  
  if(is.null(beta_init) | is.null(sigma2_init)){
    oo <- glm(unlist(Y_nc) ~ unlist(X_nc),family = 'poisson')
    
    beta_init <- oo$coefficients
    sigma2_init <- sd(oo$residuals)
  }
  
  if(is.null(rhos_init)){rhos_init <- a_s/(a_s + b_s)}
  
  if(is.null(rhot_init)){rhot_init <- a_t/(a_t + b_t)}
  
  if(jitter){
    beta_init <- jitter(beta_init)
    sigma2_init <- jitter(sigma2_init)
    rhos_init <- jitter(rhos_init)
    rhot_init <- jitter(rhot_init)
  }
  
  beta <- beta_init
  beta0 <- beta[1]
  beta1 <- beta[2]
  sigma2 <- sigma2_init
  rhos <- rhos_init
  rhot <- rhot_init
  Dtmat <- abs(outer(1:nt,1:nt,"-"))
  
  Sigs_inv <- M - rhos*Adj
  Sigs_inv_c <- chol(Sigs_inv)
  logdet_Sigs_inv <- 2*sum(log(diag(Sigs_inv_c)))
  Sigs <- solve(Sigs_inv)
  
  Sigt <- rhot^Dtmat
  diags <- c(1,rep(1+rhot^2,nt-2),1)
  # odiags <- rep(-rhot,nt-1)
  ldi <- cbind(2:nt,1:(nt-1))
  udi <- cbind(1:(nt-1),2:nt)
  Sigt_inv <- diag(diags)
  Sigt_inv[ldi] <- -rhot
  Sigt_inv[udi] <- -rhot
  Sigt_inv <- Sigt_inv/(1-(rhot^2))
  
  
  U <- t(chol(Sigt))%*%matrix(rnorm(ns*nt),nt,ns)%*%chol(Sigs)
  
  lambda <- exp(U + beta0 + X*beta1)
  
  lambda_c <- lapply(1:nt,function(i){lambda[i,Scindsl[[i]]]})
  lambda_nc <- lapply(1:nt,function(i){lambda[i,Sncindsl[[i]]]})
  
  Uquad <- sum((Sigt_inv%*%U)*(U%*%Sigs_inv))
  
  nsamp <- niter-nburn
  
  keep_beta0 <- rep(NA,nsamp)
  keep_beta1 <- rep(NA,nsamp)
  keep_sigma2 <- rep(NA,nsamp)
  keep_rhos <- rep(NA,nsamp)
  keep_rhot <- rep(NA,nsamp)
  keep_U <- array(NA,dim=c(nsamp,nt,ns))
  
  mh_beta0 <- 1
  mh_beta1 <- 1
  mh_sigma2 <- 1
  mh_rhos <- 1
  mh_rhot <- 1
  mh_U <- 1
  
  acc_beta0 <- acc_beta1 <- acc_sigma2 <- acc_rhos <- acc_rhot <- acc_U <- 0
  
  ppc <- function(x){ppois(cl,x,log.p=T)}
  
  cur_ll <- -do.call("sum",lambda_nc) + do.call("sum",Map("*",Y_nc,Map("log",lambda_nc))) + do.call("sum",Map("ppc",lambda_c))
  
  
  ## GO!
  for(i in 1:niter){
    for(ii in 1:nthin){
      
      ## Update beta0
      cur_l_beta0 <- cur_ll  - 0.5*((beta0 - beta0mn)/beta0sd)^2
      
      can_beta0 <- beta0 + mh_beta0*rnorm(1)
      
      can_lambda <- exp(can_beta0 - beta0)*lambda
      
      can_lambda_c <- lapply(1:nt,function(i){can_lambda[i,Scindsl[[i]]]})
      can_lambda_nc <- lapply(1:nt,function(i){can_lambda[i,Sncindsl[[i]]]})
      
      can_ll <- -do.call("sum",can_lambda_nc) + do.call("sum",Map("*",Y_nc,Map("log",can_lambda_nc))) + do.call("sum",Map("ppc",can_lambda_c))
      can_l_beta0 <- can_ll  - 0.5*((can_beta0 - beta0mn)/beta0sd)^2
      
      a_beta0 <- can_l_beta0 - cur_l_beta0
      if(log(runif(1)) < a_beta0){
        beta0 <- can_beta0
        lambda <- can_lambda
        lambda_nc <- can_lambda_nc
        lambda_c <- can_lambda_c
        cur_ll <- can_ll
        acc_beta0 <- acc_beta0 + 1/nthin
      }
      
      ## Update beta1
      
      cur_l_beta1 <- cur_ll - 0.5*((beta1 - beta1mn)/beta1sd)^2
      
      can_beta1 <- beta1 + mh_beta1*rnorm(1)
      
      can_lambda <- exp((can_beta1-beta1)*X)*lambda
      
      can_lambda_c <- lapply(1:nt,function(i){can_lambda[i,Scindsl[[i]]]})
      can_lambda_nc <- lapply(1:nt,function(i){can_lambda[i,Sncindsl[[i]]]})
      
      can_ll <- -do.call("sum",can_lambda_nc) + do.call("sum",Map("*",Y_nc,Map("log",can_lambda_nc))) + do.call("sum",Map("ppc",can_lambda_c))
      can_l_beta1 <- can_ll  - 0.5*((can_beta1 - beta1mn)/beta1sd)^2
      
      a_beta1 <- can_l_beta1 - cur_l_beta1
      if(log(runif(1))<a_beta1){
        beta1 <- can_beta1
        lambda <- can_lambda
        lambda_nc <- can_lambda_nc
        lambda_c <- can_lambda_c
        cur_ll <- can_ll
        acc_beta1 <- acc_beta1 + 1/nthin
      }
      
      ## Update sigma2
      
      sigma2 <- 1/rgamma(1,a_sig + ns*nt*0.5, b_sig + 0.5*Uquad)
      
      ## Update rhos
      
      cur_l_rhos <- 0.5*nt*logdet_Sigs_inv - 0.5*Uquad/sigma2 + (a_s-1)*log(rhos) + (b_s-1)*log(1-rhos)
      
      can_rhos <- 1/(1+exp(-log(rhos/(1-rhos)) - mh_rhos*rnorm(1)))
      
      can_Sigs_inv <- M - can_rhos*Adj
      can_Sigs_inv_c <- chol(can_Sigs_inv)
      can_logdet_Sigs_inv <- 2*sum(log(diag(can_Sigs_inv_c)))
      
      can_Uquad <- sum((Sigt_inv%*%U)*(U%*%can_Sigs_inv))
      
      can_l_rhos <- 0.5*nt*can_logdet_Sigs_inv - 0.5*can_Uquad/sigma2 + (a_s-1)*log(can_rhos) + (b_s-1)*log(1-can_rhos)
      
      a_rhos <- can_l_rhos - cur_l_rhos + log(can_rhos) + log(1-can_rhos) - log(rhos) - log(1-can_rhos)
      if(log(runif(1))<a_rhos){
        rhos <- can_rhos
        Sigs_inv <- can_Sigs_inv
        Sigs_inv_c <- can_Sigs_inv_c
        logdet_Sigs_inv <- can_logdet_Sigs_inv
        Uquad <- can_Uquad
        acc_rhos <- acc_rhos + 1/nthin
      }
      
      ## Update rhot
      
      cur_l_rhot <- -0.5*ns*(nt-1)*log(1-(rhot^2)) - 0.5*Uquad/sigma2 + (a_t-1)*log(rhot) + (b_t-1)*log(1-rhot)
      
      can_rhot <- 1/(1+exp(-log(rhot/(1-rhot)) - mh_rhot*rnorm(1)))
      
      can_Sigt <- can_rhot^Dtmat
      can_Sigt_inv <- diag(c(1,rep(1+can_rhot^2,nt-2),1))
      can_Sigt_inv[ldi] <- -can_rhot
      can_Sigt_inv[udi] <- -can_rhot
      can_Sigt_inv <- can_Sigt_inv/(1-(can_rhot^2))
      
      can_Uquad <- sum((can_Sigt_inv%*%U)*(U%*%Sigs_inv))
      
      can_l_rhot <- -0.5*ns*(nt-1)*log(1-(can_rhot^2)) - 0.5*can_Uquad/sigma2 + (a_t-1)*log(can_rhot) + (b_t-1)*log(1-can_rhot)
      
      a_rhot <- can_l_rhot - cur_l_rhot + log(can_rhot) + log(1-can_rhot) - log(rhot) - log(1-can_rhot)
      if(log(runif(1))<a_rhot){
        rhot <- can_rhot
        Sigt_inv <- can_Sigt_inv
        Uquad <- can_Uquad
        acc_rhot <- acc_rhot + 1/nthin
      }
      
      
      ## Update U
      
      cur_l_U <- cur_ll - 0.5*Uquad/sigma2
      
      can_U <- U + mh_U*matrix(rnorm(ns*nt),nt,ns)
      
      can_lambda <- exp(can_U - U)*lambda
      
      can_lambda_c <- lapply(1:nt,function(i){can_lambda[i,Scindsl[[i]]]})
      can_lambda_nc <- lapply(1:nt,function(i){can_lambda[i,Sncindsl[[i]]]})
      
      can_ll <- -do.call("sum",can_lambda_nc) + do.call("sum",Map("*",Y_nc,Map("log",can_lambda_nc))) + do.call("sum",Map("ppc",can_lambda_c))
      
      can_Uquad <- sum((Sigt_inv%*%can_U)*(can_U%*%Sigs_inv))
      
      can_l_U <- can_ll - 0.5*can_Uquad/sigma2
      
      a_U <- can_l_U - cur_l_U
      if(log(runif(1))<a_U){
        U <- can_U
        lambda <- can_lambda
        lambda_nc <- can_lambda_nc
        lambda_c <- can_lambda_c
        cur_ll <- can_ll
        Uquad <- can_Uquad
        acc_U <- acc_U + 1/nthin
      }
      
    }
    
    
    if(i < (nburn/2)){
      mh_beta0 <- ifelse(acc_beta0 < 0.3, mh_beta0*.8,ifelse(acc_beta0 > 0.5, mh_beta0*1.2,mh_beta0))
      mh_beta1 <- ifelse(acc_beta1 < 0.3, mh_beta1*.8,ifelse(acc_beta1 > 0.5, mh_beta1*1.2,mh_beta1))
      mh_sigma2 <- ifelse(acc_sigma2 < 0.3, mh_sigma2*.8,ifelse(acc_sigma2 > 0.5, mh_sigma2*1.2,mh_sigma2))
      mh_rhos <- ifelse(acc_rhos < 0.3, mh_rhos*.8,ifelse(acc_rhos > 0.5, mh_rhos*1.2,mh_rhos))
      mh_rhot <- ifelse(acc_rhot < 0.3, mh_rhot*.8,ifelse(acc_rhot > 0.5, mh_rhot*1.2,mh_rhot))
      mh_U <- ifelse(acc_U < 0.3, mh_U*.8,ifelse(acc_U > 0.5, mh_U*1.2,mh_U))
      
      acc_beta0 <- acc_beta1 <- acc_sigma2 <- acc_rhos <- acc_rhot <- 0
      acc_U <- 0
    }
    
    if(i > nburn){
      ind <- i - nburn
      keep_beta0[ind] <- beta0
      keep_beta1[ind] <- beta1
      keep_sigma2[ind] <- sigma2
      keep_rhos[ind] <- rhos
      keep_rhot[ind] <- rhot
      keep_U[ind,,] <- U
    }
  }
  
  ll <- list("beta0"=keep_beta0,"beta1"=keep_beta1,"sigma2"=keep_sigma2,"rhos"=keep_rhos,"rhot"=keep_rhot,"U"=keep_U)
  
  return(ll)
}

