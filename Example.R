rm(list=ls())

## Set time and space points
ns <- 114
nt <- 15

## Generate covariates
Xs <- rexp(ns,1/50000)
Xt <- exp(0.015*(0:(nt-1)))
Xm <- tcrossprod(Xs,Xt)
X <- apply(Xm,2,"jitter")

## set betas
beta0 <- 2
beta1 <- 0.08

## Mean Structure
mm <- beta0 + beta1*log(X)

## Set Spatial and temporal Structure

## Temporal Structure -- AR(1)

rho_t <- 0.8
D_t <- outer(1:nt,1:nt,"-")
Sigma_T <- rho_t^abs(D_t)
c_Sigma_T <- chol(Sigma_T)

## Spatial Structure -- CAR

rho_s <- 0.95

# Adjacency structure for MO
# Can be any other adjacency structure
AA <- read.csv("MO_Adjacency.csv")
Adj <- data.matrix(AA[,-1])

M <- diag(rowSums(Adj))

Sigma_s_inv <- M - rho_s*Adj

c_Sigma_s_inv <- chol(Sigma_s_inv)

## Generate U

U <- backsolve(c_Sigma_s_inv,(matrix(rnorm(ns*nt),ns,nt)%*%c_Sigma_T))

Lambda <- exp(mm + U)

## Generate Y 100 times
# nrep <- 100
# Yl_tr <- lapply(1:nrep,function(i){matrix(rpois(ns*nt,c(Lambda)),ns,nt)})

YY <- matrix(rpois(ns*nt,c(Lambda)),ns,nt)

## Censor Limit
Cl <- 10

## Censoring the data

C <- (YY < Cl)
cind <- which(YY < Cl, arr.ind = TRUE)
YY[cind] <- Cl
Yl_obs <- YY

## Run the function

source("MainFunc.R")

out <- CensArealCount(Yl_obs,log(X),C,Adj,cl=Cl,
                      niter=10000,nburn=50,nthin=10)

## TRACE PLOTS

plot(out$beta0,type="l",main="Intercept")
plot(out$beta1,type="l",main="Effect of log-Population")
plot(out$sigma2,type="l",main="sigma2")
plot(out$rhos,type="l",main="Spatial Association Measure")
plot(out$rhot,type="l",main="Temporal Correlation Coefficient")
