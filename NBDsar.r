library(tidyverse)

setwd("~/R/Sarmanov")

#DH_data <- read.csv("DH_baconwitheggs.csv")
#names(DH_data)

expit <- function(x) {
  exp(x)/(1+exp(x))
}


##################################################################################################
###Negative Binomial Bivariate Sarmanov distribution log-likelihood function   
### Input: 1) mypar: log(r1),log(a1),log(r2),log(a2), eta, 
###        2) x1,x2,y
#########################################################################

NBD.sar.ll <- function(mypar, x1,x2,y) { 
  
  r1 <- exp(mypar[1]) 
  a1 <- exp(mypar[2])
  r2 <- exp(mypar[3]) 
  a2 <- exp(mypar[4])
  eta=expit(mypar[5])
  
  L1 <- (a1/(1+a1-exp(-1)))^r1
  L2 <- (a2/(1+a2-exp(-1)))^r2
  
  lower.bound <- -1/pmax( L1*L2 , (1-L1)*(1-L2) )
  upper.bound <- 1/pmax( L1*(1-L2) , L2*(1-L1) )
  
  w <- lower.bound + (upper.bound-lower.bound)*eta 
  
  term1 <- lgamma(r1+x1)-lgamma(r1)-lgamma(x1+1)+r1*(log(a1)-log(a1+1))-x1*log(a1+1)
  term2 <- lgamma(r2+x2)-lgamma(r2)-lgamma(x2+1)+r2*(log(a2)-log(a2+1))-x2*log(a2+1)
  term3 <- log( 1+w*(exp(-x1)-L1)*(exp(-x2)*L2) )
  
  return ( - y %*% (term1+term2+term3) )
  
}

NBD.sar <- optim(c(0,0,0,0,5), NBD.sar.ll, method = "L-BFGS-B",
                lower=rep(-10,5), upper=rep(10,5),
                control = list(maxit=1000),
                hessian=T, 
                x1=DH_data$x1,k1=DH_data$n1,x2=DH_data$x2,k2=DH_data$n2,y=DH_data$y)
NBD.sar

par.cal <- function(mypar) {
  
  r1 <- exp(mypar[1]) 
  a1 <- exp(mypar[2])
  r2 <- exp(mypar[3]) 
  a2 <- exp(mypar[4])
  eta=expit(mypar[5])
  
  L1 <- (a1/(1+a1-exp(-1)))^r1
  L2 <- (a2/(1+a2-exp(-1)))^r2
  
  lower.bound <- -1/pmax( L1*L2 , (1-L1)*(1-L2) )
  upper.bound <- 1/pmax( L1*(1-L2) , L2*(1-L1) )
  
  w <- lower.bound + (upper.bound-lower.bound)*eta 
  
  corr <- w*(1-exp(-1))^2*sqrt(r1*r2*(1+a1)*(1+a2)/(a1*a2))*(a1/(1+a1-exp(-1)))^(r1+1)*(a2/(1+a2-exp(-1)))^(r2+1)   
  
  return(c(a1,b1,a2,b2,w, mu1, mu2, lower.bound, upper.bound, sigma1, sigma2, corr ))
}

par.cal(mypar=NBD.sar$par)
