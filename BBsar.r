library(tidyverse)

setwd("~/R/Sarmanov")

DH_data <- read.csv("DH_baconwitheggs.csv")
names(DH_data)

expit <- function(x) {
  exp(x)/(1+exp(x))
}


##################################################################################################
###Beta-Binomial Sarmanov distribution log-likelihood function  
### Input: 1) mypar: log(a1),log(b1),log(a2),log(b2), eta, 
###        2) x1,k1,x2,k2,y
#########################################################################

BB.sar.ll <- function(mypar, x1,k1,x2,k2,y) { 
  
  a1 <- exp(mypar[1]) 
  b1 <- exp(mypar[2])
  a2 <- exp(mypar[3]) 
  b2 <- exp(mypar[4])
  eta=expit(mypar[5])

  mu1 <- a1/(a1+b1) 
  mu2 <- a2/(a2+b2)  
  cc <- (a1+b1+k1)/k1*(a2+b2+k2)/k2
  lower.bound <- cc * pmax( -1/(mu1*mu2) , -1/((1-mu1)*(1-mu2)) )
  upper.bound <- cc * pmin( 1/(mu1*(1-mu2)) , 1/((1-mu1)*mu2) )

  w <- lower.bound + (upper.bound-lower.bound)*eta 
  
  term1 <- lchoose(k1,x1)+lbeta(a1+x1,b1+k1-x1)-lbeta(a1,b1)
  term2 <- lchoose(k2,x2)+lbeta(a2+x2,b2+k2-x2)-lbeta(a2,b2)
  term3 <- log( 1+w*(x1-k1*mu1)/(a1+b1+k1)*(x2-k2*mu2)/(a2+b2+k2) )
  
  return ( - y %*% (term1+term2+term3) )
  
}

BB.sar <- optim(c(0,0,0,0,5), BB.sar.ll, method = "L-BFGS-B",
                lower=rep(-10,5), upper=rep(10,5),
                control = list(maxit=1000),
                hessian=T, 
                x1=DH_data$x1,k1=DH_data$n1,x2=DH_data$x2,k2=DH_data$n2,y=DH_data$y)
BB.sar

BB.ll <- function(mypar, x1,k1,x2,k2,y) { 

  a1 <- exp(mypar[1]) 
  b1 <- exp(mypar[2])
  a2 <- exp(mypar[3]) 
  b2 <- exp(mypar[4])

  term1 <- lchoose(k1,x1)+lbeta(a1+x1,b1+k1-x1)-lbeta(a1,b1)
  term2 <- lchoose(k2,x2)+lbeta(a2+x2,b2+k2-x2)-lbeta(a2,b2)
  
  return ( - y %*% (term1+term2) )
  
}


###################################################################################################


BB <- optim(c(0,0,0,0), BB.ll, method = "L-BFGS-B",
                lower=rep(-10,4), upper=rep(10,4),
                control = list(maxit=1000),
                hessian=T, 
                x1=DH_data$x1,k1=DH_data$n1,x2=DH_data$x2,k2=DH_data$n2,y=DH_data$y)

BB

##################################################################################################
### Purpose:  Transfer the paremeter estimates in the transformed scales into the original scale
### input: parameters in the transformed scales: log(a1,b1,a2,b2), eta
###################################################################################################
par.cal <- function(mypar,k1,k2) {
  
  a1 <- exp(mypar[1]) 
  b1 <- exp(mypar[2])
  a2 <- exp(mypar[3])
  b2 <- exp(mypar[4])
  eta <- expit(mypar[5])

  mu1 <- a1/(a1+b1) 
  mu2 <- a2/(a2+b2)  
  cc <- (a1+b1+k1)/k1*(a2+b2+k2)/k2
  lower.bound <- cc * pmax( -1/(mu1*mu2) , -1/((1-mu1)*(1-mu2)) )
  upper.bound <- cc * pmin( 1/(mu1*(1-mu2)) , 1/((1-mu1)*mu2) )
  
  w <- (upper.bound-lower.bound)*eta + lower.bound
  
  sigma1 <- sqrt( mu1*(1-mu1)/ (a1+b1+1))
  sigma2 <- sqrt( mu2*(1-mu2)/ (a2+b2+1))
  corr <- w*sigma1*sigma2  
  covar <- corr*k1*k2
  
  return(c(a1,b1,a2,b2,w, mu1, mu2, lower.bound, upper.bound, sigma1, sigma2, corr,covar ))
}

par.cal(mypar=BB.sar$par,k1=max(DH_data$n1),k2=max(DH_data$n2))
