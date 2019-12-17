AR1 <- function(tau, m) {
  if(m==1) {R <- 1}
  if(m > 1) {
    R <- diag(1, m)
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        R[i,j] <- R[j,i] <- tau^(abs(i-j))
      }
    }
  }
  return(R)
}
#' function for standardize covariate
normalize = function(x){
  y = sqrt(length(x))*(x-mean(x))/sqrt(sum((x-mean(x))^2))
  return(y)
}

#' The evaluation function 
#' @param TrueBeta  the true beta
#' @param beta the result beta from the algorithm
#' @return a list of the following values
#' *FP:  false positive 
#' *FN:  false negative
#' *Se:  Senstitivity
#' *Sp:  Specificity
#' *FDP: False positive portion
#' 
#' @export
FPFNSeSpLik=function(TrueBeta=TrueBeta,beta=beta){
  FP <- length(which(TrueBeta==0 & beta!=0))
  FN <- length(which(TrueBeta!=0 & beta==0))
  Se <- length(which(TrueBeta!=0 & beta!=0))/length(which(TrueBeta!=0))
  Sp <- length(which(TrueBeta==0 & beta==0))/length(which(TrueBeta==0))
  FDP=FP/max(length(which(beta!=0)),1)
  output=c(FP, FN, Se, Sp, FDP)
  return(output)
}


#' Simulation function
#' This function is for simulating data
#' @param N sample size
#' @param p number of parameters
#' @param p_true number of true parameters
#' @param m1 default is 100
#' @import mvtnorm
#' @return a list of the following values
#' 
#' *delta: event indicator
#' *z: Covariate matrix
#' *time: the death time
#' *p: number of variables
#' *N: the sample size
#' 
#' 
#' @export
simul  <- function(N = 1000, p = 500, p_true = 5){
  m1          <- p
  Sigma_z1    <- diag(p)
  Corr1       <- AR1(0.3,m1)
  diag(Corr1) <- 1
  z           <- NULL
  j           <- 0
  
  #Simulate z
  while(j<(p/m1)){
    j    <- j+1
    z    <- cbind(z,rmvnorm(N, mean=rep(0,m1), sigma=Corr1))
  }
  z      <- apply(z,2,normalize)
  #z2     <- z**2
  
  #Simulate True beta
  TrueBeta       <- rep(0, p)
  TrueBeta_index <- sample(1:p,p_true,replace=FALSE)
  signbeta       <- sample(c(-1,1),p_true,replace=T)
  #mag            <- runif(p_true, 1,2)
  mag            <- runif(p_true, 0.5,1)
  TrueBeta[TrueBeta_index]  <- mag*signbeta
  xbeta          <- z%*%TrueBeta
  U              <- runif(N, 0, 1)
  #Simulate the time and death indicator
  pre_time       <- -log(U)/(1*exp(xbeta))
  pre_censoring  <- runif(N,1,30)
  pre_censoring  <- pre_censoring*(pre_censoring<3)+3*(pre_censoring>=3)
  tcens          <- (pre_censoring<pre_time) # censoring indicator
  delta          <- 1-tcens
  time           <- pre_time*(delta==1)+pre_censoring*(delta==0)
  #order delta, z and time by time.
  #delta          <- delta[order(time)]
  #z              <- z[order(time),]
  #time           <- time[order(time)]
  
  return( list(delta=delta,  z =z, time = time, p=p, N = N, TrueBeta = TrueBeta) )
  
}






