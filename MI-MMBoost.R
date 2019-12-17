rm(list=ls())

library(splines)
library(RcppArmadillo)
library(Rcpp)
library(mvtnorm)
library(mice)
library(SGL)

source("simul.R")
sourceCpp("ddloglik_missingdata.cpp")


nloop = 50                  #number of iterations                

time_main50 <-rep(0,nloop)    # record the time used to do varaible selection
FPFNSeSpLik1 <- NULL          # record the selection results

for (loop in 1:nloop) {
  
  set.seed(2018 + loop)
  #simulate data
  s   <- simul(N = 100, p =20, p_true = 5)
  #set data
  delta = s$delta
  z     = s$z
  time  = s$time
  p     = s$p
  N     = s$N
  true_beta = s$TrueBeta
  
  #censore rate
  sum(s$delta==1)/1000
  
  #sort time
  delta          <- delta[order(time)]
  z              <- z[order(time),]
  time           <- time[order(time)]
  s$z            <- s$z[order(time),]
  
  #get the Nelson-Aalen estimator
  data = as.data.frame(cbind(z,delta,time)) 
  HT <- nelsonaalen(data, timevar = time, statusvar = delta)
  
  #generate missing data (MCAR)############################################################################################################
  
  alpha0 <- rnorm(N,-0.3)
  Pr_R <- matrix(0,N,p)
  for (i in 1:N) {
    for (j in (p/2+1):p){
      Pr_R[i,j] = alpha0[i] + 0.5*z[i,j-(p/2)]  + 0.5*delta[i] + 0.5*HT[i]
    }
  }
  Pr_R2 <- exp(Pr_R)/(1+exp(Pr_R))
  
  Rij <- matrix(1,N,p)
  for (i in 1:N) {
    for (j in (p/2+1):p){
      Rij[i,j] = sample(c(0,1),size = 1,prob = c(Pr_R2[i,j],1-Pr_R2[i,j]))
    }
  }
  
  #missing rate (change alpha 0 to change missing rate)
  print((N*p-sum(Rij))/(N*p)*100)
  
  for (i in 1:N) {
    for (j in 1:p){
      if(Rij[i,j] == 0){
        z[i,j] = NA
      }
    }
  }
  
  #Nelson-alon
  data = cbind(z,delta,HT)
  
  #MI
  data_mi <- mice(data, method = "pmm")
  mydata <- NULL
  for(i in 1:5) {
    mydata[[i]] <- complete(data_mi,i)
  }
  zz <- NULL
  for(m in 1:5){
    zz[[m]] <- as.matrix(mydata[[m]][,1:p])
  }
  
  
  #Variable seleciton part################################################################################################################################################
  
  time_main <- proc.time()
  
  maxit = 1000
  rate = 0.001
  tol = 1e-6
  #order delta, z and time by time.
  
  
  # set initial values
  key           <- 0                     # the number of iterations 
  likelihood_all<- NULL                  # the likelihood
  j_star_all    <- NULL                  # the chosen optimal block
  converge      <- FALSE                 # record whether it converges
  beta          <- matrix(0,5,p)
  
  #set initial beta to be true
  #beta         <- s$TrueBeta
  
  #main loop
  repeat{
    key       <- key + 1
    
    GD        <- rep(0,p)
    likelihood<- rep(0,m)
    update_all<- matrix(0,5,p)
    
    #ddloglik_md2: Rcpp function. 
    #m respresents 5 cohort.  for each cohort, record the loglikelihood(scalar) and L1(p*1 vector)
    for(m in 1:5){
      result     <- ddloglik_md2(delta,zz[[m]],beta[m,])  #for imputed data set, should be ddloglik_md2(delta,zz[[m]],beta[m,])
      update_all[m,] <- result$L1                   
      likelihood[m] <- sum(result$partial_likelihood)/N
    }
    
    GD          <- apply(update_all^2, 2, sum)
    
    j_star              <- which(GD==max(GD))                # choose the maximum as the optimal block.
    j_star_all          <- c(j_star_all, j_star)
    beta[,j_star] <- beta[,j_star] + rate*sign(update_all[,j_star])  #update the beta
    likelihood_all      <- c(likelihood_all, sum(likelihood))
    
    # check the convergence
    # use relative likelihood ratio change here. 
    if (key>=10) {
      llk.diff = abs((likelihood_all[key]-likelihood_all[key-1])/(likelihood_all[key-1]))
      if(llk.diff < tol) {
        converge <- TRUE      # relative likelihood ratio change is less than the threshold, ,mark the result's convergence true.
        break
      }
    }
    
    # if the loop gets the maximum iteration, we break the loop.
    if (key==maxit){
      break
    }
  }#end loop
  
  
  which(beta[1,]!=0)
  which(true_beta!=0)
  FPFNSeSpLik1[[loop]] <- FPFNSeSpLik(TrueBeta = true_beta,beta = beta[1,] )
  
  time_main50[loop] <- (proc.time() -time_main)[3]
}



result <- c(0,0,0,0,0)

for(i in 1:nloop) {
  result <- result+ FPFNSeSpLik1[[i]]
} 

#result (FP,FN,SEN,SPE,FDP):

print(result/nloop)

#time usde for MI:
mean(time_main50)












