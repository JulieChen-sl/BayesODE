#Spline smoothing

filenames1 <- paste("output", 1:10000, sep="")
filenames13 <- paste("DATA", 1:100000, sep="")
para<-read.csv("para_unif.csv",header=TRUE)

library("splines")
# install.packages("DEoptimR")
library(DEoptimR)
f_optim_w1<-function(theta)
{
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  mu_d<-(lambda2-lambda1)*yobsv[,1]^2+(lambda1*alpha-lambda2*(1+beta))*yobsv[,1]+lambda2*beta
  return(sum((mu_d-muhat)^2))
}
f_optim_w2<-function(theta)
{
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  mu<-theta[5:17]
  N0<-1000
  delta_t<-2
  t<-1:12
  Nt<-c(N0,N0*exp((lambda1-lambda2)*0.5*(mu[-1]+mu[-13])*delta_t+lambda2*delta_t*t))
  sigma_d<-(2*(lambda1*alpha-lambda2*beta)-(lambda2+lambda1))*yobsv[,2]
  +2*(lambda2-lambda1)*mu*yobsv[,2]+mu*(lambda1-lambda2)/(2*Nt)+lambda2/(2*Nt)
  return(sum((sigma_d-sigmahat)^2))
}
f_optim<-function(theta)
{
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  N0<-1000
  delta_t<-2
  t<-1:12
  Nt<-c(N0,N0*exp((lambda1-lambda2)*0.5*(yobsv[-1,1]+yobsv[-13,1])*delta_t+lambda2*delta_t*t))
  
  mu_d<-(lambda2-lambda1)*yobsv[,1]^2+(lambda1*alpha-lambda2*(1+beta))*yobsv[,1]+lambda2*beta
  sigma_d<-(2*(lambda1*alpha-lambda2*beta)-(lambda2+lambda1))*yobsv[,2]
            +2*(lambda2-lambda1)*yobsv[,1]*yobsv[,2]+yobsv[,1]*(lambda1-lambda2)/(2*Nt)+lambda2/(2*Nt)
  return(sum((mu_d-muhat)^2/w1+(sigma_d-sigmahat)^2/w2))
}


for(i in 1:100)
{
  theta<-as.numeric(para[i,2:5])
  yobsv1<-as.matrix(read.table(paste(filenames13[i], ".csv", sep=""),header=FALSE)[,2:3])
  yobsv<-yobsv1[seq(1,37,3),]
  y_initial<-yobsv[1,]

  t<-seq(0,24,2)
  fit_mu<-splinefun(t,yobsv[,1],method="natural")
  muhat<-fit_mu(t,deriv=1)
  fit_sigma<-splinefun(t,yobsv[,2],method="natural")
  sigmahat<-fit_sigma(t,deriv=1)

  ####Calculation of weights####
  w1<-JDEoptim(rep(0,4), c(1,1,log(2),log(2)), f_optim_w1, maxiter = 1000000)[[2]]
  w2<-JDEoptim(rep(0,17), c(1,1,log(2),log(2),rep(1,13)), f_optim_w2, maxiter = 1000000)[[2]]
  w1<-w1/10  #n-p-1£¬n=13,p=2
  w2<-w2/10  #n-p-1£¬n=13,p=2
  ####
   output_optim<-JDEoptim(rep(0,4), c(1,1,log(2),log(2)), f_optim, maxiter = 3000000)
   bias<-(output_optim[[1]]-theta)^2

   write.table(c(i,bias),paste(filenames1[i], ".csv", sep=""),append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
}
