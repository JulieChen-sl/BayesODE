#Write a program to estimate parameters using 
#Differential Evolution Optimization algorithm and RK4

#Write a function that, given the values of 
#the arguments, computes the value of the objective function
f_rk4<-function(t,x,theta)
{
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  N0<-theta[5]
  
  Nt<-N0*exp((lambda1-lambda2)*x[1]*t+lambda2*t)
  dmiu<-(lambda2-lambda1)*x[1]^2+(lambda1*alpha-lambda2*(1+beta))*x[1]+lambda2*beta
  dsigma2<-(2*(lambda1*alpha-lambda2*beta)-(lambda1+lambda2))*x[2]+2*(lambda2-lambda1)*x[1]*x[2]+x[1]*(lambda1-lambda2)/(2*Nt)+lambda2/(2*Nt)
  return(list(c(dmiu,dsigma2)))
}

#The objective function, weighting, 
#mean and variance weights are 1/w1,1/w2 respectively
f_optim<-function(theta)
{
   time<-seq(0,24,2)
   theta1<-c(theta,1000)
   out<-as.data.frame(rk4(y=y_initial,times=time,func=f_rk4,parms=theta1))
   return(sum((out[[2]]-yobsv[,1])^2/w1+(out[[3]]-yobsv[,2])^2/w2))
 }

####In order to calculate the weights, write to calculate the objective function separately
f_rk4_mu<-function(t,x,theta)
{
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  N0<-theta[5]
  
  Nt<-N0*exp((lambda1-lambda2)*x[1]*t+lambda2*t)
  dmiu<-(lambda2-lambda1)*x[1]^2+(lambda1*alpha-lambda2*(1+beta))*x[1]+lambda2*beta
  #dsigma2<-(2*(lambda1*alpha-lambda2*beta)-(lambda1+lambda2))*x[2]+2*(lambda2-lambda1)*x[1]*x[2]+x[1]*(lambda1-lambda2)/(2*Nt)+lambda2/(2*Nt)
  return(list(dmiu))
}
f_rk4_sigma<-function(t,x,theta)
{
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  mu<-theta[5:17]
  N0<-theta[18]
  
  Nt<-N0*exp((lambda1-lambda2)*mu[t/2+1]*t+lambda2*t)
  #dmiu<-(lambda2-lambda1)*x[1]^2+(lambda1*alpha-lambda2*(1+beta))*x[1]+lambda2*beta
  dsigma2<-(2*(lambda1*alpha-lambda2*beta)-(lambda1+lambda2))*x+2*(lambda2-lambda1)*mu[t/2+1]*x+mu[t/2+1]*(lambda1-lambda2)/(2*Nt)+lambda2/(2*Nt)
  return(list(dsigma2))
}
f_optim_w1<-function(theta)
{
  time<-seq(0,24,2)
  theta1<-c(theta,1000)
  out<-as.data.frame(rk4(y=y_initial[1],times=time,func=f_rk4_mu,parms=theta1))
  return(sum((out[[2]]-yobsv[,1])^2))
}
f_optim_w2<-function(theta)
{
  time<-seq(0,24,2)
  theta1<-c(theta,1000)
  out<-as.data.frame(rk4(y=y_initial[2],times=time,func=f_rk4_sigma,parms=theta1))
  return(sum((out[[2]]-yobsv[,2])^2))
}


#install.packages("deSolve")
library(deSolve)
# install.packages("DEoptimR")
library(DEoptimR)

filenames1 <- paste("output", 1:10000, sep="")
filenames13 <- paste("DATA", 1:100000, sep="")
para<-read.csv("para_unif.csv",header=TRUE)
for (i in jobB:jobE)
{
  theta<-as.numeric(para[i,2:5])
  yobsv1<-as.matrix(read.table(paste(filenames13[i], ".csv", sep=""),header=FALSE)[,2:3])
  yobsv<-yobsv1[seq(1,37,3),]
  y_initial<-yobsv[1,]
  ####Calculation of weights####
  w1<-JDEoptim(rep(0,4), c(1,1,log(2),log(2)), f_optim_w1, maxiter = 3000000)[[2]]
  w2<-JDEoptim(rep(0,17), c(1,1,log(2),log(2),rep(1,13)), f_optim_w2, maxiter = 3000000)[[2]]
  w1<-w1/10  #n-p-1,n=13,p=2
  w2<-w2/10  #n-p-1,n=13,p=2
  ####
  output_optim<-JDEoptim(rep(0,4), c(1,1,log(2),log(2)), f_optim, maxiter = 3000000)
  bias<-(output_optim[[1]]-theta)^2
  
  write.table(c(i,bias,time2),paste(filenames1[i], ".csv", sep=""),append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
}
