#for server
args=commandArgs(T); 
jobB=args[1];
jobE=args[2];
DirPre=args[3];

library(coda)#for MCMC
library(MASS)
library(truncnorm)

source("logL_nomissing.R")
source("Gibbs.R")
source("check_convergence.R")
source("accept_rate_com.R")

delta_t<-2/3
n.par<-4 
n.track=5
N0<-1000
filenames1 <- paste("output", 1:10000, sep="")
filenames2 <- paste("DATA", 1:100000, sep="")

para<-read.csv("para_unif.csv",header=TRUE)

#for(i in jobB:jobE)
for(i in 1:100)  
{
  theta <- as.numeric(para[i,2:5])
  m0 <- para[i,6]
  
  yobsv <- as.matrix(read.table(paste(filenames2[i], ".csv", sep=""),
                              header=FALSE)[,2:3])
  
  f_real <- f_lambda(theta[1],theta[2],theta[3],theta[4],yobsv[,1],
                     yobsv[,2],delta_t,n.track,nrow(yobsv))
  
  output_real <- check_convergence(yobsv,n.par,delta_t,n.track)
  
  write.table(c(i,m0,theta,output_real,f_real),
              paste(filenames1[i], ".csv", sep=""),
              append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
  
}