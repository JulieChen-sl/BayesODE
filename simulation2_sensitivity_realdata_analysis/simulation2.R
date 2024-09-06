#for server
args=commandArgs(T); 
jobB=args[1];
jobE=args[2];
DirPre=args[3];

library(coda)
library(MASS)
library(truncnorm)

source("logL_missing.R")
source("add_data_new.R")
source("reset_yobsv.R")
source("Gibbs_missing_not0beta_NotEqualLambda.R")
source("check_convergence_missing_not0beta_NotEqualLambda.R")
source("accept_rate_missing.R")


delta_t<-2/3;n.data<-36
n.par<-4+2*(36-12) #the number of parameters
n.track=5;
N0<-1000#the initial total number of cell,change this parameter when do sensitivity analysis
filenames1 <- paste("output", 1:1000, sep="")
para<-read.csv("para_unif.csv",header=TRUE)

filenames13 <- paste("DATA", 1:100000, sep="")
for(i in jobB:jobE)
{
  theta<-as.numeric(para[i,2:5])
  m_0<-para[i,6]
  
  yobsv_real<-as.matrix(read.table(paste(filenames13[i], ".csv", sep=""),header=FALSE)[,2:3])
  f_real<-f_lambda(theta,yobsv_real,delta_t,n.track)
    
  yobsv_missing<-matrix(NA,nrow=13,ncol=2)
  for(a in 1:13)
  {
      yobsv_missing[a,]<-yobsv_real[3*a-2,]
  }#Getting simulated missing data

  output_real1<-check_convergence(yobsv_missing,n.par=4+2*(n.data-12),delta_t,n.track)
  write.table(c(i,m_0,theta,output_real1,f_real),paste(filenames1[i], ".csv", sep=""),
              append=TRUE,col.names=FALSE, sep = " ",eol = "\n")

}