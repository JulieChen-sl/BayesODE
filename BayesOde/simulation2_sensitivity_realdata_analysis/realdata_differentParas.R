library(coda)
library(MASS)
library(truncnorm)

source("logL_missing.R")
source("add_data_new.R")
source("reset_yobsv.R")
source("accept_rate_missing.R")
source("IC.R")
source("maxlikeli.R")

###beta!=0,lambda1!=lambda2
source("Gibbs_missing_not0beta_NotEqualLambda.R")
source("check_convergence_missing_not0beta_NotEqualLambda.R")

###beta=0,lambda1!=lambda2
source("Gibbs_missing_0beta_NotEqualLambda.R")
source("check_convergence_missing_0beta_NotEqualLambda.R")

###beta!=0,lambda1=lambda2
source("Gibbs_missing_not0beta_EqualLambda.R")
source("check_convergence_missing_not0beta_EqualLambda.R")

###beta=0,lambda1=lambda2
source("Gibbs_missing_0beta_EqualLambda.R")
source("check_convergence_missing_0beta_EqualLambda.R")



delta_t<-2/3;n.data<-36
n.par<-4+2*(36-12) #the number of parameters
n.track=5;
N0<-1000#the initial total number of cell,change this parameter when do sensitivity analysis

filenames1 <- paste("not0beta_NotEqualLambda", 1:10, sep="")
filenames2 <- paste("0beta_NotEqualLambda", 1:10, sep="")
filenames3 <- paste("not0beta_EqualLambda", 1:10, sep="")
filenames4 <- paste("0beta_EqualLambda", 1:10, sep="")

data<-read.csv("realdata.csv",header = TRUE)
yobsv_missing<-matrix(NA,nrow=13,ncol=8)
yobsv_missing[,c(1,3,5,7)]<-as.matrix(data[,c(1,3,5,7)])
yobsv_missing[,c(2,4,6,8)]<-as.matrix(data[,c(2,4,6,8)])^2

theta<-rep(0,4)
ss<-seq(1,37,by=3)
re<-rep(c(0,1,2),13)
de<-rep(0:13,each=3)

for(i in 1:4)
{
  yobsv_real<-yobsv_missing[,(2*i-1):(2*i)]
  
  output1<-check_convergence_not0beta_NotEqualLambda(yobsv_real,n.par=4+2*(n.data-12),delta_t,n.track)
  write.table(c(i,output1),paste(filenames1[i], ".csv", sep=""),
              append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
  
  output2<-check_convergence_0beta_NotEqualLambda(yobsv_real,n.par=4+2*(n.data-12),delta_t,n.track)
  write.table(c(i,output2),paste(filenames2[i], ".csv", sep=""),
              append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
  
  output3<-check_convergence_not0beta_EqualLambda(yobsv_real,n.par=4+2*(n.data-12),delta_t,n.track)
  write.table(c(i,output3),paste(filenames3[i], ".csv", sep=""),
              append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
  
  output4<-check_convergence_0beta_EqualLambda(yobsv_real,n.par=4+2*(n.data-12),delta_t,n.track)
  write.table(c(i,output4),paste(filenames4[i], ".csv", sep=""),
              append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
  
}