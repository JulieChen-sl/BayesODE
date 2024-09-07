args=commandArgs(T); 
jobB=args[1];
jobE=args[2];
DirPre=args[3];

library(truncnorm)
library(coda)
##########logL,MH############
f=function(theta,yobsv,n.track)
{
  m<-yobsv[,1]
  v<-yobsv[,2]
  T<-nrow(yobsv)
  n<-n.track
  alpha<-theta[1]
  beta<-theta[2]
  f=0
  aux1<-(1+alpha-beta)/2
  
  for(t in 2:(T-1))
  {
    f=f-n*(2*log(aux1)+log(v[t]))+(n-3)*log(v[t+1])-(n*(m[t+1]-aux1*m[t]-beta/2)^2+(n-1)*v[t+1])/(aux1^2*v[t])
  }
  
  return(0.5*f)
}

#MH RATE
R1 = function(alpha.t,alpha.star,theta.t,theta.star,yobsv,n.track)
{exp(f(theta.star,yobsv,n.track)-f(theta.t,yobsv,n.track))*dtruncnorm(alpha.t,a=0,b=1,mean=alpha.star,sd=0.001)/dtruncnorm(alpha.star,a=0,b=1,mean=alpha.t,sd=0.001)}

####Gibbs####
mcmc.sample<-function(yobsv,n.track,n,n0,mc1,mc2,mc3,mc4)
{
  chains<-NULL
  varnames <- paste("mc", 1:10, sep="")
  filenames10 <- paste("test", 1:100000, sep="")
  for(s in 1:4)
  { 
    x.val1<-matrix(data=NA,ncol=2,nrow=n/2+n0/2)
    x.val1[1:(n/2-n0/2),]<-get(varnames[s])
    for(k in (n/2-n0/2):(n/2+n0/2-1))
    {  
      para.t<-x.val1[k,1]
      theta.t<-c(para.t,x.val1[k,2])
      para.star<-rtruncnorm(1,a=0,b=1,mean=para.t,sd=0.001)
      theta.star<-c(para.star,x.val1[k,2])
      R<-R1(para.t,para.star,theta.t,theta.star,yobsv,n.track)
      p <-min(R,1)
      d <-rbinom(1,1,p) 
      x.val1[k+1,1]<-para.star*d+para.t*(1-d)
      
      para.t<-x.val1[k,2]
      theta.t<-c(x.val1[k+1,1],para.t)
      para.star<-rtruncnorm(1,a=0,b=1,mean=para.t,sd=0.001)
      theta.star<-c(x.val1[k+1,1],para.star)
      R<-R1(para.t,para.star,theta.t,theta.star,yobsv,n.track)
      p <-min(R,1)
      d <-rbinom(1,1,p)
      x.val1[k+1,2]<-para.star*d+para.t*(1-d)
    }
    chains[[s]]<-x.val1[-(1:(n0/2)),]
  }
  return(chains)
}

check_convergence<-function(yobsv,f_real,n.track,n.data)
{
  varnames <- paste("mc", 1:10, sep="")
  filenames6 <- paste("not_0betaMPSRF_lambda", 1:10000, sep="")
  filenames9 <- paste("Test", 1:10000, sep="")
  n0<-10000
  n<-n0+2
  
  chains<-NULL
  assign(varnames[1],value=runif(2,0,1))
  assign(varnames[2],value=runif(2,0,1))
  assign(varnames[3],value=runif(2,0,1))
  assign(varnames[4],value=runif(2,0,1))
  
  flag<-0
  repeat
  {
    x.val<-mcmc.sample(yobsv,n.track,n,n0,mc1,mc2,mc3,mc4)
    
    chains<-NULL
    chains[[1]]<-as.mcmc(x.val[[1]])
    chains[[2]]<-as.mcmc(x.val[[2]])
    chains[[3]]<-as.mcmc(x.val[[3]])
    chains[[4]]<-as.mcmc(x.val[[4]])
    MPSRF<-gelman.diag(chains,autoburnin=FALSE,multivariate = TRUE)[[2]]
    write.table(MPSRF,paste(filenames6[i], ".csv", sep=""),append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
    rm(chains)
    gc()
    max_iteration<-1e+07
    if(MPSRF<1.1||(n>max_iteration))
    {
        steady<-matrix(data=NA,nrow=4*n/2,ncol=2)
        steady[1:(n/2),]<-x.val[[1]]
        steady[(n/2+1):n,]<-x.val[[2]]
        steady[(n+1):(3*n/2),]<-x.val[[3]]
        steady[(3*n/2+1):(2*n),]<-x.val[[4]]
        
        num.par<-rep(0,2)
        if(theta[1]<=quantile(steady[,1],probs =c(0.025,0.975))[2]&&
           theta[1]>=quantile(steady[,1],probs =c(0.025,0.975))[1])
          num.par[1]=num.par[1]+1
        if(theta[2]<=quantile(steady[,2],probs =c(0.025,0.975))[2]&&
           theta[2]>=quantile(steady[,2],probs =c(0.025,0.975))[1])
          num.par[2]=num.par[2]+1
        
    
          confidence.int<-c(quantile(steady[,1],probs =c(0.025,0.975))[1],quantile(steady[,1],probs =c(0.025,0.975))[2],
                            quantile(steady[,2],probs =c(0.025,0.975))[1],quantile(steady[,2],probs =c(0.025,0.975))[2])
          
          x.mean<-apply(steady,2,mean)
          f_esti_ex<-f(x.mean,yobsv,n.track)     
          
          
          int.len<-c(quantile(steady[,1],probs =c(0.025,0.975))[2]-quantile(steady[,1],probs =c(0.025,0.975))[1],
                     quantile(steady[,2],probs =c(0.025,0.975))[2]-quantile(steady[,2],probs =c(0.025,0.975))[1])
          
          
          bias_ex<-c((x.mean[1]-theta[1])^2,(x.mean[2]-theta[2])^2)
          
          break
        }else
        {
          for(flag2 in 1:4)
          {
            assign(varnames[flag2],value=x.val[[flag2]])
          }
          n<-n+n0
          rm(x.val)
          gc()
        } 
  }
  
  return(c(theta[1:2],x.mean,num.par,bias_ex,int.len,confidence.int,f_esti_ex,f_real))
  
}

############main program##############
N0<-1000
n.track=5;
n.data=13;
para<-read.csv("para_unif.csv",header=TRUE)
filenames1 <- paste("output", 1:10000, sep="")
filenames13 <- paste("DATA", 1:100000, sep="")

for(i in jobB:jobE)
{
  theta<-as.numeric(para[i,2:5])
  m_0<-para[i,6]
  
  yobsv<-as.matrix(read.table(paste(filenames13[i], ".csv", sep=""),header=FALSE)[,2:3])
  f_real<-f(theta,yobsv,n.track)
  
  k1<-1;k2<-2;
  yobsv_missing1<-matrix(NA,nrow=13,ncol=2)
  m_miss<-rep(0,24); v_miss<-rep(0,24);
  for(a in 1:12)
  {
    yobsv_missing1[a,]<-yobsv[(k1+k2)*a-k2,]
    m_miss[(2*a-1):(2*a)]<-yobsv[((k1+k2)*a-k2+1):((k1+k2)*a-k2+2),1]
    v_miss[(2*a-1):(2*a)]<-yobsv[((k1+k2)*a-k2+1):((k1+k2)*a-k2+2),2]
  }
  yobsv_missing1[13,]<-yobsv[(k1+k2)*13-k2,]
  
  output_real<-check_convergence(yobsv_missing1,f_real,n.track,n.data)
  write.table(c(i,m_0,output_real),paste(filenames1[i], ".csv", sep=""),append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
}