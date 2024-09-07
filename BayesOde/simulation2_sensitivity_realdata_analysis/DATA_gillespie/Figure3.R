##Picture of the total number of cells####
library(actuar)
filenames1<- paste("NT", 1:100, sep="")

############################################################################
for(i in 1:4)
{
  if(i%%4==1)
  {
   pdf(paste(filenames1[ii%/%4+1], ".pdf", sep=""))
   par(mfrow=c(2,2))
  }

  theta<-c(runif(2,0,1),runif(2,0,0.5*log(2)))
  n.data<-36
  t=24
  delta_t<-t/n.data
  time<-seq(from=0,to=t,length.out=n.data+1)
  N0<-1000

  n<-n.data+1
  n.track=5
  m<-matrix(data=0,ncol=n,nrow=n.track);
  n.total<-matrix(data=0,ncol=n,nrow=n.track);
  v1<-c(1,0,1,0);
  v2<-c(0,1,0,1);
  s<-t/(n-1);
  x<-rep(0,n)
  y<-rep(0,n)
  m1<-runif(1,0,1)
  x[1]=round(N0*m1)
  y[1]=N0-x[1]
  
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  
  for(i in 1:n.track)
  {
    x2<-x[1]
    y2<-y[1]
    t1<-0;
    t2<-0;
    for(k in 1:(n-1))
    {
      while (t2<s*k) 
      {  
        t1<-t2
        x1<-x2
        y1<-y2
        
        a<-c(lambda1*alpha*x1,lambda1*(1-alpha)*x1,lambda2*beta*y1,lambda2*(1-beta)*y1)
        a0<-sum(a)
        r1<-runif(1,0,1)
        r2<-runif(1,0,1)
        tao<-(1/a0)*log(1/r1)
        r<-rank(c(0,a[1],sum(a[1:2]),sum(a[1:3]),sum(a[1:4]),r2*a0))
        miu<-r[6]-1;
        x2<-x1+v1[miu]
        y2<-y1+v2[miu]
        t2<-t1+tao
      }
      x[k+1]<-x1;
      y[k+1]<-y1;
      k=k+1;
    }
    m[i,]<-x/(x+y)
    n.total[i,]<-x+y
  }
  
  Nt_breed<-apply(n.total,2,mean)
  mt<-apply(m,2,mean)
  vt<-apply(m,2,var)
 
  Nt_var<-rep(0,n);
  Nt_mean2<-rep(0,n)
  Nt_var[1]<-N0
  Nt_var[2:n]<-((2*lambda1*alpha-2*lambda2*beta+lambda2-lambda1)
                *mt[-1]-lambda2+2*lambda2*beta)/(2*(lambda2-lambda1)*vt[-1])
  
  for(k in 1:n)
  {
    re<-(k-1)%%3
    if(k<4)
    {
      miu_est<-(1/2)*delta_t*(re*(mt[k-re]+mt[k]))
    }else{
      miu_est<-(1/2)*delta_t*(3*sum(mt[seq(1,k-3-re,by=3)]+mt[seq(4,k-re,by=3)])+re*(mt[k-re]+mt[k]))
    }
    Nt_mean2[k]<-N0*exp((lambda1-lambda2)*miu_est+lambda2*delta_t*(k-1))
  }
 
  lab=paste("m0=",round (m1,3),"¦Á=",round (alpha,3),"¦Â=", round (beta,3),"¦Ë1=", round (lambda1,3),"¦Ë2=", round (lambda2,3))
  mmm<-max(Nt_var,Nt_breed,Nt_mean2)
  par(mar = c(4, 4, 1, 1))
  plot(time,n.total[1,],col=3,type="l",mgp=c(2.5,1,0),xaxt='n',ylim=c(min(Nt_var),mmm),xlab="Day",ylab="Nt",
       lty=5,lwd=1,main = lab,cex.main=0.8)
  axis(1,at=seq(0,24,4),cex.axis=0.8)
  for(k in 2:n.track)
  {
    lines(time,n.total[k,],col=3,lwd=1,lty=5,type="l")
  }
  lines(time,Nt_mean2,col=2,lty=1,lwd=1,pch=18,type="o",cex=0.6)
  lines(time,Nt_var,col=4,lty=1,lwd=1,pch=15,type="o",cex=0.4)
  lines(time,Nt_breed,col=1,lty=1,lwd=1,pch=20,type="o",cex=0.6)
  legend("topleft", inset=.05,legend=c("trajectory","Nt_mean","Nt_var","Nt_average"),col=c(3,2,4,1),lty=c(5,1,1,1)
         ,pch=c(NA,18,15,20),lwd=c(1,1,1,1))
  if(i%%4==0)
  {
    dev.off()
  }

}