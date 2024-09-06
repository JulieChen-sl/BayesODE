mcmcsample<-function(yobsv,n.par,delta_t,n.track,n,n0,
                             mc1,mc2,mc3,mc4,sd_pro)#n.par是参数个数
{
  chains<-NULL
  varnames <- paste("mc", 1:10, sep="")
  T<-nrow(yobsv)
  
  for(s in 1:3)
  { 
    x.val1<-matrix(data=NA,ncol=n.par,nrow=n/2+n0/2)
    x.val1[1:(n/2-n0/2),]<-get(varnames[s])
    for(k in (n/2-n0/2):(n/2+n0/2-1))
    {  
      alpha.t<-x.val1[k,1]
      alpha.star<-rtruncnorm(1,a=0,b=1,mean=alpha.t,sd=sd_pro[1])
      R<-R1(alpha.t,alpha.star,x.val1[k,2],x.val1[k,3],x.val1[k,4],yobsv[,1],
            yobsv[,2],delta_t,n.track,T,sd_pro[1])
      p <-min(R,1)
      d <-rbinom(1,1,p)
      x.val1[k+1,1]<-alpha.star*d+alpha.t*(1-d)

      beta.t<-x.val1[k,2]
      beta.star<-rtruncnorm(1,a=0,b=1,mean=beta.t,sd=sd_pro[2])
      R<-R2(x.val1[k+1,1],beta.t,beta.star, x.val1[k,3],x.val1[k,4],yobsv[,1],
            yobsv[,2],delta_t,n.track,T,sd_pro[2])
      p <-min(R,1)
      d <-rbinom(1,1,p)
      x.val1[k+1,2]<-beta.star*d+beta.t*(1-d)
      
      #x.val1[k+1,1:2]<-theta[1:2]
      
      lambda1.t<-x.val1[k,3]
      lambda1.star<-rtruncnorm(1,a=0,b=log(2),mean=lambda1.t,sd=sd_pro[3])
      R<-R3(x.val1[k+1,1],x.val1[k+1,2],lambda1.t,lambda1.star,x.val1[k,4],yobsv[,1],
            yobsv[,2],delta_t,n.track,T,sd_pro[3])
      p <-min(R,1)
      d <-rbinom(1,1,p)
      x.val1[k+1,3]<-lambda1.star*d+lambda1.t*(1-d)
      
      #x.val1[k+1,4]<-theta[4]
      
      lambda2.t<-x.val1[k,4]
      lambda2.star<-rtruncnorm(1,a=0,b=log(2),mean=lambda2.t,sd=sd_pro[4])
      R<-R4(x.val1[k+1,1],x.val1[k+1,2],x.val1[k+1,3],lambda2.t,lambda2.star,yobsv[,1],
            yobsv[,2],delta_t,n.track,T,sd_pro[4])
      p <-min(R,1)
      d <-rbinom(1,1,p)
      x.val1[k+1,4]<-lambda2.star*d+lambda2.t*(1-d)
    }
    chains[[s]]<-x.val1[-(1:(n0/2)),]
  }
  return(chains)
}