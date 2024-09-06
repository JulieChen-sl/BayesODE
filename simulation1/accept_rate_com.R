proposal_sd<-function(yobsv)
{
  set.seed(i)
  theta_initial<-c(runif(2,0,1),runif(2,0,log(2)))
  assign("mc",value=theta_initial)
  sd_pro<-rep(0.02,4)
  T<-nrow(yobsv)
  n0<-1000
  n<-n0+2
  flag<-1
 
  while(flag==1)
  {
    d_accept<-rep(0,4)
    x.val1<-matrix(data=NA,ncol=n.par,nrow=n/2+n0/2)
    x.val1[1:(n/2-n0/2),]<-get("mc")
    
    for(k in (n/2-n0/2):(n/2+n0/2-1))
    {  
      alpha.t<-x.val1[k,1]
      alpha.star<-rtruncnorm(1,a=0,b=1,mean=alpha.t,sd=sd_pro[1])
      R<-R1(alpha.t,alpha.star,x.val1[k,2],x.val1[k,3],x.val1[k,4],yobsv[,1],
            yobsv[,2],delta_t,n.track,T,sd_pro[1])
      p <-min(R,1)
      d <-rbinom(1,1,p) 
      x.val1[k+1,1]<-alpha.star*d+alpha.t*(1-d)
      d_accept[1]<-d_accept[1]+d
      
      beta.t<-x.val1[k,2]
      beta.star<-rtruncnorm(1,a=0,b=1,mean=beta.t,sd=sd_pro[2])
      R<-R2(x.val1[k+1,1],beta.t,beta.star, x.val1[k,3],x.val1[k,4],yobsv[,1],
            yobsv[,2],delta_t,n.track,T,sd_pro[2])
      p <-min(R,1)
      d <-rbinom(1,1,p)
      x.val1[k+1,2]<-beta.star*d+beta.t*(1-d)
      d_accept[2]<-d_accept[2]+d
      
      lambda1.t<-x.val1[k,3]
      lambda1.star<-rtruncnorm(1,a=0,b=log(2),mean=lambda1.t,sd=sd_pro[3])
      R<-R3(x.val1[k+1,1],x.val1[k+1,2],lambda1.t,lambda1.star,x.val1[k,4],yobsv[,1],
            yobsv[,2],delta_t,n.track,T,sd_pro[3])
      p <-min(R,1)
      d <-rbinom(1,1,p)
      x.val1[k+1,3]<-lambda1.star*d+lambda1.t*(1-d)
      d_accept[3]<-d_accept[3]+d
      
      lambda2.t<-x.val1[k,4]
      lambda2.star<-rtruncnorm(1,a=0,b=log(2),mean=lambda2.t,sd=sd_pro[4])
      R<-R4(x.val1[k+1,1],x.val1[k+1,2],x.val1[k+1,3],lambda2.t,lambda2.star,yobsv[,1],
            yobsv[,2],delta_t,n.track,T,sd_pro[4])
      p <-min(R,1)
      d <-rbinom(1,1,p)
      x.val1[k+1,4]<-lambda2.star*d+lambda2.t*(1-d)
      d_accept[4]<-d_accept[4]+d
      
    }
    
    flag<-0
    for(m in 1:n.par)
    {
      if(((d_accept[m]/n0)<0.15)&&(sd_pro[m]>1e-3))
      {
        sd_pro[m]<-(2/3)*sd_pro[m]
        flag<-1
      }else if(((d_accept[m]/n0)>0.35)&&(sd_pro[m]<0.5))
      {
        sd_pro[m]<-1.5*sd_pro[m]
        flag<-1
      }
    }
    assign("mc",value=x.val1[-(1:(n0/2)),])
    n<-n0+n
  }
  return(sd_pro)
}
