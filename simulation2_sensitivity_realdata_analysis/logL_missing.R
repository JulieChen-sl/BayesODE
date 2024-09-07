#f_lambda####
f_lambda=function(theta,yobsv,delta_t,n.track)
{
  m<-yobsv[,1]
  v<-yobsv[,2]
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  T<-length(m)
  f=0
  
  aux1<-delta_t*(2*(lambda1*alpha-lambda2*beta)-(lambda1+lambda2))+1
  aux2<-delta_t*(lambda1*alpha-lambda2*(1+beta))+1
  aux3<-(lambda2-lambda1)*delta_t
  aux5<-lambda2*beta*delta_t
  
  ss<-seq(1,37,by=3)
  re<-rep(c(0,1,2),13)
  de<-rep(0:13,each=3)
  
  for(t in 1:(T-1))
  {
    if(t<4)
    {
      miu_est<-0.5*delta_t*(re[t]*(m[t-re[t]]+m[t]))
      miu_est2<-0.5*delta_t*((re[t]+1)*(m[t-re[t]]+m[t+1]))
    }else{
      s1<-ss[1:de[t]]
      s2<-ss[2:(de[t]+1)]
      middle<-3*sum(m[s1]+m[s2])
      miu_est<-0.5*delta_t*(middle+re[t]*(m[t-re[t]]+m[t]))
      miu_est2<-0.5*delta_t*(middle+(re[t]+1)*(m[t-re[t]]+m[t+1]))
    }
    
    Nt<-N0*exp((lambda1-lambda2)*miu_est+lambda2*delta_t*(t-1))
    Nt2<-N0*exp((lambda1-lambda2)*miu_est2+lambda2*delta_t*t)
    
    sigma<-aux1*v[t]+2*aux3*m[t]*v[t]-aux3*m[t]/(2*Nt)+lambda2*delta_t/(2*Nt)
    miu1<-aux3*m[t]^2+aux2*m[t]+aux5
    
    sigma_star<-v[t]+(aux1-1)*sigma+2*aux3*m[t+1]*sigma-aux3*m[t+1]/(2*Nt2)+lambda2*delta_t/(2*Nt2)
    miu2<-m[t]+aux3*miu1^2+(aux2-1)*miu1+aux5
    
    sigma2<-0.5*(sigma+sigma_star)
    miu<-0.5*(miu1+miu2)
    
    if(sigma2<0)
    {
      sigma2<-1e-300
      write.table("sigma2<0",paste(filenames12[i], ".csv", sep=""),
                  append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
    }
    if(v[t+1]==0)
    {
      v[t+1]<-v[t+2]
    }
    f=f-n.track*log(sigma2)+(n.track-3)*log(v[t+1])-(n.track*(m[t+1]-miu)^2+(n.track-1)*v[t+1])/sigma2
  }
  return(0.5*f)
}

f_lambda1=function(theta,yobsv,delta_t,n.track)
{
  m<-yobsv[,1]
  v<-yobsv[,2]
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  T<-length(m)
  f=0
  
  aux1<-delta_t*(2*(lambda1*alpha-lambda2*beta)-(lambda1+lambda2))+1
  aux2<-delta_t*(lambda1*alpha-lambda2*(1+beta))+1
  aux3<-(lambda2-lambda1)*delta_t
  aux5<-lambda2*beta*delta_t
  
  ss<-seq(1,37,by=3)
  re<-rep(c(0,1,2),13)
  de<-rep(0:13,each=3)
  
  for(t in 1:(T-1))
  {
    if(t<4)
    {
      miu_est<-0.5*delta_t*(re[t]*(m[t-re[t]]+m[t]))
      miu_est2<-0.5*delta_t*((re[t]+1)*(m[t-re[t]]+m[t+1]))
    }else{
      s1<-ss[1:de[t]]
      s2<-ss[2:(de[t]+1)]
      middle<-3*sum(m[s1]+m[s2])
      miu_est<-0.5*delta_t*(middle+re[t]*(m[t-re[t]]+m[t]))
      miu_est2<-0.5*delta_t*(middle+(re[t]+1)*(m[t-re[t]]+m[t+1]))
    }
    
    Nt<-N0*exp((lambda1-lambda2)*miu_est+lambda2*delta_t*(t-1))
    Nt2<-N0*exp((lambda1-lambda2)*miu_est2+lambda2*delta_t*t)
    
    sigma<-aux1*v[t]+2*aux3*m[t]*v[t]-aux3*m[t]/(2*Nt)+lambda2*delta_t/(2*Nt)
    miu1<-aux3*m[t]^2+aux2*m[t]+aux5
    
    sigma_star<-v[t]+(aux1-1)*sigma+2*aux3*m[t+1]*sigma-aux3*m[t+1]/(2*Nt2)+lambda2*delta_t/(2*Nt2)
    miu2<-m[t]+aux3*miu1^2+(aux2-1)*miu1+aux5
    
    sigma2<-0.5*(sigma+sigma_star)
    miu<-0.5*(miu1+miu2)
    
    if(sigma2<0)
    {
      sigma2<-1e-300
      write.table("sigma2<0",paste(filenames12[i], ".csv", sep=""),
                  append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
    }
    if(v[t+1]==0)
    {
      v[t+1]<-v[t+2]
    }
    f=f-n.track*log(sigma2)-(n.track*(m[t+1]-miu)^2+(n.track-1)*v[t+1])/sigma2
  }
  return(0.5*f)
}

f_lambda_m=function(theta,yobsv,delta_t,n.track,k)
{
  m<-yobsv[,1]
  v<-yobsv[,2]
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  T<-length(m)
  f=0
  
  aux1<-delta_t*(2*(lambda1*alpha-lambda2*beta)-(lambda1+lambda2))+1
  aux2<-delta_t*(lambda1*alpha-lambda2*(1+beta))+1
  aux3<-(lambda2-lambda1)*delta_t
  aux5<-lambda2*beta*delta_t
  
  ss<-seq(1,37,by=3)
  re<-rep(c(0,1,2),13)
  de<-rep(0:13,each=3)
  
  for(t in  (k-1):k)
  {
    if(t<4)
    {
      miu_est<-0.5*delta_t*(re[t]*(m[t-re[t]]+m[t]))
      miu_est2<-0.5*delta_t*((re[t]+1)*(m[t-re[t]]+m[t+1]))
    }else{
      s1<-ss[1:de[t]]
      s2<-ss[2:(de[t]+1)]
      middle<-3*sum(m[s1]+m[s2])
      miu_est<-0.5*delta_t*(middle+re[t]*(m[t-re[t]]+m[t]))
      miu_est2<-0.5*delta_t*(middle+(re[t]+1)*(m[t-re[t]]+m[t+1]))
    }
    
    Nt<-N0*exp((lambda1-lambda2)*miu_est+lambda2*delta_t*(t-1))
    Nt2<-N0*exp((lambda1-lambda2)*miu_est2+lambda2*delta_t*t)
    
    sigma<-aux1*v[t]+2*aux3*m[t]*v[t]-aux3*m[t]/(2*Nt)+lambda2*delta_t/(2*Nt)
    miu1<-aux3*m[t]^2+aux2*m[t]+aux5
    
    sigma_star<-v[t]+(aux1-1)*sigma+2*aux3*m[t+1]*sigma-aux3*m[t+1]/(2*Nt2)+lambda2*delta_t/(2*Nt2)
    miu2<-m[t]+aux3*miu1^2+(aux2-1)*miu1+aux5
    
    sigma2<-0.5*(sigma+sigma_star)
    miu<-0.5*(miu1+miu2)
    
    if(sigma2<0)
    {
      sigma2<-1e-300
      write.table("sigma2<0",paste(filenames12[i], ".csv", sep=""),
                  append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
    }
    if(v[t+1]==0)
    {
      v[t+1]<-v[t+2]
    }
    f=f-n.track*log(sigma2)-(n.track*(m[t+1]-miu)^2+(n.track-1)*v[t+1])/sigma2
  }
  return(0.5*f)
}

f_lambda_v=function(theta,yobsv,delta_t,n.track,k)
{
  m<-yobsv[,1]
  v<-yobsv[,2]
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  
  aux1<-delta_t*(2*(lambda1*alpha-lambda2*beta)-(lambda1+lambda2))+1
  aux2<-delta_t*(lambda1*alpha-lambda2*(1+beta))+1
  aux3<-(lambda2-lambda1)*delta_t
  aux5<-lambda2*beta*delta_t
  
  t<-k-1
  if(t<4)
  {
    miu_est<-0.5*delta_t*(re[t]*(m[t-re[t]]+m[t]))
    miu_est2<-0.5*delta_t*((re[t]+1)*(m[t-re[t]]+m[t+1]))
  }else{
    s1<-ss[1:de[t]]
    s2<-ss[2:(de[t]+1)]
    middle<-3*sum(m[s1]+m[s2])
    miu_est<-0.5*delta_t*(middle+re[t]*(m[t-re[t]]+m[t]))
    miu_est2<-0.5*delta_t*(middle+(re[t]+1)*(m[t-re[t]]+m[t+1]))
  }
  
  Nt<-N0*exp((lambda1-lambda2)*miu_est+lambda2*delta_t*(t-1))
  Nt2<-N0*exp((lambda1-lambda2)*miu_est2+lambda2*delta_t*t)
  
  sigma<-aux1*v[t]+2*aux3*m[t]*v[t]-aux3*m[t]/(2*Nt)+lambda2*delta_t/(2*Nt)
  miu1<-aux3*m[t]^2+aux2*m[t]+aux5
  
  sigma_star<-v[t]+(aux1-1)*sigma+2*aux3*m[t+1]*sigma-aux3*m[t+1]/(2*Nt2)+lambda2*delta_t/(2*Nt2)
  miu2<-m[t]+aux3*miu1^2+(aux2-1)*miu1+aux5
  
  sigma2<-0.5*(sigma+sigma_star)
  miu<-0.5*(miu1+miu2)
  
  if(sigma2<0)
  {
    sigma2<-1e-300
    write.table("sigma2<0",paste(filenames12[i], ".csv", sep=""),
                append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
  }
  if(v[t+1]==0)
  {
    v[t+1]<-v[t+2]
  }
  f=(n.track-3)*log(v[t+1])-(n.track*(m[t+1]-miu)^2+(n.track-1)*v[t+1])/sigma2
  
  t<-k
  if(t<4)
  {
    miu_est<-0.5*delta_t*(re[t]*(m[t-re[t]]+m[t]))
    miu_est2<-0.5*delta_t*((re[t]+1)*(m[t-re[t]]+m[t+1]))
  }else{
    s1<-ss[1:de[t]]
    s2<-ss[2:(de[t]+1)]
    middle<-3*sum(m[s1]+m[s2])
    miu_est<-0.5*delta_t*(middle+re[t]*(m[t-re[t]]+m[t]))
    miu_est2<-0.5*delta_t*(middle+(re[t]+1)*(m[t-re[t]]+m[t+1]))
  }
  
  Nt<-N0*exp((lambda1-lambda2)*miu_est+lambda2*delta_t*(t-1))
  Nt2<-N0*exp((lambda1-lambda2)*miu_est2+lambda2*delta_t*t)
  
  sigma<-aux1*v[t]+2*aux3*m[t]*v[t]-aux3*m[t]/(2*Nt)+lambda2*delta_t/(2*Nt)
  miu1<-aux3*m[t]^2+aux2*m[t]+aux5
  
  sigma_star<-v[t]+(aux1-1)*sigma+2*aux3*m[t+1]*sigma-aux3*m[t+1]/(2*Nt2)+lambda2*delta_t/(2*Nt2)
  miu2<-m[t]+aux3*miu1^2+(aux2-1)*miu1+aux5
  
  sigma2<-0.5*(sigma+sigma_star)
  miu<-0.5*(miu1+miu2)
  
  if(sigma2<0)
  {
    sigma2<-1e-300
    write.table("sigma2<0",paste(filenames12[i], ".csv", sep=""),
                append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
  }
  if(v[t+1]==0)
  {
    v[t+1]<-v[t+2]
  }
  f=f-n.track*log(sigma2)-(n.track*(m[t+1]-miu)^2+(n.track-1)*v[t+1])/sigma2
  
  return(0.5*f)
}

#alpha,beta,0~1
R1_lambda = function(para.t,para.star,theta.t,theta.star,yobsv,delta_t,sd,n.track)
{exp(f_lambda1(theta.star,yobsv,delta_t,n.track)-f_lambda1(theta.t,yobsv,delta_t,n.track)+log(dtruncnorm(para.t,a=0,b=1,mean=para.star,sd=sd))-log(dtruncnorm(para.star,a=0,b=1,mean=para.t,sd=sd)))}

#lambda,0~log(2)
R2_lambda = function(para.t,para.star,theta.t,theta.star,yobsv,delta_t,sd,n.track)
{exp(f_lambda1(theta.star,yobsv,delta_t,n.track)-f_lambda1(theta.t,yobsv,delta_t,n.track)+log(dtruncnorm(para.t,a=0,b=log(2),mean=para.star,sd=sd))-log(dtruncnorm(para.star,a=0,b=log(2),mean=para.t,sd=sd)))}

#m,0~1
R3_lambda = function(para.t,para.star,theta,yobsv.t,yobsv.star,delta_t,sd,n.track,k)
{exp(f_lambda_m(theta,yobsv.star,delta_t,n.track,k)-f_lambda_m(theta,yobsv.t,delta_t,n.track,k)+log(dtruncnorm(para.t,a=0,b=1,mean=para.star,sd=sd))-log(dtruncnorm(para.star,a=0,b=1,mean=para.t,sd=sd)))}

#v,0~0.05
R4_lambda = function(para.t,para.star,theta,yobsv.t,yobsv.star,delta_t,sd,n.track,k)
{exp(f_lambda_v(theta,yobsv.star,delta_t,n.track,k)-f_lambda_v(theta,yobsv.t,delta_t,n.track,k)+log(dtruncnorm(para.t,a=0,b=0.05,mean=para.star,sd=sd))-log(dtruncnorm(para.star,a=0,b=0.05,mean=para.t,sd=sd)))}
