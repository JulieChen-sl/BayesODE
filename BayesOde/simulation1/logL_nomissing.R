#log likelihood####
f_lambda=function(alpha,beta, lambda1, lambda2, m,v,delta_t,n.track,T)#T=length(m)
{
  f=0
  aux1<-delta_t*(2*(lambda1*alpha-lambda2*beta)-(lambda1+lambda2))+1
  aux2<-delta_t*(lambda1*alpha-lambda2*(1+beta))+1
  aux3<-(lambda2-lambda1)*delta_t
  #aux4<-(lambda1*alpha-lambda2*beta)*delta_t
  aux5<-lambda2*beta*delta_t
  Nt<-rep(0,T)
  Nt[1]<-N0
  for(i in 2:T)
  {
    Nt[i]<-N0*exp((lambda1-lambda2)*0.5*delta_t*(sum(m[1:(i-1)]+m[2:i]))+
                    lambda2*delta_t*(i-1))
  }
  
  for(t in 1:(T-1))
  {
    sigma <- aux1*v[t]+2*aux3*m[t]*v[t]-aux3*m[t]/(2*Nt[t])+lambda2*delta_t/(2*Nt[t])
    miu1 <- aux3*m[t]^2+aux2*m[t]+aux5
    
    sigma_star <- v[t]+(aux1-1)*sigma+2*aux3*miu1*sigma-aux3*miu1/(2*Nt[t+1])+
                  lambda2*delta_t/(2*Nt[t+1])
    miu2 <- m[t]+aux3*miu1^2+(aux2-1)*miu1+aux5
    
    sigma2 <- 0.5*(sigma+sigma_star)
    miu <- 0.5*(miu1+miu2)
    
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
    f=f-n.track*log(sigma2)+(n.track-3)*log(v[t+1])-
      (n.track*(m[t+1]-miu)^2+(n.track-1)*v[t+1])/sigma2
  }
  return(0.5*f)
}

R1 = function(alpha.t,alpha.star,beta, lambda1, lambda2,m,v,delta_t,n.track,T,sd)
{exp(f_lambda(alpha.star,beta, lambda1, lambda2, m,v,delta_t,n.track,T)
     -f_lambda(alpha.t,beta, lambda1, lambda2, m,v,delta_t,n.track,T)
     +log(dtruncnorm(alpha.t,a=0,b=1,mean=alpha.star,sd=sd))
     -log(dtruncnorm(alpha.star,a=0,b=1,mean=alpha.t,sd=sd)))}

R2 = function(alpha,beta.t,beta.star, lambda1, lambda2,m,v,delta_t,n.track,T,sd)
{exp(f_lambda(alpha,beta.star, lambda1, lambda2, m,v,delta_t,n.track,T)
     -f_lambda(alpha,beta.t, lambda1, lambda2, m,v,delta_t,n.track,T)
     +log(dtruncnorm(beta.t,a=0,b=1,mean=beta.star,sd=sd))
     -log(dtruncnorm(beta.star,a=0,b=1,mean=beta.t,sd=sd)))}

R3 = function(alpha,beta, lambda1.t,lambda1.star, lambda2,m,v,delta_t,n.track,T,sd)
{exp(f_lambda(alpha,beta, lambda1.star, lambda2, m,v,delta_t,n.track,T)
     -f_lambda(alpha,beta, lambda1.t, lambda2, m,v,delta_t,n.track,T)
     +log(dtruncnorm(lambda1.t,a=0,b=log(2),mean=lambda1.star,sd=sd))
     -log(dtruncnorm(lambda1.star,a=0,b=log(2),mean=lambda1.t,sd=sd)))}

R4 = function(alpha,beta, lambda1,lambda2.t, lambda2.star,m,v,delta_t,n.track,T,sd)
{exp(f_lambda(alpha,beta, lambda1, lambda2.star, m,v,delta_t,n.track,T)
     -f_lambda(alpha,beta, lambda1, lambda2.t, m,v,delta_t,n.track,T)
     +log(dtruncnorm(lambda2.t,a=0,b=log(2),mean=lambda2.star,sd=sd))
     -log(dtruncnorm(lambda2.star,a=0,b=log(2),mean=lambda2.t,sd=sd)))}

