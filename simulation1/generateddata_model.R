#generate data by our model
#(Given the true values of the parameters, 
#  generate the expectation and variance of the proportional data)

#theta: the true values of the parameters
#n.data: the length of data
#n.track: the number of trajectories
#m0: The initial mean value
#N0: The initial number of cell population

generateddata_model<-function(theta,n.data,n.track,m0,N0,delta_t)#给定参数真值，生成比例数据
{
  m<-rep(0,n.data)
  v<-rep(0,n.data)
  miu<-rep(0,n.data)
  sigma2<-rep(0,n.data)
  Nt<-rep(0,n.data)

  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  
  aux1<-delta_t*(2*(lambda1*alpha-lambda2*beta)-(lambda1+lambda2))+1
  aux2<-delta_t*(lambda1*alpha-lambda2*(1+beta))+1
  aux3<-(lambda2-lambda1)*delta_t
  #aux4<-(lambda1*alpha-lambda2*beta)*delta_t
  aux5<-lambda2*beta*delta_t
  
  m[1]<-m0
  v[1]<-runif(1,0,0.01)
  Nt[1]<-N0
  
  for(i in 1:(n.data-1))
  {
    sigma<-aux1*v[i]+2*aux3*m[i]*v[i]-aux3*m[i]/(2*Nt[i])+lambda2*delta_t/(2*Nt[i])
    miu1<-aux3*m[i]^2+aux2*m[i]+aux5
    
    Nt[i+1]<-N0*exp((lambda1-lambda2)*0.5*delta_t*(sum(m[1:i]+m[2:(i+1)]))
                  +lambda2*delta_t*i)
    
    sigma_star<-v[i]+(aux1-1)*sigma+2*aux3*miu1*sigma-aux3*miu1/(2*Nt[i+1])+lambda2*delta_t/(2*Nt[i+1])
    miu2<-m[i]+aux3*miu1^2+(aux2-1)*miu1+aux5
    
    sigma2[i+1]<-0.5*(sigma+sigma_star)
    miu[i+1]<-0.5*(miu1+miu2)
    m[i+1]<-rnorm(1,miu[i+1],sqrt(sigma2[i+1]/n.track))
    v[i+1]<-sigma2[i+1]*rchisq(1,n.track-1)/(n.track-1)
  }
  yobsv<-cbind(m,v)
  
  return(yobsv)
}

#generate 100 groups of data by our model
n.data=37
t=24
delta_t<-t/(n.data-1)
n.par<-4 
n.track=5;
N0<-1000
filenames1 <- paste("DATA", 1:10000, sep="")

para<-read.csv("para_unif.csv",header=TRUE)

for(i in 1:100)
{
  theta<-as.numeric(para[i,2:5])
  m0<-para[i,6]
  yobsv<-generateddata_model(theta,n.data,n.track,m0,N0,delta_t)
  write.table(yobsv,paste(filenames1[i], ".csv", sep=""),append=FALSE,col.names=FALSE, sep = " ",eol = "\n")
  
}
