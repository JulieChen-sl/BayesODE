# Write a function that expands the amount of data, inputs data with missing values, and outputs the expanded data.
add_data_lambda<-function(yobsv_missing,theta,k1,k2,n.track)# Insert k2 data every k1 data
{
  n<-n.track
  n_data<-nrow(yobsv_missing)
  yobsv<-matrix(NA,nrow=n_data+(n_data-1)*k2/k1,ncol=ncol(yobsv_missing))
  
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  
  aux1<-delta_t*(2*(lambda1*alpha-lambda2*beta)-(lambda1+lambda2))+1
  aux2<-delta_t*(lambda1*alpha-lambda2*(1+beta))+1
  aux3<-(lambda2-lambda1)*delta_t
  aux5<-lambda2*beta*delta_t
  
  for(kcol in 1:(ncol(yobsv_missing)/2))
  {
   for(i in 1:(n_data/k1-1))
   {
    yobsv[((k1+k2)*(i-1)+1):((k1+k2)*i-k2),]<-yobsv_missing[(k1*(i-1)+1):(k1*i),]
    sigma2<-rep(0,k2)
    miu<-rep(0,k2)
    sigma2[1]<-yobsv_missing[k1*i,2*kcol]
    miu[1]<-yobsv_missing[k1*i,2*kcol-1]
    j<-1
    while(j<=k2)
    {
      m<-yobsv[(k1+k2)*i-k2+j-1,2*kcol-1]
      v<-yobsv[(k1+k2)*i-k2+j-1,2*kcol] 
      
      re<-((k1+k2)*i-k2+j-2)%%3
      if(((k1+k2)*i-k2+j-1)<4)
      {
        miu_est<-(1/2)*delta_t*(re*(yobsv[(k1+k2)*i-k2+j-1-re,2*kcol-1]+yobsv[(k1+k2)*i-k2+j-1,2*kcol-1]))
      }else{
      miu_est<-(1/2)*delta_t*(3*sum(yobsv[seq(1,(k1+k2)*i-k2+j-4-re,by=3),2*kcol-1]+yobsv[seq(4,(k1+k2)*i-k2+j-1-re,by=3),2*kcol-1])
                              +re*(yobsv[(k1+k2)*i-k2+j-1-re,2*kcol-1]+yobsv[(k1+k2)*i-k2+j-1,2*kcol-1]))
      
      }
      Nt<-N0*exp((lambda1-lambda2)*miu_est+lambda2*delta_t*((k1+k2)*i-k2+j-2))
     
      sigma<-aux1*v+2*aux3*m*v-aux3*m/(2*Nt)+lambda2*delta_t/(2*Nt)#Ԥ??ֵ
      miu1<-aux3*m^2+aux2*m+aux5
      if(((k1+k2)*i-k2+j-1)<4)
      {
        miu_est2<-(1/2)*delta_t*((re+1)*(yobsv[(k1+k2)*i-k2+j-1-re,2*kcol-1]+miu1))
      }else{
        miu_est2<-(1/2)*delta_t*(3*sum(yobsv[seq(1,(k1+k2)*i-k2+j-4-re,by=3),2*kcol-1]+yobsv[seq(4,(k1+k2)*i-k2+j-1-re,by=3),2*kcol-1])
                                 +(re+1)*(yobsv[(k1+k2)*i-k2+j-1-re,2*kcol-1]+miu1))
      }
      Nt2<-N0*exp((lambda1-lambda2)*miu_est2+lambda2*delta_t*((k1+k2)*i-k2+j-1))
      sigma_star<-v+(aux1-1)*sigma+2*aux3*miu1*sigma-aux3*miu1/(2*Nt2)+lambda2*delta_t/(2*Nt2)
      sigma2[j+1]<-0.5*(sigma+sigma_star)
      miu2<-m+aux3*miu1^2+(aux2-1)*miu1+aux5
      miu[j+1]<-0.5*(miu1+miu2)
      
      
      if(sigma2[j+1]<=0)
      {
        j=j-1
      }
      yobsv[(k1+k2)*i-k2+j,2*kcol-1]<-rnorm(1,miu[j+1],sqrt(sigma2[j+1]/n))
      yobsv[(k1+k2)*i-k2+j,2*kcol]<-sigma2[j+1]*rchisq(1,n-1)/(n-1)  
      j<-j+1
    }
  }
    
  }
  yobsv[n_data+(n_data-1)*k2/k1,]<-yobsv_missing[n_data,]
  return(yobsv)
}
