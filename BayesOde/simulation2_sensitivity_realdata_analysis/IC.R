#k is the model parameters, n is the sample size
#x.max is the estimate of the parameter, f_esti_max posterior, pd is the log-likelihood of -2 times the sample obtained from mcmc sampling
#Information criterion for the four data sets together
IC<-function(k,n,likeli,pd)
{
  PD<-(-2)*(mean(pd)-likeli)
  DIC1<-(-2)*likeli+2*PD
  DIC2<-(-2)*likeli+var((-2)*pd)IC<-(-2)*likeli+2*k  #-2IC<-(-2)*likeli+k*log(n)#-2?n(c(DIC1,DIC2,AIC,BIC))
}
