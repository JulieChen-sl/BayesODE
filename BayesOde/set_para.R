#Generate 1000 sets of parameter truth values for simulation

k=1000
para<-matrix(0,ncol=5,nrow=2*k)
for(i in 1:k)
{
  theta<-c(runif(2,0,1),runif(2,0,log(2)))
  m_0<-runif(1,0,1)
  para[2*i-1,]<-c(theta,m_0)
  
  theta2<-c(1-theta[2],1-theta[1],theta[4],theta[3])
  para[2*i,]<-c(theta2,1-m_0)
}
write.csv(para,"para_unif.csv")