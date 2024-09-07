#generate data by Gillespie’s algorithm 
#(Given the true values of the parameters, 
#  generate the expectation and variance of the proportional data)

#theta: the true values of the parameters
#n.data: the length of data
#n.track: the number of trajectories
#m0: The initial mean value
#N0: The initial number of cell population

generateddata_gillespie<-function(theta,n.data,n.track,m0,N0)
{
  m<-matrix(data=0,ncol=n.data,nrow=n.track); #record the trajectories
  v1<-c(1,0,1,0);
  v2<-c(0,1,0,1);
  s<-t/(n.data-1);
  x<-rep(0,n.data) #record the number of cell A
  y<-rep(0,n.data) #record the number of cell B
  #The initial number of two type of cells
  x[1]=round(N0*m0) 
  y[1]=N0-x[1]
  
  alpha<-theta[1]
  beta<-theta[2]
  lambda1<-theta[3]
  lambda2<-theta[4]
  
  for(sk in 1:n.track)
  {
    x2<-x[1]
    y2<-y[1]
    t1<-0;
    t2<-0;
    for(k in 1:(n.data-1))
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
    }
    m[sk,]<-x/(x+y)
  }
  yobsv<-cbind(apply(m,2,mean),apply(m,2,var))
  return(yobsv)
}

#generate 100 groups of data by gillespie
n.data=37
t=24
delta_t<-t/(n.data-1)
n.par<-4 
n.track=5;
N0<-1000#the initial total number of cell,change this parameter when do sensitivity analysis
filenames1 <- paste("DATA", 1:10000, sep="")

para<-read.csv("para_unif.csv",header=TRUE)

for(i in 1:100)
{
  theta<-as.numeric(para[i,2:5])
  m0<-para[i,6]
  yobsv<-generateddata_gillespie(theta,n.data,n.track,m0,N0)
  write.table(yobsv,paste(filenames1[i], ".csv", sep=""),append=FALSE,col.names=FALSE, sep = " ",eol = "\n")
  
}