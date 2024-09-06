#Judging convergence
check_convergence<-function(yobsv,n.par,delta_t,n.track)
{
  varnames <- paste("mc", 1:10, sep="")
  filenames6 <- paste("not_0betaMPSRF_lambda", 1:10000, sep="")
  sd_pro<-proposal_sd(yobsv)
  
  n0<-20000
  n<-n0+2
  #Given an initial value
  assign(varnames[1],value=c(runif(2,0,1),runif(2,0,log(2))))
  assign(varnames[2],value=c(runif(2,0,1),runif(2,0,log(2))))
  assign(varnames[3],value=c(runif(2,0,1),runif(2,0,log(2))))
  
  
  flag<-0
  repeat
  {
    x.val<-mcmcsample(yobsv,n.par,delta_t,n.track,n,n0,mc1,mc2,mc3,mc4,sd_pro)
    
    chains<-NULL
    chains[[1]]<-as.mcmc(x.val[[1]])
    chains[[2]]<-as.mcmc(x.val[[2]])
    chains[[3]]<-as.mcmc(x.val[[3]])

   # MPSRF<-gelman.diag(chains,autoburnin=FALSE,multivariate = FALSE)[[1]][1]
    MPSRF<-gelman.diag(chains,autoburnin=FALSE,multivariate = TRUE)[[2]]
    write.table(MPSRF,paste(filenames6[i], ".csv", sep=""),append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
    rm(chains)
    gc()
    #if(MPSRF<1.1||(n>5e+06))
      if(MPSRF<1.1||(n>1e+07))
    {
       if(MPSRF<1.1){flag<-1}
      
        steady<-matrix(data=NA,nrow=3*n/2,ncol=n.par)
        steady[1:(n/2),]<-x.val[[1]]
        steady[(n/2+1):n,]<-x.val[[2]]
        steady[(n+1):(3*n/2),]<-x.val[[3]]
        
        num.par<-rep(0,4)
        #Determine whether the confidence interval covers the true value
        if(theta[1]<=quantile(steady[,1],probs =c(0.025,0.975))[2]&&
           theta[1]>=quantile(steady[,1],probs =c(0.025,0.975))[1])
          num.par[1]=num.par[1]+1
        if(theta[2]<=quantile(steady[,2],probs =c(0.025,0.975))[2]&&
           theta[2]>=quantile(steady[,2],probs =c(0.025,0.975))[1])
          num.par[2]=num.par[2]+1
        if(theta[3]<=quantile(steady[,3],probs =c(0.025,0.975))[2]&&
           theta[3]>=quantile(steady[,3],probs =c(0.025,0.975))[1])
          num.par[3]=num.par[3]+1
        if(theta[4]<=quantile(steady[,4],probs =c(0.025,0.975))[2]&&
           theta[4]>=quantile(steady[,4],probs =c(0.025,0.975))[1])
          num.par[4]=num.par[4]+1
        
          #confidence interval
          confidence.int<-c(quantile(steady[,1],probs =c(0.025,0.975))[1],quantile(steady[,1],probs =c(0.025,0.975))[2],
                            quantile(steady[,2],probs =c(0.025,0.975))[1],quantile(steady[,2],probs =c(0.025,0.975))[2],
                            quantile(steady[,3],probs =c(0.025,0.975))[1],quantile(steady[,3],probs =c(0.025,0.975))[2],
                            quantile(steady[,4],probs =c(0.025,0.975))[1],quantile(steady[,4],probs =c(0.025,0.975))[2])
          
          x<-apply(steady,2,mean)
          x.mean<-x[1:4]
        
          f_esti_ex<-f_lambda(x.mean[1],x.mean[2],x.mean[3],x.mean[4],yobsv[,1],
                              yobsv[,2],delta_t,n.track,nrow(yobsv))
        
        #the length of the confidence interval
        int.len<-c(quantile(steady[,1],probs =c(0.025,0.975))[2]-quantile(steady[,1],probs =c(0.025,0.975))[1],
                   quantile(steady[,2],probs =c(0.025,0.975))[2]-quantile(steady[,2],probs =c(0.025,0.975))[1],
                   quantile(steady[,3],probs =c(0.025,0.975))[2]-quantile(steady[,3],probs =c(0.025,0.975))[1],
                   quantile(steady[,4],probs =c(0.025,0.975))[2]-quantile(steady[,4],probs =c(0.025,0.975))[1])

        #square error
        bias_ex<-c((x.mean[1]-theta[1])^2,(x.mean[2]-theta[2])^2,(x.mean[3]-theta[3])^2,(x.mean[4]-theta[4])^2)
        
        break
     
    }else
    {
      for(flag2 in 1:3)
      {
        assign(varnames[flag2],value=x.val[[flag2]])
      }
      n<-n+n0 
      rm(x.val)
      gc()
    }
    
  }
  return(c(flag,x.mean,num.par,bias_ex,int.len,confidence.int,f_esti_ex))
}
