#Judging convergence,beta!=0,lambda1=lambda2
check_convergence_not0beta_EqualLambda<-function(yobsv_missing,n.par,delta_t,n.track)
{
  varnames <- paste("mc", 1:10, sep="")
  filenames6 <- paste("MPSRF", 1:100000, sep="")
  
  n0<-10000
  n<-n0+2
  sd_pro<-proposal_sd(yobsv_missing)
  sd_pro<-sd_pro[-4]
  
  #Given an initial value
  assign(varnames[1],value=runif(3,0,0.1))
  assign(varnames[2],value=c(runif(2,0.9,1),runif(1,0.6,log(2))))
  assign(varnames[3],value=c(runif(2,0,0.1),runif(1,0.6,log(2))))
  assign(varnames[4],value=c(runif(2,0.9,1),runif(1,0,0.1)))
  k1<-1;k2<-2
  for(k_mc in 1:4)
  {
    theta_initial<-get(varnames[k_mc])
    yobsv<-add_data_lambda(yobsv_missing,theta=c(theta_initial,theta_initial[3]),k1=1,k2=2,n.track)
    c1<-NULL;c2<-NULL
    for(kcol in 1:(ncol(yobsv)/2))
    {
      for(a in 1:12)
      {
        c1<-c(c1,yobsv[((k1+k2)*a-k2+1):((k1+k2)*a),2*kcol-1])
        c2<-c(c2,yobsv[((k1+k2)*a-k2+1):((k1+k2)*a),2*kcol])
      }
    }
    assign(varnames[k_mc],value=c(theta_initial,c1,c2))
  }
  
  repeat
  {
    x.val<-mcmcsample_lambdaequal(yobsv_missing,n.par,delta_t,n.track,n,n0,
                      mc1,mc2,mc3,mc4,sd_pro)
    chains<-NULL
    chains[[1]]<-as.mcmc(x.val[[1]])
    chains[[2]]<-as.mcmc(x.val[[2]])
    chains[[3]]<-as.mcmc(x.val[[3]])
    chains[[4]]<-as.mcmc(x.val[[4]])
    
    #MPSRF<-gelman.diag(chains,autoburnin=FALSE,multivariate = FALSE)[[1]][1]
    MPSRF<-gelman.diag(chains,autoburnin=FALSE,multivariate = TRUE)[[2]]
    write.table(MPSRF,paste(filenames6[i], ".csv", sep=""),append=TRUE,col.names=FALSE, sep = " ",eol = "\n")
    rm(chains)
    gc()
    if(MPSRF<1.1||(n>8e+06))
    {
      steady<-matrix(data=NA,nrow=4*n/2,ncol=n.par-1)
      steady[1:(n/2),]<-x.val[[1]]
      steady[(n/2+1):n,]<-x.val[[2]]
      steady[(n+1):(3*n/2),]<-x.val[[3]]
      steady[(3*n/2+1):(2*n),]<-x.val[[4]]
      
      num.par<-rep(0,3)
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
      
      #point estimation
      x<-apply(steady,2,mean)
      x.mean<-x[1:3]
      yobsv_ex<-reset_yobsv_lambda(c(x[1:3],x[3],x[-(1:3)]),yobsv_missing)
      f_esti_ex<-f_lambda(c(x.mean,x.mean[3]),yobsv_ex,delta_t,n.track)
      
      
      #confidence interval
      confidence.int<-c(quantile(steady[,1],probs =c(0.025,0.975))[1],quantile(steady[,1],probs =c(0.025,0.975))[2],
                        quantile(steady[,2],probs =c(0.025,0.975))[1],quantile(steady[,2],probs =c(0.025,0.975))[2],
                        quantile(steady[,3],probs =c(0.025,0.975))[1],quantile(steady[,3],probs =c(0.025,0.975))[2])
      
      #the length of the confidence interval
      int.len<-c(quantile(steady[,1],probs =c(0.025,0.975))[2]-quantile(steady[,1],probs =c(0.025,0.975))[1],
                 quantile(steady[,2],probs =c(0.025,0.975))[2]-quantile(steady[,2],probs =c(0.025,0.975))[1],
                 quantile(steady[,3],probs =c(0.025,0.975))[2]-quantile(steady[,3],probs =c(0.025,0.975))[1])
      
      #square error
      bias_ex<-c((x.mean[1]-theta[1])^2,(x.mean[2]-theta[2])^2,(x.mean[3]-theta[3])^2)
      
      #Computing IC
      max_value<-maxlikeli(cbind(steady[,1:3],steady[,3],steady[,-(1:3)]),yobsv_missing)
      pd<-max_value[(n.par+2):length(max_value)]
      IC_ex<-IC(k=n.par-1,n=13,f_esti_ex,pd)
      
      break
    }else
    {
      for(flag2 in 1:4)
      {
        assign(varnames[flag2],value=x.val[[flag2]])
      }
      n<-n+n0
      rm(x.val)
      gc()
    }
  }
  
  return(c(x.mean,num.par,bias_ex,int.len,confidence.int,f_esti_ex,IC_ex))
  
  
}