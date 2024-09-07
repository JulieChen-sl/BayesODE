mcmcsample_beta0<-function(yobsv_missing,n.par,delta_t,n.track,n,n0,
                     mc1,mc2,mc3,mc4,sd_pro)
{
  chains<-NULL
  varnames <- paste("mc", 1:10, sep="")
  
  for(s in 1:4)
  { 
    x.val1<-matrix(data=NA,ncol=n.par,nrow=n/2+n0/2)
    x.val1[1:(n/2-n0/2),]<-get(varnames[s])
    for(k in (n/2-n0/2):(n/2+n0/2-1))
    {  
      yobsv<-reset_yobsv_lambda(x.val1[k,],yobsv_missing)
      
      para.t<-x.val1[k,1]
      theta.t<-c(para.t,x.val1[k,2:4])
      para.star<-rtruncnorm(1,a=0,b=1,mean=para.t,sd=sd_pro[1])
      theta.star<-c(para.star,x.val1[k,2:4])
      R<-R1_lambda(para.t,para.star,theta.t,theta.star,yobsv,delta_t,sd_pro[1],n.track)
      p <-min(R,1)
      d <-rbinom(1,1,p) 
      x.val1[k+1,1]<-para.star*d+para.t*(1-d)
      
      x.val1[k+1,2]<-0
      
      para.t<-x.val1[k,3]
      theta.t<-c(x.val1[k+1,1:2],para.t,x.val1[k,4])
      para.star<-rtruncnorm(1,a=0,b=log(2),mean=para.t,sd=sd_pro[3])
      theta.star<-c(x.val1[k+1,1:2],para.star,x.val1[k,4])
      R<-R2_lambda(para.t,para.star,theta.t,theta.star,yobsv,delta_t,sd_pro[3],n.track)
      p <-min(R,1)
      d <-rbinom(1,1,p)  
      x.val1[k+1,3]<-para.star*d+para.t*(1-d)
      
      para.t<-x.val1[k,4]
      theta.t<-c(x.val1[k+1,1:3],para.t)
      para.star<-rtruncnorm(1,a=0,b=log(2),mean=para.t,sd=sd_pro[4])
      theta.star<-c(x.val1[k+1,1:3],para.star)
      R<-R2_lambda(para.t,para.star,theta.t,theta.star,yobsv,delta_t,sd_pro[4],n.track)
      p <-min(R,1)
      d <-rbinom(1,1,p)  
      x.val1[k+1,4]<-para.star*d+para.t*(1-d)
      
      theta.new<-x.val1[k+1,1:4]
      for(m in 5:((n.par-4)/2+4))
      {
        para.t<-x.val1[k,m]
        theta.t<-c(x.val1[k+1,1:(m-1)],x.val1[k,m:n.par])
        para.star<-rtruncnorm(1,a=0,b=1,mean=para.t,sd=sd_pro[m])
        theta.star<-c(x.val1[k+1,1:(m-1)],para.star,x.val1[k,(m+1):n.par])
        yobsv.star<-reset_yobsv_lambda(theta.star,yobsv_missing)
        yobsv.t<-reset_yobsv_lambda(theta.t,yobsv_missing)
        flag_col<-(m-5)%/%24+1
        flag_row<-(m-4-24*(flag_col-1))*1.5+0.5*(((m-4-24*(flag_col-1))%%2)==1)
        R<-R3_lambda(para.t,para.star,theta.new,yobsv.t[,(2*flag_col-1):(2*flag_col)],
                     yobsv.star[,(2*flag_col-1):(2*flag_col)],delta_t,sd_pro[m],
                     n.track,flag_row)
        p <-min(R,1)
        d <-rbinom(1,1,p)  
        x.val1[k+1,m]<-para.star*d+para.t*(1-d)
      }
      
      for(m in ((n.par-4)/2+5):(n.par-1))
      {
        para.t<-x.val1[k,m]
        theta.t<-c(x.val1[k+1,1:(m-1)],x.val1[k,m:n.par])
        para.star<-rtruncnorm(1,a=0,b=0.05,mean=para.t,sd=sd_pro[m])
        theta.star<-c(x.val1[k+1,1:(m-1)],para.star,x.val1[k,(m+1):n.par])
        yobsv.star<-reset_yobsv_lambda(theta.star,yobsv_missing)
        yobsv.t<-reset_yobsv_lambda(theta.t,yobsv_missing)
        flag_col<-(m-5-12*ncol(yobsv))%/%24+1
        flag_row<-(m-4-24*(flag_col-1+ncol(yobsv)/2))*1.5+0.5*(((m-4-24*(flag_col-1+ncol(yobsv)/2))%%2)==1)
        R<-R4_lambda(para.t,para.star,theta.new,yobsv.t[,(2*flag_col-1):(2*flag_col)],
                     yobsv.star[,(2*flag_col-1):(2*flag_col)],delta_t,sd_pro[m],
                     n.track,flag_row)
        p <-min(R,1)
        d <-rbinom(1,1,p)  
        x.val1[k+1,m]<-para.star*d+para.t*(1-d)
      }
      
      para.t<-x.val1[k,n.par]
      theta.t<-c(x.val1[k+1,1:(n.par-1)],x.val1[k,n.par])
      para.star<-rtruncnorm(1,a=0,b=0.05,mean=para.t,sd=sd_pro[n.par])
      theta.star<-c(x.val1[k+1,1:(n.par-1)],para.star)
      yobsv.star<-reset_yobsv_lambda(theta.star,yobsv_missing)
      yobsv.t<-reset_yobsv_lambda(theta.t,yobsv_missing)
      flag_row<-(n.par-4-24*(flag_col-1+ncol(yobsv)/2))*1.5+0.5*(((n.par-4-24*(flag_col-1+ncol(yobsv)/2))%%2)==1)
      R<-R4_lambda(para.t,para.star,theta.new,yobsv.t[,(ncol(yobsv)-1):ncol(yobsv)],
                   yobsv.star[,(ncol(yobsv)-1):ncol(yobsv)],delta_t,sd_pro[n.par],
                   n.track,flag_row)
      p <-min(R,1)
      d <-rbinom(1,1,p) 
      x.val1[k+1,n.par]<-para.star*d+para.t*(1-d)
      
    }
    chains[[s]]<-x.val1[-(1:(n0/2)),]
  }
  return(chains)
}