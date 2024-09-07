#Write a function that calculates the standard deviation of the proposed distribution
proposal_sd<-function(yobsv_missing)
{
  k1<-1;k2<-2
  theta_initial<-c(runif(2,0,1),runif(2,0,log(2)))
  yobsv<-add_data_lambda(yobsv_missing,theta=theta_initial,k1=k1,k2=k2,n.track)
  c1<-rep(0,k2*12);c2<-rep(0,k2*12)
  c1<-NULL;c2<-NULL
  for(kcol in 1:(ncol(yobsv)/2))
  {
    for(a in 1:12)
    {
      c1<-c(c1,yobsv[((k1+k2)*a-k2+1):((k1+k2)*a),2*kcol-1])
      c2<-c(c2,yobsv[((k1+k2)*a-k2+1):((k1+k2)*a),2*kcol])
    }
  }
  assign("mc",value=c(theta_initial,c1,c2))
  sd_proposal<-c(rep(0.02,(4+(n.par-4)/2)),rep(0.001,n.par-(4+(n.par-4)/2)))
  n0<-1000
  n<-n0+2
  flag<-1
  while(flag==1)
  {
    d_accept<-rep(0,n.par)
    x.val1<-matrix(data=NA,ncol=n.par,nrow=n/2+n0/2)
    x.val1[1:(n/2-n0/2),]<-get("mc")
    
    for(k in (n/2-n0/2):(n/2+n0/2-1))
    {  
      yobsv<-reset_yobsv_lambda(x.val1[k,],yobsv_missing)
      
      para.t<-x.val1[k,1]
      theta.t<-c(para.t,x.val1[k,2:4])
      para.star<-rtruncnorm(1,a=0,b=1,mean=para.t,sd=sd_proposal[1])
      theta.star<-c(para.star,x.val1[k,2:4])
      R<-R1_lambda(para.t,para.star,theta.t,theta.star,yobsv,delta_t,sd_proposal[1],n.track)
      p <-min(R,1)
      d <-rbinom(1,1,p)  
      x.val1[k+1,1]<-para.star*d+para.t*(1-d)
      d_accept[1]<-d_accept[1]+d
      
      para.t<-x.val1[k,2]
      theta.t<-c(x.val1[k+1,1],para.t,x.val1[k,3:4])
      para.star<-rtruncnorm(1,a=0,b=1,mean=para.t,sd=sd_proposal[2])
      theta.star<-c(x.val1[k+1,1],para.star,x.val1[k,3:4])
      R<-R1_lambda(para.t,para.star,theta.t,theta.star,yobsv,delta_t,sd_proposal[2],n.track)
      p <-min(R,1)
      d <-rbinom(1,1,p)  
      x.val1[k+1,2]<-para.star*d+para.t*(1-d)
      d_accept[2]<-d_accept[2]+d
      
      para.t<-x.val1[k,3]
      theta.t<-c(x.val1[k+1,1:2],para.t,x.val1[k,4])
      para.star<-rtruncnorm(1,a=0,b=log(2),mean=para.t,sd=sd_proposal[3])
      theta.star<-c(x.val1[k+1,1:2],para.star,x.val1[k,4])
      R<-R2_lambda(para.t,para.star,theta.t,theta.star,yobsv,delta_t,sd_proposal[3],n.track)
      p <-min(R,1)
      d <-rbinom(1,1,p) 
      x.val1[k+1,3]<-para.star*d+para.t*(1-d)
      d_accept[3]<-d_accept[3]+d
      
      para.t<-x.val1[k,4]
      theta.t<-c(x.val1[k+1,1:3],para.t)
      para.star<-rtruncnorm(1,a=0,b=log(2),mean=para.t,sd=sd_proposal[4])
      theta.star<-c(x.val1[k+1,1:3],para.star)
      R<-R2_lambda(para.t,para.star,theta.t,theta.star,yobsv,delta_t,sd_proposal[4],n.track)
      p <-min(R,1)
      d <-rbinom(1,1,p)  
      x.val1[k+1,4]<-para.star*d+para.t*(1-d)
      d_accept[4]<-d_accept[4]+d
      
      theta.new<-x.val1[k+1,1:4]
      for(m in 5:((n.par-4)/2+4))
      {
        para.t<-x.val1[k,m]
        theta.t<-c(x.val1[k+1,1:(m-1)],x.val1[k,m:n.par])
        para.star<-rtruncnorm(1,a=0,b=1,mean=para.t,sd=sd_proposal[m])
        theta.star<-c(x.val1[k+1,1:(m-1)],para.star,x.val1[k,(m+1):n.par])
        yobsv.star<-reset_yobsv_lambda(theta.star,yobsv_missing)
        yobsv.t<-reset_yobsv_lambda(theta.t,yobsv_missing)
        flag_col<-(m-5)%/%24+1
        flag_row<-(m-4-24*(flag_col-1))*1.5+0.5*(((m-4-24*(flag_col-1))%%2)==1)
        R<-R3_lambda(para.t,para.star,theta.new,yobsv.t[,(2*flag_col-1):(2*flag_col)],
                     yobsv.star[,(2*flag_col-1):(2*flag_col)],delta_t,sd_proposal[m],
                     n.track,flag_row)
        p <-min(R,1)
        d <-rbinom(1,1,p) 
        x.val1[k+1,m]<-para.star*d+para.t*(1-d)
        d_accept[m]<-d_accept[m]+d
      }
      for(m in ((n.par-4)/2+5):(n.par-1))
      {
        para.t<-x.val1[k,m]
        theta.t<-c(x.val1[k+1,1:(m-1)],x.val1[k,m:n.par])
        para.star<-rtruncnorm(1,a=0,b=0.05,mean=para.t,sd=sd_proposal[m])
        theta.star<-c(x.val1[k+1,1:(m-1)],para.star,x.val1[k,(m+1):n.par])
        yobsv.star<-reset_yobsv_lambda(theta.star,yobsv_missing)
        yobsv.t<-reset_yobsv_lambda(theta.t,yobsv_missing)
        flag_col<-(m-5-12*ncol(yobsv))%/%24+1
        flag_row<-(m-4-24*(flag_col-1+ncol(yobsv)/2))*1.5+0.5*(((m-4-24*(flag_col-1+ncol(yobsv)/2))%%2)==1)
        R<-R4_lambda(para.t,para.star,theta.new,yobsv.t[,(2*flag_col-1):(2*flag_col)],
                     yobsv.star[,(2*flag_col-1):(2*flag_col)],delta_t,sd_proposal[m],
                     n.track,flag_row)
        p <-min(R,1)
        d <-rbinom(1,1,p)  
        x.val1[k+1,m]<-para.star*d+para.t*(1-d)
        d_accept[m]<-d_accept[m]+d
      }
      
      para.t<-x.val1[k,n.par]
      theta.t<-c(x.val1[k+1,1:(n.par-1)],x.val1[k,n.par])
      para.star<-rtruncnorm(1,a=0,b=0.05,mean=para.t,sd=sd_proposal[n.par])
      theta.star<-c(x.val1[k+1,1:(n.par-1)],para.star)
      yobsv.star<-reset_yobsv_lambda(theta.star,yobsv_missing)
      yobsv.t<-reset_yobsv_lambda(theta.t,yobsv_missing)
      flag_row<-(n.par-4-24*(flag_col-1+ncol(yobsv)/2))*1.5+0.5*(((n.par-4-24*(flag_col-1+ncol(yobsv)/2))%%2)==1)
      R<-R4_lambda(para.t,para.star,theta.new,yobsv.t[,(ncol(yobsv)-1):ncol(yobsv)],
                   yobsv.star[,(ncol(yobsv)-1):ncol(yobsv)],delta_t,sd_proposal[n.par],
                   n.track,flag_row)
      p <-min(R,1)
      d <-rbinom(1,1,p)  
      x.val1[k+1,n.par]<-para.star*d+para.t*(1-d)
      d_accept[n.par]<-d_accept[n.par]+d
    }
    
      flag<-0
      for(m in 1:(4+(n.par-4)/2))
      {
        if(((d_accept[m]/n0)<0.15)&&(sd_proposal[m]>1e-3))
        {
          sd_proposal[m]<-(2/3)*sd_proposal[m]
          flag<-1
        }else if(((d_accept[m]/n0)>0.35)&&(sd_proposal[m]<0.5))
        {
          sd_proposal[m]<-1.5*sd_proposal[m]
          flag<-1
        }
      }
      for(m in (5+(n.par-4)/2):n.par)
      {
        if(((d_accept[m]/n0)<0.15)&&(sd_proposal[m]>1e-7))
        {
          sd_proposal[m]<-(2/3)*sd_proposal[m]
          flag<-1
        }else if(((d_accept[m]/n0)>0.35)&&(sd_proposal[m]<0.01))
        {
          sd_proposal[m]<-1.5*sd_proposal[m]
          flag<-1
        }
      }
    assign("mc",value=x.val1[-(1:(n0/2)),])
    n<-n0+n
  }
  return(sd_proposal)
}
