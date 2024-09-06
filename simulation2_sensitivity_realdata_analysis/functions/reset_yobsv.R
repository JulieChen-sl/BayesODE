#Write a function that arranges the inserted observations back into their corresponding positions
reset_yobsv_lambda<-function(x,yobsv_missing)
{
  k1<-1;k2<-2
  yobsv<-matrix(NA,nrow=n.data+1,ncol=ncol(yobsv_missing))
  for(kcol in 1:(ncol(yobsv_missing)/2))
  {
    c1<-x[(5+24*(kcol-1)):(4+24*kcol)]
    c2<-x[(5+24*(kcol+ncol(yobsv_missing)/2-1)):(4+24*(kcol+ncol(yobsv_missing)/2))]
    for(a in 1:12)
    {
      yobsv[(k1+k2)*a-(k1+k2-1),(2*kcol-1):(2*kcol)]<-yobsv_missing[a,(2*kcol-1):(2*kcol)]
      yobsv[((k1+k2)*a-k2+1):((k1+k2)*a),2*kcol-1]<-c1[(k2*a-k2+1):(k2*a)]
      yobsv[((k1+k2)*a-k2+1):((k1+k2)*a),2*kcol]<-c2[(k2*a-k2+1):(k2*a)]
    }
  }
  yobsv[n.data+1,]<-yobsv_missing[13,]
  return(yobsv)
}
