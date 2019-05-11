bandwidth1<-function()
{
  bb1<-as.numeric(readline("Input Batas Bawah hpred 1 = "))
  ba1<-as.numeric(readline("Input Batas atas hpred 1 = "))
  inc1<-1
  h1<-seq(bb1,ba1,inc1)
  nh1<-length(h1)
  mlcv1<-rep(0,nh1)
  for(a in 1:nh1)
  {
    mlcv1[a]<-MLCV(h1[a])
  }
  cat("  TABEL BANDWIDTH OPTIMAL PREDIKTOR 2 \n")
  #cat("  TABEL BANDWIDTH OPTIMAL PREDIKTOR 2 \n")
  cat("h1	\t MLCV \n")
  for(i in 1:nh1)
  {
    cat(h1[i],"\t",mlcv1[i],"\n")
  }
  hnmlcv1<-matrix(c(h1,mlcv1),nh1,2)
  MLCVmax1<-max(mlcv1)
  hopt1<-hnmlcv1[hnmlcv1[,2]==MLCVmax1,1]
  cat("bandwidth optimal prediktor 1 =",hopt1,"\n")
  plot(h1,mlcv1,xlab="h",ylab="mlcv",col="blue",type="l")
  title(main="PLOT BANDWIDTH PADA PREDIKTOR 2 THD MLCV",col=2)
  #title(main="PLOT BANDWIDTH PADA PREDIKTOR 2 THD MLCV",col=2)
}
