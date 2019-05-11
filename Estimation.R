library("MASS")
library(foreign)
Kgauss<-function(x)
{
  Kgaus<-(1/sqrt(2*pi))*exp(-(x^2)/2)
  return(Kgaus)
}
estimasi<-function(data)
{
  y<-data[,1]
  #respon diganti 1 dan 2
  x11<-data[,2]
  x12<-data[,3]
  n<-nrow(data)
  datanew<-data.frame(y,x11,x12)
  nbGLM<-glm.nb(y~x11+x12, data=datanew)
  h1<-12927
  h2<-5
  cat("ESTIMASI BETA PARAMETRIK=\n")
  print(nbGLM)
  beta<-c(coef(nbGLM)[1],coef(nbGLM)[2],coef(nbGLM)[3])
  beta0<-coef(nbGLM)[1]
  beta11<-coef(nbGLM)[2]
  beta12<-coef(nbGLM)[3]
  ytopipar<-exp(beta0+beta11*x11+beta12*x12)
  alfa<-nbGLM$theta
  mparam=matrix(0,n,4)
  ytopi=matrix(0,n,1)
  kuadrat=matrix(0,n,1)
  param<-as.numeric(c(beta,alfa))
  for(a in 1:n)
  {
    wparam<-rep(0,4)
    x01<-x11[a]
    x02<-x12[a]
    x1<-x11-x01
    x2<-x12-x02
    X12<-cbind(x1,x2)
    p<-ncol(X12)
    K1<-(1/h1)*Kgauss(x1/h1)
    K2<-(1/h2)*Kgauss(x2/h2)
    K<-cbind(K1,K2)
    X<-cbind((rep(1,n)),X12)
    likelihood<-function(wparam)
    {
      betaa<-wparam[1:(p+1)]
      alf<-abs(wparam[p+2])
      m=exp(X%*%betaa)
      k=X%*%betaa
      lnf=rep(0,n)
      prodK=rep(0,n)
      A=rep(0,n)
      for(i in 1:n)
      {
        lnf[i]<-lgamma(y[i]+(1/alf))-lfactorial(y[i])-lgamma(1/alf)-(1/alf)*log(1+alf*m[i])+(y[i]*log(alf))+(y[i]*k[i])-y[i]*log(1+alf*m[i])
        prodK[i]=prod(K[i,])
        A[i]<-lnf[i]*prodK[i]
      }
      hasil=sum(A[1:n])
    }
    #maxit yg optimal 122
    newrap=optim(param,likelihood,control=list(fnscale=-1,maxit=122),hessian=TRUE)
    mparam[a,]<-as.numeric(newrap$par)
    ytopi[a]<-exp(mparam[a,1]+mparam[a,2]*(x1[a])+mparam[a,3]*(x2[a]))
    hessian<-newrap$hess
    s<-c(1,1,1,1)
    kuadrat[a]<-t(s)%*%hessian%*%s
  }
  for(i in 1:n)
  {
    if (kuadrat[i]<0) cat("data ke-",i,"=",kuadrat[i],"(definit negatif)\n")
    else cat("data ke-",i,"=",kuadrat[i],"(tidak definit negatif)\n")
  }
  cat("MATRIKS PARAMETER=\n")
  alfab=mparam[,4]
  errorpar<-y-ytopipar
  error<-y-ytopi
  cat("\t yasli \t yparametrik \t error \t\t nonpar \t error \n")
  for(s in 1:n)
  {					cat("\t",y[s],"\t",ytopipar[s],"\t",errorpar[s],"\t",ytopi[s],"\t\t",error[s],"\n")
  }
  MSEpar<-sum((errorpar)^2)/n
  MSE<-sum((error)^2)/n
  cat("MSE PAR =",MSEpar,"\n")
  cat("MSE NONPAR =",MSE,"\n")
  JKT1<-t(y-(mean(y)))%*%(y-(mean(y)))
  JKG1<-t(y-ytopipar)%*%(y-ytopipar)
  JKT2<-t(y-(mean(y)))%*%(y-(mean(y)))
  JKG2<-t(y-ytopi)%*%(y-ytopi)
  R2par=1-(JKG1/JKT1)
  R2nonpar=1-(JKG2/JKT2)
  cat("R2 PAR =",R2par,"\n")
  cat("R2 NONPAR =",R2nonpar,"\n")
  devpar=nbGLM$deviance
  edev1=y*log(y/ytopi)
  edev2=y+(1/alfab)
  edev3=log((y+(1/alfab))/(ytopi+(1/alfab)))
  dev=sum(edev1-(edev2*edev3))
  chisq<-qchisq(0.01,n-p)
  cat(" DEV_PAR\t\tDEV_NONPAR\t\tCHISQUARE =\n")
  cat(devpar,"\t\t",dev,"\t\t",chisq,	"\n")
  if(dev<chisq)
    cat ("KARENA NILAI DEVIANS < NILAI CHISQUARE TABEL, \n MAKA HO DITERIMA, ARTINYA MODEL SESUAI \n")
  else cat("KARENA NILAI DEVIANS> NILAI CHISQUARE TABEL, MAKA HO DITOLAK, ARTINYA MODEL TIDAK SESUAI \n")
  i<-rep(1:n,1)
  win.graph()
  plot(i,ytopi,xlab="City/District",ylab="Number of AIDS Patient",type="n")
  plot(i,y,xlab="City/District",ylab="Number of AIDS Patient",col="red",type="p")
  lines(i,ytopi,xlab="City/District",ylab="Number of AIDS Patient",col="green",lwd=3)
  lines(i,ytopipar,xlab="City/District",ylab="Number of AIDS Patient",col="blue",lwd=3)
  
  title(main="Plot Observation and Parametric Estimation (blue) and Nonparametric Estimation (green)")
}
