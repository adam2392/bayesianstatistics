########################################################################
##sensitivity analysis
pdf("plot/fig1.pdf",family="Times",height=4.4,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0),oma=c(0,0,.5,0))
par(mfrow=c(2,2))
theta<-seq(0,1,length=100)
a<-1; b<-1
n<-4 ; y<-1
plot(theta,dbeta(theta,a+y,b+n-y),type="l",ylab=
       expression(paste(italic("p("),theta,"|y)",sep="")), xlab=expression(theta), 
     lwd=2)
mtext(expression(paste("beta(1,1) prior,  ", italic("n"),"=4  ",italic(sum(y[i])),"=1",sep="")), side=3,line=.1)
abline(v=c((a+y-1)/(a+b+n-2), (a+y)/(a+b+n)) , col=c("green","red") )
lines(theta,dbeta(theta,a,b),type="l",col="gray",lwd=2)
legend(.45,2.4,legend=c("prior","posterior"),lwd=c(2,2),col=c("gray","black"), bty="n")

a<-3; b<-2
n<-4 ; y<-1
plot(theta,dbeta(theta,a+y,b+n-y),type="l",ylab=
       expression(paste(italic("p("),theta,"|y)",sep="")), xlab=expression(theta), 
     lwd=2)
   expression(italic(paste("p(",theta,"|y)",sep=""))), xlab=expression(theta),lwd=2)
mtext(expression(paste("beta(3,2) prior,  ", italic("n"),"=4  ",italic(sum(y[i])),"=1",sep="")), side=3,line=.1)
abline(v=c((a+y-1)/(a+b+n-2), (a+y)/(a+b+n)) , col=c("green","red") )
lines(theta,dbeta(theta,a,b),type="l",col="gray",lwd=2)

a<-1 ; b<-1
n<-100; y<-25
plot(theta,dbeta(theta,a+y,b+n-y),type="l",ylab=
       expression(paste(italic("p("),theta,"|y)",sep="")), xlab=expression(theta), 
     lwd=2)
 expression(italic(paste("p(",theta,"|y)",sep=""))),    xlab=expression(theta),lwd=2)
mtext(expression(paste("beta(1,1) prior,  ", italic("n"),"=100  ",italic(sum(y[i])),"=25",sep="")), side=3,line=.1)
abline(v=c((a+y-1)/(a+b+n-2), (a+y)/(a+b+n)) , col=c("green","red") )
lines(theta,dbeta(theta,a,b),type="l",col="gray",lwd=2)

a<-3 ; b<-2
n<-100; y<-25
plot(theta,dbeta(theta,a+y,b+n-y),type="l",ylab=
       expression(paste(italic("p("),theta,"|y)",sep="")), xlab=expression(theta), 
     lwd=2)
   expression(italic(paste("p(",theta,"|y)",sep=""))),xlab=expression(theta),
   lwd=2)
mtext(expression(paste("beta(3,2) prior,  ", italic("n"),"=100  ",italic(sum(y[i])),"=25",sep="")), side=3,line=.1)
abline(v=c((a+y-1)/(a+b+n-2), (a+y)/(a+b+n)) , col=c("green","red") )
lines(theta,dbeta(theta,a,b),type="l",col="gray",lwd=2)
dev.off()

############################################################################
######HPD
pdf("plot/HPD.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
a = 1; b=1; n=20; y=15;
theta.support<-seq(0,1,length=5000)
plot(theta.support, dbeta(theta.support, a+y, b+n-y), type="l",
     xlab=expression(theta),ylab=expression(paste(italic("p("),theta,"|y)"))) 
pth<-dbeta(theta.support, a+y, b+n-y)
pth<-pth
ord<- order(-pth)
xpx<-cbind(theta.support[ord], pth[ord])
xpx<-cbind(xpx,cumsum(xpx[,2])/sum(xpx[,2]))

hpd<-function(x,dx,p){
  md<-x[dx==max(dx)]
  px<-dx/sum(dx)
  pxs<--sort(-px)
  ct<-min(pxs[cumsum(pxs)< p])
  list(hpdr=range(x[px>=ct]),mode=md) }

tmp<-hpd(xpx[,1],xpx[,2],.5)$hpdr
lines( x=c(tmp[1],tmp[1],tmp[2],tmp[2]),
       y=dbeta(c(0,tmp[1],tmp[2],0),a+y,b+n-y)  ,col=gray(.75),lwd=2   )
tmp<-hpd(xpx[,1],xpx[,2],.75)$hpdr
lines( x=c(tmp[1],tmp[1],tmp[2],tmp[2]),
       y=dbeta(c(0,tmp[1],tmp[2],0),a+y,b+n-y)  ,col=gray(.5),lwd=2   )
tmp<-hpd(xpx[,1],xpx[,2],.95)$hpdr
lines( x=c(tmp[1],tmp[1],tmp[2],tmp[2]),
       y=dbeta(c(0,tmp[1],tmp[2],0),a+y,b+n-y)  ,col=gray(0),lwd=2   )

tmp<-qbeta( c(.025,.975), a+y,b+n-y)
lines( x=c(tmp[1],tmp[1],tmp[2],tmp[2]),
       y=dbeta(c(0,tmp[1],tmp[2],0),a+y,b+n-y)  ,col=gray(0),lwd=2 ,lty=2  )


legend(.5, 2.75, c("50% HPD","75% HPD","95% HPD","95% quantile-based"), 
       col=c(gray(.75),gray(.5),
             gray(0),gray(0)),lty=c(1,1,1,2),lwd=c(2,2,2,2),
       bty="n")

dev.off()
########################################################################################
##Gamma prior
pdf("plot/gamma.pdf",family="Times",height=4,width=6)  
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(2,3))

a<-1 ; b<-1
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",1," b=",1,sep=""))),side=3,line=.12,cex=.8)

a<-2 ; b<-2
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",2," b=",2,sep=""))),side=3,line=.12,cex=.8)

a<-4 ; b<-4
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",4," b=",4,sep=""))),side=3,line=.12,cex=.8)

a<-2 ; b<-1
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",2," b=",1,sep=""))),side=3,line=.12,cex=.8)

a<-8 ; b<-4
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",8," b=",4,sep=""))),side=3,line=.12,cex=.8)

a<-32 ; b<-16
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",32," b=",16,sep=""))),side=3,line=.12,cex=.8)

dev.off()
#############################################
##Birth example
pdf("plot/birth.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
a<-2
b<-1
s1<-217
n1<-111
s2<-66
  n2<-44
xtheta<-seq(0,5,length=1000)
plot(xtheta,dgamma(xtheta,a+s1,b+n1),type="l",col=gray(.5),xlab=expression(theta),
     ylab=expression(paste(italic("p("),theta,"|",y[1],"...",y[n],")",sep="")))
lines(xtheta,dgamma(xtheta,a+s2,b+n2),col=gray(0),lwd=2)
lines(xtheta,dgamma(xtheta,a,b),type="l",lty=2,lwd=2)
abline(h=0,col="black")
y<-(0:12)
plot(y-.1, dnbinom(y, size=(a+s1), mu=(a+s1)/(b+n1)) , col=gray(.5) ,type="h",
     ylab=expression(paste(italic("p("),y[n+1],"|",y[1],"...",y[n],")",sep="")), 
     xlab=expression(italic(y[n+1])),ylim=c(0,.35),lwd=3)
points(y+.1, dnbinom(y, size=(a+s2), mu=(a+s2)/(b+n2)) , col=gray(0) ,type="h",lwd=3)
legend(1,.375,legend=c("Less than bachelor's","Bachelor's or higher"),bty="n",
       lwd=c(3,3),col=c(gray(.5),gray(0)))
dev.off()
######