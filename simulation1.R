rm(list=ls(all=TRUE))
library(parallel)   #Library for PARALLEL Processing

simFunc <- function(N){
  B<- 500           #Number of Estimates.
  set.seed(2016)

  #Set the covariates.
  #Set the first parameter scenario.
  beta1 = 0.5 
  beta2 = 0.5 
  beta3 = 1.5
  beta4 =   2
  beta5 =  -3
  beta6 =   1
  beta7 =  -2
  beta8 = 0.75
  varp=c(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8)
  lime = 10
  
  # Declaring vectors and counters.
  emle  <-matrix(nrow=B,ncol=8)
  old  =  options(digits=10)
  o<-1
  ite<-1
  pc1<-rep(0,times=B)
  pc2<-rep(0,times=B)
  pc3<-rep(0,times=B)
  pc4<-rep(0,times=B)
  pc5<-rep(0,times=B)
  pc6<-rep(0,times=B)
  pc7<-rep(0,times=B)
  pc8<-rep(0,times=B)
  censura = matrix(numeric(B),B,1)
  
  model <- function(N,varp){
    a       =  exp(varp[1])*exp(y*varp[2])
    lam     =  exp(varp[3])*exp(y*varp[4])
    
    # Parameters - WEIGHT IN ZEROS AND INFS: GAMMA ZERO AND GAMMA ONES WITH A COVARIABLE
    gamm_zer = (exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))  
    gamm_inf = (exp(varp[7])*exp(y*varp[8]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))
    
    U=runif(N)
    V=runif(N,gamm_zer,1-gamm_inf)
    y=numeric(N)
    for(i in 1:N){
      if(U[i]<=gamm_zer[i]){ 
        y[i] = 0
      }else 
        if(U[i]>gamm_zer[i] & U[i] <= 1-gamm_inf[i]){
          #y[i]=qweibull(V[i], a, lam)
          y[i] = lam[i]*((-log(1 - (V[i]-gamm_zer[i])/(1-gamm_zer[i]-gamm_inf[i])))^(1/a[i]))
        }else {      
          y[i]=+Inf
        }
    }
    return(y)
  }
  
  

  GPE=function(varp)
  {
    a       =   exp(varp[1])*exp(y*varp[2])
    lam     =   exp(varp[3])*exp(y*varp[4])
    gamm_zer = (exp(varp[5])*exp(y*varp[6]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))  
    gamm_inf = (exp(varp[7])*exp(y*varp[8]))/(1+exp(varp[5])*exp(y*varp[6])+exp(varp[7])*exp(y*varp[8]))
    f0   = gamm_zer
    fw = (exp(log(a)-log(lam)))*((exp(log(time)-log(lam)))^(a-1))*(exp(-(exp(log(time)-log(lam)))^a))
    Sw = (exp(-(exp(log(time)-log(lam)))^a))
    fpop = (1-gamm_zer-gamm_inf)*fw
    Spop = gamm_inf + (1-gamm_zer-gamm_inf)*Sw
    f    = (fpop^delta)*(Spop^(1-delta))
    g   <- ifelse(time==0, f0, f)
    adFunc   = sum(log(g))
    return(adFunc)
  }
  
  
  while(o<=B){
    y=rbinom(N,1,0.5)
    tp = model(N,varp)
    tM = max(tp[tp<+Inf])
    Z=runif(N,0,tM)
    time=numeric(N)
    for (i in 1:N)
    {time[i] = min(tp[i],Z[i])
    }
    delta=numeric(N)
    for (i in 1:N) {
      if (time[i] < Z[i]){
        delta[i]=1
      }
    }
    
    t<-time
    d<-delta
    x<-y 
    
    fit=try(optim(varp,GPE,method="BFGS",hessian=TRUE,control=list(fnscale=-1)))
    estc = try(fit$par)
    Hc=try(fit$hessian); varic=try(-solve(Hc))
    
    if ( is.finite(estc[1]) & is.finite(estc[2]) & is.finite(estc[3]) & is.finite(estc[4]) & 
         is.finite(estc[5]) & is.finite(estc[6]) & is.finite(estc[7]) & is.finite(estc[8])) {
      if ( is.finite(varic[1,1]) & is.finite(varic[2,2]) & is.finite(varic[3,3]) & is.finite(varic[4,4]) & 
           is.finite(varic[5,5]) & is.finite(varic[6,6]) & is.finite(varic[7,7]) & is.finite(varic[8,8])){
        if ( varic[1,1] > 0 & varic[2,2] > 0 & varic[3,3] > 0 & varic[4,4] > 0 & 
             varic[5,5] > 0 & varic[6,6] > 0 & varic[7,7] > 0 & varic[8,8] > 0 &
             estc[1] < lime & estc[2] < lime & estc[3] < lime & estc[4] < lime &
             estc[5] < lime & estc[6] < lime & estc[7] < lime & estc[8] < lime ){
         
            var.varp1 =varic[1,1]; Lvarp1 =estc[1]-1.96*sqrt(var.varp1); Uvarp1=estc[1]+1.96*sqrt(var.varp1)
            var.varp2 =varic[2,2]; Lvarp2 =estc[2]-1.96*sqrt(var.varp2); Uvarp2=estc[2]+1.96*sqrt(var.varp2)
            var.varp3 =varic[3,3]; Lvarp3 =estc[3]-1.96*sqrt(var.varp3); Uvarp3=estc[3]+1.96*sqrt(var.varp3)
            var.varp4 =varic[4,4]; Lvarp4 =estc[4]-1.96*sqrt(var.varp4); Uvarp4=estc[4]+1.96*sqrt(var.varp4)
            var.varp5 =varic[5,5]; Lvarp5 =estc[5]-1.96*sqrt(var.varp5); Uvarp5=estc[5]+1.96*sqrt(var.varp5)
            var.varp6 =varic[6,6]; Lvarp6 =estc[6]-1.96*sqrt(var.varp6); Uvarp6=estc[6]+1.96*sqrt(var.varp6)
            var.varp7 =varic[7,7]; Lvarp7 =estc[7]-1.96*sqrt(var.varp7); Uvarp7=estc[7]+1.96*sqrt(var.varp7)
            var.varp8 =varic[8,8]; Lvarp8 =estc[8]-1.96*sqrt(var.varp8); Uvarp8=estc[8]+1.96*sqrt(var.varp8)
            
            emle[o,1] <- estc[1]; emle[o,2] <- estc[2]; emle[o,3] <- estc[3]; emle[o,4] <- estc[4]; 
            emle[o,5] <- estc[5]; emle[o,6] <- estc[6]; emle[o,7] <- estc[7]; emle[o,8] <- estc[8]; 
                {
                  if (beta1 >= Lvarp1 && beta1 <= Uvarp1) pc1[o]<-1
                  if (beta2 >= Lvarp2 && beta2 <= Uvarp2) pc2[o]<-1
                  if (beta3 >= Lvarp3 && beta3 <= Uvarp3) pc3[o]<-1
                  if (beta4 >= Lvarp4 && beta4 <= Uvarp4) pc4[o]<-1
                  if (beta5 >= Lvarp5 && beta5 <= Uvarp5) pc5[o]<-1
                  if (beta6 >= Lvarp6 && beta6 <= Uvarp6) pc6[o]<-1
                  if (beta7 >= Lvarp7 && beta7 <= Uvarp7) pc7[o]<-1
                  if (beta8 >= Lvarp8 && beta8 <= Uvarp8) pc8[o]<-1
                  
                  censura[o,] = table(d)[1]/N 
                }
                cat(o,"     ",ite,"  ",round(sum(pc1)/o,3),"  ",round(sum(pc2)/o,3),"  ",round(sum(pc3)/o,3),"  ",round(sum(pc4)/o,3),"  ",round(sum(pc5)/o,3),"  ",round(sum(pc6)/o,3),"  ",round(sum(pc7)/o,3),"  ",round(sum(pc8)/o,3),"\n"); o<-(o+1);
              }
            
      }}
    ite<-ite+1
    }
  
  qtite  <-c(ite)
  s_cens = mean(censura[,1])
  
  ESTI1c <-c(mean(emle[,1]),mean(emle[,2]),mean(emle[,3]),mean(emle[,4]),mean(emle[,5]),mean(emle[,6]),mean(emle[,7]),mean(emle[,8]))
  BIAS1c <-c(mean(emle[,1])-beta1, mean(emle[,2])-beta2, mean(emle[,3])-beta3, mean(emle[,4])-beta4, mean(emle[,5])-beta5, mean(emle[,6])-beta6, mean(emle[,7])-beta7, mean(emle[,8])-beta8)
  RMSE1c <-c(sqrt(mean((emle[,1]-beta1)^2)),sqrt(mean((emle[,2]-beta2)^2)),sqrt(mean((emle[,3]-beta3)^2)),sqrt(mean((emle[,4]-beta4)^2)),sqrt(mean((emle[,5]-beta5)^2)),sqrt(mean((emle[,6]-beta6)^2)),sqrt(mean((emle[,7]-beta7)^2)),sqrt(mean((emle[,8]-beta8)^2)))
  PCc    <-c(sum(pc1)/B,sum(pc2)/B,sum(pc3)/B,sum(pc4)/B,sum(pc5)/B,sum(pc6)/B,sum(pc7)/B,sum(pc8)/B)
  
  return(list(qtite   =  qtite  ,
              s_cens  =  s_cens ,
              ESTI1c  =  ESTI1c ,
              BIAS1c  =  BIAS1c ,
              RMSE1c  =  RMSE1c ,
              PCc     =  PCc   ) )
}


numCores <- 5
cl <- makeCluster(numCores)
finalResults1 <- parLapply(cl, c(100,250,500,750,1000), simFunc)
stopCluster(cl)

#Plot

beta1 = 0.5 
beta2 = 0.5 
beta3 = 1.5
beta4 =   2
beta5 =  -3
beta6 =   1
beta7 =  -2
beta8 = 0.75
beta1c=c(beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8)

size = c(100,250,500,750,1000)

AV1_c = rbind( finalResults1[[1]]$ESTI1c,
               finalResults1[[2]]$ESTI1c,
               finalResults1[[3]]$ESTI1c,
               finalResults1[[4]]$ESTI1c,
               finalResults1[[5]]$ESTI1c)

RMSE1_c = rbind(finalResults1[[1]]$RMSE1c,
                finalResults1[[2]]$RMSE1c,
                finalResults1[[3]]$RMSE1c,
                finalResults1[[4]]$RMSE1c,
                finalResults1[[5]]$RMSE1c)

B1_c = rbind(    finalResults1[[1]]$BIAS1c,
                 finalResults1[[2]]$BIAS1c,
                 finalResults1[[3]]$BIAS1c,
                 finalResults1[[4]]$BIAS1c,
                 finalResults1[[5]]$BIAS1c)


cob1_c = rbind(finalResults1[[1]]$PCc,
               finalResults1[[2]]$PCc,
               finalResults1[[3]]$PCc,
               finalResults1[[4]]$PCc,
               finalResults1[[5]]$PCc)


par(mfrow=c(4,3))

## Param[5]
a5   =    0# min(RMSE1_c[,5]) - 0.25
b5   =       max(RMSE1_c[,5]) + 0.25
aa5   =      min(B1_c[,5]) - 0.25
bb5   =      max(B1_c[,5]) + 0.25

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,B1_c[,5],xlim=c(min(size),max(size)), xlab="",ylim=c(aa5,bb5),ylab = expression(paste("BIAS ( ", hat(beta)[10], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,RMSE1_c[,5],xlim=c(min(size),max(size)), xlab="",ylim=c(a5,b5),ylab = expression(paste("RMSE ( ", hat(beta)[10], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,cob1_c[,5],xlim=c(min(size),max(size)), xlab="",ylim=c(0.80,1),ylab=expression(paste("CP"," (",hat(beta)[10],")")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0.95, v = 0.95, col = "gray60")


## Param[6]
a6   =    0# min(RMSE1_c[,6]) - 0.25
b6   =       max(RMSE1_c[,6]) + 0.25
aa6   =      min(B1_c[,6]) - 0.25
bb6   =      max(B1_c[,6]) + 0.25

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,B1_c[,6],xlim=c(min(size),max(size)), xlab="",ylim=c(aa6,bb6),ylab = expression(paste("BIAS ( ", hat(beta)[11], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,RMSE1_c[,6],xlim=c(min(size),max(size)), xlab="",ylim=c(a6,b6),ylab = expression(paste("RMSE ( ", hat(beta)[11], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,cob1_c[,6],xlim=c(min(size),max(size)), xlab="",ylim=c(0.80,1),ylab=expression(paste("CP"," (",hat(beta)[11],")")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0.95, v = 0.95, col = "gray60")

## Param[7]
a7   =    0# min(RMSE1_c[,7])- 0.25
b7   =       max(RMSE1_c[,7])+ 0.25
aa7   =      min(B1_c[,7])- 0.25
bb7   =      max(B1_c[,7])+ 0.25


par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,B1_c[,7],xlim=c(min(size),max(size)), xlab="",ylim=c(aa7,bb7),ylab = expression(paste("BIAS ( ", hat(beta)[20], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,RMSE1_c[,7],xlim=c(min(size),max(size)), xlab="",ylim=c(a7,b7),ylab = expression(paste("RMSE ( ", hat(beta)[20], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,cob1_c[,7],xlim=c(min(size),max(size)), xlab="",ylim=c(0.80,1),ylab=expression(paste("CP"," (",hat(beta)[20],")")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0.95, v = 0.95, col = "gray60")

## Param[8]
a8   =    0# min(RMSE1_c[,8])- 0.25
b8   =       max(RMSE1_c[,8])+ 0.25
aa8   =      min(B1_c[,8])- 0.25
bb8   =      max(B1_c[,8])+ 0.25

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,B1_c[,8],xlim=c(min(size),max(size)), xlab="",ylim=c(aa8,bb8),ylab = expression(paste("BIAS ( ", hat(beta)[21], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,RMSE1_c[,8],xlim=c(min(size),max(size)), xlab="",ylim=c(a8,b8),ylab = expression(paste("RMSE ( ", hat(beta)[21], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,cob1_c[,8],xlim=c(min(size),max(size)), xlab="",ylim=c(0.80,1),ylab=expression(paste("CP"," (",hat(beta)[21],")")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0.95, v = 0.95, col = "gray60")
dev.off()

par(mfrow=c(4,3))
## Param[1]
a1   =     0#min(RMSE1_c[,1])- 0.25
b1   =       max(RMSE1_c[,1])+ 0.25
aa1   =      min(B1_c[,1])- 0.25
bb1   =      max(B1_c[,1])+ 0.25

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,B1_c[,1],xlim=c(min(size),max(size)), xlab="",ylim=c(aa1,bb1),ylab = expression(paste("BIAS ( ", hat(beta)[30], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,RMSE1_c[,1],xlim=c(min(size),max(size)), xlab="",ylim=c(a1,b1),ylab = expression(paste("RMSE ( ", hat(beta)[30], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,cob1_c[,1],xlim=c(min(size),max(size)), xlab="",ylim=c(0.80,1),ylab=expression(paste("CP"," (",hat(beta)[30],")")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0.95, v = 0.95, col = "gray60")

## Param[2]
a2   =    0# min(RMSE1_c[,2])- 0.25
b2   =       max(RMSE1_c[,2])+ 0.25
aa2   =      min(B1_c[,2])- 0.25
bb2   =      max(B1_c[,2])+ 0.25

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,B1_c[,2],xlim=c(min(size),max(size)), xlab="",ylim=c(aa2,bb2),ylab = expression(paste("BIAS ( ", hat(beta)[31], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,RMSE1_c[,2],xlim=c(min(size),max(size)), xlab="",ylim=c(a2,b2),ylab = expression(paste("RMSE ( ", hat(beta)[31], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,cob1_c[,2],xlim=c(min(size),max(size)), xlab="",ylim=c(0.80,1),ylab=expression(paste("CP"," (",hat(beta)[31],")")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0.95, v = 0.95, col = "gray60")

## Param[3]
a3   =    0# min(RMSE1_c[,3])- 0.25
b3   =       max(RMSE1_c[,3])+ 0.25
aa3   =      min(B1_c[,3])- 0.25
bb3   =      max(B1_c[,3])+ 0.25

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,B1_c[,3],xlim=c(min(size),max(size)), xlab="",ylim=c(aa3,bb3),ylab = expression(paste("BIAS ( ", hat(beta)[40], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,RMSE1_c[,3],xlim=c(min(size),max(size)), xlab="",ylim=c(a3,b3),ylab = expression(paste("RMSE ( ", hat(beta)[40], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,cob1_c[,3],xlim=c(min(size),max(size)), xlab="",ylim=c(0.80,1),ylab=expression(paste("CP"," (",hat(beta)[40],")")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0.95, v = 0.95, col = "gray60")

## Param[4]
a4   =    0# min(RMSE1_c[,4])- 0.25
b4   =       max(RMSE1_c[,4])+ 0.25
aa4   =      min(B1_c[,4])- 0.25
bb4   =      max(B1_c[,4])+ 0.25

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,B1_c[,4],xlim=c(min(size),max(size)), xlab="",ylim=c(aa4,bb4),ylab = expression(paste("BIAS ( ", hat(beta)[41], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,RMSE1_c[,4],xlim=c(min(size),max(size)), xlab="",ylim=c(a4,b4),ylab = expression(paste("RMSE ( ", hat(beta)[41], " )")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0, v = 1, col = "gray60")

par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,cob1_c[,4],xlim=c(min(size),max(size)), xlab="",ylim=c(0.80,1),ylab=expression(paste("CP"," (",hat(beta)[41],")")),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = 0.95, v = 0.95, col = "gray60")
dev.off()

par(mfrow=c(4,2))
inf5   =       min(AV1_c[,5]) - 0.25
sup5   =       max(AV1_c[,5]) + 0.25
par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,  AV1_c[,5],xlim=c(min(size),max(size)), ylim=c(inf5,sup5), xlab="",ylab = expression(paste(hat(beta)[10])),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = beta1c[5], v = 1, col = "gray60")

inf6   =       min(AV1_c[,6]) - 0.25
sup6   =       max(AV1_c[,6]) + 0.25
par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,  AV1_c[,6],xlim=c(min(size),max(size)), ylim=c(inf6,sup6), xlab="",ylab = expression(paste(hat(beta)[11])),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = beta1c[6], v = 1, col = "gray60")

inf7   =       min(AV1_c[,7]) - 0.25
sup7   =       max(AV1_c[,7]) + 0.25
par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,  AV1_c[,7],xlim=c(min(size),max(size)), ylim=c(inf7,sup7), xlab="",ylab = expression(paste(hat(beta)[20])),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = beta1c[7], v = 1, col = "gray60")

inf8   =       min(AV1_c[,8]) - 0.25
sup8   =       max(AV1_c[,8]) + 0.25
par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,  AV1_c[,8],xlim=c(min(size),max(size)), ylim=c(inf8,sup8), xlab="",ylab = expression(paste(hat(beta)[21])),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = beta1c[8], v = 1, col = "gray60")

inf1   =       min(AV1_c[,1]) - 0.25
sup1   =       max(AV1_c[,1]) + 0.25
par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,  AV1_c[,1],xlim=c(min(size),max(size)), ylim=c(inf1,sup1), xlab="",ylab = expression(paste(hat(beta)[30])),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = beta1c[1], v = 1, col = "gray60")

inf2   =       min(AV1_c[,2]) - 0.25
sup2   =       max(AV1_c[,2]) + 0.25
par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,  AV1_c[,2],xlim=c(min(size),max(size)), ylim=c(inf2,sup2), xlab="",ylab = expression(paste(hat(beta)[31])),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = beta1c[2], v = 1, col = "gray60")

inf3   =       min(AV1_c[,3]) - 0.25
sup3   =       max(AV1_c[,3]) + 0.25
par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,  AV1_c[,3],xlim=c(min(size),max(size)), ylim=c(inf3,sup3), xlab="",ylab = expression(paste(hat(beta)[40])),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = beta1c[3], v = 1, col = "gray60")

inf4   =       min(AV1_c[,4]) - 0.25
sup4   =       max(AV1_c[,4]) + 0.25
par(mai=c(0.45,0.9,0.1, 0.1)) 
plot(size,  AV1_c[,4],xlim=c(min(size),max(size)), ylim=c(inf4,sup4), xlab="",ylab = expression(paste(hat(beta)[41])),type="o", pch = 49,col="1",lwd=1,lty=2,main=NULL)
abline(h = beta1c[4], v = 1, col = "gray60")
dev.off()

save.image("model1.RData") 