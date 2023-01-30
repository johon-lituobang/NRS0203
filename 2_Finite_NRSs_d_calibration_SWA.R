
#require foreach and doparallel for parallel processing of bootstrap (not available for some types of computers)
if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)

#require randtoolbox for random number generations
if (!require("randtoolbox")) install.packages("randtoolbox")
library(randtoolbox)
if (!require("Rfast")) install.packages("Rfast")
library(Rfast)

numCores <- detectCores()
#registering clusters, can set a smaller number using numCores-1 

registerDoParallel(numCores)


#bootsize for bootstrap approximation of the distributions of the kernel of U-statistics.
n <- 1.8*10^4
(n%%10)==0
# maximum order of moments
morder <- 4

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                     mixed = FALSE, method = "C", start = 1)
quasiuni<-rbind(quasiunisobol)

quasiunisobol<-c()

quasiuni_sorted2 <- na.omit(rowSort(quasiuni[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3 <- na.omit(rowSort(quasiuni[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4 <- na.omit(rowSort(quasiuni, descend = FALSE, stable = FALSE, parallel = TRUE))
# Forever...

roundunique<-function(orderlist1,dimension,size){
  roundedlist1<-Round(orderlist1*(size-1),digit=0,na.rm = FALSE)+1
  find1<-unique(roundedlist1)
  if (length(find1)<dimension){
    return(rep(NA, times=dimension))
  }else{
    return(roundedlist1)
  }
}

#simulate the order list for later bootstrap and distribution simulation(SE) (this is the most important part, because the order list ensure the performance of deterministic simulation.)
removelist<-function(orderlist1){
  orderlist1<-orderlist1[!duplicated(orderlist1), ]
  remended1<-(length(orderlist1[,1])%%8)
  if (remended1==0){
    return(orderlist1)
  }else{
    orderlist1<-orderlist1[-c(1:remended1),]
    return(orderlist1)
  }
}

extract1<-function(orderlist,sortedx,dimension){
  if (dimension==2){
    ss1<-sortedx[orderlist[1]]
    ss2<-sortedx[orderlist[2]]
    return(c(ss1,ss2))
  }else if (dimension==3){
    ss1<-sortedx[orderlist[1]]
    ss2<-sortedx[orderlist[2]]
    ss3<-sortedx[orderlist[3]]
    return(c(ss1,ss2,ss3))
  }else if (dimension==4){
    ss1<-sortedx[orderlist[1]]
    ss2<-sortedx[orderlist[2]]
    ss3<-sortedx[orderlist[3]]
    ss4<-sortedx[orderlist[4]]
    return(c(ss1,ss2,ss3,ss4))
  }
}
factdivide<-function(n1,n2){
  decin1<-n1-floor(n1)
  if(decin1==0){decin1=1}
  decin2<-n2-floor(n2)
  if(decin2==0){decin2=1}
  n1seq<-seq(decin1,n1, by=1)
  n2seq<-seq(decin2,n2,by=1)
  all<-list(n1seq,n2seq)
  maxlen <- max(lengths(all))
  all2 <- as.data.frame(lapply(all, function(lst) c(lst, rep(1, maxlen - length(lst)))))
  division<-all2[,1]/all2[,2]
  answer<-exp(sum(log(division)))*(gamma(decin1)/gamma(decin2))
  return(answer)
}
unbiasedsd<-function (x){
  n<-length(x)
  if(n==1){
    return(1000000)
  }
  sd1<-sd(x)
  if(n==2){
    c1<-0.7978845608
  }else if (n==3){
    c1<-0.8862269255 
  }else if (n==4){
    c1<-0.9213177319 
  }else{
    c1<-sqrt(2/(n-1))*factdivide(n1=((n/2)-1),n2=(((n-1)/2)-1))
  }
  listall<-sd1*c1
  (listall)
}

#load the deterministic simulation functions of 9 common unimodal distributions
dsexp<-function (uni,scale=1) {
  sample1<-qexp(uni,rate=scale)
  sample1
}
dsRayleigh<-function (uni,scale=1) {
  sample1 <- scale * sqrt(-2 * log((uni)))
  sample1[scale <= 0] <- NaN
  rev(sample1)
}
dsnorm<-function (uni,location=0,scale=1) {
  sample1<-qnorm(uni,mean =location,sd=scale)
  sample1
}
dsLaplace<-function (uni,location=0,scale=1) {
  sample1<-location - sign(uni - 0.5) * scale * (log(2) + ifelse(uni < 0.5, log(uni), log1p(-uni)))
  sample1
}
dslogis<-function (uni,location=0,scale=1) {
  sample1<-qlogis(uni,location=location,scale=scale)
  sample1
}
dsPareto<-function (uni,shape,scale=1) {
  sample1 <- scale*((uni))^(-1/shape)
  sample1[scale <= 0] <- NaN
  sample1[shape <= 0] <- NaN
  rev(sample1)
}
dslnorm<-function (uni,location=0,scale) {
  sample1 <- qlnorm(uni,meanlog=location,sdlog = scale)
  sample1
}
dsgamma<-function (uni,shape,scale = 1) {
  sample1<-qgamma(uni,shape=shape,scale=scale)
  sample1
}
dsWeibull<-function (uni,shape, scale = 1){
  sample1<-qweibull(uni,shape=shape, scale = scale)
  sample1
}


moments<-function (x){
  n<-length(x)
  m1<-mean(x)
  var1<-(sum((x - m1)^2)/n)
  tm1<-(sum((x - m1)^3)/n)
  fm1<-(sum((x - m1)^4)/(n))
  listall<-c(mean=m1,variance=var1,tm=tm1,fm=fm1)
  (listall)
}
unbiasedmoments<-function (x){
  n<-length(x)
  m1<-mean(x)
  var1<-sd(x)^2
  var2<-(sum((x - m1)^2)/n)
  tm1<-(sum((x - m1)^3)/n)*(n^2/((n-1)*(n-2)))
  fm1<-(sum((x - m1)^4)/n)
  ufm1<--3*var2^2*(2*n-3)*n/((n-1)*(n-2)*(n-3))+(n^2-2*n+3)*fm1*n/((n-1)*(n-2)*(n-3))
  listall<-c(mean=m1,variance=var1,tm=tm1,fm=ufm1)
  (listall)
}
se_mean<-function (x){
  n<-length(x)
  usd<-unbiasedsd(x)
  usd/sqrt(n)
}
se_sd<-function (x){
  n<-length(x)
  m1<-mean(x)
  var1<-sd(x)^2
  var2<-(sum((x - m1)^2)/n)
  fm1<-(sum((x - m1)^4)/n)
  ufm1<--3*var2^2*(2*n-3)*n/((n-1)*(n-2)*(n-3))+(n^2-2*n+3)*fm1*n/((n-1)*(n-2)*(n-3))
  sqrt((ufm1/(4*n*var1))-((n-3)/(4*n*(n-1)))*var1)
}

standardizedmoments<-function (x){
  n<-length(x)
  m1<-mean(x)
  var1<-sd(x)^2
  var2<-(sum((x - m1)^2)/n)
  tm1<-(sum((x - m1)^3)/n)*(n^2/((n-1)*(n-2)))
  fm1<-(sum((x - m1)^4)/n)
  ufm1<--3*var2^2*(2*n-3)*n/((n-1)*(n-2)*(n-3))+(n^2-2*n+3)*fm1*n/((n-1)*(n-2)*(n-3))
  listall<-c(mean=m1,variance=var1,skewness=tm1/((var1)^(3/2)),kurtosis=ufm1/((var1)^(2)))
  (listall)
}
skewness<-function (x){
  n<-length(x)
  m1<-mean(x)
  sd1<-sd(x)
  tm1<-(sum((x - m1)^3)/n)*(n^2/((n-1)*(n-2)))
  listall<-c(skewness=tm1/(sd1)^3)
  (listall)
}
kurtosis<-function (x){
  n<-length(x)
  m1<-mean(x)
  sd1<-sd(x)
  fm1<-(sum((x - m1)^4)/n)
  listall<-c(kurtosis=fm1/sd1^4)
  (listall)
}
greatest_common_divisor<- function(a, b) {
  if (b == 0) a else Recall(b, a %% b)
}
least_common_multiple<-function(a,b){
  g<-greatest_common_divisor(a, b)
  return(a/g * b)
}
data_augmentation<-function (x,targetsize){
  lengthx<-length(x)
  meanx<-mean(x)
  disn<-least_common_multiple(lengthx,targetsize)
  orderedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  disori<-rep(orderedx, each = disn/lengthx)
  group<-rep(1:targetsize, each =disn/targetsize)
  data_augmentationresult<-sapply(split(disori, group), mean)
  return(data_augmentationresult)
}
mediansorted<-function(sortedx,lengthx){
  if (lengthx%%2==0){
    return((sortedx[lengthx/2]+sortedx[(lengthx/2)+1])/2)
  }
  else{return((sortedx[(lengthx+1)/2]))}
}


SWA<-function (x,interval=8,fast=TRUE,batch="auto"){
  lengthx<-length(x)
  if (batch=="auto" ){
    batch<-ceiling(500000/lengthx)+1
  }
  if (interval!=8 ){
    return("interval must be 8 ")
  }
  Ksamples<-lengthx/interval
  IntKsamples<-ceiling(Ksamples)
  target1<-IntKsamples*interval
  Bmass<-function(x1,IntKsamples,target1,sorted=FALSE){
    if(sorted){
      x_ordered<-x1
    }else{
      x_ordered<-Sort(x1,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    }
    if ((target1/IntKsamples)==8){
      partiallist<-(c(x_ordered[(IntKsamples+1):(2*(IntKsamples))],x_ordered[(6*IntKsamples+1):(7*IntKsamples)]))
      alllist<-c(x_ordered[(3*(IntKsamples)+1):(5*(IntKsamples))])
      #comlist1<-(c(x_ordered[(3*(IntKsamples)+1):(4*(IntKsamples))],x_ordered[(4*(IntKsamples)+1):(5*IntKsamples)]))
      comlist<-(c(x_ordered[(2*(IntKsamples)+1):(3*(IntKsamples))],x_ordered[(5*(IntKsamples)+1):(6*IntKsamples)]))
      tlist<-(c(x_ordered[1:(1*(IntKsamples))],x_ordered[(7*IntKsamples+1):(target1)]))
      
      sumpartiallist<-sum(partiallist)
      sumalllist<-sum(alllist)
      #sumcomlist1<-sum(comlist1)
      sumcomlist<-sum(comlist)
      sumtlist<-sum(tlist)
      
      median1<-mediansorted(sortedx=x_ordered,lengthx=target1)
      
      tm3<-(sumalllist)/(IntKsamples*(2))
      
      tm2<-(sumalllist+sumcomlist)/(IntKsamples*(4))
      
      tm1<-(sumpartiallist+sumalllist+sumcomlist)/(IntKsamples*(6))
      
      winsor<-(c(x_ordered[(IntKsamples+1)],x_ordered[(((target1-IntKsamples)))]))
      winsor1<-sum(winsor)*IntKsamples
      wm1<-(tm1*IntKsamples*(6)+winsor1)/(IntKsamples*(8))
      
      BM1<-(sumpartiallist*4-2*sumcomlist+2*sumalllist)/(IntKsamples*(8))
      
      mean1<-(sumpartiallist+sumalllist+sumcomlist+sumtlist)/(IntKsamples*(8))
      
      sq8<-(c(x_ordered[(IntKsamples):(IntKsamples+1)],x_ordered[(3*IntKsamples):(3*IntKsamples+1)],x_ordered[(5*IntKsamples):(5*IntKsamples+1)],x_ordered[(7*IntKsamples):(7*IntKsamples+1)]))
      
      sqm1<-sum(sq8)/8
      
      wm2<-(tm1*IntKsamples*(6)+sumpartiallist)/(IntKsamples*(8))
      
      Groupmean<-c(mean1=mean1,BM1=BM1,sqm1=sqm1,wm1=wm1,wm2=wm2,tm1=tm1,tm2=tm2,tm3=tm3,median1=median1)
      
      return(Groupmean)
    }
    else{
      return("Not supported yet.")
    }
  }
  if (Ksamples%%1!=0 ){
    if (fast==TRUE & lengthx<10000){
      x_ordered<-data_augmentation(x,target1)
      Groupmean<-Bmass(x_ordered,IntKsamples,target1,sorted=TRUE)
      return(Groupmean)
    }
    else{
      addedx<-matrix(sample(x,size=(target1-lengthx)*batch,replace=TRUE),nrow=batch)
      xt<-t(as.data.frame(x))
      xmatrix<-as.data.frame(lapply(xt, rep, batch))
      allmatrix<-cbind(addedx,xmatrix)
      batchresults<-apply(allmatrix,1,Bmass,IntKsamples=IntKsamples,target1=target1,sorted=FALSE)
      return(mean((batchresults)))}
  }
  else{
    Groupmean<-Bmass(x,IntKsamples,target1,sorted=FALSE)
    return(Groupmean)}
} 



CDF<-function(x,xevaluated,sorted=FALSE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  lengthx<-length(x)
  quantile2<-min(which(sortedx>(xevaluated)))
  if(is.infinite(quantile2)){
    quantile2=length(sortedx)
  }
  sortedquantile2<-sortedx[quantile2]
  sortedquantile21<-sortedx[quantile2-1]
  if(quantile2==1){
    sortedquantile21<-sortedx[quantile2]
  }
  if(sortedquantile21==sortedquantile2){
    result<-(quantile2-1)/lengthx
  }else{
    result<-((xevaluated-sortedquantile21)/(sortedquantile2-sortedquantile21)+quantile2-1)/lengthx
  }
  if(result>1){
    return(1)
  }else if (result<0){
    return(0)
  }else{
    return(result)
  }
}

finddmmm<-function(expectanalytic,expecttrue,x,SWA,median,isboot=TRUE,interval=8,fast=TRUE,batch=1000,sorted=FALSE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  lengthx<-length(x)
  
  quatileexpectanalytic<-CDF(x=sortedx,xevaluated=expectanalytic,sorted=TRUE)
  quatileexpecttrue<-CDF(x=sortedx,xevaluated=expecttrue,sorted=TRUE)
  if(isboot){
    expectboot=mean(x)
    quatileexpectboot<-CDF(x=sortedx,xevaluated=expectboot,sorted=TRUE)
  }else{
    expectboot=expecttrue
    quatileexpectboot<-quatileexpecttrue
  }
  mx1<-CDF(x=sortedx,xevaluated=SWA,sorted=TRUE)
  listd<-c(expectanalytic=expectanalytic,expecttrue=expecttrue,expectboot=expectboot,SWA=SWA,median=median,mx1=mx1,quatileexpectanalytic=quatileexpectanalytic,quatileexpecttrue=quatileexpecttrue,quatileexpectboot=quatileexpectboot)
  return(listd)
}

findd<-function(expectanalytic=NULL,expecttrue=NULL,expectboot=NULL,SWA1=NULL,median1=NULL,mx1=NULL,quatileexpectanalytic=NULL,
                quatileexpecttrue=NULL,quatileexpectboot=NULL){
  
  if (mx1>0.5){
    dqm1boot<-(quatileexpectboot-mx1)/(mx1-0.5)
    dqm1analytic<-(quatileexpectanalytic-mx1)/(mx1-0.5)
    dqm1true<-(quatileexpecttrue-mx1)/(mx1-0.5)
  }else if(mx1==0.5){
    dqm1boot<-0
    dqm1analytic<-0
    dqm1true<-0
  }else{
    quatileexpectboot<-1-quatileexpectboot
    quatileexpectanalytic<-1-quatileexpectanalytic
    quatileexpecttrue<-1-quatileexpecttrue
    mx1<-1-mx1
    dqm1boot<-(quatileexpectboot-mx1)/(mx1-0.5)
    dqm1analytic<-(quatileexpectanalytic-mx1)/(mx1-0.5)
    dqm1true<-(quatileexpecttrue-mx1)/(mx1-0.5)
  }
  
  drm1boot<-(expectboot-SWA1)/(SWA1-median1)
  drm1analytic<-(expectanalytic-SWA1)/(SWA1-median1)
  drm1true<-(expecttrue-SWA1)/(SWA1-median1)
  listd<-c(drm1analytic=drm1analytic,drm1true=drm1true,drm1boot=drm1boot,dqm1analytic=dqm1analytic,dqm1true=dqm1true,dqm1boot=dqm1boot)
  return(listd)
}


finddall<-function (x,targetm,targetvar,targettm,targetfm,orderlist1_sorted2=NULL,orderlist1_sorted3=NULL,orderlist1_sorted4=NULL,interval=8,fast=TRUE,batch="auto",boot=TRUE){
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  lengthx<-length(sortedx)
  
  samplemoments<-unbiasedmoments(sortedx)
  SWAmmm<-SWA(x=sortedx,interval=8,fast=TRUE,batch="auto")
  mmm1BM<-finddmmm(expectanalytic=targetm,expecttrue=samplemoments[1],x=sortedx,SWA=SWAmmm[2],median=SWAmmm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
  
  finddmmm1BM<-findd(expectanalytic=mmm1BM[1],expecttrue=mmm1BM[2],expectboot=mmm1BM[3],SWA1=mmm1BM[4],median1=mmm1BM[5],mx1=mmm1BM[6],quatileexpectanalytic=mmm1BM[7],
                   quatileexpecttrue=mmm1BM[8],quatileexpectboot=mmm1BM[9])
  
  mmm1sqm<-finddmmm(expectanalytic=targetm,expecttrue=samplemoments[1],x=sortedx,SWA=SWAmmm[3],median=SWAmmm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
  
  finddmmm1sqm<-findd(expectanalytic=mmm1sqm[1],expecttrue=mmm1sqm[2],expectboot=mmm1sqm[3],SWA1=mmm1sqm[4],median1=mmm1sqm[5],mx1=mmm1sqm[6],quatileexpectanalytic=mmm1sqm[7],
                     quatileexpecttrue=mmm1sqm[8],quatileexpectboot=mmm1sqm[9])
  
  mmm1wm1<-finddmmm(expectanalytic=targetm,expecttrue=samplemoments[1],x=sortedx,SWA=SWAmmm[4],median=SWAmmm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
  
  finddmmm1wm1<-findd(expectanalytic=mmm1wm1[1],expecttrue=mmm1wm1[2],expectboot=mmm1wm1[3],SWA1=mmm1wm1[4],median1=mmm1wm1[5],mx1=mmm1wm1[6],quatileexpectanalytic=mmm1wm1[7],
                      quatileexpecttrue=mmm1wm1[8],quatileexpectboot=mmm1wm1[9])
  
  mmm1wm2<-finddmmm(expectanalytic=targetm,expecttrue=samplemoments[1],x=sortedx,SWA=SWAmmm[5],median=SWAmmm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
  
  finddmmm1wm2<-findd(expectanalytic=mmm1wm2[1],expecttrue=mmm1wm2[2],expectboot=mmm1wm2[3],SWA1=mmm1wm2[4],median1=mmm1wm2[5],mx1=mmm1wm2[6],quatileexpectanalytic=mmm1wm2[7],
                      quatileexpecttrue=mmm1wm2[8],quatileexpectboot=mmm1wm2[9])
  
  mmm1tm1<-finddmmm(expectanalytic=targetm,expecttrue=samplemoments[1],x=sortedx,SWA=SWAmmm[6],median=SWAmmm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
  
  finddmmm1tm1<-findd(expectanalytic=mmm1tm1[1],expecttrue=mmm1tm1[2],expectboot=mmm1tm1[3],SWA1=mmm1tm1[4],median1=mmm1tm1[5],mx1=mmm1tm1[6],quatileexpectanalytic=mmm1tm1[7],
                      quatileexpecttrue=mmm1tm1[8],quatileexpectboot=mmm1tm1[9])
  
  mmm1tm2<-finddmmm(expectanalytic=targetm,expecttrue=samplemoments[1],x=sortedx,SWA=SWAmmm[7],median=SWAmmm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
  
  finddmmm1tm2<-findd(expectanalytic=mmm1tm2[1],expecttrue=mmm1tm2[2],expectboot=mmm1tm2[3],SWA1=mmm1tm2[4],median1=mmm1tm2[5],mx1=mmm1tm2[6],quatileexpectanalytic=mmm1tm2[7],
                      quatileexpecttrue=mmm1tm2[8],quatileexpectboot=mmm1tm2[9])
  
  mmm1tm3<-finddmmm(expectanalytic=targetm,expecttrue=samplemoments[1],x=sortedx,SWA=SWAmmm[8],median=SWAmmm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
  
  finddmmm1tm3<-findd(expectanalytic=mmm1tm3[1],expecttrue=mmm1tm3[2],expectboot=mmm1tm3[3],SWA1=mmm1tm3[4],median1=mmm1tm3[5],mx1=mmm1tm3[6],quatileexpectanalytic=mmm1tm3[7],
                      quatileexpecttrue=mmm1tm3[8],quatileexpectboot=mmm1tm3[9])
  
  if (boot){
    
    bootstrappedsample2<-t(as.data.frame(apply(orderlist1_sorted2,MARGIN=1,FUN=extract1,sortedx=sortedx,dimension=2)))
    
    getvar<-function(vector){ 
      ((vector[1]-vector[2])^2)/2
    }
    
    dp2varx<-apply(bootstrappedsample2,MARGIN=1,FUN=getvar)
    dp2varx<-Sort(x=dp2varx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    bootstrappedsample2<-c()
    SWAvar<-SWA(x=dp2varx,interval=8,fast=TRUE,batch="auto")
    
    varmoBM<-finddmmm(expectanalytic=targetvar,expecttrue=samplemoments[2],x=dp2varx,SWA=SWAvar[2],median=SWAvar[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddvarmoBM<-findd(expectanalytic=varmoBM[1],expecttrue=varmoBM[2],expectboot=varmoBM[3],SWA1=varmoBM[4],median1=varmoBM[5],mx1=varmoBM[6],quatileexpectanalytic=varmoBM[7],
                       quatileexpecttrue=varmoBM[8],quatileexpectboot=varmoBM[9])
    
    varmosqm<-finddmmm(expectanalytic=targetvar,expecttrue=samplemoments[2],x=dp2varx,SWA=SWAvar[3],median=SWAvar[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddvarmosqm<-findd(expectanalytic=varmosqm[1],expecttrue=varmosqm[2],expectboot=varmosqm[3],SWA1=varmosqm[4],median1=varmosqm[5],mx1=varmosqm[6],quatileexpectanalytic=varmosqm[7],
                        quatileexpecttrue=varmosqm[8],quatileexpectboot=varmosqm[9])
    
    varmowm1<-finddmmm(expectanalytic=targetvar,expecttrue=samplemoments[2],x=dp2varx,SWA=SWAvar[4],median=SWAvar[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddvarmowm1<-findd(expectanalytic=varmowm1[1],expecttrue=varmowm1[2],expectboot=varmowm1[3],SWA1=varmowm1[4],median1=varmowm1[5],mx1=varmowm1[6],quatileexpectanalytic=varmowm1[7],
                        quatileexpecttrue=varmowm1[8],quatileexpectboot=varmowm1[9])
    
    varmowm2<-finddmmm(expectanalytic=targetvar,expecttrue=samplemoments[2],x=dp2varx,SWA=SWAvar[5],median=SWAvar[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddvarmowm2<-findd(expectanalytic=varmowm2[1],expecttrue=varmowm2[2],expectboot=varmowm2[3],SWA1=varmowm2[4],median1=varmowm2[5],mx1=varmowm2[6],quatileexpectanalytic=varmowm2[7],
                        quatileexpecttrue=varmowm2[8],quatileexpectboot=varmowm2[9])
    
    varmotm1<-finddmmm(expectanalytic=targetvar,expecttrue=samplemoments[2],x=dp2varx,SWA=SWAvar[6],median=SWAvar[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddvarmotm1<-findd(expectanalytic=varmotm1[1],expecttrue=varmotm1[2],expectboot=varmotm1[3],SWA1=varmotm1[4],median1=varmotm1[5],mx1=varmotm1[6],quatileexpectanalytic=varmotm1[7],
                        quatileexpecttrue=varmotm1[8],quatileexpectboot=varmotm1[9])
    
    varmotm2<-finddmmm(expectanalytic=targetvar,expecttrue=samplemoments[2],x=dp2varx,SWA=SWAvar[7],median=SWAvar[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddvarmotm2<-findd(expectanalytic=varmotm2[1],expecttrue=varmotm2[2],expectboot=varmotm2[3],SWA1=varmotm2[4],median1=varmotm2[5],mx1=varmotm2[6],quatileexpectanalytic=varmotm2[7],
                        quatileexpecttrue=varmotm2[8],quatileexpectboot=varmotm2[9])
    
    varmotm3<-finddmmm(expectanalytic=targetvar,expecttrue=samplemoments[2],x=dp2varx,SWA=SWAvar[8],median=SWAvar[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddvarmotm3<-findd(expectanalytic=varmotm3[1],expecttrue=varmotm3[2],expectboot=varmotm3[3],SWA1=varmotm3[4],median1=varmotm3[5],mx1=varmotm3[6],quatileexpectanalytic=varmotm3[7],
                        quatileexpecttrue=varmotm3[8],quatileexpectboot=varmotm3[9])
    
    
    dp2varx<-c()
    
    bootstrappedsample3<-t(as.data.frame(apply(orderlist1_sorted3,MARGIN=1,FUN=extract1,sortedx=sortedx,dimension=3)))
    
    gettm<-function(vector){ 
      ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
    }
    
    dp3tmx<-apply(bootstrappedsample3,MARGIN=1,FUN=gettm)
    dp3tmx<-Sort(x=dp3tmx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    bootstrappedsample3<-c()
    SWAtm<-SWA(x=dp3tmx,interval=8,fast=TRUE,batch="auto")
    
    tmmoBM<-finddmmm(expectanalytic=targettm,expecttrue=samplemoments[3],x=dp3tmx,SWA=SWAtm[2],median=SWAtm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddtmmoBM<-findd(expectanalytic=tmmoBM[1],expecttrue=tmmoBM[2],expectboot=tmmoBM[3],SWA1=tmmoBM[4],median1=tmmoBM[5],mx1=tmmoBM[6],quatileexpectanalytic=tmmoBM[7],
                        quatileexpecttrue=tmmoBM[8],quatileexpectboot=tmmoBM[9])
    
    tmmosqm<-finddmmm(expectanalytic=targettm,expecttrue=samplemoments[3],x=dp3tmx,SWA=SWAtm[3],median=SWAtm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddtmmosqm<-findd(expectanalytic=tmmosqm[1],expecttrue=tmmosqm[2],expectboot=tmmosqm[3],SWA1=tmmosqm[4],median1=tmmosqm[5],mx1=tmmosqm[6],quatileexpectanalytic=tmmosqm[7],
                         quatileexpecttrue=tmmosqm[8],quatileexpectboot=tmmosqm[9])
    
    tmmowm1<-finddmmm(expectanalytic=targettm,expecttrue=samplemoments[3],x=dp3tmx,SWA=SWAtm[4],median=SWAtm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddtmmowm1<-findd(expectanalytic=tmmowm1[1],expecttrue=tmmowm1[2],expectboot=tmmowm1[3],SWA1=tmmowm1[4],median1=tmmowm1[5],mx1=tmmowm1[6],quatileexpectanalytic=tmmowm1[7],
                         quatileexpecttrue=tmmowm1[8],quatileexpectboot=tmmowm1[9])
    
    tmmowm2<-finddmmm(expectanalytic=targettm,expecttrue=samplemoments[3],x=dp3tmx,SWA=SWAtm[5],median=SWAtm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddtmmowm2<-findd(expectanalytic=tmmowm2[1],expecttrue=tmmowm2[2],expectboot=tmmowm2[3],SWA1=tmmowm2[4],median1=tmmowm2[5],mx1=tmmowm2[6],quatileexpectanalytic=tmmowm2[7],
                         quatileexpecttrue=tmmowm2[8],quatileexpectboot=tmmowm2[9])
    
    tmmotm1<-finddmmm(expectanalytic=targettm,expecttrue=samplemoments[3],x=dp3tmx,SWA=SWAtm[6],median=SWAtm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddtmmotm1<-findd(expectanalytic=tmmotm1[1],expecttrue=tmmotm1[2],expectboot=tmmotm1[3],SWA1=tmmotm1[4],median1=tmmotm1[5],mx1=tmmotm1[6],quatileexpectanalytic=tmmotm1[7],
                         quatileexpecttrue=tmmotm1[8],quatileexpectboot=tmmotm1[9])
    
    tmmotm2<-finddmmm(expectanalytic=targettm,expecttrue=samplemoments[3],x=dp3tmx,SWA=SWAtm[7],median=SWAtm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddtmmotm2<-findd(expectanalytic=tmmotm2[1],expecttrue=tmmotm2[2],expectboot=tmmotm2[3],SWA1=tmmotm2[4],median1=tmmotm2[5],mx1=tmmotm2[6],quatileexpectanalytic=tmmotm2[7],
                         quatileexpecttrue=tmmotm2[8],quatileexpectboot=tmmotm2[9])
    
    tmmotm3<-finddmmm(expectanalytic=targettm,expecttrue=samplemoments[3],x=dp3tmx,SWA=SWAtm[8],median=SWAtm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddtmmotm3<-findd(expectanalytic=tmmotm3[1],expecttrue=tmmotm3[2],expectboot=tmmotm3[3],SWA1=tmmotm3[4],median1=tmmotm3[5],mx1=tmmotm3[6],quatileexpectanalytic=tmmotm3[7],
                         quatileexpecttrue=tmmotm3[8],quatileexpectboot=tmmotm3[9])
    
    dp3tmx<-c()
    
    bootstrappedsample4<-t(as.data.frame(apply(orderlist1_sorted4,MARGIN=1,FUN=extract1,sortedx=sortedx,dimension=4)))
    
    getfm<-function(vector){ 
      resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                    4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                    vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                    6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                    vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                                 6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
      return(resd)
    }
    
    dp4fmx<-apply(bootstrappedsample4,MARGIN=1,FUN=getfm)
    dp4fmx<-Sort(x=dp4fmx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    bootstrappedsample4<-c()
    SWAfm<-SWA(x=dp4fmx,interval=8,fast=TRUE,batch="auto")
    
    fmmoBM<-finddmmm(expectanalytic=targetfm,expecttrue=samplemoments[4],x=dp4fmx,SWA=SWAfm[2],median=SWAfm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddfmmoBM<-findd(expectanalytic=fmmoBM[1],expecttrue=fmmoBM[2],expectboot=fmmoBM[3],SWA1=fmmoBM[4],median1=fmmoBM[5],mx1=fmmoBM[6],quatileexpectanalytic=fmmoBM[7],
                       quatileexpecttrue=fmmoBM[8],quatileexpectboot=fmmoBM[9])
    
    fmmosqm<-finddmmm(expectanalytic=targetfm,expecttrue=samplemoments[4],x=dp4fmx,SWA=SWAfm[3],median=SWAfm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddfmmosqm<-findd(expectanalytic=fmmosqm[1],expecttrue=fmmosqm[2],expectboot=fmmosqm[3],SWA1=fmmosqm[4],median1=fmmosqm[5],mx1=fmmosqm[6],quatileexpectanalytic=fmmosqm[7],
                        quatileexpecttrue=fmmosqm[8],quatileexpectboot=fmmosqm[9])
    
    fmmowm1<-finddmmm(expectanalytic=targetfm,expecttrue=samplemoments[4],x=dp4fmx,SWA=SWAfm[4],median=SWAfm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddfmmowm1<-findd(expectanalytic=fmmowm1[1],expecttrue=fmmowm1[2],expectboot=fmmowm1[3],SWA1=fmmowm1[4],median1=fmmowm1[5],mx1=fmmowm1[6],quatileexpectanalytic=fmmowm1[7],
                        quatileexpecttrue=fmmowm1[8],quatileexpectboot=fmmowm1[9])
    
    fmmowm2<-finddmmm(expectanalytic=targetfm,expecttrue=samplemoments[4],x=dp4fmx,SWA=SWAfm[5],median=SWAfm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddfmmowm2<-findd(expectanalytic=fmmowm2[1],expecttrue=fmmowm2[2],expectboot=fmmowm2[3],SWA1=fmmowm2[4],median1=fmmowm2[5],mx1=fmmowm2[6],quatileexpectanalytic=fmmowm2[7],
                        quatileexpecttrue=fmmowm2[8],quatileexpectboot=fmmowm2[9])
    
    fmmotm1<-finddmmm(expectanalytic=targetfm,expecttrue=samplemoments[4],x=dp4fmx,SWA=SWAfm[6],median=SWAfm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddfmmotm1<-findd(expectanalytic=fmmotm1[1],expecttrue=fmmotm1[2],expectboot=fmmotm1[3],SWA1=fmmotm1[4],median1=fmmotm1[5],mx1=fmmotm1[6],quatileexpectanalytic=fmmotm1[7],
                        quatileexpecttrue=fmmotm1[8],quatileexpectboot=fmmotm1[9])
    
    fmmotm2<-finddmmm(expectanalytic=targetfm,expecttrue=samplemoments[4],x=dp4fmx,SWA=SWAfm[7],median=SWAfm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddfmmotm2<-findd(expectanalytic=fmmotm2[1],expecttrue=fmmotm2[2],expectboot=fmmotm2[3],SWA1=fmmotm2[4],median1=fmmotm2[5],mx1=fmmotm2[6],quatileexpectanalytic=fmmotm2[7],
                        quatileexpecttrue=fmmotm2[8],quatileexpectboot=fmmotm2[9])
    
    fmmotm3<-finddmmm(expectanalytic=targetfm,expecttrue=samplemoments[4],x=dp4fmx,SWA=SWAfm[8],median=SWAfm[9],isboot=FALSE,interval=interval,fast=fast,batch=batch,sorted=TRUE)
    
    finddfmmotm3<-findd(expectanalytic=fmmotm3[1],expecttrue=fmmotm3[2],expectboot=fmmotm3[3],SWA1=fmmotm3[4],median1=fmmotm3[5],mx1=fmmotm3[6],quatileexpectanalytic=fmmotm3[7],
                        quatileexpecttrue=fmmotm3[8],quatileexpectboot=fmmotm3[9])
    
    dp4fmx<-c()
    
  }else{
    if (lengthn>5000){
      print("Warning: The computational complexity is (e*n/2)^2, bootstrap is recommended")
    }
    subtract2x<-sapply(sortedx, "-", sortedx)
    subtract2x[lower.tri(subtract2x)] <- NA
    diag(subtract2x)=NA
    subtract2x<-na.omit(as.vector(subtract2x))
    dp22x<-subtract2x[subtract2x>0]
    subtract2x<-c()
    dp2varx<-(dp22x^2)/2
    dp22x<-c()
    
    varmoraw<-finddmmm(expectanalytic=targetvar,expecttrue=samplemoments[2],x=dp2varx,isboot=TRUE,interval=interval,fast=fast,batch=batch,sorted=FALSE)
    
    finddvar<-findd(expectanalytic=varmoraw[1],expecttrue=varmoraw[2],expectboot=varmoraw[3],SWA1=varmoraw[4],median1=varmoraw[5],mx1=varmoraw[6],quatileexpectanalytic=varmoraw[7],
                    quatileexpecttrue=varmoraw[8],quatileexpectboot=varmoraw[9])
    
    dp2varx<-c()
    if (lengthn>300){
      print("Warning: The computational complexity is n^3, bootstrap is recommended.")
    }
    subtract3x<-t(combn(sortedx, 3))
    
    gettm<-function(vector){ 
      ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
    }
    
    dp3tmx<-apply(subtract3x,MARGIN=1,FUN=gettm)
    
    subtract3x<-c()
    
    tmmoraw<-finddmmm(expectanalytic=targettm,expecttrue=samplemoments[3],x=dp3tmx,isboot=TRUE,interval=interval,fast=fast,batch=batch,sorted=FALSE)
    
    finddtm<-findd(expectanalytic=tmmoraw[1],expecttrue=tmmoraw[2],expectboot=tmmoraw[3],SWA1=tmmoraw[4],median1=tmmoraw[5],mx1=tmmoraw[6],quatileexpectanalytic=tmmoraw[7],
                   quatileexpecttrue=tmmoraw[8],quatileexpectboot=tmmoraw[9])
    
    dp3tmx<-c()
    if (lengthn>300){
      print("Warning: The computational complexity is n^4, bootstrap is recommended.")
    }
    subtract4x<-t(combn(sortedx, 4))
    
    getfm<-function(vector){ 
      resd<-1/12*(3*vector[1]^4 + 3*vector[2]^4 + 3*vector[3]^4 + 6*(vector[2]^2)*vector[3]*vector[4] - 4*(vector[3]^3)*vector[4] - 
                    4*vector[3]*(vector[4]^3) + 3*(vector[4]^4) - 4*(vector[2]^3)*(vector[3] + vector[4]) - 4*(vector[1]^3)*(vector[2]+vector[3]+vector[4])+ 
                    vector[2]*(-4*(vector[3]^3)+6*(vector[3]^2)*vector[4]+6*(vector[3])*(vector[4]^2) - 4*(vector[4]^3)) + 
                    6*(vector[1]^2)*(vector[3]*vector[4] + vector[2]*(vector[3] + vector[4])) + 
                    vector[1]*(-4*(vector[2]^3) - 4*(vector[3]^3) + 6*(vector[3]^2)*vector[4] + 6*vector[3]*(vector[4]^2) - 4*(vector[4]^3) + 
                                 6*(vector[2]^2)*(vector[3] + vector[4]) + 6*vector[2]*((vector[3]^2) - 6*vector[3]*vector[4] + vector[4]^2)))
      return(resd)
    }
    
    dp4fmx<-apply(subtract4x,MARGIN=1,FUN=getfm)
    
    subtract4x<-c()
    fmmoraw<-finddmmm(expectanalytic=targetfm,expecttrue=samplemoments[4],x=dp4fmx,isboot=TRUE,interval=interval,fast=fast,batch=batch,sorted=FALSE)
    
    finddfm<-findd(expectanalytic=fmmoraw[1],expecttrue=fmmoraw[2],expectboot=fmmoraw[3],SWA1=fmmoraw[4],median1=fmmoraw[5],mx1=fmmoraw[6],quatileexpectanalytic=fmmoraw[7],
                   quatileexpecttrue=fmmoraw[8],quatileexpectboot=fmmoraw[9])
    
    dp4fmx<-c()
  }
  
  all<-c(mmm1BM=mmm1BM,mmm1sqm=mmm1sqm,mmm1wm1=mmm1wm1,mmm1wm2=mmm1wm2,mmm1tm1=mmm1tm1,mmm1tm2=mmm1tm2,mmm1tm3=mmm1tm3,
         varmoBM=varmoBM,varmosqm=varmosqm,varmowm1=varmowm1,varmowm2=varmowm2,varmotm1=varmotm1,varmotm2=varmotm2,varmotm3=varmotm3,
         tmmoBM=tmmoBM,tmmosqm=tmmosqm,tmmowm1=tmmowm1,tmmowm2=tmmowm2,tmmotm1=tmmotm1,tmmotm2=tmmotm2,tmmotm3=tmmotm3,
         fmmoBM=fmmoBM,fmmosqm=fmmosqm,fmmowm1=fmmowm1,fmmowm2=fmmowm2,fmmotm1=fmmotm1,fmmotm2=fmmotm2,fmmotm3=fmmotm3
  )
  findall<-c(finddmmm1BM=finddmmm1BM,finddmmm1sqm=finddmmm1sqm,finddmmm1wm1=finddmmm1wm1,finddmmm1wm2=finddmmm1wm2,finddmmm1tm1=finddmmm1tm1,finddmmm1tm2=finddmmm1tm2,finddmmm1tm3=finddmmm1tm3,
             finddvarmoBM=finddvarmoBM,finddvarmosqm=finddvarmosqm,finddvarmowm1=finddvarmowm1,finddvarmowm2=finddvarmowm2,finddvarmotm1=finddvarmotm1,finddvarmotm2=finddvarmotm2,finddvarmotm3=finddvarmotm3,
             finddtmmoBM=finddtmmoBM,finddtmmosqm=finddtmmosqm,finddtmmowm1=finddtmmowm1,finddtmmowm2=finddtmmowm2,finddtmmotm1=finddtmmotm1,finddtmmotm2=finddtmmotm2,finddtmmotm3=finddtmmotm3,
             finddfmmoBM=finddfmmoBM,finddfmmosqm=finddfmmosqm,finddfmmowm1=finddfmmowm1,finddfmmowm2=finddfmmowm2,finddfmmotm1=finddfmmotm1,finddfmmotm2=finddfmmotm2,finddfmmotm3=finddfmmotm3)
  
  allcombine<-c(findall,all)
  return(allcombine)
  
}

findderror<-function(expectanalytic_error=NULL,expecttrue_error=NULL,expectboot_error=NULL,SWA1_error=NULL,median1_error=NULL,mx1_error=NULL,quatileexpectanalytic_error=NULL,
                     quatileexpecttrue_error=NULL,quatileexpectboot_error=NULL,expectanalytic=NULL,expecttrue=NULL,expectboot=NULL,SWA1=NULL,median1=NULL,mx1=NULL,quatileexpectanalytic=NULL,
                     quatileexpecttrue=NULL,quatileexpectboot=NULL,CovRobust=NULL,CovQuantile=NULL){
  mx2_error=0
  mx2=1/2
  
  dqm1boot<-(abs(1/(mx1-mx2))^2)*(quatileexpectboot_error^2)+(mx1_error^2)*((abs(-((quatileexpectboot-mx1)/((mx1-mx2)^2))-(1/(mx1-mx2))))^2)
  dqm1analytic<-(abs(1/(mx1-mx2))^2)*(quatileexpectanalytic_error^2)+(mx1_error^2)*((abs(-((quatileexpectanalytic-mx1)/((mx1-mx2)^2))-(1/(mx1-mx2))))^2)
  dqm1true<-(abs(1/(mx1-mx2))^2)*(quatileexpecttrue_error^2)+(mx1_error^2)*((abs(-((quatileexpecttrue-mx1)/((mx1-mx2)^2))-(1/(mx1-mx2))))^2)
  
  drm1boot<-(abs(-(expectboot-SWA1)/((SWA1-median1)^2)-1/(SWA1-median1))^2)*(SWA1_error^2)+((abs((expectboot-SWA1)/((SWA1-median1)^2)))^2)*(median1_error^2)
  drm1analytic<-(abs(-(expectanalytic-SWA1)/((SWA1-median1)^2)-1/(SWA1-median1))^2)*(SWA1_error^2)+((abs((expectanalytic-SWA1)/((SWA1-median1)^2)))^2)*(median1_error^2)
  drm1true<-(abs(-(expecttrue-SWA1)/((SWA1-median1)^2)-1/(SWA1-median1))^2)*(SWA1_error^2)+((abs((expecttrue-SWA1)/((SWA1-median1)^2)))^2)*(median1_error^2)
  
  listd<-c(drm1analytic=sqrt(drm1analytic),drm1true=sqrt(drm1true),drm1boot=sqrt(drm1boot),dqm1analytic=sqrt(dqm1analytic),dqm1true=sqrt(dqm1true),dqm1boot=sqrt(dqm1boot))
  return(listd)
}

kurtWeibull<- read.csv(("kurtWeibull_28260.csv"))

allkurtWeibull<-unlist(kurtWeibull)

samplesize=5400

orderlist1_AB2<-removelist(na.omit(t(apply(quasiuni_sorted2,MARGIN=1,FUN=roundunique,dimension=2,size=samplesize))))
orderlist1_AB3<-removelist(na.omit(t(apply(quasiuni_sorted3,MARGIN=1,FUN=roundunique,dimension=3,size=samplesize))))
orderlist1_AB4<-removelist(na.omit(t(apply(quasiuni_sorted4,MARGIN=1,FUN=roundunique,dimension=4,size=samplesize))))

batchsizebase=2000

batchsize=(batchsizebase)

n <- samplesize

setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)

#Then, start the Monte Simulation
simulatedbatchWeibull_bias_Monte<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
  library(Rfast)
  if (!require("foreach")) install.packages("foreach")
  library(foreach)
  if (!require("doParallel")) install.packages("doParallel")
  library(doParallel)
  #registering clusters, can set a smaller number using numCores-1 
  
  #require randtoolbox for random number generations
  if (!require("randtoolbox")) install.packages("randtoolbox")
  library(randtoolbox)
  #require Rfast for faster computation
  if (!require("Rfast")) install.packages("Rfast")
  library(Rfast)
  if (!require("gnorm")) install.packages("gnorm")
  library(gnorm)
  
  a=allkurtWeibull[batchnumber]
  
  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
      
    x<-c(dsWeibull(uni=unibatch[,batch1], shape=a/1, scale = 1))
    
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    x<-c()
    dall<-finddall(x=sortedx,targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,fast=TRUE,batch="auto")
    sortedx<-c()
    all1<-t(c(samplesize,kurtx,dall))
    SEbataches<-rbind(SEbataches,all1)
  }
    
  write.csv(SEbataches,paste("Weibull_raw_d_calibration",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
    
  SEbataches
}
resultcolnames1<- read.csv(paste("Weibull_raw_d_calibration",5400,9,".csv", sep = ","))

colnames(simulatedbatchWeibull_bias_Monte)<-colnames(resultcolnames1)

Monte_Weibull<-simulatedbatchWeibull_bias_Monte[,c(1,2,171:ncol(simulatedbatchWeibull_bias_Monte))]

simulatedbatchWeibull_bias_Monte<-c()

write.csv(Monte_Weibull,paste("Weibull_ds_raw_dcalibration_SWA_Finite_Monte.csv", sep = ","), row.names = FALSE)

meanall<-foreach(estimators = c(1:ncol(Monte_Weibull)), .combine = 'cbind') %dopar% {
  Monte_sample_Estimator1 <- data.frame(Size = ( Monte_Weibull[,1]),
                                        Kurtosis = (Monte_Weibull[,2]),
                                        Estimator = Monte_Weibull[,estimators])
  
  rownames(Monte_sample_Estimator1)<-c()
  
  Finite_1<-tapply(Monte_sample_Estimator1$Estimator, Monte_sample_Estimator1$Kurtosis, mean)
  
  Finite_1
}
colnames(meanall)<-colnames(Monte_Weibull)
write.csv(meanall,paste("Weibull_ds_raw_dcalibration_SWA_Finite_Monte.csv", sep = ","), row.names = FALSE)

musall1<-meanall[,3:ncol(meanall)]
finddall<-c()
for (estimators in (1:28)){
  mmm1<-musall1[,(9*estimators-8):(9*estimators)]
  mmall<-c()
  for (allestimators in (1:nrow(meanall))){
    mmm1one<-as.numeric(mmm1[allestimators,])
    finddmmm1one1<-findd(expectanalytic=mmm1one[1],expecttrue=mmm1one[2],expectboot=mmm1one[3],SWA1=mmm1one[4],median1=mmm1one[5],mx1=mmm1one[6],quatileexpectanalytic=mmm1one[7],
                         quatileexpecttrue=mmm1one[8],quatileexpectboot=mmm1one[9])
    mmall<-rbind(mmall,finddmmm1one1)
  }
  finddall<-cbind(finddall,mmall)
}
Monte_d_Weibull<-cbind(meanall[,1:2],finddall[,c(seq(from=1, to=166, by=3))])

Label_Weibull1<- read.csv(("finite_d_label.csv"))

colnames(Monte_d_Weibull)<-colnames(Label_Weibull1)

write.csv(Monte_d_Weibull,paste("d_value_finite",samplesize,".csv", sep = ","), row.names = FALSE)

Asymptotic_Weibull<- read.csv(("asymptotic_d_Weibull_SWA.csv"))

Asymptotic_Weibull<- cbind(Size=rep(1800000,nrow(Asymptotic_Weibull)),Asymptotic_Weibull)

colnames(Asymptotic_Weibull)<-colnames(Label_Weibull1)

AllFinal_Weibull<-rbind(Monte_d_Weibull,Asymptotic_Weibull)

write.csv(AllFinal_Weibull,paste("d_value_Weibull.csv", sep = ","), row.names = FALSE)

simulatedbatchWeibull_bias_Monte_SE<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
  library(Rfast)
  if (!require("foreach")) install.packages("foreach")
  library(foreach)
  if (!require("doParallel")) install.packages("doParallel")
  library(doParallel)
  #registering clusters, can set a smaller number using numCores-1 
  
  #require randtoolbox for random number generations
  if (!require("randtoolbox")) install.packages("randtoolbox")
  library(randtoolbox)
  #require Rfast for faster computation
  if (!require("Rfast")) install.packages("Rfast")
  library(Rfast)
  if (!require("gnorm")) install.packages("gnorm")
  library(gnorm)
  
  a=allkurtWeibull[batchnumber]
    
  n <- samplesize
  
  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  batchsize=batchsizebase
  
  SEbataches<- read.csv(paste("Weibull_raw_d_calibration",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  standarderrors1<-apply(SEbataches,2,se_mean)
  standarderrors1
}
Allstandarderror_each<-cbind(meanall[,1:2],simulatedbatchWeibull_bias_Monte_SE[,3:ncol(simulatedbatchWeibull_bias_Monte_SE)])

All_each1<-meanall[,3:ncol(meanall)]
Allstandarderror_each1<-Allstandarderror_each[,171:ncol(Allstandarderror_each)]

finddall_error<-c()
for (estimators in (1:28)){
  mmm1_error<-Allstandarderror_each1[,(9*estimators-8):(9*estimators)]
  mmm1_11<-All_each1[,(9*estimators-8):(9*estimators)]
  
  mmall_error<-c()
  for (allestimators in (1:nrow(Allstandarderror_each1))){
    mmm1one_error<-as.numeric(mmm1_error[allestimators,])
    mmm1one<-as.numeric(mmm1_11[allestimators,])
    
    finddmmm1one1_error<-findderror(expectanalytic_error=mmm1one_error[1],expecttrue_error=mmm1one_error[2],expectboot_error=mmm1one_error[3],SWA1_error=mmm1one_error[4],median1_error=mmm1one_error[5],mx1_error=mmm1one_error[6],quatileexpectanalytic_error=mmm1one_error[7],
                                    quatileexpecttrue_error=mmm1one_error[8],quatileexpectboot_error=mmm1one_error[9],expectanalytic=mmm1one[1],expecttrue=mmm1one[2],expectboot=mmm1one[3],SWA1=mmm1one[4],median1=mmm1one[5],mx1=mmm1one[6],quatileexpectanalytic=mmm1one[7],
                                    quatileexpecttrue=mmm1one[8],quatileexpectboot=mmm1one[9])
    mmall_error<-rbind(mmall_error,finddmmm1one1_error)
  }
  finddall_error<-cbind(finddall_error,mmall_error)
}
Monte_d_Weibull_error<-cbind(Allstandarderror_each[,1:2],finddall_error[,c(seq(from=1, to=166, by=3))])

colnames(Monte_d_Weibull_error)<-colnames(Label_Weibull1)

write.csv(Monte_d_Weibull_error,paste("finite_d_error.csv", sep = ","), row.names = FALSE)


registerDoSEQ()

