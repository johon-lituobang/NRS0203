
#require foreach and doparallel for parallel processing
if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)
#registering clusters, can set a smaller number using numCores-1 
numCores <- detectCores()

#require randtoolbox for random number generations
if (!require("randtoolbox")) install.packages("randtoolbox")
library(randtoolbox)
#require Rfast for faster computation
if (!require("Rfast")) install.packages("Rfast")
library(Rfast)
if (!require("gnorm")) install.packages("gnorm")
library(gnorm)
if (!require("stats")) install.packages("stats")
library(stats)

registerDoParallel(numCores)


asymptotic_n <- 1.8*10^6
(asymptotic_n%%10)==0
# maximum order of moments
morder <- 5
#large sample size (asymptotic bias)
largesize<-1.8*10^6

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol_asymptotic<-sobol(n=asymptotic_n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                mixed = FALSE, method = "C", start = 1)

quasiuni_asymptotic<-rbind(quasiunisobol_asymptotic)

quasiunisobol_asymptotic<-c()

quasiuni_sorted2_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:4], descend = FALSE, stable = FALSE, parallel = TRUE))

quasiuni_sorted5_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:5], descend = FALSE, stable = FALSE, parallel = TRUE))
# Forever...

quasiuni_asymptotic<-Sort(quasiuni_asymptotic[,1])

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
  }else if (dimension==5){
    ss1<-sortedx[orderlist[1]]
    ss2<-sortedx[orderlist[2]]
    ss3<-sortedx[orderlist[3]]
    ss4<-sortedx[orderlist[4]]
    ss5<-sortedx[orderlist[5]]
    return(c(ss1,ss2,ss3,ss4,ss5))
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
orderlist1_AB2_asymptotic<-removelist(na.omit(t(apply(quasiuni_sorted2_asymptotic,MARGIN=1,FUN=roundunique,dimension=2,size=largesize))))
orderlist1_AB3_asymptotic<-removelist(na.omit(t(apply(quasiuni_sorted3_asymptotic,MARGIN=1,FUN=roundunique,dimension=3,size=largesize))))
orderlist1_AB4_asymptotic<-removelist(na.omit(t(apply(quasiuni_sorted4_asymptotic,MARGIN=1,FUN=roundunique,dimension=4,size=largesize))))
orderlist1_AB5_asymptotic<-removelist(na.omit(t(apply(quasiuni_sorted5_asymptotic,MARGIN=1,FUN=roundunique,dimension=5,size=largesize))))


quasiuni_sorted2_asymptotic<-c()
quasiuni_sorted3_asymptotic<-c()
quasiuni_sorted4_asymptotic<-c()
quasiuni_sorted5_asymptotic<-c()
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
dsgnorm<-function (uni,shape, scale = 1){
  sample1<-qgnorm(p=uni, mu = 0, alpha = scale, beta = shape)
  sample1
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

sm<-function (x,interval=9,fast=TRUE,batch="auto"){
  lengthx<-length(x)
  if (batch=="auto" ){
    batch<-ceiling(500000/lengthx)+1
  }
  if (interval!=9 ){
    return("interval must be 9 ")
  }
  Ksamples<-lengthx/interval
  IntKsamples<-ceiling(Ksamples)
  target1<-IntKsamples*interval
  smass<-function(x1,IntKsamples,target1,sorted=FALSE){
    if(sorted){
      x_ordered<-x1
    }else{
      x_ordered<-Sort(x1,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    }
    if ((target1/IntKsamples)==9){
      onelist<-(c(x_ordered[1:(1*(IntKsamples))],x_ordered[(target1-1*IntKsamples+1):(target1)]))
      twolist<-(c(x_ordered[(IntKsamples+1):(2*(IntKsamples))],x_ordered[(target1-2*IntKsamples+1):(target1-IntKsamples)]))
      threelist<-(c(x_ordered[(2*(IntKsamples)+1):(3*(IntKsamples))],x_ordered[(6*(IntKsamples)+1):(target1-2*IntKsamples)]))
      fourlist<-(c(x_ordered[(3*(IntKsamples)+1):(4*(IntKsamples))],x_ordered[(5*(IntKsamples)+1):(target1-3*IntKsamples)]))
      fivelist<-(x_ordered[(4*(IntKsamples)+1):(5*(IntKsamples))])
      
      sonelist<-sum(onelist)
      stwolist<-sum(twolist)
      sthreelist<-sum(threelist)
      sfourlist<-sum(fourlist)
      sfivelist<-sum(fivelist)
      
      winsor<-(c(x_ordered[(IntKsamples+1)],x_ordered[((target1-IntKsamples))]))
      winsor1<-sum(winsor)*IntKsamples
      wm1<-(stwolist+sthreelist+sfourlist+sfivelist+winsor1)/(IntKsamples*(9))
      
      Groupmean<-c(sm=(sfivelist+stwolist)/(IntKsamples*3),wm=wm1)
      
      return(Groupmean)
    }
    else{
      return("Not supported yet.")
    }
  }
  if (Ksamples%%1!=0 ){
    if (fast==TRUE & lengthx<10000){
      x_ordered<-data_augmentation(x,target1)
      Groupmean<-smass(x_ordered,IntKsamples,target1,sorted=TRUE)
      return(Groupmean)
    }
    else{
      addedx<-matrix(sample(x,size=(target1-lengthx)*batch,replace=TRUE),nrow=batch)
      xt<-t(as.data.frame(x))
      xmatrix<-as.data.frame(lapply(xt, rep, batch))
      allmatrix<-cbind(addedx,xmatrix)
      batchresults<-apply(allmatrix,1,smass,IntKsamples=IntKsamples,target1=target1,sorted=FALSE)
      return(mean((batchresults)))}
  }
  else{
    Groupmean<-smass(x,IntKsamples,target1,sorted=FALSE)
    return(Groupmean)}
} 


#moments for checking the accuracy of bootstrap 

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

mom<-function (x, bend = 1.172,medianx=NULL,madx=NULL) {
  indicator1<-(x>medianx+bend*madx)
  indicator2<-(x<medianx-bend*madx)
  indicator <- rep(TRUE,length(x))
  indicator[indicator1]<-FALSE
  indicator[indicator2]<-FALSE
  result <- mean(x[indicator])
  return(result)
}
HuberPsi<-function(x,bend=1.172){
  result<-ifelse(abs(x)<=bend,x,bend*sign(x))
  return(result)
}
onestep<-function (x, bend = 1.172) {
  medianx<-median(x)
  madx<-mad(x)
  firstloc = mom(x, bend = bend,medianx=medianx,madx=madx)
  indicator <- (x - firstloc)/mad(x,constant=1)
  result<-medianx + madx * (sum(HuberPsi(x=indicator, bend=bend)))/(length(x[abs(indicator) <= bend]))
  return(c(onestep=result))
}

mmmraw<-function(x,interval=8,fast=TRUE,batch="auto"){
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  SWA1<-SWA(x=sortedx,interval=interval,fast=fast,batch=batch)
  if(SWA1[2]==Inf){
    return(print("BM is infinity, due to the double precision floating point limits. Usually, the solution is transforming your original data."))
  }
  lengthx<-length(x)
  mx1BM<-CDF(x=sortedx,xevaluated=SWA1[2],sorted=TRUE)
  mx1sqm<-CDF(x=sortedx,xevaluated=SWA1[3],sorted=TRUE)
  mx1wm1<-CDF(x=sortedx,xevaluated=SWA1[4],sorted=TRUE)
  mx1wm2<-CDF(x=sortedx,xevaluated=SWA1[5],sorted=TRUE)
  mx1tm1<-CDF(x=sortedx,xevaluated=SWA1[6],sorted=TRUE)
  mx1tm2<-CDF(x=sortedx,xevaluated=SWA1[7],sorted=TRUE)
  mx1tm3<-CDF(x=sortedx,xevaluated=SWA1[8],sorted=TRUE)
  output1<-c(mean=SWA1[1],BM=SWA1[2],SQM=SWA1[3],wm1=SWA1[4],wm2=SWA1[5],tm1=SWA1[6],tm2=SWA1[7],tm3=SWA1[8],median=SWA1[9],mx1BM=mx1BM,mx1sqm=mx1sqm,mx1wm1=mx1wm1,mx1wm2=mx1wm2,mx1tm1=mx1tm1,mx1tm2=mx1tm2,mx1tm3=mx1tm3)
  return(output1)
}

mmmprocessrm<-function(x,interval=8,SWA,median,mx1,drm=0.375){
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  lengthx<-length(x)
  rm1<--drm*median+SWA+drm*SWA
  names(rm1)<-NULL
  output1<-c(rm=rm1)
  return(output1)
}
mmmprocessqm<-function(x,interval=8,SWA,median,mx1,dqm=0.567){
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  lengthx<-length(x)
  if (mx1==0.5){
    quatiletarget<-0.5
  }else if (mx1>0.5){
    quatiletarget<-mx1+(mx1-0.5)*dqm
  }else{
    mx1<-1-mx1
    quatiletarget<-mx1+(mx1-0.5)*dqm
    quatiletarget<-1-quatiletarget
  }
  upper1<-(1-1/interval)
  lower1<-1/interval
  if (!is.na(quatiletarget) & quatiletarget>(upper1)){
    print(paste("Warning: the percentile exceeds ",as.character(upper1*interval),"/",as.character(interval),", the robustness shrinks."))
  }else if(!is.na(quatiletarget) & quatiletarget<(lower1)){
    print(paste("Warning: the percentile exceeds ",as.character(lower1*interval),"/",as.character(interval),", the robustness shrinks."))
  }
  if(quatiletarget>upper1){
    quatiletarget=upper1
  }else if(quatiletarget<lower1){
    quatiletarget=lower1
  }
  qm1<-quantile(sortedx,quatiletarget,type=8)
  output1<-c(qm=qm1)
  return(output1)
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

Balltest<-function (x,orderlist1_sorted2=NULL,orderlist1_sorted3=NULL,orderlist1_sorted4=NULL,orderlist1_sorted5=NULL,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=(1/10)){
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  lengthx<-length(sortedx)
  
  mmm1raw<-mmmraw(x=sortedx,interval=interval,fast=fast,batch=batch)
  
  mmm1_BM_rm_exp1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[2],median=mmm1raw[9],mx1=mmm1raw[10],drm=0.37523)
  
  mmm1_BM_qm_exp1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[2],median=mmm1raw[9],mx1=mmm1raw[10],dqm=0.32128)
  
  bootstrappedsample2<-t(as.data.frame(apply(orderlist1_sorted2,MARGIN=1,FUN=extract1,sortedx=sortedx,dimension=2)))
  
  hlkernel<-function(vector){ 
    return(mean(vector))
  }
  
  dp2x<-apply(bootstrappedsample2,MARGIN=1,FUN=hlkernel)
  
  HL2<-SWA(dp2x,interval=8,fast=TRUE,batch="auto")
  
  bootstrappedsample2<-c()
  
  bootstrappedsample3<-t(as.data.frame(apply(orderlist1_sorted3,MARGIN=1,FUN=extract1,sortedx=sortedx,dimension=3)))
  
  dp3x<-apply(bootstrappedsample3,MARGIN=1,FUN=hlkernel)
  HL3<-SWA(dp3x,interval=8,fast=TRUE,batch="auto")
  
  bootstrappedsample3<-c()
  
  bootstrappedsample4<-t(as.data.frame(apply(orderlist1_sorted4,MARGIN=1,FUN=extract1,sortedx=sortedx,dimension=4)))
  
  dp4x<-apply(bootstrappedsample4,MARGIN=1,FUN=hlkernel)
  
  HL4<-SWA(dp4x,interval=8,fast=TRUE,batch="auto")
  
  bootstrappedsample4<-c()
  
  bootstrappedsample5<-t(as.data.frame(apply(orderlist1_sorted5,MARGIN=1,FUN=extract1,sortedx=sortedx,dimension=5)))
  
  dp5x<-apply(bootstrappedsample5,MARGIN=1,FUN=hlkernel)
  
  HL5<-SWA(dp5x,interval=8,fast=TRUE,batch="auto")
  
  bootstrappedsample5<-c()
  
  finallall<-c(mmm1raw=mmm1raw,mmm1_BM_rm_exp1=mmm1_BM_rm_exp1,mmm1_BM_qm_exp1=mmm1_BM_qm_exp1,HL2=HL2,HL3=HL3,HL4=HL4,HL5=HL5)
  #finallall<-c(HL1=HL1,mmm1raw=mmm1raw,mmm1exp1=mmm1exp1,mmm1exp2=mmm1exp2,mmm1Weibull1=mmm1Weibull1,mmm1Weibull2=mmm1Weibull2,varmoraw=varmoraw,varmoexp1=varmoexp1,varmoexp2=varmoexp2,varmoWeibull1=varmoWeibull1,varmoWeibull2=varmoWeibull2,tmmoraw=tmmoraw,tmmoexp1=tmmoexp1,tmmoexp2=tmmoexp2,tmmoWeibull1=tmmoWeibull1,tmmoWeibull2=tmmoWeibull2,fmmoraw=fmmoraw,fmmoexp1=fmmoexp1,fmmoexp2=fmmoexp2,fmmoWeibull1=fmmoWeibull1,fmmoWeibull2=fmmoWeibull2,sdall=sdall)
  return(finallall)
}


#set the convergence criterion
criterionset=1/20

#Weibull

kurtWeibull<- read.csv(("kurtWeibull_31150.csv"))
allkurtWeibull<-unlist(kurtWeibull)

simulatedbatchWeibull_asymptoticbias<-foreach(batchnumber = (1:length(allkurtWeibull)), .combine = 'rbind') %dopar% {
  library(Rfast)
  a=allkurtWeibull[batchnumber]
  x<-c(dsWeibull(uni=quasiuni_asymptotic, shape=a/1, scale = 1))
  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  
  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  onestepx<-onestep(x=sortedx, bend = 1.172)
  SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
  dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,orderlist1_sorted5=orderlist1_AB5_asymptotic,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
  
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  
  allrawmoBias<-c(
    firstbias=abs(c(onestepx,SMWM9,dall)-dall[1])/(sqrt(momentsx[2]))
    )
  
  all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,allrawmoBias))
}

write.csv(simulatedbatchWeibull_asymptoticbias,paste("asymptotic_Weibull_SRM_Process",largesize,".csv", sep = ","), row.names = FALSE)

#gamma

kurtgamma<- read.csv(("kurtgamma_31150.csv"))
allkurtgamma<-unlist(kurtgamma)

simulatedbatchgamma_asymptoticbias<-foreach(batchnumber = (1:length(allkurtgamma)), .combine = 'rbind') %dopar% {
  library(Rfast)
  a=allkurtgamma[batchnumber]
  x<-c(dsgamma(uni=quasiuni_asymptotic, shape=a, scale = 1))
  targetm<-a
  targetvar<-(a)
  targettm<-((sqrt(a))^3)*2/sqrt(a)
  targetfm<-((sqrt(a))^4)*((6/(a))+3)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  
  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  onestepx<-onestep(x=sortedx, bend = 1.172)
  SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
  dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,orderlist1_sorted5=orderlist1_AB5_asymptotic,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
  
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  
  allrawmoBias<-c(
    firstbias=abs(c(onestepx,SMWM9,dall)-dall[1])/(sqrt(momentsx[2]))
  )
  
  all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,allrawmoBias))
}

write.csv(simulatedbatchgamma_asymptoticbias,paste("asymptotic_gamma_SRM_Process",largesize,".csv", sep = ","), row.names = FALSE)


#Pareto
kurtPareto<- read.csv(("kurtPareto_91210.csv"))
allkurtPareto<-unlist(kurtPareto)

simulatedbatchPareto_asymptoticbias<-foreach(batchnumber = (1:length(allkurtPareto)), .combine = 'rbind') %dopar% {
  library(Rfast)
  a=allkurtPareto[batchnumber]
  x<-c(dsPareto(uni=quasiuni_asymptotic, shape=a, scale = 1))
  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  
  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  onestepx<-onestep(x=sortedx, bend = 1.172)
  SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
  dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,orderlist1_sorted5=orderlist1_AB5_asymptotic,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
  
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  
  allrawmoBias<-c(
    firstbias=abs(c(onestepx,SMWM9,dall)-dall[1])/(sqrt(momentsx[2]))
  )
  
  all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,allrawmoBias))
}

write.csv(simulatedbatchPareto_asymptoticbias,paste("asymptotic_Pareto_SRM_Process",largesize,".csv", sep = ","), row.names = FALSE)



#lognorm

kurtlognorm<- read.csv(("kurtlognorm_31150.csv"))
allkurtlognorm<-unlist(kurtlognorm)

simulatedbatchlognorm_asymptoticbias<-foreach(batchnumber = (1:length(allkurtlognorm)), .combine = 'rbind') %dopar% {
  library(Rfast)
  a=allkurtlognorm[batchnumber]
  x<-c(dslnorm(uni=quasiuni_asymptotic, location=0, scale = a/1))
  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  
  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  onestepx<-onestep(x=sortedx, bend = 1.172)
  SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
  dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,orderlist1_sorted5=orderlist1_AB5_asymptotic,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
  
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  
  allrawmoBias<-c(
    firstbias=abs(c(onestepx,SMWM9,dall)-dall[1])/(sqrt(momentsx[2]))
  )
  
  all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,allrawmoBias))
}

write.csv(simulatedbatchlognorm_asymptoticbias,paste("asymptotic_lognorm_SRM_Process",largesize,".csv", sep = ","), row.names = FALSE)



registerDoSEQ()

