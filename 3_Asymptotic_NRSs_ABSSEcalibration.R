
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


#bootsize for bootstrap approximation of the distributions of the kernel of U-statistics.
n <- 1.8*10^6
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

#load asymptotic d for two parameter distributions
Weibull_d<- read.csv(("d_value_Weibull.csv"))
# #load asymptotic d for two parameter distributions
# gamma_d<- read.csv(("asymptotic_d_gamma319.csv"))
# #load asymptotic d for two parameter distributions
# lognormal_d<- read.csv(("asymptotic_d_lognorm319.csv"))
# # #load asymptotic d for two parameter distributions
# Pareto_d<- read.csv(("asymptotic_d_Pareto919.csv"))

d_adjust<-function(size,kurt,dlist,type){
  dlist[,2]<-round(dlist[,2],digits=1)
  indextypelist<-c("mean_BM_drm","mean_BM_dqm","mean_sqm_drm","mean_sqm_dqm","mean_wm1_drm","mean_wm1_dqm","mean_wm2_drm","mean_wm2_dqm","mean_tm1_drm","mean_tm1_dqm","mean_tm2_drm","mean_tm2_dqm","mean_tm3_drm","mean_tm3_dqm","var_BM_drm","var_BM_dqm","var_sqm_drm","var_sqm_dqm","var_wm1_drm","var_wm1_dqm","var_wm2_drm","var_wm2_dqm","var_tm1_drm","var_tm1_dqm","var_tm2_drm","var_tm2_dqm","var_tm3_drm","var_tm3_dqm","tm_BM_drm","tm_BM_dqm","tm_sqm_drm","tm_sqm_dqm","tm_wm1_drm","tm_wm1_dqm","tm_wm2_drm","tm_wm2_dqm","tm_tm1_drm","tm_tm1_dqm","tm_tm2_drm","tm_tm2_dqm","tm_tm3_drm","tm_tm3_dqm","fm_BM_drm","fm_BM_dqm","fm_sqm_drm","fm_sqm_dqm","fm_wm1_drm","fm_wm1_dqm","fm_wm2_drm","fm_wm2_dqm","fm_tm1_drm","fm_tm1_dqm","fm_tm2_drm","fm_tm2_dqm","fm_tm3_drm","fm_tm3_dqm")
  indextype<-which(indextypelist==(type))+2
  if(size%in% dlist[,1]){
    if (kurt%in% dlist[,2]){
      result1<-dlist[dlist[,1] == size & dlist[,2]==kurt,indextype]
    }else{
      rown2<-as.numeric(dlist[,2])
      infn2<-max(rown2[rown2 < kurt])
      if(is.infinite(infn2)){
        infn2<-rown2[2]
        #print(c("The kurtosis is out of range supported.",kurt))
      }
      supn2<-min(rown2[rown2 > kurt])
      if(is.infinite(supn2)){
        supn2<-max(rown2)
        #print(c("The kurtosis is out of range supported.",kurt))
      }
      d1<-dlist[dlist[,1]==size & dlist[,2]==infn2,indextype]
      d2<-dlist[dlist[,1]==size & dlist[,2]==supn2,indextype]
      if(d1==d2){
        result1<-d1
      }else{
        result1<-(((d2-d1)*((kurt-infn2)/(supn2-infn2)))+d1)
      }
    }
  }else{
    rown<-as.numeric(dlist[,1])
    infn<-max(rown[rown < size])
    if(is.infinite(infn)){
      infn<-rown[1]
    }
    supn<-min(rown[rown > size])
    if(is.infinite(supn)){
      supn<-max(rown)
    }
    
    rown2<-as.numeric(dlist[,2])
    infn2<-max(rown2[rown2 < kurt])
    if(is.infinite(infn2)){
      infn2<-rown2[2]
    }
    supn2<-min(rown2[rown2 > kurt])
    if(is.infinite(supn2)){
      supn2<-max(rown2)
    }
    
    d1<-dlist[dlist[,1]==infn & dlist[,2]==infn2,indextype]
    d2<-dlist[dlist[,1]==supn & dlist[,2]==supn2,indextype]
    
    
    if(d1==d2){
      result1<-d1
    }else if (supn!=infn){
      result1<-(((d2-d1)*((size-infn)/(supn-infn)))+d1)
    }else{
      result1<-(((d2-d1)*((kurt-infn2)/(supn2-infn2)))+d1)
    }
  }
  return(result1)
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


Balltest<-function (x,orderlist1_sorted2=NULL,orderlist1_sorted3=NULL,orderlist1_sorted4=NULL,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=(1/10)){
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  lengthx<-length(sortedx)
  
  bootstrappedsample2<-t(as.data.frame(apply(orderlist1_sorted2,MARGIN=1,FUN=extract1,sortedx=sortedx,dimension=2)))
  
  getvar<-function(vector){ 
    ((vector[1]-vector[2])^2)/2
  }
  
  dp2varx<-apply(bootstrappedsample2,MARGIN=1,FUN=getvar)
  
  getHL<-function(vector){ 
    ((vector[1]+vector[2]))/2
  }
  
  dp2HLx<-apply(bootstrappedsample2,MARGIN=1,FUN=getHL)
  
  HL1<-median(dp2HLx)
  
  bootstrappedsample2<-c()
  
  bootstrappedsample3<-t(as.data.frame(apply(orderlist1_sorted3,MARGIN=1,FUN=extract1,sortedx=sortedx,dimension=3)))
  
  gettm<-function(vector){ 
    ((1/6)*(2*vector[1]-vector[2]-vector[3])*(-1*vector[1]+2*vector[2]-vector[3])*(-vector[1]-vector[2]+2*vector[3]))
  }
  
  dp3tmx<-apply(bootstrappedsample3,MARGIN=1,FUN=gettm)
  
  bootstrappedsample3<-c()
  
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
  
  bootstrappedsample4<-c()
  
  mmm1raw<-mmmraw(x=sortedx,interval=interval,fast=fast,batch=batch)
  
  varmoraw<-mmmraw(x=dp2varx,interval=interval,fast=fast,batch=batch)
  
  tmmoraw<-mmmraw(x=dp3tmx,interval=interval,fast=fast,batch=batch)
  
  fmmoraw<-mmmraw(x=dp4fmx,interval=interval,fast=fast,batch=batch)
  
  standist_d=Weibull_d
  
  #exponential
  startpoint=9
  
  mean_BM_drm1<-d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_BM_drm")
  var_BM_drm1<-d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_BM_drm")
  tm_BM_drm1<-d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_BM_drm")
  fm_BM_drm1<-d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_BM_drm")
  
  mmm1_BM_rm_exp1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[2],median=mmm1raw[9],mx1=mmm1raw[10],drm=mean_BM_drm1)
  
  varmo_BM_rm_exp1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[2],median=varmoraw[9],mx1=varmoraw[10],drm=var_BM_drm1)
  
  tmmo_BM_rm_exp1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[2],median=tmmoraw[9],mx1=tmmoraw[10],drm=tm_BM_drm1)
  
  fmmo_BM_rm_exp1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[2],median=fmmoraw[9],mx1=fmmoraw[10],drm=fm_BM_drm1)
  
  
  rkurt_BM_rm<-((fmmo_BM_rm_exp1))/(varmo_BM_rm_exp1^2)
  
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_BM_rm && rkurt_BM_rm<(startpoint*(1)/(1-criterion)))){
    mmm1_BM_rm_Weibull1<-mmm1_BM_rm_exp1
    varmo_BM_rm_Weibull1<-varmo_BM_rm_exp1
    tmmo_BM_rm_Weibull1<-tmmo_BM_rm_exp1
    fmmo_BM_rm_Weibull1<-fmmo_BM_rm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_BM_drm1<-d_adjust(size=lengthx,kurt=rkurt_BM_rm,dlist=standist_d,type="var_BM_drm")
      
      fm_BM_drm1<-d_adjust(size=lengthx,kurt=rkurt_BM_rm,dlist=standist_d,type="fm_BM_drm")
      
      varmo_BM_rm_Weibull1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[2],median=varmoraw[9],mx1=varmoraw[10],drm=var_BM_drm1)
      
      fmmo_BM_rm_Weibull1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[2],median=fmmoraw[9],mx1=fmmoraw[10],drm=fm_BM_drm1)
      
      newrrkurt_BM_rm<-((fmmo_BM_rm_Weibull1))/(varmo_BM_rm_Weibull1^2)
      
      if ((((rkurt_BM_rm*(1-criterion))<newrrkurt_BM_rm && newrrkurt_BM_rm<(rkurt_BM_rm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_BM_drm1<-d_adjust(size=lengthx,kurt=rkurt_BM_rm,dlist=standist_d,type="mean_BM_drm")
        
        tm_BM_drm1<-d_adjust(size=lengthx,kurt=rkurt_BM_rm,dlist=standist_d,type="tm_BM_drm")
        
        mmm1_BM_rm_Weibull1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[2],median=mmm1raw[9],mx1=mmm1raw[10],drm=mean_BM_drm1)
        
        tmmo_BM_rm_Weibull1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[2],median=tmmoraw[9],mx1=tmmoraw[10],drm=tm_BM_drm1)
        
        break
      }
      rkurt_BM_rm<-newrrkurt_BM_rm
      
    }
    
  }
  
  mean_BM_dqm1<-d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_BM_dqm")
  var_BM_dqm1<-d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_BM_dqm")
  tm_BM_dqm1<-d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_BM_dqm")
  fm_BM_dqm1<-d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_BM_dqm")
  
  mmm1_BM_qm_exp1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[2],median=mmm1raw[9],mx1=mmm1raw[10],dqm=mean_BM_dqm1)
  
  varmo_BM_qm_exp1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[2],median=varmoraw[9],mx1=varmoraw[10],dqm=var_BM_dqm1)
  
  tmmo_BM_qm_exp1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[2],median=tmmoraw[9],mx1=tmmoraw[10],dqm=tm_BM_dqm1)
  
  fmmo_BM_qm_exp1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[2],median=fmmoraw[9],mx1=fmmoraw[10],dqm=fm_BM_dqm1)
  
  rkurt_BM_qm<-((fmmo_BM_qm_exp1))/(varmo_BM_qm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_BM_qm && rkurt_BM_qm<(startpoint*(1)/(1-criterion)))){
    mmm1_BM_qm_Weibull1<-mmm1_BM_qm_exp1
    varmo_BM_qm_Weibull1<-varmo_BM_qm_exp1
    tmmo_BM_qm_Weibull1<-tmmo_BM_qm_exp1
    fmmo_BM_qm_Weibull1<-fmmo_BM_qm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_BM_dqm1<-d_adjust(size=lengthx,kurt=rkurt_BM_qm,dlist=standist_d,type="var_BM_dqm")
      
      fm_BM_dqm1<-d_adjust(size=lengthx,kurt=rkurt_BM_qm,dlist=standist_d,type="fm_BM_dqm")
      
      varmo_BM_qm_Weibull1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[2],median=varmoraw[9],mx1=varmoraw[10],dqm=var_BM_dqm1)
      
      fmmo_BM_qm_Weibull1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[2],median=fmmoraw[9],mx1=fmmoraw[10],dqm=fm_BM_dqm1)
      
      newrrkurt_BM_qm<-((fmmo_BM_qm_Weibull1))/(varmo_BM_qm_Weibull1^2)
      
      if ((((rkurt_BM_qm*(1-criterion))<newrrkurt_BM_qm && newrrkurt_BM_qm<(rkurt_BM_qm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_BM_dqm1<-d_adjust(size=lengthx,kurt=rkurt_BM_qm,dlist=standist_d,type="mean_BM_dqm")
        
        tm_BM_dqm1<-d_adjust(size=lengthx,kurt=rkurt_BM_qm,dlist=standist_d,type="tm_BM_dqm")
        
        mmm1_BM_qm_Weibull1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[2],median=mmm1raw[9],mx1=mmm1raw[10],dqm=mean_BM_dqm1)
        
        tmmo_BM_qm_Weibull1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[2],median=tmmoraw[9],mx1=tmmoraw[10],dqm=tm_BM_dqm1)
        
        break
      }
      rkurt_BM_qm<-newrrkurt_BM_qm
      
    }
    
  }
  
  mean_sqm_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_sqm_drm")
  var_sqm_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_sqm_drm")
  tm_sqm_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_sqm_drm")
  fm_sqm_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_sqm_drm")
  
  
  mmm1_sqm_rm_exp1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[3],median=mmm1raw[9],mx1=mmm1raw[11],drm=mean_sqm_drm1)
  
  varmo_sqm_rm_exp1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[3],median=varmoraw[9],mx1=varmoraw[11],drm=var_sqm_drm1)
  
  tmmo_sqm_rm_exp1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[3],median=tmmoraw[9],mx1=tmmoraw[11],drm=tm_sqm_drm1)
  
  fmmo_sqm_rm_exp1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[3],median=fmmoraw[9],mx1=fmmoraw[11],drm=fm_sqm_drm1)
  
  
  rkurt_sqm_rm<-((fmmo_sqm_rm_exp1))/(varmo_sqm_rm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_sqm_rm && rkurt_sqm_rm<(startpoint*(1)/(1-criterion)))){
    mmm1_sqm_rm_Weibull1<-mmm1_sqm_rm_exp1
    varmo_sqm_rm_Weibull1<-varmo_sqm_rm_exp1
    tmmo_sqm_rm_Weibull1<-tmmo_sqm_rm_exp1
    fmmo_sqm_rm_Weibull1<-fmmo_sqm_rm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_sqm_drm1<-d_adjust(size=lengthx,kurt=rkurt_sqm_rm,dlist=standist_d,type="var_sqm_drm")
      
      fm_sqm_drm1<-d_adjust(size=lengthx,kurt=rkurt_sqm_rm,dlist=standist_d,type="fm_sqm_drm")
      
      varmo_sqm_rm_Weibull1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[3],median=varmoraw[9],mx1=varmoraw[11],drm=var_sqm_drm1)
      
      fmmo_sqm_rm_Weibull1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[3],median=fmmoraw[9],mx1=fmmoraw[11],drm=fm_sqm_drm1)
      
      newrrkurt_sqm_rm<-((fmmo_sqm_rm_Weibull1))/(varmo_sqm_rm_Weibull1^2)
      
      if ((((rkurt_sqm_rm*(1-criterion))<newrrkurt_sqm_rm && newrrkurt_sqm_rm<(rkurt_sqm_rm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_sqm_drm1<-d_adjust(size=lengthx,kurt=rkurt_sqm_rm,dlist=standist_d,type="mean_sqm_drm")
        
        tm_sqm_drm1<-d_adjust(size=lengthx,kurt=rkurt_sqm_rm,dlist=standist_d,type="tm_sqm_drm")
        
        mmm1_sqm_rm_Weibull1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[3],median=mmm1raw[9],mx1=mmm1raw[11],drm=mean_sqm_drm1)
        
        tmmo_sqm_rm_Weibull1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[3],median=tmmoraw[9],mx1=tmmoraw[11],drm=tm_sqm_drm1)
        
        break
      }
      rkurt_sqm_rm<-newrrkurt_sqm_rm
      
    }
    
  }
  
  
  mean_sqm_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_sqm_dqm")
  var_sqm_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_sqm_dqm")
  tm_sqm_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_sqm_dqm")
  fm_sqm_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_sqm_dqm")
  
  mmm1_sqm_qm_exp1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[3],median=mmm1raw[9],mx1=mmm1raw[11],dqm=mean_sqm_dqm1)
  
  varmo_sqm_qm_exp1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[3],median=varmoraw[9],mx1=varmoraw[11],dqm=var_sqm_dqm1)
  
  tmmo_sqm_qm_exp1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[3],median=tmmoraw[9],mx1=tmmoraw[11],dqm=tm_sqm_dqm1)
  
  fmmo_sqm_qm_exp1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[3],median=fmmoraw[9],mx1=fmmoraw[11],dqm=fm_sqm_dqm1)
  
  rkurt_sqm_qm<-((fmmo_sqm_qm_exp1))/(varmo_sqm_qm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_sqm_qm && rkurt_sqm_qm<(startpoint*(1)/(1-criterion)))){
    mmm1_sqm_qm_Weibull1<-mmm1_sqm_qm_exp1
    varmo_sqm_qm_Weibull1<-varmo_sqm_qm_exp1
    tmmo_sqm_qm_Weibull1<-tmmo_sqm_qm_exp1
    fmmo_sqm_qm_Weibull1<-fmmo_sqm_qm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_sqm_dqm1<-d_adjust(size=lengthx,kurt=rkurt_sqm_qm,dlist=standist_d,type="var_sqm_dqm")
      
      fm_sqm_dqm1<-d_adjust(size=lengthx,kurt=rkurt_sqm_qm,dlist=standist_d,type="fm_sqm_dqm")
      
      varmo_sqm_qm_Weibull1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[3],median=varmoraw[9],mx1=varmoraw[11],dqm=var_sqm_dqm1)
      
      fmmo_sqm_qm_Weibull1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[3],median=fmmoraw[9],mx1=fmmoraw[11],dqm=fm_sqm_dqm1)
      
      newrrkurt_sqm_qm<-((fmmo_sqm_qm_Weibull1))/(varmo_sqm_qm_Weibull1^2)
      
      if ((((rkurt_sqm_qm*(1-criterion))<newrrkurt_sqm_qm && newrrkurt_sqm_qm<(rkurt_sqm_qm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_sqm_dqm1<-d_adjust(size=lengthx,kurt=rkurt_sqm_qm,dlist=standist_d,type="mean_sqm_dqm")
        
        tm_sqm_dqm1<-d_adjust(size=lengthx,kurt=rkurt_sqm_qm,dlist=standist_d,type="tm_sqm_dqm")
        
        mmm1_sqm_qm_Weibull1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[3],median=mmm1raw[9],mx1=mmm1raw[11],dqm=mean_sqm_dqm1)
        
        tmmo_sqm_qm_Weibull1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[3],median=tmmoraw[9],mx1=tmmoraw[11],dqm=tm_sqm_dqm1)
        
        break
      }
      rkurt_sqm_qm<-newrrkurt_sqm_qm
      
    }
    
  }
  
  mean_wm1_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_wm1_drm")
  var_wm1_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_wm1_drm")
  tm_wm1_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_wm1_drm")
  fm_wm1_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_wm1_drm")
  
  
  mmm1_wm1_rm_exp1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[4],median=mmm1raw[9],mx1=mmm1raw[12],drm=mean_wm1_drm1)
  
  varmo_wm1_rm_exp1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[4],median=varmoraw[9],mx1=varmoraw[12],drm=var_wm1_drm1)
  
  tmmo_wm1_rm_exp1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[4],median=tmmoraw[9],mx1=tmmoraw[12],drm=tm_wm1_drm1)
  
  fmmo_wm1_rm_exp1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[4],median=fmmoraw[9],mx1=fmmoraw[12],drm=fm_wm1_drm1)
  
  
  rkurt_wm1_rm<-((fmmo_wm1_rm_exp1))/(varmo_wm1_rm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_wm1_rm && rkurt_wm1_rm<(startpoint*(1)/(1-criterion)))){
    mmm1_wm1_rm_Weibull1<-mmm1_wm1_rm_exp1
    varmo_wm1_rm_Weibull1<-varmo_wm1_rm_exp1
    tmmo_wm1_rm_Weibull1<-tmmo_wm1_rm_exp1
    fmmo_wm1_rm_Weibull1<-fmmo_wm1_rm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_wm1_drm1<-d_adjust(size=lengthx,kurt=rkurt_wm1_rm,dlist=standist_d,type="var_wm1_drm")
      
      fm_wm1_drm1<-d_adjust(size=lengthx,kurt=rkurt_wm1_rm,dlist=standist_d,type="fm_wm1_drm")
      
      varmo_wm1_rm_Weibull1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[4],median=varmoraw[9],mx1=varmoraw[12],drm=var_wm1_drm1)
      
      fmmo_wm1_rm_Weibull1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[4],median=fmmoraw[9],mx1=fmmoraw[12],drm=fm_wm1_drm1)
      
      newrrkurt_wm1_rm<-((fmmo_wm1_rm_Weibull1))/(varmo_wm1_rm_Weibull1^2)
      
      if ((((rkurt_wm1_rm*(1-criterion))<newrrkurt_wm1_rm && newrrkurt_wm1_rm<(rkurt_wm1_rm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_wm1_drm1<-d_adjust(size=lengthx,kurt=rkurt_wm1_rm,dlist=standist_d,type="mean_wm1_drm")
        
        tm_wm1_drm1<-d_adjust(size=lengthx,kurt=rkurt_wm1_rm,dlist=standist_d,type="tm_wm1_drm")
        
        mmm1_wm1_rm_Weibull1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[4],median=mmm1raw[9],mx1=mmm1raw[12],drm=mean_wm1_drm1)
        
        tmmo_wm1_rm_Weibull1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[4],median=tmmoraw[9],mx1=tmmoraw[12],drm=tm_wm1_drm1)
        
        break
      }
      rkurt_wm1_rm<-newrrkurt_wm1_rm
      
    }
    
  }
  
  
  mean_wm1_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_wm1_dqm")
  var_wm1_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_wm1_dqm")
  tm_wm1_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_wm1_dqm")
  fm_wm1_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_wm1_dqm")
  
  
  mmm1_wm1_qm_exp1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[4],median=mmm1raw[9],mx1=mmm1raw[12],dqm=mean_wm1_dqm1)
  
  varmo_wm1_qm_exp1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[4],median=varmoraw[9],mx1=varmoraw[12],dqm=var_wm1_dqm1)
  
  tmmo_wm1_qm_exp1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[4],median=tmmoraw[9],mx1=tmmoraw[12],dqm=tm_wm1_dqm1)
  
  fmmo_wm1_qm_exp1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[4],median=fmmoraw[9],mx1=fmmoraw[12],dqm=fm_wm1_dqm1)
  
  rkurt_wm1_qm<-((fmmo_wm1_qm_exp1))/(varmo_wm1_qm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_wm1_qm && rkurt_wm1_qm<(startpoint*(1)/(1-criterion)))){
    mmm1_wm1_qm_Weibull1<-mmm1_wm1_qm_exp1
    varmo_wm1_qm_Weibull1<-varmo_wm1_qm_exp1
    tmmo_wm1_qm_Weibull1<-tmmo_wm1_qm_exp1
    fmmo_wm1_qm_Weibull1<-fmmo_wm1_qm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_wm1_dqm1<-d_adjust(size=lengthx,kurt=rkurt_wm1_qm,dlist=standist_d,type="var_wm1_dqm")
      
      fm_wm1_dqm1<-d_adjust(size=lengthx,kurt=rkurt_wm1_qm,dlist=standist_d,type="fm_wm1_dqm")
      
      varmo_wm1_qm_Weibull1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[4],median=varmoraw[9],mx1=varmoraw[12],dqm=var_wm1_dqm1)
      
      fmmo_wm1_qm_Weibull1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[4],median=fmmoraw[9],mx1=fmmoraw[12],dqm=fm_wm1_dqm1)
      
      newrrkurt_wm1_qm<-((fmmo_wm1_qm_Weibull1))/(varmo_wm1_qm_Weibull1^2)
      
      if ((((rkurt_wm1_qm*(1-criterion))<newrrkurt_wm1_qm && newrrkurt_wm1_qm<(rkurt_wm1_qm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_wm1_dqm1<-d_adjust(size=lengthx,kurt=rkurt_wm1_qm,dlist=standist_d,type="mean_wm1_dqm")
        
        tm_wm1_dqm1<-d_adjust(size=lengthx,kurt=rkurt_wm1_qm,dlist=standist_d,type="tm_wm1_dqm")
        
        mmm1_wm1_qm_Weibull1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[4],median=mmm1raw[9],mx1=mmm1raw[12],dqm=mean_wm1_dqm1)
        
        tmmo_wm1_qm_Weibull1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[4],median=tmmoraw[9],mx1=tmmoraw[12],dqm=tm_wm1_dqm1)
        
        break
      }
      rkurt_wm1_qm<-newrrkurt_wm1_qm
      
    }
    
  }
  
  
  mean_wm2_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_wm2_drm")
  var_wm2_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_wm2_drm")
  tm_wm2_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_wm2_drm")
  fm_wm2_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_wm2_drm")
  
  mmm1_wm2_rm_exp1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[5],median=mmm1raw[9],mx1=mmm1raw[13],drm=mean_wm2_drm1)
  
  varmo_wm2_rm_exp1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[5],median=varmoraw[9],mx1=varmoraw[13],drm=var_wm2_drm1)
  
  tmmo_wm2_rm_exp1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[5],median=tmmoraw[9],mx1=tmmoraw[13],drm=tm_wm2_drm1)
  
  fmmo_wm2_rm_exp1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[5],median=fmmoraw[9],mx1=fmmoraw[13],drm=fm_wm2_drm1)
  
  
  rkurt_wm2_rm<-((fmmo_wm2_rm_exp1))/(varmo_wm2_rm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_wm2_rm && rkurt_wm2_rm<(startpoint*(1)/(1-criterion)))){
    mmm1_wm2_rm_Weibull1<-mmm1_wm2_rm_exp1
    varmo_wm2_rm_Weibull1<-varmo_wm2_rm_exp1
    tmmo_wm2_rm_Weibull1<-tmmo_wm2_rm_exp1
    fmmo_wm2_rm_Weibull1<-fmmo_wm2_rm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_wm2_drm1<-d_adjust(size=lengthx,kurt=rkurt_wm2_rm,dlist=standist_d,type="var_wm2_drm")
      
      fm_wm2_drm1<-d_adjust(size=lengthx,kurt=rkurt_wm2_rm,dlist=standist_d,type="fm_wm2_drm")
      
      varmo_wm2_rm_Weibull1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[5],median=varmoraw[9],mx1=varmoraw[13],drm=var_wm2_drm1)
      
      fmmo_wm2_rm_Weibull1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[5],median=fmmoraw[9],mx1=fmmoraw[13],drm=fm_wm2_drm1)
      
      newrrkurt_wm2_rm<-((fmmo_wm2_rm_Weibull1))/(varmo_wm2_rm_Weibull1^2)
      
      if ((((rkurt_wm2_rm*(1-criterion))<newrrkurt_wm2_rm && newrrkurt_wm2_rm<(rkurt_wm2_rm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_wm2_drm1<-d_adjust(size=lengthx,kurt=rkurt_wm2_rm,dlist=standist_d,type="mean_wm2_drm")
        
        tm_wm2_drm1<-d_adjust(size=lengthx,kurt=rkurt_wm2_rm,dlist=standist_d,type="tm_wm2_drm")
        
        mmm1_wm2_rm_Weibull1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[5],median=mmm1raw[9],mx1=mmm1raw[13],drm=mean_wm2_drm1)
        
        tmmo_wm2_rm_Weibull1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[5],median=tmmoraw[9],mx1=tmmoraw[13],drm=tm_wm2_drm1)
        
        break
      }
      rkurt_wm2_rm<-newrrkurt_wm2_rm
      
    }
    
  }
  
  
  mean_wm2_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_wm2_dqm")
  var_wm2_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_wm2_dqm")
  tm_wm2_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_wm2_dqm")
  fm_wm2_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_wm2_dqm")
  
  
  mmm1_wm2_qm_exp1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[5],median=mmm1raw[9],mx1=mmm1raw[13],dqm=mean_wm2_dqm1)
  
  varmo_wm2_qm_exp1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[5],median=varmoraw[9],mx1=varmoraw[13],dqm=var_wm2_dqm1)
  
  tmmo_wm2_qm_exp1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[5],median=tmmoraw[9],mx1=tmmoraw[13],dqm=tm_wm2_dqm1)
  
  fmmo_wm2_qm_exp1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[5],median=fmmoraw[9],mx1=fmmoraw[13],dqm=fm_wm2_dqm1)
  
  rkurt_wm2_qm<-((fmmo_wm2_qm_exp1))/(varmo_wm2_qm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_wm2_qm && rkurt_wm2_qm<(startpoint*(1)/(1-criterion)))){
    mmm1_wm2_qm_Weibull1<-mmm1_wm2_qm_exp1
    varmo_wm2_qm_Weibull1<-varmo_wm2_qm_exp1
    tmmo_wm2_qm_Weibull1<-tmmo_wm2_qm_exp1
    fmmo_wm2_qm_Weibull1<-fmmo_wm2_qm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_wm2_dqm1<-d_adjust(size=lengthx,kurt=rkurt_wm2_qm,dlist=standist_d,type="var_wm2_dqm")
      
      fm_wm2_dqm1<-d_adjust(size=lengthx,kurt=rkurt_wm2_qm,dlist=standist_d,type="fm_wm2_dqm")
      
      varmo_wm2_qm_Weibull1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[5],median=varmoraw[9],mx1=varmoraw[13],dqm=var_wm2_dqm1)
      
      fmmo_wm2_qm_Weibull1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[5],median=fmmoraw[9],mx1=fmmoraw[13],dqm=fm_wm2_dqm1)
      
      newrrkurt_wm2_qm<-((fmmo_wm2_qm_Weibull1))/(varmo_wm2_qm_Weibull1^2)
      
      if ((((rkurt_wm2_qm*(1-criterion))<newrrkurt_wm2_qm && newrrkurt_wm2_qm<(rkurt_wm2_qm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_wm2_dqm1<-d_adjust(size=lengthx,kurt=rkurt_wm2_qm,dlist=standist_d,type="mean_wm2_dqm")
        
        tm_wm2_dqm1<-d_adjust(size=lengthx,kurt=rkurt_wm2_qm,dlist=standist_d,type="tm_wm2_dqm")
        
        mmm1_wm2_qm_Weibull1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[5],median=mmm1raw[9],mx1=mmm1raw[13],dqm=mean_wm2_dqm1)
        
        tmmo_wm2_qm_Weibull1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[5],median=tmmoraw[9],mx1=tmmoraw[13],dqm=tm_wm2_dqm1)
        
        break
      }
      rkurt_wm2_qm<-newrrkurt_wm2_qm
      
    }
    
  }
  
  mean_tm1_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_tm1_drm")
  var_tm1_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_tm1_drm")
  tm_tm1_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_tm1_drm")
  fm_tm1_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_tm1_drm")
  
  mmm1_tm1_rm_exp1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[6],median=mmm1raw[9],mx1=mmm1raw[14],drm=mean_tm1_drm1)
  
  varmo_tm1_rm_exp1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[6],median=varmoraw[9],mx1=varmoraw[14],drm=var_tm1_drm1)
  
  tmmo_tm1_rm_exp1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[6],median=tmmoraw[9],mx1=tmmoraw[14],drm=tm_tm1_drm1)
  
  fmmo_tm1_rm_exp1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[6],median=fmmoraw[9],mx1=fmmoraw[14],drm=fm_tm1_drm1)
  
  
  rkurt_tm1_rm<-((fmmo_tm1_rm_exp1))/(varmo_tm1_rm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_tm1_rm && rkurt_tm1_rm<(startpoint*(1)/(1-criterion)))){
    mmm1_tm1_rm_Weibull1<-mmm1_tm1_rm_exp1
    varmo_tm1_rm_Weibull1<-varmo_tm1_rm_exp1
    tmmo_tm1_rm_Weibull1<-tmmo_tm1_rm_exp1
    fmmo_tm1_rm_Weibull1<-fmmo_tm1_rm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_tm1_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm1_rm,dlist=standist_d,type="var_tm1_drm")
      
      fm_tm1_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm1_rm,dlist=standist_d,type="fm_tm1_drm")
      
      varmo_tm1_rm_Weibull1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[6],median=varmoraw[9],mx1=varmoraw[14],drm=var_tm1_drm1)
      
      fmmo_tm1_rm_Weibull1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[6],median=fmmoraw[9],mx1=fmmoraw[14],drm=fm_tm1_drm1)
      
      newrrkurt_tm1_rm<-((fmmo_tm1_rm_Weibull1))/(varmo_tm1_rm_Weibull1^2)
      
      if ((((rkurt_tm1_rm*(1-criterion))<newrrkurt_tm1_rm && newrrkurt_tm1_rm<(rkurt_tm1_rm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_tm1_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm1_rm,dlist=standist_d,type="mean_tm1_drm")
        
        tm_tm1_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm1_rm,dlist=standist_d,type="tm_tm1_drm")
        
        mmm1_tm1_rm_Weibull1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[6],median=mmm1raw[9],mx1=mmm1raw[14],drm=mean_tm1_drm1)
        
        tmmo_tm1_rm_Weibull1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[6],median=tmmoraw[9],mx1=tmmoraw[14],drm=tm_tm1_drm1)
        
        break
      }
      rkurt_tm1_rm<-newrrkurt_tm1_rm
      
    }
    
  }
  
  
  mean_tm1_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_tm1_dqm")
  var_tm1_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_tm1_dqm")
  tm_tm1_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_tm1_dqm")
  fm_tm1_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_tm1_dqm")
  
  mmm1_tm1_qm_exp1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[6],median=mmm1raw[9],mx1=mmm1raw[14],dqm=mean_tm1_dqm1)
  
  varmo_tm1_qm_exp1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[6],median=varmoraw[9],mx1=varmoraw[14],dqm=var_tm1_dqm1)
  
  tmmo_tm1_qm_exp1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[6],median=tmmoraw[9],mx1=tmmoraw[14],dqm=tm_tm1_dqm1)
  
  fmmo_tm1_qm_exp1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[6],median=fmmoraw[9],mx1=fmmoraw[14],dqm=fm_tm1_dqm1)
  
  rkurt_tm1_qm<-((fmmo_tm1_qm_exp1))/(varmo_tm1_qm_exp1^2)
  
  #Weibull
  
  
  if (((startpoint*(1-criterion))<rkurt_tm1_qm && rkurt_tm1_qm<(startpoint*(1)/(1-criterion)))){
    mmm1_tm1_qm_Weibull1<-mmm1_tm1_qm_exp1
    varmo_tm1_qm_Weibull1<-varmo_tm1_qm_exp1
    tmmo_tm1_qm_Weibull1<-tmmo_tm1_qm_exp1
    fmmo_tm1_qm_Weibull1<-fmmo_tm1_qm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_tm1_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm1_qm,dlist=standist_d,type="var_tm1_dqm")
      
      fm_tm1_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm1_qm,dlist=standist_d,type="fm_tm1_dqm")
      
      varmo_tm1_qm_Weibull1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[6],median=varmoraw[9],mx1=varmoraw[14],dqm=var_tm1_dqm1)
      
      fmmo_tm1_qm_Weibull1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[6],median=fmmoraw[9],mx1=fmmoraw[14],dqm=fm_tm1_dqm1)
      
      newrrkurt_tm1_qm<-((fmmo_tm1_qm_Weibull1))/(varmo_tm1_qm_Weibull1^2)
      
      if ((((rkurt_tm1_qm*(1-criterion))<newrrkurt_tm1_qm && newrrkurt_tm1_qm<(rkurt_tm1_qm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_tm1_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm1_qm,dlist=standist_d,type="mean_tm1_dqm")
        
        tm_tm1_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm1_qm,dlist=standist_d,type="tm_tm1_dqm")
        
        mmm1_tm1_qm_Weibull1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[6],median=mmm1raw[9],mx1=mmm1raw[14],dqm=mean_tm1_dqm1)
        
        tmmo_tm1_qm_Weibull1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[6],median=tmmoraw[9],mx1=tmmoraw[14],dqm=tm_tm1_dqm1)
        
        break
      }
      rkurt_tm1_qm<-newrrkurt_tm1_qm
      
    }
    
  }
  
  mean_tm2_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_tm2_drm")
  var_tm2_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_tm2_drm")
  tm_tm2_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_tm2_drm")
  fm_tm2_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_tm2_drm")
  
  
  mmm1_tm2_rm_exp1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[7],median=mmm1raw[9],mx1=mmm1raw[15],drm=mean_tm2_drm1)
  
  varmo_tm2_rm_exp1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[7],median=varmoraw[9],mx1=varmoraw[15],drm=var_tm2_drm1)
  
  tmmo_tm2_rm_exp1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[7],median=tmmoraw[9],mx1=tmmoraw[15],drm=tm_tm2_drm1)
  
  fmmo_tm2_rm_exp1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[7],median=fmmoraw[9],mx1=fmmoraw[15],drm=fm_tm2_drm1)
  
  
  rkurt_tm2_rm<-((fmmo_tm2_rm_exp1))/(varmo_tm2_rm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_tm2_rm && rkurt_tm2_rm<(startpoint*(1)/(1-criterion)))){
    mmm1_tm2_rm_Weibull1<-mmm1_tm2_rm_exp1
    varmo_tm2_rm_Weibull1<-varmo_tm2_rm_exp1
    tmmo_tm2_rm_Weibull1<-tmmo_tm2_rm_exp1
    fmmo_tm2_rm_Weibull1<-fmmo_tm2_rm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_tm2_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm2_rm,dlist=standist_d,type="var_tm2_drm")
      
      fm_tm2_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm2_rm,dlist=standist_d,type="fm_tm2_drm")
      
      varmo_tm2_rm_Weibull1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[7],median=varmoraw[9],mx1=varmoraw[15],drm=var_tm2_drm1)
      
      fmmo_tm2_rm_Weibull1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[7],median=fmmoraw[9],mx1=fmmoraw[15],drm=fm_tm2_drm1)
      
      newrrkurt_tm2_rm<-((fmmo_tm2_rm_Weibull1))/(varmo_tm2_rm_Weibull1^2)
      
      if ((((rkurt_tm2_rm*(1-criterion))<newrrkurt_tm2_rm && newrrkurt_tm2_rm<(rkurt_tm2_rm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_tm2_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm2_rm,dlist=standist_d,type="mean_tm2_drm")
        
        tm_tm2_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm2_rm,dlist=standist_d,type="tm_tm2_drm")
        
        mmm1_tm2_rm_Weibull1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[7],median=mmm1raw[9],mx1=mmm1raw[15],drm=mean_tm2_drm1)
        
        tmmo_tm2_rm_Weibull1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[7],median=tmmoraw[9],mx1=tmmoraw[15],drm=tm_tm2_drm1)
        
        break
      }
      rkurt_tm2_rm<-newrrkurt_tm2_rm
      
    }
    
  }
  
  
  mean_tm2_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_tm2_dqm")
  var_tm2_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_tm2_dqm")
  tm_tm2_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_tm2_dqm")
  fm_tm2_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_tm2_dqm")
  
  
  mmm1_tm2_qm_exp1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[7],median=mmm1raw[9],mx1=mmm1raw[15],dqm=mean_tm2_dqm1)
  
  varmo_tm2_qm_exp1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[7],median=varmoraw[9],mx1=varmoraw[15],dqm=var_tm2_dqm1)
  
  tmmo_tm2_qm_exp1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[7],median=tmmoraw[9],mx1=tmmoraw[15],dqm=tm_tm2_dqm1)
  
  fmmo_tm2_qm_exp1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[7],median=fmmoraw[9],mx1=fmmoraw[15],dqm=fm_tm2_dqm1)
  
  rkurt_tm2_qm<-((fmmo_tm2_qm_exp1))/(varmo_tm2_qm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_tm2_qm && rkurt_tm2_qm<(startpoint*(1)/(1-criterion)))){
    mmm1_tm2_qm_Weibull1<-mmm1_tm2_qm_exp1
    varmo_tm2_qm_Weibull1<-varmo_tm2_qm_exp1
    tmmo_tm2_qm_Weibull1<-tmmo_tm2_qm_exp1
    fmmo_tm2_qm_Weibull1<-fmmo_tm2_qm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_tm2_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm2_qm,dlist=standist_d,type="var_tm2_dqm")
      
      fm_tm2_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm2_qm,dlist=standist_d,type="fm_tm2_dqm")
      
      varmo_tm2_qm_Weibull1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[7],median=varmoraw[9],mx1=varmoraw[15],dqm=var_tm2_dqm1)
      
      fmmo_tm2_qm_Weibull1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[7],median=fmmoraw[9],mx1=fmmoraw[15],dqm=fm_tm2_dqm1)
      
      newrrkurt_tm2_qm<-((fmmo_tm2_qm_Weibull1))/(varmo_tm2_qm_Weibull1^2)
      
      if ((((rkurt_tm2_qm*(1-criterion))<newrrkurt_tm2_qm && newrrkurt_tm2_qm<(rkurt_tm2_qm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_tm2_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm2_qm,dlist=standist_d,type="mean_tm2_dqm")
        
        tm_tm2_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm2_qm,dlist=standist_d,type="tm_tm2_dqm")
        
        mmm1_tm2_qm_Weibull1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[7],median=mmm1raw[9],mx1=mmm1raw[15],dqm=mean_tm2_dqm1)
        
        tmmo_tm2_qm_Weibull1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[7],median=tmmoraw[9],mx1=tmmoraw[15],dqm=tm_tm2_dqm1)
        
        break
      }
      rkurt_tm2_qm<-newrrkurt_tm2_qm
      
    }
    
  }
  
  mean_tm3_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_tm3_drm")
  var_tm3_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_tm3_drm")
  tm_tm3_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_tm3_drm")
  fm_tm3_drm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_tm3_drm")
  
  mmm1_tm3_rm_exp1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[8],median=mmm1raw[9],mx1=mmm1raw[16],drm=mean_tm3_drm1)
  
  varmo_tm3_rm_exp1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[8],median=varmoraw[9],mx1=varmoraw[16],drm=var_tm3_drm1)
  
  tmmo_tm3_rm_exp1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[8],median=tmmoraw[9],mx1=tmmoraw[16],drm=tm_tm3_drm1)
  
  fmmo_tm3_rm_exp1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[8],median=fmmoraw[9],mx1=fmmoraw[16],drm=fm_tm3_drm1)
  
  
  rkurt_tm3_rm<-((fmmo_tm3_rm_exp1))/(varmo_tm3_rm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_tm3_rm && rkurt_tm3_rm<(startpoint*(1)/(1-criterion)))){
    mmm1_tm3_rm_Weibull1<-mmm1_tm3_rm_exp1
    varmo_tm3_rm_Weibull1<-varmo_tm3_rm_exp1
    tmmo_tm3_rm_Weibull1<-tmmo_tm3_rm_exp1
    fmmo_tm3_rm_Weibull1<-fmmo_tm3_rm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_tm3_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm3_rm,dlist=standist_d,type="var_tm3_drm")
      
      fm_tm3_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm3_rm,dlist=standist_d,type="fm_tm3_drm")
      
      varmo_tm3_rm_Weibull1<-mmmprocessrm(x=dp2varx,interval=interval,SWA=varmoraw[8],median=varmoraw[9],mx1=varmoraw[16],drm=var_tm3_drm1)
      
      fmmo_tm3_rm_Weibull1<-mmmprocessrm(x=dp4fmx,interval=interval,SWA=fmmoraw[8],median=fmmoraw[9],mx1=fmmoraw[16],drm=fm_tm3_drm1)
      
      newrrkurt_tm3_rm<-((fmmo_tm3_rm_Weibull1))/(varmo_tm3_rm_Weibull1^2)
      
      if ((((rkurt_tm3_rm*(1-criterion))<newrrkurt_tm3_rm && newrrkurt_tm3_rm<(rkurt_tm3_rm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_tm3_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm3_rm,dlist=standist_d,type="mean_tm3_drm")
        
        tm_tm3_drm1<-d_adjust(size=lengthx,kurt=rkurt_tm3_rm,dlist=standist_d,type="tm_tm3_drm")
        
        mmm1_tm3_rm_Weibull1<-mmmprocessrm(x=sortedx,interval=interval,SWA=mmm1raw[8],median=mmm1raw[9],mx1=mmm1raw[16],drm=mean_tm3_drm1)
        
        tmmo_tm3_rm_Weibull1<-mmmprocessrm(x=dp3tmx,interval=interval,SWA=tmmoraw[8],median=tmmoraw[9],mx1=tmmoraw[16],drm=tm_tm3_drm1)
        
        break
      }
      rkurt_tm3_rm<-newrrkurt_tm3_rm
      
    }
    
  }
  
  
  mean_tm3_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="mean_tm3_dqm")
  var_tm3_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="var_tm3_dqm")
  tm_tm3_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="tm_tm3_dqm")
  fm_tm3_dqm1=d_adjust(size=lengthx,kurt=startpoint,dlist=standist_d,type="fm_tm3_dqm")
  
  mmm1_tm3_qm_exp1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[8],median=mmm1raw[9],mx1=mmm1raw[16],dqm=mean_tm3_dqm1)
  
  varmo_tm3_qm_exp1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[8],median=varmoraw[9],mx1=varmoraw[16],dqm=var_tm3_dqm1)
  
  tmmo_tm3_qm_exp1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[8],median=tmmoraw[9],mx1=tmmoraw[16],dqm=tm_tm3_dqm1)
  
  fmmo_tm3_qm_exp1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[8],median=fmmoraw[9],mx1=fmmoraw[16],dqm=fm_tm3_dqm1)
  
  rkurt_tm3_qm<-((fmmo_tm3_qm_exp1))/(varmo_tm3_qm_exp1^2)
  
  #Weibull
  
  if (((startpoint*(1-criterion))<rkurt_tm3_qm && rkurt_tm3_qm<(startpoint*(1)/(1-criterion)))){
    mmm1_tm3_qm_Weibull1<-mmm1_tm3_qm_exp1
    varmo_tm3_qm_Weibull1<-varmo_tm3_qm_exp1
    tmmo_tm3_qm_Weibull1<-tmmo_tm3_qm_exp1
    fmmo_tm3_qm_Weibull1<-fmmo_tm3_qm_exp1
  }else{
    step1 <- 1
    repeat {
      step1 <-step1 + 1
      
      var_tm3_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm3_qm,dlist=standist_d,type="var_tm3_dqm")
      
      fm_tm3_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm3_qm,dlist=standist_d,type="fm_tm3_dqm")
      
      varmo_tm3_qm_Weibull1<-mmmprocessqm(x=dp2varx,interval=interval,SWA=varmoraw[8],median=varmoraw[9],mx1=varmoraw[16],dqm=var_tm3_dqm1)
      
      fmmo_tm3_qm_Weibull1<-mmmprocessqm(x=dp4fmx,interval=interval,SWA=fmmoraw[8],median=fmmoraw[9],mx1=fmmoraw[16],dqm=fm_tm3_dqm1)
      
      newrrkurt_tm3_qm<-((fmmo_tm3_qm_Weibull1))/(varmo_tm3_qm_Weibull1^2)
      
      if ((((rkurt_tm3_qm*(1-criterion))<newrrkurt_tm3_qm && newrrkurt_tm3_qm<(rkurt_tm3_qm*(1)/(1-criterion)))) || (step1 == stepsize)){
        
        mean_tm3_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm3_qm,dlist=standist_d,type="mean_tm3_dqm")
        
        tm_tm3_dqm1<-d_adjust(size=lengthx,kurt=rkurt_tm3_qm,dlist=standist_d,type="tm_tm3_dqm")
        
        mmm1_tm3_qm_Weibull1<-mmmprocessqm(x=sortedx,interval=interval,SWA=mmm1raw[8],median=mmm1raw[9],mx1=mmm1raw[16],dqm=mean_tm3_dqm1)
        
        tmmo_tm3_qm_Weibull1<-mmmprocessqm(x=dp3tmx,interval=interval,SWA=tmmoraw[8],median=tmmoraw[9],mx1=tmmoraw[16],dqm=tm_tm3_dqm1)
        
        break
      }
      rkurt_tm3_qm<-newrrkurt_tm3_qm
      
    }
    
  }
  
  
  sdall<-c(sdvar=sd(dp2varx),sdtm=sd(dp3tmx),sdfm=sd(dp4fmx))
  
  
  finallall<-c(HL1=HL1,mmm1raw=mmm1raw,mmm1_BM_rm_exp1=mmm1_BM_rm_exp1,mmm1_BM_rm_Weibull1=mmm1_BM_rm_Weibull1,mmm1_BM_qm_exp1=mmm1_BM_qm_exp1,mmm1_BM_qm_Weibull1=mmm1_BM_qm_Weibull1,
               mmm1_sqm_rm_exp1=mmm1_sqm_rm_exp1,mmm1_sqm_rm_Weibull1=mmm1_sqm_rm_Weibull1,mmm1_sqm_qm_exp1=mmm1_sqm_qm_exp1,mmm1_sqm_qm_Weibull1=mmm1_sqm_qm_Weibull1,
               mmm1_wm1_rm_exp1=mmm1_wm1_rm_exp1,mmm1_wm1_rm_Weibull1=mmm1_wm1_rm_Weibull1,mmm1_wm1_qm_exp1=mmm1_wm1_qm_exp1,mmm1_wm1_qm_Weibull1=mmm1_wm1_qm_Weibull1,
               mmm1_wm2_rm_exp1=mmm1_wm2_rm_exp1,mmm1_wm2_rm_Weibull1=mmm1_wm2_rm_Weibull1,mmm1_wm2_qm_exp1=mmm1_wm2_qm_exp1,mmm1_wm2_qm_Weibull1=mmm1_wm2_qm_Weibull1,
               mmm1_tm1_rm_exp1=mmm1_tm1_rm_exp1,mmm1_tm1_rm_Weibull1=mmm1_tm1_rm_Weibull1,mmm1_tm1_qm_exp1=mmm1_tm1_qm_exp1,mmm1_tm1_qm_Weibull1=mmm1_tm1_qm_Weibull1,
               mmm1_tm2_rm_exp1=mmm1_tm2_rm_exp1,mmm1_tm2_rm_Weibull1=mmm1_tm2_rm_Weibull1,mmm1_tm2_qm_exp1=mmm1_tm2_qm_exp1,mmm1_tm2_qm_Weibull1=mmm1_tm2_qm_Weibull1,
               mmm1_tm3_rm_exp1=mmm1_tm3_rm_exp1,mmm1_tm3_rm_Weibull1=mmm1_tm3_rm_Weibull1,mmm1_tm3_qm_exp1=mmm1_tm3_qm_exp1,mmm1_tm3_qm_Weibull1=mmm1_tm3_qm_Weibull1,
               
               varmoraw=varmoraw,varmo_BM_rm_exp1=varmo_BM_rm_exp1,varmo_BM_rm_Weibull1=varmo_BM_rm_Weibull1,varmo_BM_qm_exp1=varmo_BM_qm_exp1,varmo_BM_qm_Weibull1=varmo_BM_qm_Weibull1,
               varmo_sqm_rm_exp1=varmo_sqm_rm_exp1,varmo_sqm_rm_Weibull1=varmo_sqm_rm_Weibull1,varmo_sqm_qm_exp1=varmo_sqm_qm_exp1,varmo_sqm_qm_Weibull1=varmo_sqm_qm_Weibull1,
               varmo_wm1_rm_exp1=varmo_wm1_rm_exp1,varmo_wm1_rm_Weibull1=varmo_wm1_rm_Weibull1,varmo_wm1_qm_exp1=varmo_wm1_qm_exp1,varmo_wm1_qm_Weibull1=varmo_wm1_qm_Weibull1,
               varmo_wm2_rm_exp1=varmo_wm2_rm_exp1,varmo_wm2_rm_Weibull1=varmo_wm2_rm_Weibull1,varmo_wm2_qm_exp1=varmo_wm2_qm_exp1,varmo_wm2_qm_Weibull1=varmo_wm2_qm_Weibull1,
               varmo_tm1_rm_exp1=varmo_tm1_rm_exp1,varmo_tm1_rm_Weibull1=varmo_tm1_rm_Weibull1,varmo_tm1_qm_exp1=varmo_tm1_qm_exp1,varmo_tm1_qm_Weibull1=varmo_tm1_qm_Weibull1,
               varmo_tm2_rm_exp1=varmo_tm2_rm_exp1,varmo_tm2_rm_Weibull1=varmo_tm2_rm_Weibull1,varmo_tm2_qm_exp1=varmo_tm2_qm_exp1,varmo_tm2_qm_Weibull1=varmo_tm2_qm_Weibull1,
               varmo_tm3_rm_exp1=varmo_tm3_rm_exp1,varmo_tm3_rm_Weibull1=varmo_tm3_rm_Weibull1,varmo_tm3_qm_exp1=varmo_tm3_qm_exp1,varmo_tm3_qm_Weibull1=varmo_tm3_qm_Weibull1,
               
               tmmoraw=tmmoraw,tmmo_BM_rm_exp1=tmmo_BM_rm_exp1,tmmo_BM_rm_Weibull1=tmmo_BM_rm_Weibull1,tmmo_BM_qm_exp1=tmmo_BM_qm_exp1,tmmo_BM_qm_Weibull1=tmmo_BM_qm_Weibull1,
               tmmo_sqm_rm_exp1=tmmo_sqm_rm_exp1,tmmo_sqm_rm_Weibull1=tmmo_sqm_rm_Weibull1,tmmo_sqm_qm_exp1=tmmo_sqm_qm_exp1,tmmo_sqm_qm_Weibull1=tmmo_sqm_qm_Weibull1,
               tmmo_wm1_rm_exp1=tmmo_wm1_rm_exp1,tmmo_wm1_rm_Weibull1=tmmo_wm1_rm_Weibull1,tmmo_wm1_qm_exp1=tmmo_wm1_qm_exp1,tmmo_wm1_qm_Weibull1=tmmo_wm1_qm_Weibull1,
               tmmo_wm2_rm_exp1=tmmo_wm2_rm_exp1,tmmo_wm2_rm_Weibull1=tmmo_wm2_rm_Weibull1,tmmo_wm2_qm_exp1=tmmo_wm2_qm_exp1,tmmo_wm2_qm_Weibull1=tmmo_wm2_qm_Weibull1,
               tmmo_tm1_rm_exp1=tmmo_tm1_rm_exp1,tmmo_tm1_rm_Weibull1=tmmo_tm1_rm_Weibull1,tmmo_tm1_qm_exp1=tmmo_tm1_qm_exp1,tmmo_tm1_qm_Weibull1=tmmo_tm1_qm_Weibull1,
               tmmo_tm2_rm_exp1=tmmo_tm2_rm_exp1,tmmo_tm2_rm_Weibull1=tmmo_tm2_rm_Weibull1,tmmo_tm2_qm_exp1=tmmo_tm2_qm_exp1,tmmo_tm2_qm_Weibull1=tmmo_tm2_qm_Weibull1,
               tmmo_tm3_rm_exp1=tmmo_tm3_rm_exp1,tmmo_tm3_rm_Weibull1=tmmo_tm3_rm_Weibull1,tmmo_tm3_qm_exp1=tmmo_tm3_qm_exp1,tmmo_tm3_qm_Weibull1=tmmo_tm3_qm_Weibull1,
               
               fmmoraw=fmmoraw,fmmo_BM_rm_exp1=fmmo_BM_rm_exp1,fmmo_BM_rm_Weibull1=fmmo_BM_rm_Weibull1,fmmo_BM_qm_exp1=fmmo_BM_qm_exp1,fmmo_BM_qm_Weibull1=fmmo_BM_qm_Weibull1,
               fmmo_sqm_rm_exp1=fmmo_sqm_rm_exp1,fmmo_sqm_rm_Weibull1=fmmo_sqm_rm_Weibull1,fmmo_sqm_qm_exp1=fmmo_sqm_qm_exp1,fmmo_sqm_qm_Weibull1=fmmo_sqm_qm_Weibull1,
               fmmo_wm1_rm_exp1=fmmo_wm1_rm_exp1,fmmo_wm1_rm_Weibull1=fmmo_wm1_rm_Weibull1,fmmo_wm1_qm_exp1=fmmo_wm1_qm_exp1,fmmo_wm1_qm_Weibull1=fmmo_wm1_qm_Weibull1,
               fmmo_wm2_rm_exp1=fmmo_wm2_rm_exp1,fmmo_wm2_rm_Weibull1=fmmo_wm2_rm_Weibull1,fmmo_wm2_qm_exp1=fmmo_wm2_qm_exp1,fmmo_wm2_qm_Weibull1=fmmo_wm2_qm_Weibull1,
               fmmo_tm1_rm_exp1=fmmo_tm1_rm_exp1,fmmo_tm1_rm_Weibull1=fmmo_tm1_rm_Weibull1,fmmo_tm1_qm_exp1=fmmo_tm1_qm_exp1,fmmo_tm1_qm_Weibull1=fmmo_tm1_qm_Weibull1,
               fmmo_tm2_rm_exp1=fmmo_tm2_rm_exp1,fmmo_tm2_rm_Weibull1=fmmo_tm2_rm_Weibull1,fmmo_tm2_qm_exp1=fmmo_tm2_qm_exp1,fmmo_tm2_qm_Weibull1=fmmo_tm2_qm_Weibull1,
               fmmo_tm3_rm_exp1=fmmo_tm3_rm_exp1,fmmo_tm3_rm_Weibull1=fmmo_tm3_rm_Weibull1,fmmo_tm3_qm_exp1=fmmo_tm3_qm_exp1,fmmo_tm3_qm_Weibull1=fmmo_tm3_qm_Weibull1,
               
               sdall=sdall)
  #finallall<-c(HL1=HL1,mmm1raw=mmm1raw,mmm1exp1=mmm1exp1,mmm1exp2=mmm1exp2,mmm1Weibull1=mmm1Weibull1,mmm1Weibull2=mmm1Weibull2,varmoraw=varmoraw,varmoexp1=varmoexp1,varmoexp2=varmoexp2,varmoWeibull1=varmoWeibull1,varmoWeibull2=varmoWeibull2,tmmoraw=tmmoraw,tmmoexp1=tmmoexp1,tmmoexp2=tmmoexp2,tmmoWeibull1=tmmoWeibull1,tmmoWeibull2=tmmoWeibull2,fmmoraw=fmmoraw,fmmoexp1=fmmoexp1,fmmoexp2=fmmoexp2,fmmoWeibull1=fmmoWeibull1,fmmoWeibull2=fmmoWeibull2,sdall=sdall)
  return(finallall)
}

#set the convergence criterion
criterionset=1/20

kurtWeibull<- read.csv(("kurtWeibull_28260.csv"))

allkurtWeibull<-unlist(kurtWeibull)

samplesize=1.8*10^6
batchsizebase=30
orderlist1_AB2<-removelist(na.omit(t(apply(quasiuni_sorted2,MARGIN=1,FUN=roundunique,dimension=2,size=samplesize))))
orderlist1_AB3<-removelist(na.omit(t(apply(quasiuni_sorted3,MARGIN=1,FUN=roundunique,dimension=3,size=samplesize))))
orderlist1_AB4<-removelist(na.omit(t(apply(quasiuni_sorted4,MARGIN=1,FUN=roundunique,dimension=4,size=samplesize))))

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
  
  batchsize=30#ceiling(batchsizebase/(kurtx^(1/4)))
  
  n <- samplesize
  
  unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)
  
  unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
    
    x<-c(dsWeibull(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    onestepx<-onestep(x=sortedx, bend = 1.172)
    SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
    dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
    momentsx<-unbiasedmoments(x=sortedx)
    sortedx<-c()
    momentssd<-c(sqrt(momentsx[2]),dall[178],dall[179],dall[180])
    allrawmoBias<-c(
      firstbias=(c(onestepx,SMWM9,dall[1:45])-targetm),
      secondbias=(c(dall[46:89])-targetvar),
      thirdbias=(c(dall[90:133])-targettm),
      fourbias=(c(dall[134:177])-targetfm))
    
    all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,momentssd,allrawmoBias))
  
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("Weibull_raw_w_calibration_asymptotic",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  SEbatachesmean <- colMeans(SEbataches)
  
  AB_mexp<-(which(abs(SEbatachesmean[seq(from=196+21, to=196+35, by=2)])==(min(abs(SEbatachesmean[seq(from=196+21, to=196+35, by=2)])))))[1]
  AB_mWeibull<-(which(abs(SEbatachesmean[seq(from=196+22, to=196+36, by=2)])==(min(abs(SEbatachesmean[seq(from=196+22, to=196+36, by=2)])))))[1]
  
  AB_varexp<-(which(abs(SEbatachesmean[seq(from=196+65, to=196+80, by=2)])==(min(abs(SEbatachesmean[seq(from=196+65, to=196+80, by=2)])))))[1]
  AB_varWeibull<-(which(abs(SEbatachesmean[seq(from=196+66, to=196+80, by=2)])==(min(abs(SEbatachesmean[seq(from=196+66, to=196+80, by=2)])))))[1]
  
  AB_tmexp<-(which(abs(SEbatachesmean[seq(from=196+109, to=196+124, by=2)])==(min(abs(SEbatachesmean[seq(from=196+109, to=196+124, by=2)])))))[1]
  AB_tmWeibull<-(which(abs(SEbatachesmean[seq(from=196+110, to=196+124, by=2)])==(min(abs(SEbatachesmean[seq(from=196+110, to=196+124, by=2)])))))[1]
  
  AB_fmexp<-(which(abs(SEbatachesmean[seq(from=196+153, to=196+168, by=2)])==(min(abs(SEbatachesmean[seq(from=196+153, to=196+168, by=2)])))))[1]
  AB_fmWeibull<-(which(abs(SEbatachesmean[seq(from=196+154, to=196+168, by=2)])==(min(abs(SEbatachesmean[seq(from=196+154, to=196+168, by=2)])))))[1]
  
  ABrank1<-c(AB_mexp=AB_mexp,AB_mWeibull=AB_mWeibull,AB_varexp=AB_varexp,AB_varWeibull=AB_varWeibull,AB_tmexp=AB_tmexp,AB_tmWeibull=AB_tmWeibull,AB_fmexp=AB_fmexp,AB_fmWeibull=AB_fmWeibull)
  
  ratiomean1<-c(SEbatachesmean[2:49])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:49]), 2, unbiasedsd)
  
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:49])/ratiomean1)
  
  meansd1<-apply(mean_SEbatachesmeanprocess, 2, unbiasedsd)
  ratiovar1<-SEbatachesmean[50:93]/SEbatachesmean[50]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,50:93]), 2, unbiasedsd)
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,50:93])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  ratiotm1<-SEbatachesmean[94:137]/SEbatachesmean[94]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,94:137]), 2, unbiasedsd)
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,94:137])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  
  ratiofm1<-SEbatachesmean[138:181]/SEbatachesmean[138]
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,138:181]), 2, unbiasedsd)
  
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,138:181])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  
  allSE_unstan<-c(meansd_unscaled1=meansd_unscaled1,varsd_unscaled1=varsd_unscaled1,tmsd_unscaled1=tmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE_unstand<-c(meansd1=meansd1,varsd1=varsd1,tmsd1=tmsd1,fmsd1=fmsd1
  )
  
  mexp<-which(meansd1[seq(from=21, to=36, by=2)]==(min(meansd1[seq(from=21, to=36, by=2)])))[1]
  mWeibull<-which(meansd1[seq(from=22, to=37, by=2)]==(min(meansd1[seq(from=22, to=37, by=2)])))[1]
  varexp<-which(varsd1[seq(from=17, to=32, by=2)]==(min(varsd1[seq(from=17, to=32, by=2)])))[1]
  varWeibull<-which(varsd1[seq(from=18, to=33, by=2)]==(min(varsd1[seq(from=18, to=33, by=2)])))[1]
  tmexp<-which(tmsd1[seq(from=17, to=32, by=2)]==(min(tmsd1[seq(from=17, to=32, by=2)])))[1]
  tmWeibull<-which(tmsd1[seq(from=18, to=33, by=2)]==(min(tmsd1[seq(from=18, to=33, by=2)])))[1]
  fmexp<-which(fmsd1[seq(from=17, to=32, by=2)]==(min(fmsd1[seq(from=17, to=32, by=2)])))[1]
  fmWeibull<-which(fmsd1[seq(from=18, to=33, by=2)]==(min(fmsd1[seq(from=18, to=33, by=2)])))[1]
  
  rank1<-c(mexp=mexp,mWeibull=mWeibull,varexp=varexp,varWeibull=varWeibull,tmexp=tmexp,tmWeibull=tmWeibull,fmexp=fmexp,fmWeibull=fmWeibull)
  
  allresultsSE<-c(samplesize,SEbatachesmean,allSE_unstan,allSSE_unstand,ABrank1,rank1)
  
  allresultsSE
}

write.csv(simulatedbatchWeibull_bias_Monte,paste("Weibull_ABSE_Monte_asymptotic.csv", sep = ","), row.names = FALSE)

Label_ABSE_Weibull1<- read.csv(("ABSSE_w_Weibull_SWA_finite_Label.csv"))

Optimum_ABSE<-cbind(simulatedbatchWeibull_bias_Monte[,1:2],simulatedbatchWeibull_bias_Monte[,738:ncol(simulatedbatchWeibull_bias_Monte)])

colnames(Optimum_ABSE)<-colnames(Label_ABSE_Weibull1)

write.csv(Optimum_ABSE,paste("ABSSE_w_Weibull_SWA_asymptotic.csv", sep = ","), row.names = FALSE)

