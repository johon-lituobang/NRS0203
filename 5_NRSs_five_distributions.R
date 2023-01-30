
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

asymptotic_n <- 1.8*10^6
(asymptotic_n%%10)==0
# maximum order of moments
morder <- 4
#large sample size (asymptotic bias)
largesize<-1.8*10^6

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol_asymptotic<-sobol(n=asymptotic_n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                mixed = FALSE, method = "C", start = 1)

quasiuni_asymptotic<-rbind(quasiunisobol_asymptotic)

quasiunisobol_asymptotic<-c()

quasiuni_sorted2_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4_asymptotic <- na.omit(rowSort(quasiuni_asymptotic, descend = FALSE, stable = FALSE, parallel = TRUE))
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


quasiuni_sorted2_asymptotic<-c()
quasiuni_sorted3_asymptotic<-c()
quasiuni_sorted4_asymptotic<-c()

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
#load asymptotic d for two parameter distributions
Weibull_d<- read.csv(("d_value_Weibull.csv"))
Weibull_w_ABSSE<- read.csv(("ABSSE_w_Weibull_SWA.csv"))
# #load asymptotic d for two parameter distributions
# gamma_d<- read.csv(("asymptotic_d_gamma319.csv"))
# #load asymptotic d for two parameter distributions
# lognormal_d<- read.csv(("asymptotic_d_lognorm319.csv"))
# # #load asymptotic d for two parameter distributions
# Pareto_d<- read.csv(("asymptotic_d_Pareto919.csv"))


w_adjust<-function(size,kurt,wlist,type){
  wlist[,2]<-round(wlist[,2],digits=1)
  indextypelist<-c("mean_AB_exp","var_AB_exp","tm_AB_exp","fm_AB_exp","mean_AB_Weibull","var_AB_Weibull","tm_AB_Weibull","fm_AB_Weibull","mean_SSE_exp","var_SSE_exp","tm_SSE_exp","fm_SSE_exp","mean_SSE_Weibull","var_SSE_Weibull","tm_SSE_Weibull","fm_SSE_Weibull")
  indextype<-which(indextypelist==(type))+2
  if(size%in% wlist[,1]){
    if (kurt%in% wlist[,2]){
      result1<-round(wlist[wlist[,1] == size & wlist[,2]==kurt,indextype],digits = 0)
    }else{
      rown2<-as.numeric(wlist[,2])
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
      
      d1<-wlist[wlist[,1]==size & wlist[,2]==infn2,indextype]
      d2<-wlist[wlist[,1]==size & wlist[,2]==supn2,indextype]
      
      
      if(d1==d2){
        result1<-round(d1,digits = 0)
      }else{
        result1<-round((((d2-d1)*((kurt-infn2)/(supn2-infn2)))+d1),digits = 0)
      }
    }
  }else{
    rown<-as.numeric(wlist[,1])
    infn<-max(rown[rown < size])
    if(is.infinite(infn)){
      infn<-rown[1]
    }
    supn<-min(rown[rown > size])
    if(is.infinite(supn)){
      supn<-max(rown)
    }
    
    rown2<-as.numeric(wlist[,2])
    infn2<-max(rown2[rown2 < kurt])
    if(is.infinite(infn2)){
      infn2<-rown2[2]
    }
    supn2<-min(rown2[rown2 > kurt])
    if(is.infinite(supn2)){
      supn2<-max(rown2)
    }
    
    d1<-wlist[wlist[,1]==infn & wlist[,2]==infn2,indextype]
    d2<-wlist[wlist[,1]==supn & wlist[,2]==supn2,indextype]
    
    
    if(d1==d2){
      result1<-round(d1,digits = 0)
    }else if (supn!=infn){
      result1<-round((((d2-d1)*((size-infn)/(supn-infn)))+d1),digits = 0)
    }else{
      result1<-round((((d2-d1)*((kurt-infn2)/(supn2-infn2)))+d1),digits = 0)
    }
    
  }
  if(result1<1){
    result1<-1
  }
  if(result1>8){
    result1<-8
  }
  return(result1)
}

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
  meanallselectexp<-c(mmm1_BM_rm_exp1=mmm1_BM_rm_exp1,mmm1_BM_qm_exp1=mmm1_BM_qm_exp1,
                      mmm1_sqm_rm_exp1=mmm1_sqm_rm_exp1,mmm1_sqm_qm_exp1=mmm1_sqm_qm_exp1,
                      mmm1_wm1_rm_exp1=mmm1_wm1_rm_exp1,mmm1_wm1_qm_exp1=mmm1_wm1_qm_exp1,
                      mmm1_wm2_rm_exp1=mmm1_wm2_rm_exp1,mmm1_wm2_qm_exp1=mmm1_wm2_qm_exp1#,
                      #mmm1_tm1_rm_exp1=mmm1_tm1_rm_exp1,mmm1_tm1_qm_exp1=mmm1_tm1_qm_exp1
  )
  varallselectexp<-c(
    varmo_BM_rm_exp1=varmo_BM_rm_exp1,varmo_BM_qm_exp1=varmo_BM_qm_exp1,
    varmo_sqm_rm_exp1=varmo_sqm_rm_exp1,varmo_sqm_qm_exp1=varmo_sqm_qm_exp1,
    varmo_wm1_rm_exp1=varmo_wm1_rm_exp1,varmo_wm1_qm_exp1=varmo_wm1_qm_exp1,
    varmo_wm2_rm_exp1=varmo_wm2_rm_exp1,varmo_wm2_qm_exp1=varmo_wm2_qm_exp1#,
    #varmo_tm1_rm_exp1=varmo_tm1_rm_exp1,varmo_tm1_qm_exp1=varmo_tm1_qm_exp1
  )
  
  tmallselectexp<-c(             
    tmmo_BM_rm_exp1=tmmo_BM_rm_exp1,tmmo_BM_qm_exp1=tmmo_BM_qm_exp1,
    tmmo_sqm_rm_exp1=tmmo_sqm_rm_exp1,tmmo_sqm_qm_exp1=tmmo_sqm_qm_exp1,
    tmmo_wm1_rm_exp1=tmmo_wm1_rm_exp1,tmmo_wm1_qm_exp1=tmmo_wm1_qm_exp1,
    tmmo_wm2_rm_exp1=tmmo_wm2_rm_exp1,tmmo_wm2_qm_exp1=tmmo_wm2_qm_exp1#,
    #tmmo_tm1_rm_exp1=tmmo_tm1_rm_exp1,tmmo_tm1_qm_exp1,fmmo_tm1_qm_exp1
  )
  
  fmallselectexp<-c(                          
    fmmo_BM_rm_exp1=fmmo_BM_rm_exp1,fmmo_BM_qm_exp1=fmmo_BM_qm_exp1,
    fmmo_sqm_rm_exp1=fmmo_sqm_rm_exp1,fmmo_sqm_qm_exp1=fmmo_sqm_qm_exp1,
    fmmo_wm1_rm_exp1=fmmo_wm1_rm_exp1,fmmo_wm1_qm_exp1=fmmo_wm1_qm_exp1,
    fmmo_wm2_rm_exp1=fmmo_wm2_rm_exp1,fmmo_wm2_qm_exp1=fmmo_wm2_qm_exp1#,
    #fmmo_tm1_rm_exp1=fmmo_tm1_rm_exp1,fmmo_tm1_qm_exp1=fmmo_tm1_qm_exp1
  )
  
  
  meanallselectWeibull<-c(mmm1_BM_rm_Weibull1=mmm1_BM_rm_Weibull1,mmm1_BM_qm_Weibull1=mmm1_BM_qm_Weibull1,
                          mmm1_sqm_rm_Weibull1=mmm1_sqm_rm_Weibull1,mmm1_sqm_qm_Weibull1=mmm1_sqm_qm_Weibull1,
                          mmm1_wm1_rm_Weibull1=mmm1_wm1_rm_Weibull1,mmm1_wm1_qm_Weibull1=mmm1_wm1_qm_Weibull1,
                          mmm1_wm2_rm_Weibull1=mmm1_wm2_rm_Weibull1,mmm1_wm2_qm_Weibull1=mmm1_wm2_qm_Weibull1#,
                          #mmm1_tm1_rm_Weibull1=mmm1_tm1_rm_Weibull1,mmm1_tm1_qm_Weibull1=mmm1_tm1_qm_Weibull1
  )
  varallselectWeibull<-c(
    
    varmo_BM_rm_Weibull1=varmo_BM_rm_Weibull1,varmo_BM_qm_Weibull1=varmo_BM_qm_Weibull1,
    varmo_sqm_rm_Weibull1=varmo_sqm_rm_Weibull1,varmo_sqm_qm_Weibull1=varmo_sqm_qm_Weibull1,
    varmo_wm1_rm_Weibull1=varmo_wm1_rm_Weibull1,varmo_wm1_qm_Weibull1=varmo_wm1_qm_Weibull1,
    varmo_wm2_rm_Weibull1=varmo_wm2_rm_Weibull1,varmo_wm2_qm_Weibull1=varmo_wm2_qm_Weibull1#,
    #varmo_tm1_rm_Weibull1=varmo_tm1_rm_Weibull1,varmo_tm1_qm_Weibull1=varmo_tm1_qm_Weibull1
  )
  
  tmallselectWeibull<-c(             
    tmmo_BM_rm_Weibull1=tmmo_BM_rm_Weibull1,tmmo_BM_qm_Weibull1=tmmo_BM_qm_Weibull1,
    tmmo_sqm_rm_Weibull1=tmmo_sqm_rm_Weibull1,tmmo_sqm_qm_Weibull1=tmmo_sqm_qm_Weibull1,
    tmmo_wm1_rm_Weibull1=tmmo_wm1_rm_Weibull1,tmmo_wm1_qm_Weibull1=tmmo_wm1_qm_Weibull1,
    tmmo_wm2_rm_Weibull1=tmmo_wm2_rm_Weibull1,tmmo_wm2_qm_Weibull1=tmmo_wm2_qm_Weibull1#,
    #tmmo_tm1_rm_Weibull1=tmmo_tm1_rm_Weibull1,tmmo_tm1_qm_Weibull1=tmmo_tm1_qm_Weibull1
  )
  
  fmallselectWeibull<-c(                          
    fmmo_BM_rm_Weibull1=fmmo_BM_rm_Weibull1,fmmo_BM_qm_Weibull1=fmmo_BM_qm_Weibull1,
    fmmo_sqm_rm_Weibull1=fmmo_sqm_rm_Weibull1,fmmo_sqm_qm_Weibull1=fmmo_sqm_qm_Weibull1,
    fmmo_wm1_rm_Weibull1=fmmo_wm1_rm_Weibull1,fmmo_wm1_qm_Weibull1=fmmo_wm1_qm_Weibull1,
    fmmo_wm2_rm_Weibull1=fmmo_wm2_rm_Weibull1,fmmo_wm2_qm_Weibull1=fmmo_wm2_qm_Weibull1#,
    #fmmo_tm1_rm_Weibull1=fmmo_tm1_rm_Weibull1,fmmo_tm1_qm_Weibull1=fmmo_tm1_qm_Weibull1
  )
  
  standist_w_ABSSE=Weibull_w_ABSSE
  
  #exp
  kurt_adaptive_AB_exp<-(((fmmo_BM_rm_exp1))/(varmo_BM_rm_exp1^2))
  step1 <- 1
  repeat {
    step1 <-step1 + 1
    
    var_adaptive_AB_exp<-w_adjust(size=lengthx,kurt=kurt_adaptive_AB_exp,wlist=standist_w_ABSSE,type="var_AB_exp")
    
    fm_adaptive_AB_exp<-w_adjust(size=lengthx,kurt=kurt_adaptive_AB_exp,wlist=standist_w_ABSSE,type="fm_AB_exp")
    
    varmo_adaptive_AB_exp<-varallselectexp[var_adaptive_AB_exp]
    
    fmmo_adaptive_AB_exp<-fmallselectexp[fm_adaptive_AB_exp]
    
    newkurt_adaptive_AB_exp<-((fmmo_adaptive_AB_exp))/(varmo_adaptive_AB_exp^2)
    
    if ((((kurt_adaptive_AB_exp*(1-criterion))<newkurt_adaptive_AB_exp && newkurt_adaptive_AB_exp<(kurt_adaptive_AB_exp*(1)/(1-criterion)))) || (step1 == stepsize)){
      
      mean_adaptive_AB_exp<-w_adjust(size=lengthx,kurt=kurt_adaptive_AB_exp,wlist=standist_w_ABSSE,type="mean_AB_exp")
      
      tm_adaptive_AB_exp<-w_adjust(size=lengthx,kurt=kurt_adaptive_AB_exp,wlist=standist_w_ABSSE,type="tm_AB_exp")
      
      meanmo_adaptive_AB_exp<-meanallselectexp[mean_adaptive_AB_exp]
      
      tmmo_adaptive_AB_exp<-tmallselectexp[tm_adaptive_AB_exp]
      
      break
    }
    kurt_adaptive_AB_exp<-newkurt_adaptive_AB_exp
    
  }
  
  #Weibull
  kurt_adaptive_AB_Weibull<-(((fmmo_BM_rm_exp1))/(varmo_BM_rm_exp1^2))
  step1 <- 1
  repeat {
    step1 <-step1 + 1
    
    var_adaptive_AB_Weibull<-w_adjust(size=lengthx,kurt=kurt_adaptive_AB_Weibull,wlist=standist_w_ABSSE,type="var_AB_Weibull")
    
    fm_adaptive_AB_Weibull<-w_adjust(size=lengthx,kurt=kurt_adaptive_AB_Weibull,wlist=standist_w_ABSSE,type="fm_AB_Weibull")
    
    varmo_adaptive_AB_Weibull<-varallselectWeibull[var_adaptive_AB_Weibull]
    
    fmmo_adaptive_AB_Weibull<-fmallselectWeibull[fm_adaptive_AB_Weibull]
    
    newkurt_adaptive_AB_Weibull<-((fmmo_adaptive_AB_Weibull))/(varmo_adaptive_AB_Weibull^2)
    
    if ((((kurt_adaptive_AB_Weibull*(1-criterion))<newkurt_adaptive_AB_Weibull && newkurt_adaptive_AB_Weibull<(kurt_adaptive_AB_Weibull*(1)/(1-criterion)))) || (step1 == stepsize)){
      
      mean_adaptive_AB_Weibull<-w_adjust(size=lengthx,kurt=kurt_adaptive_AB_Weibull,wlist=standist_w_ABSSE,type="mean_AB_Weibull")
      
      tm_adaptive_AB_Weibull<-w_adjust(size=lengthx,kurt=kurt_adaptive_AB_Weibull,wlist=standist_w_ABSSE,type="tm_AB_Weibull")
      
      meanmo_adaptive_AB_Weibull<-meanallselectWeibull[mean_adaptive_AB_Weibull]
      
      tmmo_adaptive_AB_Weibull<-tmallselectWeibull[tm_adaptive_AB_Weibull]
      
      break
    }
    kurt_adaptive_AB_Weibull<-newkurt_adaptive_AB_Weibull
    
  }
  
  
  #exp
  kurt_adaptive_SSE_exp<-(((fmmo_BM_rm_exp1))/(varmo_BM_rm_exp1^2))
  step1 <- 1
  repeat {
    step1 <-step1 + 1
    
    var_adaptive_SSE_exp<-w_adjust(size=lengthx,kurt=kurt_adaptive_SSE_exp,wlist=standist_w_ABSSE,type="var_SSE_exp")
    
    fm_adaptive_SSE_exp<-w_adjust(size=lengthx,kurt=kurt_adaptive_SSE_exp,wlist=standist_w_ABSSE,type="fm_SSE_exp")
    
    varmo_adaptive_SSE_exp<-varallselectexp[var_adaptive_SSE_exp]
    
    fmmo_adaptive_SSE_exp<-fmallselectexp[fm_adaptive_SSE_exp]
    
    newkurt_adaptive_SSE_exp<-((fmmo_adaptive_SSE_exp))/(varmo_adaptive_SSE_exp^2)
    
    if ((((kurt_adaptive_SSE_exp*(1-criterion))<newkurt_adaptive_SSE_exp && newkurt_adaptive_SSE_exp<(kurt_adaptive_SSE_exp*(1)/(1-criterion)))) || (step1 == stepsize)){
      
      mean_adaptive_SSE_exp<-w_adjust(size=lengthx,kurt=kurt_adaptive_SSE_exp,wlist=standist_w_ABSSE,type="mean_SSE_exp")
      
      tm_adaptive_SSE_exp<-w_adjust(size=lengthx,kurt=kurt_adaptive_SSE_exp,wlist=standist_w_ABSSE,type="tm_SSE_exp")
      
      meanmo_adaptive_SSE_exp<-meanallselectexp[mean_adaptive_SSE_exp]
      
      tmmo_adaptive_SSE_exp<-tmallselectexp[tm_adaptive_SSE_exp]
      
      break
    }
    kurt_adaptive_SSE_exp<-newkurt_adaptive_SSE_exp
    
  }
  
  #Weibull
  kurt_adaptive_SSE_Weibull<-(((fmmo_BM_rm_exp1))/(varmo_BM_rm_exp1^2))
  step1 <- 1
  repeat {
    step1 <-step1 + 1
    
    var_adaptive_SSE_Weibull<-w_adjust(size=lengthx,kurt=kurt_adaptive_SSE_Weibull,wlist=standist_w_ABSSE,type="var_SSE_Weibull")
    
    fm_adaptive_SSE_Weibull<-w_adjust(size=lengthx,kurt=kurt_adaptive_SSE_Weibull,wlist=standist_w_ABSSE,type="fm_SSE_Weibull")
    
    varmo_adaptive_SSE_Weibull<-varallselectWeibull[var_adaptive_SSE_Weibull]
    
    fmmo_adaptive_SSE_Weibull<-fmallselectWeibull[fm_adaptive_SSE_Weibull]
    
    newkurt_adaptive_SSE_Weibull<-((fmmo_adaptive_SSE_Weibull))/(varmo_adaptive_SSE_Weibull^2)
    
    if ((((kurt_adaptive_SSE_Weibull*(1-criterion))<newkurt_adaptive_SSE_Weibull && newkurt_adaptive_SSE_Weibull<(kurt_adaptive_SSE_Weibull*(1)/(1-criterion)))) || (step1 == stepsize)){
      
      mean_adaptive_SSE_Weibull<-w_adjust(size=lengthx,kurt=kurt_adaptive_SSE_Weibull,wlist=standist_w_ABSSE,type="mean_SSE_Weibull")
      
      tm_adaptive_SSE_Weibull<-w_adjust(size=lengthx,kurt=kurt_adaptive_SSE_Weibull,wlist=standist_w_ABSSE,type="tm_SSE_Weibull")
      
      meanmo_adaptive_SSE_Weibull<-meanallselectWeibull[mean_adaptive_SSE_Weibull]
      
      tmmo_adaptive_SSE_Weibull<-tmallselectWeibull[tm_adaptive_SSE_Weibull]
      
      break
    }
    kurt_adaptive_SSE_Weibull<-newkurt_adaptive_SSE_Weibull
    
  }
  finallall<-c(HL1=HL1,mmm1raw=mmm1raw,mmm1_BM_rm_exp1=mmm1_BM_rm_exp1,mmm1_BM_rm_Weibull1=mmm1_BM_rm_Weibull1,mmm1_BM_qm_exp1=mmm1_BM_qm_exp1,mmm1_BM_qm_Weibull1=mmm1_BM_qm_Weibull1,
               mmm1_sqm_rm_exp1=mmm1_sqm_rm_exp1,mmm1_sqm_rm_Weibull1=mmm1_sqm_rm_Weibull1,mmm1_sqm_qm_exp1=mmm1_sqm_qm_exp1,mmm1_sqm_qm_Weibull1=mmm1_sqm_qm_Weibull1,
               mmm1_wm1_rm_exp1=mmm1_wm1_rm_exp1,mmm1_wm1_rm_Weibull1=mmm1_wm1_rm_Weibull1,mmm1_wm1_qm_exp1=mmm1_wm1_qm_exp1,mmm1_wm1_qm_Weibull1=mmm1_wm1_qm_Weibull1,
               mmm1_wm2_rm_exp1=mmm1_wm2_rm_exp1,mmm1_wm2_rm_Weibull1=mmm1_wm2_rm_Weibull1,mmm1_wm2_qm_exp1=mmm1_wm2_qm_exp1,mmm1_wm2_qm_Weibull1=mmm1_wm2_qm_Weibull1,
               mmm1_tm1_rm_exp1=mmm1_tm1_rm_exp1,mmm1_tm1_rm_Weibull1=mmm1_tm1_rm_Weibull1,mmm1_tm1_qm_exp1=mmm1_tm1_qm_exp1,mmm1_tm1_qm_Weibull1=mmm1_tm1_qm_Weibull1,
               mmm1_tm2_rm_exp1=mmm1_tm2_rm_exp1,mmm1_tm2_rm_Weibull1=mmm1_tm2_rm_Weibull1,mmm1_tm2_qm_exp1=mmm1_tm2_qm_exp1,mmm1_tm2_qm_Weibull1=mmm1_tm2_qm_Weibull1,
               mmm1_tm3_rm_exp1=mmm1_tm3_rm_exp1,mmm1_tm3_rm_Weibull1=mmm1_tm3_rm_Weibull1,mmm1_tm3_qm_exp1=mmm1_tm3_qm_exp1,mmm1_tm3_qm_Weibull1=mmm1_tm3_qm_Weibull1,
               meanmo_adaptive_AB_exp=meanmo_adaptive_AB_exp,meanmo_adaptive_AB_Weibull=meanmo_adaptive_AB_Weibull,meanmo_adaptive_SSE_exp=meanmo_adaptive_SSE_exp,meanmo_adaptive_SSE_Weibull=meanmo_adaptive_SSE_Weibull,
               varmoraw=varmoraw,varmo_BM_rm_exp1=varmo_BM_rm_exp1,varmo_BM_rm_Weibull1=varmo_BM_rm_Weibull1,varmo_BM_qm_exp1=varmo_BM_qm_exp1,varmo_BM_qm_Weibull1=varmo_BM_qm_Weibull1,
               varmo_sqm_rm_exp1=varmo_sqm_rm_exp1,varmo_sqm_rm_Weibull1=varmo_sqm_rm_Weibull1,varmo_sqm_qm_exp1=varmo_sqm_qm_exp1,varmo_sqm_qm_Weibull1=varmo_sqm_qm_Weibull1,
               varmo_wm1_rm_exp1=varmo_wm1_rm_exp1,varmo_wm1_rm_Weibull1=varmo_wm1_rm_Weibull1,varmo_wm1_qm_exp1=varmo_wm1_qm_exp1,varmo_wm1_qm_Weibull1=varmo_wm1_qm_Weibull1,
               varmo_wm2_rm_exp1=varmo_wm2_rm_exp1,varmo_wm2_rm_Weibull1=varmo_wm2_rm_Weibull1,varmo_wm2_qm_exp1=varmo_wm2_qm_exp1,varmo_wm2_qm_Weibull1=varmo_wm2_qm_Weibull1,
               varmo_tm1_rm_exp1=varmo_tm1_rm_exp1,varmo_tm1_rm_Weibull1=varmo_tm1_rm_Weibull1,varmo_tm1_qm_exp1=varmo_tm1_qm_exp1,varmo_tm1_qm_Weibull1=varmo_tm1_qm_Weibull1,
               varmo_tm2_rm_exp1=varmo_tm2_rm_exp1,varmo_tm2_rm_Weibull1=varmo_tm2_rm_Weibull1,varmo_tm2_qm_exp1=varmo_tm2_qm_exp1,varmo_tm2_qm_Weibull1=varmo_tm2_qm_Weibull1,
               varmo_tm3_rm_exp1=varmo_tm3_rm_exp1,varmo_tm3_rm_Weibull1=varmo_tm3_rm_Weibull1,varmo_tm3_qm_exp1=varmo_tm3_qm_exp1,varmo_tm3_qm_Weibull1=varmo_tm3_qm_Weibull1,
               varmo_adaptive_AB_exp=varmo_adaptive_AB_exp,
               varmo_adaptive_AB_Weibull=varmo_adaptive_AB_Weibull,
               varmo_adaptive_SSE_exp=varmo_adaptive_SSE_exp,
               varmo_adaptive_SSE_Weibull=varmo_adaptive_SSE_Weibull,
               
               tmmoraw=tmmoraw,tmmo_BM_rm_exp1=tmmo_BM_rm_exp1,tmmo_BM_rm_Weibull1=tmmo_BM_rm_Weibull1,tmmo_BM_qm_exp1=tmmo_BM_qm_exp1,tmmo_BM_qm_Weibull1=tmmo_BM_qm_Weibull1,
               tmmo_sqm_rm_exp1=tmmo_sqm_rm_exp1,tmmo_sqm_rm_Weibull1=tmmo_sqm_rm_Weibull1,tmmo_sqm_qm_exp1=tmmo_sqm_qm_exp1,tmmo_sqm_qm_Weibull1=tmmo_sqm_qm_Weibull1,
               tmmo_wm1_rm_exp1=tmmo_wm1_rm_exp1,tmmo_wm1_rm_Weibull1=tmmo_wm1_rm_Weibull1,tmmo_wm1_qm_exp1=tmmo_wm1_qm_exp1,tmmo_wm1_qm_Weibull1=tmmo_wm1_qm_Weibull1,
               tmmo_wm2_rm_exp1=tmmo_wm2_rm_exp1,tmmo_wm2_rm_Weibull1=tmmo_wm2_rm_Weibull1,tmmo_wm2_qm_exp1=tmmo_wm2_qm_exp1,tmmo_wm2_qm_Weibull1=tmmo_wm2_qm_Weibull1,
               tmmo_tm1_rm_exp1=tmmo_tm1_rm_exp1,tmmo_tm1_rm_Weibull1=tmmo_tm1_rm_Weibull1,tmmo_tm1_qm_exp1=tmmo_tm1_qm_exp1,tmmo_tm1_qm_Weibull1=tmmo_tm1_qm_Weibull1,
               tmmo_tm2_rm_exp1=tmmo_tm2_rm_exp1,tmmo_tm2_rm_Weibull1=tmmo_tm2_rm_Weibull1,tmmo_tm2_qm_exp1=tmmo_tm2_qm_exp1,tmmo_tm2_qm_Weibull1=tmmo_tm2_qm_Weibull1,
               tmmo_tm3_rm_exp1=tmmo_tm3_rm_exp1,tmmo_tm3_rm_Weibull1=tmmo_tm3_rm_Weibull1,tmmo_tm3_qm_exp1=tmmo_tm3_qm_exp1,tmmo_tm3_qm_Weibull1=tmmo_tm3_qm_Weibull1,
               tmmo_adaptive_AB_exp=tmmo_adaptive_AB_exp,
               tmmo_adaptive_AB_Weibull=tmmo_adaptive_AB_Weibull,
               tmmo_adaptive_SSE_exp=tmmo_adaptive_SSE_exp,
               tmmo_adaptive_SSE_Weibull=tmmo_adaptive_SSE_Weibull,
               
               fmmoraw=fmmoraw,fmmo_BM_rm_exp1=fmmo_BM_rm_exp1,fmmo_BM_rm_Weibull1=fmmo_BM_rm_Weibull1,fmmo_BM_qm_exp1=fmmo_BM_qm_exp1,fmmo_BM_qm_Weibull1=fmmo_BM_qm_Weibull1,
               fmmo_sqm_rm_exp1=fmmo_sqm_rm_exp1,fmmo_sqm_rm_Weibull1=fmmo_sqm_rm_Weibull1,fmmo_sqm_qm_exp1=fmmo_sqm_qm_exp1,fmmo_sqm_qm_Weibull1=fmmo_sqm_qm_Weibull1,
               fmmo_wm1_rm_exp1=fmmo_wm1_rm_exp1,fmmo_wm1_rm_Weibull1=fmmo_wm1_rm_Weibull1,fmmo_wm1_qm_exp1=fmmo_wm1_qm_exp1,fmmo_wm1_qm_Weibull1=fmmo_wm1_qm_Weibull1,
               fmmo_wm2_rm_exp1=fmmo_wm2_rm_exp1,fmmo_wm2_rm_Weibull1=fmmo_wm2_rm_Weibull1,fmmo_wm2_qm_exp1=fmmo_wm2_qm_exp1,fmmo_wm2_qm_Weibull1=fmmo_wm2_qm_Weibull1,
               fmmo_tm1_rm_exp1=fmmo_tm1_rm_exp1,fmmo_tm1_rm_Weibull1=fmmo_tm1_rm_Weibull1,fmmo_tm1_qm_exp1=fmmo_tm1_qm_exp1,fmmo_tm1_qm_Weibull1=fmmo_tm1_qm_Weibull1,
               fmmo_tm2_rm_exp1=fmmo_tm2_rm_exp1,fmmo_tm2_rm_Weibull1=fmmo_tm2_rm_Weibull1,fmmo_tm2_qm_exp1=fmmo_tm2_qm_exp1,fmmo_tm2_qm_Weibull1=fmmo_tm2_qm_Weibull1,
               fmmo_tm3_rm_exp1=fmmo_tm3_rm_exp1,fmmo_tm3_rm_Weibull1=fmmo_tm3_rm_Weibull1,fmmo_tm3_qm_exp1=fmmo_tm3_qm_exp1,fmmo_tm3_qm_Weibull1=fmmo_tm3_qm_Weibull1,
               fmmo_adaptive_AB_exp=fmmo_adaptive_AB_exp,
               fmmo_adaptive_AB_Weibull=fmmo_adaptive_AB_Weibull,
               fmmo_adaptive_SSE_exp=fmmo_adaptive_SSE_exp,
               fmmo_adaptive_SSE_Weibull=fmmo_adaptive_SSE_Weibull,
               sdall=sdall)
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
  dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  
  momentssd<-c(sqrt(momentsx[2]),dall[194],dall[195],dall[196])
  
  allrawmoBias<-c(
    firstbias=abs(c(onestepx,SMWM9,dall[1:49])-dall[2])/((momentssd[1])),
    secondbias=abs(c(dall[50:97])-dall[50])/((momentssd[2])),
    thirdbias=abs(c(dall[98:145])-dall[98])/((momentssd[3])),
    fourbias=abs(c(dall[146:193])-dall[146])/((momentssd[4])))
  
  all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,allrawmoBias,momentssd))
}

write.csv(simulatedbatchWeibull_asymptoticbias,paste("asymptotic_Weibull_raw_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatchWeibull_asymptoticbias[1:length(allkurtWeibull),1],simulatedbatchWeibull_asymptoticbias[1:length(allkurtWeibull),209:408]),paste("asymptotic_Weibull",largesize,".csv", sep = ","), row.names = FALSE)

samplesize=5400
batchsizebase=2000
orderlist1_AB2<-removelist(na.omit(t(apply(quasiuni_sorted2,MARGIN=1,FUN=roundunique,dimension=2,size=samplesize))))
orderlist1_AB3<-removelist(na.omit(t(apply(quasiuni_sorted3,MARGIN=1,FUN=roundunique,dimension=3,size=samplesize))))
orderlist1_AB4<-removelist(na.omit(t(apply(quasiuni_sorted4,MARGIN=1,FUN=roundunique,dimension=4,size=samplesize))))

batchsize=batchsizebase

n <- samplesize
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)

simulatedbatchWeibull_ABSE<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
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
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    onestepx<-onestep(x=sortedx, bend = 1.172)
    SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
    dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
    
    momentsx<-unbiasedmoments(x=sortedx)
    sortedx<-c()
    momentssd<-c(sqrt(momentsx[2]),dall[194],dall[195],dall[196])
    
    allrawmoBias<-c(
      firstbias=(c(onestepx,SMWM9,dall[1:49],momentsx[1])-targetm),
      secondbias=(c(dall[50:97],momentsx[2])-targetvar),
      thirdbias=(c(dall[98:145],momentsx[3])-targettm),
      fourbias=(c(dall[146:193],momentsx[4])-targetfm))
    
    all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,momentssd,allrawmoBias))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("Weibull_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,213:265])^2))/simulatedbatchWeibull_asymptoticbias[batchnumber,405]
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,266:314])^2))/simulatedbatchWeibull_asymptoticbias[batchnumber,406]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,315:363])^2))/simulatedbatchWeibull_asymptoticbias[batchnumber,407]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,364:412])^2))/simulatedbatchWeibull_asymptoticbias[batchnumber,408]
  
  AB1_mean<-abs(colMeans((SEbataches[,213:265])))/simulatedbatchWeibull_asymptoticbias[batchnumber,405]
  
  AB1_var<-abs(colMeans((SEbataches[,266:314])))/simulatedbatchWeibull_asymptoticbias[batchnumber,406]
  
  AB1_tm<-abs(colMeans((SEbataches[,315:363])))/simulatedbatchWeibull_asymptoticbias[batchnumber,407]
  
  AB1_fm<-abs(colMeans((SEbataches[,364:412])))/simulatedbatchWeibull_asymptoticbias[batchnumber,408]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,205])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,405]
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,206])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,406]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,207])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,407]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,208])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,408]
  
  ratiosamplemean1<-c(SEbatachesmean[205])/SEbatachesmean[6]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,205])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/simulatedbatchWeibull_asymptoticbias[batchnumber,405]
  
  ratiosamplevar1<-c(SEbatachesmean[206])/SEbatachesmean[54]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,206])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,406]
  
  ratiosampletm1<-c(SEbatachesmean[207])/SEbatachesmean[102]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,207])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,407]
  
  ratiosamplefm1<-c(SEbatachesmean[208])/SEbatachesmean[150]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,208])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  samplefm_SSE1<-samplefmsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,408]
  
  ratiomean1<-c(SEbatachesmean[2:53])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:53]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,405]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:53])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/simulatedbatchWeibull_asymptoticbias[batchnumber,405]
  
  ratiovar1<-SEbatachesmean[54:101]/SEbatachesmean[54]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,54:101]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,406]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,54:101])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,406]
  ratiotm1<-SEbatachesmean[102:149]/SEbatachesmean[102]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,102:149]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,407]
  
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,102:149])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,407]
  
  ratiofm1<-SEbatachesmean[150:197]/SEbatachesmean[150]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,150:197]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,408]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,150:197])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,408]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatchWeibull_ABSE,paste("Weibull_ABSSE.csv", sep = ","), row.names = FALSE)



simulatedbatchWeibull_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtWeibull))), .combine = 'rbind') %dopar% {
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
  
  SEbataches<- read.csv(paste("Weibull_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,205])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,405]
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,206])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,406]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,207])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,407]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,208])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,408]
  
  ratiosamplemean1<-c(SEbatachesmean[205])/SEbatachesmean[6]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,205])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/simulatedbatchWeibull_asymptoticbias[batchnumber,405]
  
  ratiosamplevar1<-c(SEbatachesmean[206])/SEbatachesmean[54]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,206])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,406]
  
  ratiosampletm1<-c(SEbatachesmean[207])/SEbatachesmean[102]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,207])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,407]
  
  ratiosamplefm1<-c(SEbatachesmean[208])/SEbatachesmean[150]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,208])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  samplefm_SSE1<-samplefmsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,408]
  
  ratiomean1<-c(SEbatachesmean[2:53])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:53]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,405]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:53])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/simulatedbatchWeibull_asymptoticbias[batchnumber,405]
  
  ratiovar1<-SEbatachesmean[54:101]/SEbatachesmean[54]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,54:101]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,406]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,54:101])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,406]
  ratiotm1<-SEbatachesmean[102:149]/SEbatachesmean[102]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,102:149]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,407]
  
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,102:149])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,407]
  
  ratiofm1<-SEbatachesmean[150:197]/SEbatachesmean[150]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,150:197]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatchWeibull_asymptoticbias[batchnumber,408]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,150:197])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatchWeibull_asymptoticbias[batchnumber,408]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean_all1=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}

write.csv(simulatedbatchWeibull_ABSE_SE,paste("Weibull_ABSSE_error.csv", sep = ","), row.names = FALSE)


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
  dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  
  momentssd<-c(sqrt(momentsx[2]),dall[194],dall[195],dall[196])
  
  allrawmoBias<-c(
    firstbias=abs(c(onestepx,SMWM9,dall[1:49])-dall[2])/((momentssd[1])),
    secondbias=abs(c(dall[50:97])-dall[50])/((momentssd[2])),
    thirdbias=abs(c(dall[98:145])-dall[98])/((momentssd[3])),
    fourbias=abs(c(dall[146:193])-dall[146])/((momentssd[4])))
  
  all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,allrawmoBias,momentssd))
}

write.csv(simulatedbatchgamma_asymptoticbias,paste("asymptotic_gamma_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatchgamma_asymptoticbias[1:length(allkurtgamma),1],simulatedbatchgamma_asymptoticbias[1:length(allkurtgamma),209:408]),paste("asymptotic_gamma",largesize,".csv", sep = ","), row.names = FALSE)


simulatedbatchgamma_ABSE<-foreach(batchnumber =c((1:length(allkurtgamma))), .combine = 'rbind') %dopar% {
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
  
  
  a=allkurtgamma[batchnumber]
  targetm<-a
  targetvar<-(a)
  targettm<-((sqrt(a))^3)*2/sqrt(a)
  targetfm<-((sqrt(a))^4)*((6/(a))+3)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
    x<-c(dsgamma(uni=unibatch[,batch1], shape=a, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    onestepx<-onestep(x=sortedx, bend = 1.172)
    SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
    dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
    
    momentsx<-unbiasedmoments(x=sortedx)
    sortedx<-c()
    momentssd<-c(sqrt(momentsx[2]),dall[194],dall[195],dall[196])
    
    allrawmoBias<-c(
      firstbias=(c(onestepx,SMWM9,dall[1:49],momentsx[1])-targetm),
      secondbias=(c(dall[50:97],momentsx[2])-targetvar),
      thirdbias=(c(dall[98:145],momentsx[3])-targettm),
      fourbias=(c(dall[146:193],momentsx[4])-targetfm))
    
    all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,momentssd,allrawmoBias))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("gamma_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,213:265])^2))/simulatedbatchgamma_asymptoticbias[batchnumber,405]
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,266:314])^2))/simulatedbatchgamma_asymptoticbias[batchnumber,406]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,315:363])^2))/simulatedbatchgamma_asymptoticbias[batchnumber,407]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,364:412])^2))/simulatedbatchgamma_asymptoticbias[batchnumber,408]
  
  AB1_mean<-abs(colMeans((SEbataches[,213:265])))/simulatedbatchgamma_asymptoticbias[batchnumber,405]
  
  AB1_var<-abs(colMeans((SEbataches[,266:314])))/simulatedbatchgamma_asymptoticbias[batchnumber,406]
  
  AB1_tm<-abs(colMeans((SEbataches[,315:363])))/simulatedbatchgamma_asymptoticbias[batchnumber,407]
  
  AB1_fm<-abs(colMeans((SEbataches[,364:412])))/simulatedbatchgamma_asymptoticbias[batchnumber,408]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,205])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,405]
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,206])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,406]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,207])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,407]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,208])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,408]
  
  ratiosamplemean1<-c(SEbatachesmean[205])/SEbatachesmean[6]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,205])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/simulatedbatchgamma_asymptoticbias[batchnumber,405]
  
  ratiosamplevar1<-c(SEbatachesmean[206])/SEbatachesmean[54]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,206])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatchgamma_asymptoticbias[batchnumber,406]
  
  ratiosampletm1<-c(SEbatachesmean[207])/SEbatachesmean[102]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,207])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatchgamma_asymptoticbias[batchnumber,407]
  
  ratiosamplefm1<-c(SEbatachesmean[208])/SEbatachesmean[150]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,208])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  samplefm_SSE1<-samplefmsd1/simulatedbatchgamma_asymptoticbias[batchnumber,408]
  
  ratiomean1<-c(SEbatachesmean[2:53])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:53]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,405]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:53])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/simulatedbatchgamma_asymptoticbias[batchnumber,405]
  
  ratiovar1<-SEbatachesmean[54:101]/SEbatachesmean[54]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,54:101]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,406]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,54:101])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatchgamma_asymptoticbias[batchnumber,406]
  ratiotm1<-SEbatachesmean[102:149]/SEbatachesmean[102]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,102:149]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,407]
  
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,102:149])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatchgamma_asymptoticbias[batchnumber,407]
  
  ratiofm1<-SEbatachesmean[150:197]/SEbatachesmean[150]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,150:197]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,408]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,150:197])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatchgamma_asymptoticbias[batchnumber,408]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}

write.csv(simulatedbatchgamma_ABSE,paste("gamma_ABSSE.csv", sep = ","), row.names = FALSE)


simulatedbatchgamma_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtgamma))), .combine = 'rbind') %dopar% {
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
  
  a=allkurtgamma[batchnumber]
  targetm<-a
  targetvar<-(a)
  targettm<-((sqrt(a))^3)*2/sqrt(a)
  targetfm<-((sqrt(a))^4)*((6/(a))+3)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  SEbataches<- read.csv(paste("gamma_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,205])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,405]
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,206])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,406]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,207])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,407]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,208])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,408]
  
  ratiosamplemean1<-c(SEbatachesmean[205])/SEbatachesmean[6]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,205])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/simulatedbatchgamma_asymptoticbias[batchnumber,405]
  
  ratiosamplevar1<-c(SEbatachesmean[206])/SEbatachesmean[54]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,206])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatchgamma_asymptoticbias[batchnumber,406]
  
  ratiosampletm1<-c(SEbatachesmean[207])/SEbatachesmean[102]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,207])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatchgamma_asymptoticbias[batchnumber,407]
  
  ratiosamplefm1<-c(SEbatachesmean[208])/SEbatachesmean[150]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,208])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  samplefm_SSE1<-samplefmsd1/simulatedbatchgamma_asymptoticbias[batchnumber,408]
  
  ratiomean1<-c(SEbatachesmean[2:53])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:53]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,405]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:53])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/simulatedbatchgamma_asymptoticbias[batchnumber,405]
  
  ratiovar1<-SEbatachesmean[54:101]/SEbatachesmean[54]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,54:101]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,406]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,54:101])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatchgamma_asymptoticbias[batchnumber,406]
  ratiotm1<-SEbatachesmean[102:149]/SEbatachesmean[102]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,102:149]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,407]
  
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,102:149])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatchgamma_asymptoticbias[batchnumber,407]
  
  ratiofm1<-SEbatachesmean[150:197]/SEbatachesmean[150]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,150:197]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatchgamma_asymptoticbias[batchnumber,408]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,150:197])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatchgamma_asymptoticbias[batchnumber,408]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean_all1=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}

write.csv(simulatedbatchgamma_ABSE_SE,paste("gamma_ABSSE_error.csv", sep = ","), row.names = FALSE)


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
  dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  
  momentssd<-c(sqrt(momentsx[2]),dall[194],dall[195],dall[196])
  
  allrawmoBias<-c(
    firstbias=abs(c(onestepx,SMWM9,dall[1:49])-dall[2])/((momentssd[1])),
    secondbias=abs(c(dall[50:97])-dall[50])/((momentssd[2])),
    thirdbias=abs(c(dall[98:145])-dall[98])/((momentssd[3])),
    fourbias=abs(c(dall[146:193])-dall[146])/((momentssd[4])))
  
  all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,allrawmoBias,momentssd))
}

write.csv(simulatedbatchPareto_asymptoticbias,paste("asymptotic_Pareto_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatchPareto_asymptoticbias[1:length(allkurtPareto),1],simulatedbatchPareto_asymptoticbias[1:length(allkurtPareto),209:408]),paste("asymptotic_Pareto",largesize,".csv", sep = ","), row.names = FALSE)


simulatedbatchPareto_ABSE<-foreach(batchnumber =c((1:length(allkurtPareto))), .combine = 'rbind') %dopar% {
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
  
  
  a=allkurtPareto[batchnumber]
  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
    x<-c(dsPareto(uni=unibatch[,batch1], shape=a, scale = 1))
    
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    onestepx<-onestep(x=sortedx, bend = 1.172)
    SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
    dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
    
    momentsx<-unbiasedmoments(x=sortedx)
    sortedx<-c()
    momentssd<-c(sqrt(momentsx[2]),dall[194],dall[195],dall[196])
    
    allrawmoBias<-c(
      firstbias=(c(onestepx,SMWM9,dall[1:49],momentsx[1])-targetm),
      secondbias=(c(dall[50:97],momentsx[2])-targetvar),
      thirdbias=(c(dall[98:145],momentsx[3])-targettm),
      fourbias=(c(dall[146:193],momentsx[4])-targetfm))
    
    all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,momentssd,allrawmoBias))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("Pareto_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,213:265])^2))/simulatedbatchPareto_asymptoticbias[batchnumber,405]
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,266:314])^2))/simulatedbatchPareto_asymptoticbias[batchnumber,406]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,315:363])^2))/simulatedbatchPareto_asymptoticbias[batchnumber,407]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,364:412])^2))/simulatedbatchPareto_asymptoticbias[batchnumber,408]
  
  AB1_mean<-abs(colMeans((SEbataches[,213:265])))/simulatedbatchPareto_asymptoticbias[batchnumber,405]
  
  AB1_var<-abs(colMeans((SEbataches[,266:314])))/simulatedbatchPareto_asymptoticbias[batchnumber,406]
  
  AB1_tm<-abs(colMeans((SEbataches[,315:363])))/simulatedbatchPareto_asymptoticbias[batchnumber,407]
  
  AB1_fm<-abs(colMeans((SEbataches[,364:412])))/simulatedbatchPareto_asymptoticbias[batchnumber,408]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,205])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,405]
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,206])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,406]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,207])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,407]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,208])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,408]
  
  ratiosamplemean1<-c(SEbatachesmean[205])/SEbatachesmean[6]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,205])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/simulatedbatchPareto_asymptoticbias[batchnumber,405]
  
  ratiosamplevar1<-c(SEbatachesmean[206])/SEbatachesmean[54]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,206])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatchPareto_asymptoticbias[batchnumber,406]
  
  ratiosampletm1<-c(SEbatachesmean[207])/SEbatachesmean[102]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,207])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatchPareto_asymptoticbias[batchnumber,407]
  
  ratiosamplefm1<-c(SEbatachesmean[208])/SEbatachesmean[150]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,208])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  samplefm_SSE1<-samplefmsd1/simulatedbatchPareto_asymptoticbias[batchnumber,408]
  
  ratiomean1<-c(SEbatachesmean[2:53])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:53]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,405]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:53])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/simulatedbatchPareto_asymptoticbias[batchnumber,405]
  
  ratiovar1<-SEbatachesmean[54:101]/SEbatachesmean[54]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,54:101]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,406]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,54:101])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatchPareto_asymptoticbias[batchnumber,406]
  ratiotm1<-SEbatachesmean[102:149]/SEbatachesmean[102]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,102:149]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,407]
  
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,102:149])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatchPareto_asymptoticbias[batchnumber,407]
  
  ratiofm1<-SEbatachesmean[150:197]/SEbatachesmean[150]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,150:197]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,408]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,150:197])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatchPareto_asymptoticbias[batchnumber,408]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}


write.csv(simulatedbatchPareto_ABSE,paste("Pareto_ABSSE.csv", sep = ","), row.names = FALSE)


simulatedbatchPareto_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtPareto))), .combine = 'rbind') %dopar% {
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
  
  
  a=allkurtPareto[batchnumber]
  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  SEbataches<- read.csv(paste("Pareto_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,205])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,405]
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,206])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,406]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,207])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,407]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,208])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,408]
  
  ratiosamplemean1<-c(SEbatachesmean[205])/SEbatachesmean[6]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,205])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/simulatedbatchPareto_asymptoticbias[batchnumber,405]
  
  ratiosamplevar1<-c(SEbatachesmean[206])/SEbatachesmean[54]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,206])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatchPareto_asymptoticbias[batchnumber,406]
  
  ratiosampletm1<-c(SEbatachesmean[207])/SEbatachesmean[102]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,207])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatchPareto_asymptoticbias[batchnumber,407]
  
  ratiosamplefm1<-c(SEbatachesmean[208])/SEbatachesmean[150]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,208])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  samplefm_SSE1<-samplefmsd1/simulatedbatchPareto_asymptoticbias[batchnumber,408]
  
  ratiomean1<-c(SEbatachesmean[2:53])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:53]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,405]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:53])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/simulatedbatchPareto_asymptoticbias[batchnumber,405]
  
  ratiovar1<-SEbatachesmean[54:101]/SEbatachesmean[54]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,54:101]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,406]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,54:101])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatchPareto_asymptoticbias[batchnumber,406]
  ratiotm1<-SEbatachesmean[102:149]/SEbatachesmean[102]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,102:149]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,407]
  
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,102:149])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatchPareto_asymptoticbias[batchnumber,407]
  
  ratiofm1<-SEbatachesmean[150:197]/SEbatachesmean[150]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,150:197]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatchPareto_asymptoticbias[batchnumber,408]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,150:197])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatchPareto_asymptoticbias[batchnumber,408]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean_all1=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}


write.csv(simulatedbatchPareto_ABSE_SE,paste("Pareto_ABSSE_error.csv", sep = ","), row.names = FALSE)


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
  dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  
  momentssd<-c(sqrt(momentsx[2]),dall[194],dall[195],dall[196])
  
  allrawmoBias<-c(
    firstbias=abs(c(onestepx,SMWM9,dall[1:49])-dall[2])/((momentssd[1])),
    secondbias=abs(c(dall[50:97])-dall[50])/((momentssd[2])),
    thirdbias=abs(c(dall[98:145])-dall[98])/((momentssd[3])),
    fourbias=abs(c(dall[146:193])-dall[146])/((momentssd[4])))
  
  all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,allrawmoBias,momentssd))
}

write.csv(simulatedbatchlognorm_asymptoticbias,paste("asymptotic_lognorm_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatchlognorm_asymptoticbias[1:length(allkurtlognorm),1],simulatedbatchlognorm_asymptoticbias[1:length(allkurtlognorm),209:408]),paste("asymptotic_lognorm",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatchlognorm_ABSE<-foreach(batchnumber =c((1:length(allkurtlognorm))), .combine = 'rbind') %dopar% {
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
  
  
  a=allkurtlognorm[batchnumber]
  
  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
    x<-c(dslnorm(uni=unibatch[,batch1], location=0, scale = a/1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    onestepx<-onestep(x=sortedx, bend = 1.172)
    SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
    dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
    
    momentsx<-unbiasedmoments(x=sortedx)
    sortedx<-c()
    momentssd<-c(sqrt(momentsx[2]),dall[194],dall[195],dall[196])
    
    allrawmoBias<-c(
      firstbias=(c(onestepx,SMWM9,dall[1:49],momentsx[1])-targetm),
      secondbias=(c(dall[50:97],momentsx[2])-targetvar),
      thirdbias=(c(dall[98:145],momentsx[3])-targettm),
      fourbias=(c(dall[146:193],momentsx[4])-targetfm))
    
    all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,momentssd,allrawmoBias))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("lognorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,213:265])^2))/simulatedbatchlognorm_asymptoticbias[batchnumber,405]
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,266:314])^2))/simulatedbatchlognorm_asymptoticbias[batchnumber,406]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,315:363])^2))/simulatedbatchlognorm_asymptoticbias[batchnumber,407]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,364:412])^2))/simulatedbatchlognorm_asymptoticbias[batchnumber,408]
  
  AB1_mean<-abs(colMeans((SEbataches[,213:265])))/simulatedbatchlognorm_asymptoticbias[batchnumber,405]
  
  AB1_var<-abs(colMeans((SEbataches[,266:314])))/simulatedbatchlognorm_asymptoticbias[batchnumber,406]
  
  AB1_tm<-abs(colMeans((SEbataches[,315:363])))/simulatedbatchlognorm_asymptoticbias[batchnumber,407]
  
  AB1_fm<-abs(colMeans((SEbataches[,364:412])))/simulatedbatchlognorm_asymptoticbias[batchnumber,408]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,205])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,405]
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,206])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,406]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,207])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,407]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,208])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,408]
  
  ratiosamplemean1<-c(SEbatachesmean[205])/SEbatachesmean[6]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,205])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/simulatedbatchlognorm_asymptoticbias[batchnumber,405]
  
  ratiosamplevar1<-c(SEbatachesmean[206])/SEbatachesmean[54]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,206])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,406]
  
  ratiosampletm1<-c(SEbatachesmean[207])/SEbatachesmean[102]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,207])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,407]
  
  ratiosamplefm1<-c(SEbatachesmean[208])/SEbatachesmean[150]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,208])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  samplefm_SSE1<-samplefmsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,408]
  
  ratiomean1<-c(SEbatachesmean[2:53])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:53]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,405]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:53])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/simulatedbatchlognorm_asymptoticbias[batchnumber,405]
  
  ratiovar1<-SEbatachesmean[54:101]/SEbatachesmean[54]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,54:101]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,406]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,54:101])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,406]
  ratiotm1<-SEbatachesmean[102:149]/SEbatachesmean[102]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,102:149]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,407]
  
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,102:149])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,407]
  
  ratiofm1<-SEbatachesmean[150:197]/SEbatachesmean[150]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,150:197]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,408]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,150:197])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,408]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}


write.csv(simulatedbatchlognorm_ABSE,paste("lognorm_ABSSE.csv", sep = ","), row.names = FALSE)


simulatedbatchlognorm_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtlognorm))), .combine = 'rbind') %dopar% {
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
  
  
  a=allkurtlognorm[batchnumber]
  
  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  SEbataches<- read.csv(paste("lognorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,205])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,405]
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,206])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,406]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,207])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,407]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,208])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,408]
  
  ratiosamplemean1<-c(SEbatachesmean[205])/SEbatachesmean[6]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,205])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/simulatedbatchlognorm_asymptoticbias[batchnumber,405]
  
  ratiosamplevar1<-c(SEbatachesmean[206])/SEbatachesmean[54]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,206])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,406]
  
  ratiosampletm1<-c(SEbatachesmean[207])/SEbatachesmean[102]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,207])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,407]
  
  ratiosamplefm1<-c(SEbatachesmean[208])/SEbatachesmean[150]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,208])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  samplefm_SSE1<-samplefmsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,408]
  
  ratiomean1<-c(SEbatachesmean[2:53])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:53]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,405]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:53])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/simulatedbatchlognorm_asymptoticbias[batchnumber,405]
  
  ratiovar1<-SEbatachesmean[54:101]/SEbatachesmean[54]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,54:101]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,406]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,54:101])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,406]
  ratiotm1<-SEbatachesmean[102:149]/SEbatachesmean[102]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,102:149]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,407]
  
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,102:149])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,407]
  
  ratiofm1<-SEbatachesmean[150:197]/SEbatachesmean[150]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,150:197]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatchlognorm_asymptoticbias[batchnumber,408]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,150:197])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatchlognorm_asymptoticbias[batchnumber,408]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean_all1=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}


write.csv(simulatedbatchlognorm_ABSE_SE,paste("lognorm_ABSSE_error.csv", sep = ","), row.names = FALSE)


#gnorm
kurtgnorm<- read.csv(("kurtgnorm_31150.csv"))
allkurtgnorm<-unlist(kurtgnorm)

simulatedbatchgnorm_asymptoticbias<-foreach(batchnumber = (1:length(allkurtgnorm)), .combine = 'rbind') %dopar% {
  library(Rfast)
  a=allkurtgnorm[batchnumber]
  x<-c(dsgnorm(uni=quasiuni_asymptotic, shape=a, scale = 1))
  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  
  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  onestepx<-onestep(x=sortedx, bend = 1.172)
  SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
  dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  
  momentssd<-c(sqrt(momentsx[2]),dall[194],dall[195],dall[196])
  
  allrawmoBias<-c(
    firstbias=abs(c(onestepx,SMWM9,dall[1:49])-dall[2])/((momentssd[1])),
    secondbias=abs(c(dall[50:97])-dall[50])/((momentssd[2])),
    thirdbias=abs(c(dall[98:145])-dall[98])/((momentssd[3])),
    fourbias=abs(c(dall[146:193])-dall[146])/((momentssd[4])))
  
  all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,allrawmoBias,momentssd))
}

write.csv(simulatedbatchgnorm_asymptoticbias,paste("asymptotic_gnorm_Process",largesize,".csv", sep = ","), row.names = FALSE)

write.csv(cbind(simulatedbatchgnorm_asymptoticbias[1:length(allkurtgnorm),1],simulatedbatchgnorm_asymptoticbias[1:length(allkurtgnorm),209:408]),paste("asymptotic_gnorm",largesize,".csv", sep = ","), row.names = FALSE)

simulatedbatchgnorm_ABSE<-foreach(batchnumber =c((1:length(allkurtgnorm))), .combine = 'rbind') %dopar% {
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
  
  
  a=allkurtgnorm[batchnumber]
  
  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
    x<-c(dsgnorm(uni=unibatch[,batch1], shape=a, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    onestepx<-onestep(x=sortedx, bend = 1.172)
    SMWM9<-sm(x=sortedx,interval=9,fast=TRUE,batch="auto")
    dall<-Balltest(x=sortedx,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,interval=8,fast=TRUE,batch="auto",startpoint=9,stepsize=1000,criterion=criterionset)
    
    momentsx<-unbiasedmoments(x=sortedx)
    sortedx<-c()
    momentssd<-c(sqrt(momentsx[2]),dall[194],dall[195],dall[196])
    
    allrawmoBias<-c(
      firstbias=(c(onestepx,SMWM9,dall[1:49],momentsx[1])-targetm),
      secondbias=(c(dall[50:97],momentsx[2])-targetvar),
      thirdbias=(c(dall[98:145],momentsx[3])-targettm),
      fourbias=(c(dall[146:193],momentsx[4])-targetfm))
    
    all1<-t(c(kurtx,onestepx,SMWM9,dall,targetall,momentsx,momentssd,allrawmoBias))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("gnorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,213:265])^2))/simulatedbatchgnorm_asymptoticbias[batchnumber,405]
  
  RMSE1_var<-sqrt(colMeans((SEbataches[,266:314])^2))/simulatedbatchgnorm_asymptoticbias[batchnumber,406]
  
  RMSE1_tm<-sqrt(colMeans((SEbataches[,315:363])^2))/simulatedbatchgnorm_asymptoticbias[batchnumber,407]
  
  RMSE1_fm<-sqrt(colMeans((SEbataches[,364:412])^2))/simulatedbatchgnorm_asymptoticbias[batchnumber,408]
  
  AB1_mean<-abs(colMeans((SEbataches[,213:265])))/simulatedbatchgnorm_asymptoticbias[batchnumber,405]
  
  AB1_var<-abs(colMeans((SEbataches[,266:314])))/simulatedbatchgnorm_asymptoticbias[batchnumber,406]
  
  AB1_tm<-abs(colMeans((SEbataches[,315:363])))/simulatedbatchgnorm_asymptoticbias[batchnumber,407]
  
  AB1_fm<-abs(colMeans((SEbataches[,364:412])))/simulatedbatchgnorm_asymptoticbias[batchnumber,408]
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,205])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,405]
  
  samplevarsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,206])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,406]
  
  sampletmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,207])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,407]
  
  samplefmsd_unscaled1<-unbiasedsd(x=SEbataches[1:batchsize,208])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,408]
  
  ratiosamplemean1<-c(SEbatachesmean[205])/SEbatachesmean[6]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,205])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, unbiasedsd)
  samplemean_SSE1<-samplemeansd1/simulatedbatchgnorm_asymptoticbias[batchnumber,405]
  
  ratiosamplevar1<-c(SEbatachesmean[206])/SEbatachesmean[54]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,206])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, unbiasedsd)
  samplevar_SSE1<-samplevarsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,406]
  
  ratiosampletm1<-c(SEbatachesmean[207])/SEbatachesmean[102]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,207])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, unbiasedsd)
  sampletm_SSE1<-sampletmsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,407]
  
  ratiosamplefm1<-c(SEbatachesmean[208])/SEbatachesmean[150]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,208])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, unbiasedsd)
  samplefm_SSE1<-samplefmsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,408]
  
  ratiomean1<-c(SEbatachesmean[2:53])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:53]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,405]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:53])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  mean_SSE1<-meansd1/simulatedbatchgnorm_asymptoticbias[batchnumber,405]
  
  ratiovar1<-SEbatachesmean[54:101]/SEbatachesmean[54]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,54:101]), 2, unbiasedsd)
  
  var_SE1<-varsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,406]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,54:101])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, unbiasedsd)
  
  var_SSE1<-varsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,406]
  ratiotm1<-SEbatachesmean[102:149]/SEbatachesmean[102]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,102:149]), 2, unbiasedsd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,407]
  
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,102:149])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, unbiasedsd)
  tm_SSE1<-tmsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,407]
  
  ratiofm1<-SEbatachesmean[150:197]/SEbatachesmean[150]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,150:197]), 2, unbiasedsd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,408]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,150:197])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, unbiasedsd)
  fm_SSE1<-fmsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,408]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],RMSE1_mean=RMSE1_mean,RMSE1_var=RMSE1_var,RMSE1_tm=RMSE1_tm,RMSE1_fm=RMSE1_fm,AB1_mean=AB1_mean,AB1_var=AB1_var,AB1_tm=AB1_tm,AB1_fm=AB1_fm,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}


write.csv(simulatedbatchgnorm_ABSE,paste("gnorm_ABSSE.csv", sep = ","), row.names = FALSE)



simulatedbatchgnorm_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtgnorm))), .combine = 'rbind') %dopar% {
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
  
  
  a=allkurtgnorm[batchnumber]
  
  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  SEbataches<- read.csv(paste("gnorm_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  SEbatachesmean <- colMeans(SEbataches)
  
  samplemeansd_unscaled1<-se_sd(x=SEbataches[1:batchsize,205])
  
  samplemean_SE1<-samplemeansd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,405]
  
  samplevarsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,206])
  
  samplevar_SE1<-samplevarsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,406]
  
  sampletmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,207])
  
  sampletm_SE1<-sampletmsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,407]
  
  samplefmsd_unscaled1<-se_sd(x=SEbataches[1:batchsize,208])
  
  samplefm_SE1<-samplefmsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,408]
  
  ratiosamplemean1<-c(SEbatachesmean[205])/SEbatachesmean[6]
  
  samplemean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,205])/ratiosamplemean1)
  
  samplemeansd1<-apply((samplemean_SEbatachesmeanprocess), 2, se_sd)
  samplemean_SSE1<-samplemeansd1/simulatedbatchgnorm_asymptoticbias[batchnumber,405]
  
  ratiosamplevar1<-c(SEbatachesmean[206])/SEbatachesmean[54]
  
  samplevar_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,206])/ratiosamplevar1)
  
  samplevarsd1<-apply((samplevar_SEbatachesmeanprocess), 2, se_sd)
  samplevar_SSE1<-samplevarsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,406]
  
  ratiosampletm1<-c(SEbatachesmean[207])/SEbatachesmean[102]
  
  sampletm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,207])/ratiosampletm1)
  
  sampletmsd1<-apply((sampletm_SEbatachesmeanprocess), 2, se_sd)
  sampletm_SSE1<-sampletmsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,407]
  
  ratiosamplefm1<-c(SEbatachesmean[208])/SEbatachesmean[150]
  
  samplefm_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,208])/ratiosamplefm1)
  
  samplefmsd1<-apply((samplefm_SEbatachesmeanprocess), 2, se_sd)
  samplefm_SSE1<-samplefmsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,408]
  
  ratiomean1<-c(SEbatachesmean[2:53])/SEbatachesmean[6]
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,2:53]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,405]
  mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,2:53])/ratiomean1)
  
  meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  mean_SSE1<-meansd1/simulatedbatchgnorm_asymptoticbias[batchnumber,405]
  
  ratiovar1<-SEbatachesmean[54:101]/SEbatachesmean[54]
  
  varsd_unscaled1<-apply((SEbataches[1:batchsize,54:101]), 2, se_sd)
  
  var_SE1<-varsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,406]
  
  var_SEbatachesvarprocess<-(t(t(SEbataches[1:batchsize,54:101])/ratiovar1))
  
  varsd1<-apply(var_SEbatachesvarprocess, 2, se_sd)
  
  var_SSE1<-varsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,406]
  ratiotm1<-SEbatachesmean[102:149]/SEbatachesmean[102]
  
  tmsd_unscaled1<-apply((SEbataches[1:batchsize,102:149]), 2, se_sd)
  
  tm_SE1<-tmsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,407]
  
  
  tm_SEbatachestmprocess<-(t(t(SEbataches[1:batchsize,102:149])/ratiotm1))
  tmsd1<-apply(tm_SEbatachestmprocess, 2, se_sd)
  tm_SSE1<-tmsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,407]
  
  ratiofm1<-SEbatachesmean[150:197]/SEbatachesmean[150]
  
  fmsd_unscaled1<-apply((SEbataches[1:batchsize,150:197]), 2, se_sd)
  
  fm_SE1<-fmsd_unscaled1/simulatedbatchgnorm_asymptoticbias[batchnumber,408]
  
  fm_SEbatachesfmprocess<-(t(t(SEbataches[1:batchsize,150:197])/ratiofm1))
  fmsd1<-apply(fm_SEbatachesfmprocess, 2, se_sd)
  fm_SSE1<-fmsd1/simulatedbatchgnorm_asymptoticbias[batchnumber,408]
  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1],samplevar_SE1=samplevar_SE1,var_SE1=var_SE1,SEbatachesmean[1],sampletm_SE1=sampletm_SE1,tm_SE1=tm_SE1,SEbatachesmean[1],samplefm_SE1=samplefm_SE1,fm_SE1=fm_SE1
  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1,SEbatachesmean[1],samplevarsd_unscaled1=samplevarsd_unscaled1,varsd_unscaled1=varsd_unscaled1,SEbatachesmean[1],
                  sampletmsd_unscaled1=sampletmsd_unscaled1,
                  tmsd_unscaled1=tmsd_unscaled1,SEbatachesmean[1],samplefmsd_unscaled1=samplefmsd_unscaled1,fmsd_unscaled1=fmsd_unscaled1
  )
  allSSE<-c(SEbatachesmean[1],mean_SSE1=mean_SSE1,SEbatachesmean[1],samplevar_SSE1=samplevar_SSE1,var_SSE1=var_SSE1,SEbatachesmean[1],
            sampletm_SSE1=sampletm_SSE1,tm_SSE1=tm_SSE1,SEbatachesmean[1],samplefm_SSE1=samplefm_SSE1,fm_SSE1=fm_SSE1
  )
  allSSE_unstand<-c(SEbatachesmean[1],meansd1=meansd1,SEbatachesmean[1],samplevarsd1=samplevarsd1,varsd1=varsd1,SEbatachesmean[1],
                    sampletmsd1=sampletmsd1,tmsd1=tmsd1,SEbatachesmean[1],samplefmsd1=samplefmsd1,fmsd1=fmsd1
  )
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean_all1=se_mean_all1,allSE=allSE,allSSE=allSSE,allSE_unstan=allSE_unstan,allSSE_unstand=allSSE_unstand,SEbatachesmean=SEbatachesmean)
  
  allErrors
}


write.csv(simulatedbatchgnorm_ABSE_SE,paste("gnorm_ABSSE_error.csv", sep = ","), row.names = FALSE)


registerDoSEQ()






