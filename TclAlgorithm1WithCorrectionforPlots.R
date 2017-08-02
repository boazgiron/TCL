#Ver6
#library(stringr)
#library("mclust", lib.loc="~/R/win-library/3.2") 
#library(ggplot2)
#library(MASS)
#library("kernlab", lib.loc="~/R/win-library/3.2")
#library("numDeriv", lib.loc="~/R/win-library/3.2")
#col=rev(rainbow(255))
#col=(heat.colors(20))
#col  = terrain.colors(255)
#source("tools.R")

#Variable----------------

require(signal)
#Addition

library("kernlab", lib.loc="~/R/win-library/3.2")
library("e1071")

col  = topo.colors(255)
colbr = c("black","red","green","blue")
ty  = c(paste0("Area",1:8),"Peak9")
bf <- butter(2, 0.2, type="low")
bf4 <- butter(4, 0.2, type="low")
args<-commandArgs(TRUE)
PLOTS = TRUE
PLOTS2 = FALSE
UseValidationMode = TRUE
WRITE_RESUTLS = FALSE
LUDA = FALSE

#ACCEPTANCE CRITERIA
USE_ACCEPTANCE_CRITERIA = FALSE
bAddmedianInformation = TRUE
medianDat = NULL
QCResults = list(ErrorNum = 0,QCarra = NULL, ErrorString = NULL )



#FUNCTIONS------------------
#new020817

subformcol <-function(df,v){
  
  ou = df
  
  for(i in 1:dim(df)[2]){
    
    ou[,i] <- ou[,i] - v[i] 
  }
  
  ou
  
}

divformcol <-function(df,v){
  
  ou = df
  
  for(i in 1:dim(df)[2]){
    
    ou[,i] <- ou[,i]/v[i] 
  }
  
  ou
  
}


geG <- function(inputA,v1,v2,fscale){
  
  
  #inputA = filedatal_sel45
  #v1 = "Area6"
  #v2 = "Area1"
  #fscale = FALSE
  
  filedatalsSample = inputA[sample(nrow(inputA),1e3),]
  df = cbind(filedatalsSample[,v1],filedatalsSample[,v2])
  mes = colMeans(df,na.rm = T)
  mds = apply(df,2,sd,na.rm = T)
  if(fscale){
    sp = scale(df)
  }else{
    sp = df
  }

  sc <- specc(sp,kernel = "laplacedot",kpar = list(sigma= 0.4), centers = 2 )
  df = data.frame(sp,cl = factor(sc@.Data))
  #with(df,plot(X1,X2,col = cl,pch = 19 ,cex = 0.3,ylim = c(-1,7)))
  svm_model <- svm(cl ~ ., data = df)# ,probability = TRUE)
  sp1 = (cbind(inputA[,v1],inputA[,v2]))
  dfa = data.frame(sp1)
  if(fscale){
    
    dff = subformcol(dfa,mes)
    dff = divformcol(dff,mds)
  }else{
    
    dff = dfa
  }
  
  dff[is.infinite(dff[,1]),1] = NA
  dff[is.infinite(dff[,2]),2] = NA
  dff = dff[complete.cases(dff),]
  pre = predict(svm_model,dff)#,probability = TRUE)
  
  #sp1[is.infinite(sp1[,1]),1] = NA
  #sp1 = sp1[complete.cases(sp1),]
  
  if(fscale){
    
    dff = divformcol(dff,1/mds)
    dff = subformcol(dff,-mes)
  }
  
  g1x = mean(dff[pre == 1,1])# + mes[1] 
  g1y = mean(dff[pre == 1,2])# + mes[2] 
  
  g2x = mean(dff[pre == 2,1])# + mes[1] 
  g2y = mean(dff[pre == 2,2])# + mes[2] 
  
  
  # p = attr(pre, "probabilities")
  # 
  # 
  # selg2 = pre == 2
  # selg1 = pre == 1
  # 
  # g1x = sum(p[selg1,1] * dff[selg1,1])/sum(p[selg1,1])
  # g1y = sum(p[selg1,2] * dff[selg1,2])/sum(p[selg1,2])
  # 
  # g2x = sum(p[selg2,1] * dff[selg2,1])/sum(p[selg2,1])
  # g2y = sum(p[selg2,2] * dff[selg2,2])/sum(p[selg2,2])
  # 
  # points(c(g1x,g2x),c(g1y,g2y),pch = 19 ,cex =3 ,col =3)
  # points(c(g1xx,g2xx),c(g1yy,g2yy),pch = 19 ,cex =5 ,col = 3 )
  # 
  
  if(g1x > g2x ){
    
    pre = ifelse(pre == 1,2,1)
    
  }
  
  if(PLOTS){
    
    png(filename= paste0(plotpath,v1,"_",v2,"_SVM.png"))
    plot(dff,col = pre,pch = 19 ,cex = 0.3)
    points(c(g1x,g2x),c(g1y,g2y),col = 3, cex = 2,pch = 19)
    dev.off()
  }
  
  pre
  
}


GetLabel <- function(m){
  
  se = dim(m)
  
  lables = matrix(rep(0,se[1]*se[2]),se[1],se[2])
  
  components = NULL
  
  for(i in 1:se[1] ){
    
    for(j in 1:se[2] ) {
      
      if ( m[i,j] != 0 ){
        
        up = ifelse(i > 1,lables[i-1,j],0) 
        
        be = ifelse(j > 1,lables[i,j-1],0)
        
        if ( (up == 0) & (be == 0 )) {
          
          lables[i,j]  = max(unique(lables)) + 1
          
        }else if(be == up ) {
          
          lables[i,j] = be
          
        }else if (( be != 0)&(up == 0 ))  {
          
          lables[i,j] = be
          
        }else if(( up != 0)&(be == 0 ))  {
          
          lables[i,j] = up
          
        }else if( (up != 0) & (be != 0 ) & (be != up ) ) {
          
          lables[i,j]  = max(unique(lables)) + 1
          
          lables[(lables == up ) | (lables == be )] = lables[i,j]
          
          #components = c(components,list(c(lables[i,j], be , up)))
        }
      }
    }  
  }
  
  for( j in  1:length(components) ){
    
    v = unlist(components[j])
    lables[lables %in% v[2:length(v)]] = v[1]
    
  }
  
  lables
}

#linear desciriminate analysis
flda<- function(u1,u2,da1,da2){
  
  
  n1 = dim(da1)[1]
  n2 = dim(da2)[1]
  
  C1 = cov(da1)
  C2 = cov(da2)
  
  C = (1/(n1+n2)) * ( n1*C1 + n2*C2 )
  
  C_1 = solve(C)
  
  betaout = C_1 %*% ( u1 - u2 )
  
  betaout
  
}

fldaProjection <- function(x,u1,u2,betain,d1,d2){
  
  u = (u1+u2)/2
  
  re = apply(x,1,function(y){ (y -u) %*% betain }) 
  
  re 
  
}  


findbycut <- function(maxrv,h,alfa,type1){
  
  #  maxrv = h_vpArea1FCS_max$maxl
  # h = h_vpArea1FCS
  #alpa = 5e-2
  #type1 = 1
  
  vmax = h$counts[ maxrv ]
  if (type1 == 0 ) {
    
    if (maxrv > 1){
      
      vmin  = min( h$counts[1:maxrv] )
      vcut = vmin  + (vmax - vmin) * alfa
      mincutind = min( which( h$counts[1:maxrv] >= vcut ) ) 
      out = 2 * maxrv - mincutind  
    }else{
      
      out = 1
      
    }
    
    
    
  }else {
    
    if (length(h$counts) > maxrv ){
      
      vmin  = min( h$counts[maxrv:length(h$counts)] )
      vcut = vmin  + (vmax - vmin) * alfa
      wv = which( h$counts[maxrv:length(h$counts)] <= vcut )
    
      if(length(wv) > 0 ){
        
        mincutind = min( wv )
        
      }else{
        
        mincutind = length( maxrv:length(h$counts) )
        
      }
      
      
      out = maxrv - mincutind
    
    }else{
      
      out  = length(h$counts)
    }
  }
  
  out
  #h_vpArea1FCS$mids[maxrv + mincutind]
}  

#Add 210617DD

SimpleLeftCut <- function(maxrv,h,alfa){
  
  vcut = h$counts[ maxrv ] * alfa
  #min(which(h[maxrv:1]$counts < vcut))
  w = which(h$counts[ maxrv:1 ] < vcut)
  
  if(length(w) > 0 ){
    backindex = min(w)
    out = max(maxrv - backindex,1)
  }
  else{
    
    out = 1
  }
  
  out 
}


SimpleRightCut <- function(maxrv,h,alfa){
  
  maxrv = h_lm68hv_max$maxl$indx
  h = h_lm68hv
  alfa = 4e-1
  
  vcut = h$counts[ maxrv ] * alfa
  #min(which(h[maxrv:1]$counts < vcut))
  w = which(h$counts[maxrv:length(h$mids)] < vcut)
  
  if(length(w) > 0 ){ 
   
     backindex = min(w)
     out = min(maxrv + backindex,length(h$mids))
     
  }else{
    
    out = length(h$mids)
  }
  
  out 
}


getCartrigeNameAdnDevice <-function(filePath){
  
  filename = strsplit(filePath,"/")[[1]]
  filepart = strsplit(filePath,"_")[[1]]
  selDeviceS = grep("AX",filepart)
  DeviceNum <- as.numeric((strsplit(filepart[selDeviceS], "[^[:digit:]]")[[1]])[3])
  gc = grep("C",filepart)
  Cartnum = filepart[gc[!(gc %in% grep("TC",filepart))]]
  Cartnum = as.numeric((strsplit(Cartnum, "[^[:digit:]]")[[1]])[2])
  out = list(DeviceNum = DeviceNum,Cartnum = Cartnum)
  out
}

average2 <- function(x){ 
  
  buf = NULL
  
  for(i in 1:(length(x)-1)){
    
    buf = c(buf,mean(x[i:(i+1)]))
  }
  
  buf
}


filter1 <- function(flitercoff,delay,x,y){
  out  = list(x = x[(delay+1):length(x)]  ,y = filter( flitercoff, y )[(delay+1):length(x)])
  out 
}

getDensity <- function(x, si , fc ){
  #x = sumCh
  #si =200
  #fc = bf
  h = hist(x,si,plot = F )
  maxXhY = max(h$counts)
  hd = density(x)
  hd$y = filter(fc,hd$y)
  hd$y = maxXhY * hd$y/max(hd$y)
  hd
} 



maxdetcetion <- function(sqin,maxpre,maxDiff,minpre,minDiff,plotFlag,xin,...){
  
  sq = sqin - min(sqin)+ 0.1 
  maxb = TRUE
  raxt  = sq[1]
  raxtIndex = 1 
  maxindx = NULL
  minindx = NULL
  
  i =0 
  for(p in sq ){
    i =i+1 
    #p = sq[i]
    if (maxb == TRUE ) {      
      if ( p > raxt ) {
        
        raxtIndex = i
        raxt = p
        
      }else if( ( p/raxt < maxpre ) & ( abs( p - raxt ) > maxDiff  ) ) {
        maxindx =  c(maxindx ,raxtIndex )
        maxb = FALSE
      }
    }else{ 
      
      if ( p < raxt ) {
        
        raxtIndex = i
        raxt = p
        
      }else if( ( ( p/ raxt ) > minpre ) & ( abs( p - raxt ) > minDiff  ) )  {
        minindx =  c(minindx ,raxtIndex )
        maxb = TRUE
      }   
      
    }
    
  }
  
  
  
  
  #Add last min
  if( length(maxindx) > length(minindx) ){
    
    xs = maxindx[length(maxindx)]:length(sqin)
    minlast  = min(sqin[xs])[1]
    minlastindex =  which(minlast == sqin[xs])[1]
    minindx = c(minindx, minlastindex + maxindx[length(maxindx)]  - 1 )
  }
  
  #Add first min
  xs = 1:maxindx[1]
  minfirst  = min(sqin[xs])[1]
  minfirstindex =  which(minfirst == sqin[xs])[1]
  minindx = c(minindx, minfirstindex )
  
  if(length(minindx) > 0 ){
    
    minindx = minindx[order(minindx)]
  }
  
  
  if( plotFlag == TRUE ){
    
    plot(xin,sqin,...)#,ylim =c(-10,15))
    points(xin[maxindx],sqin[maxindx],pch =2,col = 2,cex=2)
    points(xin[minindx],sqin[minindx],pch =3,col = 3,cex=2)
    
  }
  
  
  maxl  = data.frame(indx = maxindx,x = xin[maxindx],y = sqin[maxindx])
  mins  = data.frame(indx = minindx,x = xin[minindx],y = sqin[minindx])
  
  
  list(maxl = maxl,mins = mins)
}

#oldErrorVersion
errorH1 <- function(UseValidationMode1,ErrorString,plotpath1){}

errorH <- function(UseValidationMode1,ErrorString,plotpath1,QcResultsV){
  
  if(UseValidationMode1 == TRUE) {
    
    Variable = c(
      "Device",
      "Cartrige",
      "CD45_total",
      "CD45_live",
      "CD45_dead",
      "CD3_total",
      "CD3_live",
      "CD3_dead",
      "CD4_live",
      "CD8_live",
      "Cd45liveRatio",
      "Cd3liveRatio",
      "Cd3Cd45Ratio",
      #"CD4CD8Ratio",
      "CD4CD3Ratio",
      "CD8CD3Ratio",
      "CD4CD8Ratio",
      "ErrorNum",
      "WarningNums")
    
    Values  = rep(0,length(Variable))
    Values[length(Values)-1]  = QcResultsV$ErrorNum
    Values  = as.character( Values )
    
    write(ErrorString,paste0(plotpath1,"ErrorText.txt"))  
    
    df = data.frame(Values,Variable)
    df
  }
  else{
    Variable = c(
      #"Device",
      "Cartrige",
      "CD45_total",
      "CD45_live",
      "CD45_dead",
      "CD3_total",
      "CD3_live",
      "CD3_dead",
      "CD4_live",
      "CD8_live",
      "Cd45liveRatio",
      "Cd3liveRatio",
      "Cd3Cd45Ratio",
      #"CD4CD8Ratio",
      "CD4CD3Ratio",
      "CD8CD3Ratio",
      "CD4CD8Ratio",
      "ErrorNum",
      "WarningNums") 
    
    Values  = rep(0,length(Variable))
    Values[length(Values)-1]  = QcResultsV$ErrorNum
    Values  = as.character( Values )
    
    df = data.frame(Values,Variable)
    
    paste0(df[,c("Values")],collapse = ",")
    
  }
  
}

maxiNu = function(x,y,PLOTS1,plotpath,filename,DEVOFF,minmax ,filter_length  ){
  
  #x= xx
  #y = yy
  #Debug
  #x= h_lm68hv$mids
  #y = h_lm68hv$counts
  #minmax = 30
  #plot(x,y,"l")
  #plot(x[5:(length(x)-6)] ,buff,"l")
  #abline(h = 0)
  #filter_length = 6
  filter_length_1 = filter_length - 1
  filter_middle = filter_length/2
  buff = NULL
  xs  = NULL
  for( i in 1:(length(x) - filter_length)){
    
    x1 = x[i:(i+filter_length_1)]
    y1 = y[i:(i+filter_length_1)]
    xs = c(xs,mean(y1))
    model  = lm(y1~ x1)
    buff = c(buff, model$coefficients[2])
    
  }
  
  le = length(buff)
  
  selup  = which(buff[1:(le-1)] >= 0  & buff[2:le] <= 0 )
  
  if ( length(selup) > 0 ){
    
    selup  = selup[ y[selup+filter_middle] > minmax ]
    
  }
  
  if ( length(selup) > 0 ){
    
    xinterval = findInterval(x,x[selup+filter_middle])
    seldown = NULL
    ii  = unique(xinterval)[1]
    for(ii in unique(xinterval) ) {
      
      ind = which(ii == xinterval)
      yii  = y[ind]
      miny = min(yii)
      ind  = ind[which(miny == yii)[1]]
      seldown = c(seldown,ind)
      
    }
    
  }
  
  #remove low max
  if ( length(selup) > 0 ){
    
    bufiin = NULL
    
    for (iin in 1:length(selup) ){
      
      upv = y[selup[iin]+filter_middle]
      dov1 = y[seldown[iin]]
      dov2 = y[seldown[iin + 1]]
      if(  ( (upv  - dov2) < minmax ) & ( (upv  - dov1) < minmax ) ){
        
        bufiin = c(bufiin , iin )
      }
    }
  
    if( !is.null(bufiin)) {
      
      selup = selup[ -bufiin ]
      
      if ( length(selup) > 0 ){
        
        xinterval = findInterval(x,x[selup+filter_middle])
        seldown = NULL
        ii  = unique(xinterval)[1]
        for(ii in unique(xinterval) ) {
          
          ind = which(ii == xinterval)
          yii  = y[ind]
          miny = min(yii)
          ind  = ind[which(miny == yii)[1]]
          seldown = c(seldown,ind)
          
        }
        
      }
      
    }
  }
  
  #no max find simple max
  if ( length(selup) == 0 ){
    
    selup = which(max(y) == y)[1]
    selup = selup[y[selup] >  0.8 * minmax ]
    
    if ( length(selup) > 0 ){
      
      xinterval = findInterval(x,x[min(selup+filter_middle,length(x))])
      seldown = NULL
      ii  = unique(xinterval)[1]
      for(ii in unique(xinterval) ) {
        
        ind = which(ii == xinterval)
        yii  = y[ind]
        miny = min(yii)
        ind  = ind[which(miny == yii)[1]]
        seldown = c(seldown,ind)
        
      }
      
      selup = selup - filter_middle
    }

  }  
  
  if(PLOTS1) {
    
    png(filename= paste0(plotpath,filename,".png"))
    
    plot(x,y,main = filename ,"l" )
    
    if ( length(selup) > 0 ){
      
      points(x[selup+filter_middle],y[selup+filter_middle],col = 4, pch =19 , cex = 2)
      points(x[seldown],y[seldown],col = 2, pch =19 ,cex =2)
      
    }
    
    if( DEVOFF ){
      
      dev.off()
    }
    
  }else{
    
    if(PLOTS2){
    
      #plot(x,y,"l" )
      
      if ( length(selup) > 0 ){
        
        #points(x[selup+filter_middle],y[selup+filter_middle],col = 4, pch =19 , cex = 2)
        #points(x[seldown],y[seldown],col = 2, pch =19 ,cex =2)
      }
    }
    
  }
  
  
  if ( length(selup) > 0 ){
    maxl  = data.frame(indx = selup+filter_middle,x = x[selup+filter_middle],y = y[selup+filter_middle])
    mins  = data.frame(indx = seldown,x = x[seldown],y = y[seldown])
  
  }else{
    
    maxl  = data.frame(indx = NULL,x = NULL,y = NULL)
    mins  = data.frame(indx = NULL,x = NULL,y = NULL)
    
  }
  
  list(maxl = maxl,mins = mins)
  
}

#Error Handle function ----------------------------------------------------------

checkType  = c("QCNumberOf45","QC45","QC3","QCDraq7","QC4","QC8")
RatioMaxlimits = c(1e8,1e8,1e8,1e8,1e8,1e8)
RatioMinlimits = c(1.0,1.6,1.8,1.5,1.0,4.0)
sizelimits = c(1000,300,300,300,300,300)
Ratiolimits =  data.frame(checkType,RatioMinlimits,RatioMaxlimits)

sizeLimits = data.frame(Type = c("NumberOf45","NumOfCD3Live","NumOfCD3LiveNegative","NumCD3Dead","NumCD4","NumCD8"),Value = c(1000,300,300,300,300,300) ,stringsAsFactors = F)

#main function
main <- function(wrkingFilepath){
  
#In Function --------------------
  #QCNumberOf45!
QCumberOf45f <- function(QCResultin){
  
  resultsString = NULL
  
  if ( dim(filedatal45pl)[1] <= sizeLimits[sizeLimits$Type == "NumberOf45","Value"] ) {
    
    #Error low number of CD45
    QCResultin$ErrorNum = 1
    QCResultin$ErrorString = paste0("NumberofCD45 = ",dim(filedatal45pl)[1]," < ",sizeLimits[sizeLimits$Type == "NumberOf45","Value"])
  }
  
  QCResultin
}

#QC45!
QC45f<- function(QCResultin){
  
  resultsString = NULL
  #Not enough data unable to perform  QC
  if( NumOfCD3Live <= sizeLimits[sizeLimits$Type == "NumOfCD3Live","Value"]  ) {
    
    QCResultin$QCarra = c(QCResultin$QCarra,1)
    
    # Check color CD45 
  }else{
    
    CD3Area1 = medianDat[medianDat$Type == "CD3","Area1"]
    CD3Area6 = medianDat[medianDat$Type == "CD3","Area6"]
    
    ratio = CD3Area1/CD3Area6
    
    if ( ratio < Ratiolimits[Ratiolimits == "QC45","RatioMinlimits"] | ratio > Ratiolimits[Ratiolimits == "QC45","RatioMaxlimits"]  ){
      
      #Error color CD45
      QCResultin$ErrorNum = 2
      QCResultin$ErrorString = paste0("CD45 = ",round(ratio,2)," < ",Ratiolimits[Ratiolimits == "QC45","RatioMinlimits"] )
    }
    
  }
  
  QCResultin
  
}

#QC3!
QC3f<- function(QCResultin){
  
  #Not enough data unable to perform  QC
  if( NumOfCD3LiveNegative < sizeLimits[sizeLimits$Type == "NumOfCD3LiveNegative","Value"] |  NumOfCD3Live < sizeLimits[sizeLimits$Type == "NumOfCD3Live","Value"] ) {
    
    QCResultin$QCarra = c(QCResultin$QCarra,2)
    
    # Check color CD3   
  }else {
    CD3Area6 = medianDat[medianDat$Type == "CD3","Area6"]
    NoCD3Area6 = medianDat[medianDat$Type == "NoCD3","Area6"]
    ratio = CD3Area6 / NoCD3Area6
    if ( ratio < Ratiolimits[Ratiolimits == "QC3","RatioMinlimits"] | ratio > Ratiolimits[Ratiolimits == "QC3","RatioMaxlimits"]  ){
      
      #Error color CD3  
      QCResultin$ErrorNum = 3
      QCResultin$ErrorString = paste0("QC3 ratio = ",round(ratio,2)," < ",Ratiolimits[Ratiolimits == "QC3","RatioMinlimits"] )
    }
  }
  
  QCResultin
  
}


#QCDraq7!
QCDraq7f<- function(QCResultin){
  
  #Not enough data unable to perform  QC
  if( NumCD3Dead <  sizeLimits[sizeLimits$Type == "NumCD3Dead","Value"] | NumOfCD3Live < sizeLimits[sizeLimits$Type == "NumOfCD3Live","Value"] ) {
    
    QCResultin$QCarra = c(QCResultin$QCarra,3)
    
    # Check color QCDraq7   
  }else{
    
    CD3LiveArea8 = medianDat[medianDat$Type == "CD3Live","Area8"]
    CD3DeadArea8 = medianDat[medianDat$Type == "CD3Dead","Area8"]
    ratio = CD3DeadArea8 / CD3LiveArea8
    if ( ratio < Ratiolimits[Ratiolimits == "QCDraq7","RatioMinlimits"] | ratio > Ratiolimits[Ratiolimits == "QCDraq7","RatioMaxlimits"]  ){
      
      #Error color Draq7
      QCResultin$ErrorNum = 4
      QCResultin$ErrorString = paste0("QC Draq7 ratio = ",round(ratio,2)," < ",Ratiolimits[Ratiolimits == "QCDraq7","RatioMinlimits"] )
    }
  }
  
  QCResultin
  
}

#QC4!
QC4f<- function(QCResultin){
  
  #Not enough data unable to perform  QC
  if( NumCD4 < sizeLimits[sizeLimits$Type == "NumCD4","Value"] | NumOfCD3LiveNegative < sizeLimits[sizeLimits$Type == "NumOfCD3LiveNegative","Value"] ) {
    
    QCResultin$QCarra = c(QCResultin$QCarra,5)
    
    # Check color CD4 
  }else{
    
    NoCD3Area3 = medianDat[medianDat$Type == "NoCD3","Area3"]
    CD4Area3 = medianDat[medianDat$Type == "CD4","Area3"]
    ratio = CD4Area3 / NoCD3Area3
    if ( ratio < Ratiolimits[Ratiolimits == "QC4","RatioMinlimits"] | ratio > Ratiolimits[Ratiolimits == "QC4","RatioMaxlimits"]  ){
      
      #Error color CD4
      QCResultin$ErrorNum = 6
      QCResultin$ErrorString = paste0("QC4 ratio = ",round(ratio,2)," < ",Ratiolimits[Ratiolimits == "QC4","RatioMinlimits"] )
    }
  }
  
  QCResultin
  
}


#QC8!
QC8f<- function(QCResultin){
  
  #Not enough data unable to perform  QC
  if(NumCD8 < sizeLimits[sizeLimits$Type == "NumCD8","Value"] | NumOfCD3LiveNegative < sizeLimits[sizeLimits$Type == "NumOfCD3LiveNegative","Value"] ) {
    
    QCResultin$QCarra = c(QCResultin$QCarra,7)
    
    # Check color CD8   
  }else{
    
    NoCD3Area4 = medianDat[medianDat$Type == "NoCD3","Area4"]
    CD8Area4 = medianDat[medianDat$Type == "CD8","Area4"]
    ratio = CD8Area4 / NoCD3Area4
    if ( ratio < Ratiolimits[Ratiolimits == "QC8","RatioMinlimits"] | ratio > Ratiolimits[Ratiolimits == "QC8","RatioMaxlimits"]  ){
      
      #Error color CD8
      QCResultin$ErrorNum = 8
      QCResultin$ErrorString =  paste0( "QC8 ratio = ",round(ratio,2)," < ",Ratiolimits[Ratiolimits == "QC8","RatioMinlimits"] )
    }
  }
  
  QCResultin
}

QCDumyf<- function(QCResultin){
  
  QCResultin
  
}

QcCheck <- function(checktype1,QCResultin){
  
  QCResult1 = switch( checktype1,
                      "QCNumberOf45"=  QCumberOf45f(QCResultin),
                      "QC45" = QC45f(QCResultin),
                      "QC3" = QC3f(QCResultin),
                      "QCDraq7" = QCDraq7f(QCResultin),
                      "QC4" = QC4f(QCResultin),
                      "QC8" = QC8f(QCResultin),
                      QCDumyf(QCResultin)
  )
  
  QCResult1
}

RemoveLowSumCH <- function(){
  
  SumCH_MinMAX = 10
  sumCh = filedatals$sumCh 
  hsumch  = hist(sumCh,200,plot = F )
  fr = filter1(bf,5,hsumch$mids, hsumch$counts)
  hsumch$mids = fr$x
  hsumch$counts = fr$y
  fmax_hsumch = maxiNu(hsumch$mids,hsumch$counts,PLOTS,plotpath,"sumCh",FALSE,30,10)
  
  #Try 190717DD 
  if(FALSE){
    
    hsumch = getDensity(sumCh,200,bf)
    plot(hsumch$x,hsumch$y,"l")
    
    
   
    
    reu = NULL
    alfa =  0.05
    r = hsumch$y[1]
    for( i in 1:length(hsumch$y) ){
      r = alfa* hsumch$y[i]  + (1-alfa)*r  
      reu = c(reu,r)
      
    }
    
    which(max(reu) == reu)
    which(max(hsumch$y) == hsumch$y)
    
    
    points(hsumch$x,reu,"l",col = 2)
    
    hsumch$y = reu
    
    fmax_hsumch = maxiNu(hsumch$x,hsumch$y,PLOTS,plotpath,"sumCh1",TRUE,30,10)
    
    d = diff( fmax_hsumch$maxl$x ) / diff( range(fmax_hsumch$maxl$x) )
    ch = which( d  < 0.2 )
    
    if( length(ch) > 0 ){
      
      for( i in ch ){
        
        
        
      }    
      
    }
    
    dev.off()
    
  }
  
  
  if(length(fmax_hsumch$maxl) ==  0 ){
    dev.off()
    fmax_hsumch = maxiNu(hsumch$mids,hsumch$counts,PLOTS,plotpath,"sumCh",FALSE,30,10)
  }
  
  
  bHighMiddle = FALSE
  
  if(length(fmax_hsumch$maxl) > 0 ){
    
    if( length(fmax_hsumch$maxl$x) > 1  ) {
      
      bHighMiddle = fmax_hsumch$mins$y[2] > SumCH_MinMAX 
      
      if ( FALSE ){
        
        h  = hist((sumCh[sumCh < bondaryv & fmax_hsumch$maxl$x[1] < sumCh]),100,plot = F)
        df = data.frame(x = h$mids,y = h$counts)
        le = length(df$x)
        bu = NULL
        i=2
        for(i in 2:(le-2)){
          
          bu  = c(bu ,sd(df$y[(i+1):le]) - sd(df$y[1:i]) )
          
        }
        
        #plot(df$x[2:(le-2)],bu)
        sumcut85  = 0.85 * max(bu,na.rm = T)
        maxlocation = which(max(bu) == bu)
        selcut85  = min(maxlocation + which( bu[maxlocation:length(bu)] < sumcut85))
        bondaryv1 = fmax_hsumch$maxl$x[1] + df$x[ selcut85 ]
        
      }else{
        
        bondaryv = fmax_hsumch$mins$x[2]  
        arr = sumCh[ sumCh < bondaryv]
        bondaryv1  = quantile(arr,0.75)
      }
      
    }else{
      
      bondaryv = hsumch$mids[SimpleLeftCut(fmax_hsumch$maxl$indx[1],hsumch,5e-2)]
      bondaryv1 = bondaryv
    }
    
  }else{
    
    bondaryv = 1 
    bondaryv1 = 1
    #Add in the  070617DD
    #NO max case
    
  }
  
  #Old Code---------
  if(FALSE)
  {
    
    
    #--
    
    #delete 210617DD
    #le = dim(filedatal)[1]
    #samp  = sample(1:le,min(le,2e4))
    #with(filedatal[samp,],plot(sumCh,Area7 , pch = 19, cex= 0.2, 
    #                           xlab = "sumCh",ylab = "Area7"))
    
    qi =  quantile(filedatal$sumCh,c(0.1,0.9))
    x1 = median(filedatal$sumCh[filedatal$sumCh < qi[1]])
    x2 = median(filedatal$sumCh[filedatal$sumCh > qi[2]])
    y1 = median(filedatal$Area7[filedatal$sumCh < qi[1]])
    y2 = median(filedatal$Area7[filedatal$sumCh > qi[2]])
    x2 = 1/(4*(y2 - y1)/(x2 - x1))
    x2 = min(10,x2)
    vp = c(1,x2)/sqrt(1 + x2^2 )
    vv  = c(1,-1/x2)/sqrt(1 +  ( 1/x2 )^2 )
    sumCh = cbind( filedatal$sumCh,filedatal$Area7 ) %*% (vp)
    
    #delete--
    
    #orignal-- 
    # qi =  quantile(filedatal$sumCh,c(0.1,0.9))
    # x1 = median(filedatal$sumCh[filedatal$sumCh < qi[1]])
    # x2 = median(filedatal$sumCh[filedatal$sumCh > qi[2]])
    # y1 = median(filedatal$FCS[filedatal$sumCh < qi[1]])
    # y2 = median(filedatal$FCS[filedatal$sumCh > qi[2]])
    # x2 = 1/(4*(y2 - y1)/(x2 - x1))
    # x2 = min(10,x2)
    # 
    # #x2 = 10  # 3/20
    # vp = c(1,x2)/sqrt(1 + x2^2 )
    # vv  = c(1,-1/x2)/sqrt(1 +  ( 1/x2 )^2 )
    # 
    # sumCh = cbind( filedatal$sumCh,filedatal$FCS ) %*% (vp)
    #--
    #plot(filedatal$sumCh,filedatal$FCS,pch = 19, cex = 0.2)
    
    #Remove low sumCH events----
    if(PLOTS){
      
      hsumch  =  hist(sumCh,200,,plot = F )
      
    }else{
      
      hsumch  =  hist(sumCh,200,plot = F )
    }
    
    
    fr = filter1(bf,5,hsumch$mids, hsumch$counts)
    hsumch$mids = fr$x
    hsumch$counts = fr$y
    
    fmax_hsumch = maxiNu(hsumch$mids,hsumch$counts,PLOTS,plotpath,"sumCh",FALSE,30,10)
    
    #Ver2.1
    #min of at  least 15%
    #sellowPre = fmax_hsumch$maxl$pre < 0.85
    #fmax_hsumch$maxl  = fmax_hsumch$maxl[sellowPre,]
    #fmax_hsumch$mins = fmax_hsumch$mins[sellowPre,]
    #
    #Change if 210617DD 
    #if ( length(fmax_hsumch$maxl$x)  %in% c(2,3) ) {
    if ( length(fmax_hsumch$maxl$x)  > 1 ) {
      
      #middle12 = hsumch$counts[fmax_hsumch$maxl[1]:fmax_hsumch$maxl[2]]
      #minind =  which(min(middle12) == middle12 ) + fmax_hsumch$maxl[1] - 1
      #bondaryv = hsumch$mids[minind[1]]
      bondaryv = fmax_hsumch$mins$x[2] 
      
    }else{
      
      #Change 210617DD
      #cutv = fmax_hsumch$maxl$y[1]* 0.05
      #20.6 from + to -
      #bondaryv  = hsumch$mids[fmax_hsumch$maxl$indx[1] + min(which(hsumch$counts[fmax_hsumch$maxl$indx[1]:length(hsumch$counts)] < cutv)) - 1]
      #bondaryv  = hsumch$mids[fmax_hsumch$maxl$indx[1] - min(which(hsumch$counts[fmax_hsumch$maxl$indx[1]:length(hsumch$counts)] < cutv)) - 1]
      bondaryv = hsumch$mids[SimpleLeftCut(fmax_hsumch$maxl$indx[1],hsumch,5e-2)]
      #--
      
    }
    
    
  }
  #--------------------
  
  if(PLOTS) {
    
    abline(v = bondaryv1 , col = 2 ,lty = 2 ,lwd =2)
    dev.off()
  }
  
  out = sumCh > bondaryv1
  #filedatal = filedatals[sel_sumCh,]
  
  if(PLOTS) {
    
    png(filename= paste0(plotpath,"SumCh_Area7 boundary bondaryv1",".png"))
    le = dim(filedatals)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatals[samp,],plot(sumCh,Area7 , pch = 19, cex= 0.2, 
                                xlab = "sum FL",ylab = "Area7",
                                col = colbr[1+ out[samp]],
                                ylim = c(0,5),
                                main = paste0("SumCh_Area7 boundary bondaryv1 ", CarName)))
    dev.off()
  }
  
  out 
  
}

Sel45f <- function(){
  
  qi =  quantile(filedatal$Area1,c(0.05,0.95))
  px1 = median( filedatal$Area1[filedatal$Area1 < qi[1]] )
  px2 = median( filedatal$Area1[filedatal$Area1 > qi[2]] )
  py1 = median( filedatal$FCS[filedatal$Area1 < qi[1]] )
  py2 = median( filedatal$FCS[filedatal$Area1 > qi[2]] )
  x2 = 1/(4*(py2 - py1)/(px2 - px1))
  x2 = min(10,x2)
  vp = c(1,x2)/sqrt(1 + x2^2 )
  vv  = c(1,-1/x2)/sqrt(1 +  ( 1/x2 )^2 )
  
  
  Area1_fcs = cbind( filedatal$Area1,log10(filedatal$FCS) ) %*% (vp)
  h_Area1_fcs  = hist(Area1_fcs,200,plot = F )
  fr = filter1(bf,5,h_Area1_fcs$mids, h_Area1_fcs$counts)
  h_Area1_fcs$mids = fr$x
  h_Area1_fcs$counts = fr$y
  fmax_h_Area1_fcs = maxiNu(h_Area1_fcs$mids,h_Area1_fcs$counts,PLOTS,plotpath,"h_Area1_fcs",TRUE,30,10)
  if(length(fmax_h_Area1_fcs$maxl) > 0 ){
    
    if( length(fmax_h_Area1_fcs$maxl$x) > 1 ){
      
      bounadry_h_Area1_fcs = fmax_h_Area1_fcs$mins$x[2]
    }else{
      
      bounadry_h_Area1_fcs  = h_Area1_fcs$mids[SimpleLeftCut(fmax_h_Area1_fcs$maxl$indx[1],h_Area1_fcs,5e-2)]
      
    }
  }
  
  
  
  cl =  data.frame(cluster = (1 + (Area1_fcs > bounadry_h_Area1_fcs))) 
  if(PLOTS) {
    
    png(filename= paste0(plotpath,"pre 45 Area1_FCS ",".png"))
    le = dim(filedatal)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal[samp,],plot(Area1,FCS , pch = 19, cex= 0.2, 
                               xlab = "Area1",ylab = "FCS",
                               col = colbr[cl$cluster[samp]],
                               #ylim = c(0,5) ,
                               main = paste0("Area1_FCS ", CarName)))
    points(c(px1,px2),c(py1,py2),pch  = 19,cex = 2,col = 5) 
    dev.off()
  }
  
  
  m1 = mean(filedatal$Area1[cl$cluster == 1])
  m2 = mean(filedatal$Area1[cl$cluster == 2])
  if(m1 > m2){
    
    out = cl$cluster == 1
    
  }else{
    
    out = cl$cluster == 2
  }
  
  out
  
}

Sel6f<- function(){
  
  h_Area6_sel45  = hist(filedatal_sel45$Area6,200,plot = F )
  fr = filter1(bf,5,h_Area6_sel45$mids, h_Area6_sel45$counts)
  h_Area6_sel45$mids = fr$x
  h_Area6_sel45$counts = fr$y
  fmax_h_Area6_sel45 = maxiNu(h_Area6_sel45$mids,h_Area6_sel45$counts,PLOTS,plotpath,"Cd3 h_Area6",FALSE,30,10)
  if(length(fmax_h_Area6_sel45$maxl) > 0 ){
    
    if( length(fmax_h_Area6_sel45$maxl$x) > 1 ){
      
      bounadry_h_Area1_fcs = fmax_h_Area6_sel45$mins$x[2]
    }else{
      
      bounadry_h_Area1_fcs  = h_Area6_sel45$mids[SimpleLeftCut(fmax_h_Area6_sel45$maxl$indx[1],h_Area6_sel45,5e-2)]
      
    }
  }
  
  if(PLOTS) {
    
    abline(v = bounadry_h_Area1_fcs , col = 2 ,lty = 2 ,lwd =2)
    dev.off()
  }
  
  
  out = filedatal$Area6  > bounadry_h_Area1_fcs
  
  
  if(PLOTS) {
    
    png(filename= paste0(plotpath,"cd3 Area6_Area1 ",".png"))
    sel6_ = out[sel45]  
    le = dim(filedatal_sel45)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal_sel45[samp,],plot(Area6,Area1 , pch = 19, cex= 0.2, 
                                     xlab = "Area6",ylab = "Area1",
                                     col = colbr[1 + sel6_[samp]] ,
                                     #ylim = c(0,5) ,
                                     main = paste0("CD3 Area6_Area1 ", CarName)) )
    dev.off()
  }
  
  out
}

SelPreDeadf <- function(){
  
  lm68 = with(filedatal45pl[filedatal45pl[,"Area6"] > 1,],lm(Area8 ~ Area6 ))
  x2 = lm68$coefficients[2]
  vp = c(1,x2)/sqrt(1 + x2^2 )
  vv  = c(1,-1/x2)/sqrt(1 +  ( 1/x2 )^2 )
  
  # create parallel for cutting low dots
  vpArea6Area8 =  cbind( filedatal45pl$Area6,filedatal45pl$Area8 ) %*% (vp)
  
  if(PLOTS){
    h_vpArea6Area8 = hist(vpArea6Area8,100,plot = F)
  }else{
    h_vpArea6Area8 = hist(vpArea6Area8,100,plot = F)
  }
  
  fr = filter1(bf,5,h_vpArea6Area8$mids, h_vpArea6Area8$counts)
  h_vpArea6Area8$mids = fr$x
  h_vpArea6Area8$counts = fr$y
  h_vpArea6Area8_max = maxiNu(h_vpArea6Area8$mids,h_vpArea6Area8$counts,PLOTS,plotpath,"h_vpArea6Area8 SelPreDead",FALSE,30,10)
  
  le = length(h_vpArea6Area8_max$maxl$indx)
  
  if(le >= 2){
    
    h_lm68hpmiddsleMin = h_vpArea6Area8_max$mins$x[2] 
    
  }else if ( le == 1 ) {
    
    hmaxs  = SimpleLeftCut(h_vpArea6Area8_max$maxl$indx[1],h_vpArea6Area8,2e-2)
    
    h_lm68hpmiddsleMin =  h_vpArea6Area8$mids[hmaxs]
    
  }else{
    
    print(paste0( "le = length(h_vpArea6Area8_max$maxl$indx )le =",le ))
    #Problem
  }
  
  
  if(PLOTS){
    abline(v = h_lm68hpmiddsleMin,col  = 3,lwd = 2,lty =2)
    dev.off()
  }
  
  vpArea6Area8 =  cbind( filedatal$Area6,filedatal$Area8 ) %*% (vp)
  out  = vpArea6Area8 > h_lm68hpmiddsleMin
  
  if(PLOTS){
    
    png(filename= paste0(plotpath,"SelPreDead",".png"))#plot11
    le = dim(filedatal45pl)[1]
    samp  = sample(1:le,min(le,2e4))
    
    with(filedatal[samp,],plot(Area6,Area8 , pch = 19, cex= 0.2, 
                               xlab = "Area6",ylab = "Area8",
                               main = paste0("SelPreDead Cut tail ", CarName),
                               col = colbr[1+out[samp,]]))
    
    abline(lm68,col = "blue",lwd = 2)
    dev.off()
  }
  
  out
}

SelDeadf <- function(){
  
  selpredead = SelPreDeadf()
  
  #Cut tail select up Area6 and Are8-------------------------- 
  filedatal45pl_predead = filedatal[selpredead & sel45,]
  
  #FIRST iteration
  le = dim(filedatal45pl_predead)[1]
  samp  = sample(1:le,min(le,4e4))
  Area6 = c(rep(-2,1e4),filedatal45pl_predead[samp,"Area6"])
  Area8 = c(rep(-2,1e4),filedatal45pl_predead[samp,"Area8"])
  lm68s = lm(Area8 ~ Area6 )
  
  
  #Detect Dead------------------------------------------------
  
  x2 = lm68s$coefficients[2]
  vp = c(1,x2)/sqrt(1 + x2^2)
  vv  = c(1,-1/x2)/sqrt(1 + (1/x2)^2)
  vvs = cbind( filedatal45pl_predead$Area6,filedatal45pl_predead$Area8) %*% vv
  
  if(PLOTS){
    h_lm68hv = hist(vvs,200,plot = F )
  }else{
    h_lm68hv = hist(vvs,200,plot = F )
  }  
  
  fr =  filter1(bf,1,h_lm68hv$mids,h_lm68hv$counts)
  h_lm68hv$mids = fr$x
  h_lm68hv$counts = fr$y
  
  h_lm68hv_max = maxiNu(h_lm68hv$mids,h_lm68hv$counts,PLOTS,plotpath,"h_lm68hv Detect Deads",FALSE,30,10)
  
  if( length(h_lm68hv_max$maxl) == 0 ){
    
    dev.off()
    h_lm68hv_max = maxiNu(h_lm68hv$mids,h_lm68hv$counts,PLOTS,plotpath,"h_lm68hv Detect Deads second",FALSE,20,10)
    
  }
  
  #No max
  if( length(h_lm68hv_max$maxl) != 0 ){
    
    if(length(h_lm68hv_max$maxl$x) == 1){
      
      cut0.02ind = max(1,SimpleLeftCut(h_lm68hv_max$maxl$indx,h_lm68hv,2e-1) - 1)
      
    }else{
      
      cut0.02ind = h_lm68hv_max$mins$indx[length(h_lm68hv_max$mins$indx) - 1]
    }
    
  }else{
    
    cut0.02ind = 1
    
  }
  
  
  
  seldeadv = h_lm68hv$mids[cut0.02ind]
  
  if(PLOTS){
    
    abline(v = seldeadv ,col = 2,lwd = 2,lty = 2)
    dev.off()
  }
  
  
  vvs = cbind( filedatal$Area6,filedatal$Area8) %*% vv
  seldead = vvs < seldeadv
  
  vvsAll = cbind( filedatal$Area6,filedatal$Area8) %*% vv
  out = ( vvsAll < seldeadv ) & selpredead
  
  #DeadNum  = sum(seldeadAll & sel45 )
  
  if(PLOTS){
    
    png(filename= paste0(plotpath,"Dead",".png"))
    le = dim(filedatal[selpredead,])[1]
    samp  = sample(1:le,min(le,2e4))
    filedatalp = filedatal[selpredead,]
    
    with( filedatalp[samp,],plot( Area6 , Area8 , pch = 19, cex= 0.2, 
                                  xlab = "Area6",ylab = "Area8",
                                  xlim = c( 0,quantile(Area6,0.999,na.rm = T)) , 
                                  ylim = c( 0.1,max(Area8,na.rm = T) ) ,
                                  main = paste0("Area8 vs Area6 ", CarName),
                                  col = colbr[1+(seldead[selpredead])[samp]]))
    
    legend("topleft",c("live","Dead"),col = c(1,2),pch=19 )
    #abline(lm68live,col  = 5, lwd = 2 )
    dev.off()
    
  }
  
  out
}

SelNonTailCD48f <- function(){
  
  if(PLOTS){  
    
    abline(lm34,col = 3,lwd = 2)
    dev.off()
    
  }
  
  x2 = lm34$coefficients[2]
  vp = c(1,x2)/sqrt(1 + x2^2)
  vv  = c(1,-1/x2)/sqrt(1 + (1/x2)^2)
  vps34 = with(filedatalCD3,cbind(Area4,Area3) %*% vp)
  
  if(PLOTS){
    
    h_vps34  = hist(vps34,200,plot = F)
  }else{
    
    h_vps34  = hist(vps34,200,plot = F)
  }
  
  fr =  filter1(bf,1,h_vps34$mids,h_vps34$counts)
  h_vps34$mids = fr$x
  h_vps34$counts = fr$y
  
  h_lm68hv_max = maxiNu(h_vps34$mids,h_vps34$counts,PLOTS,plotpath,"h_vps34 Cut Tail",FALSE,30,10)
  if( length(h_lm68hv_max$maxl) == 0 ){
    
    dev.off()
    h_lm68hv_max = maxiNu(h_vps34$mids,h_vps34$counts,PLOTS,plotpath,"maxNe_h_vps34",FALSE,20,10)
  }
  
  leMax = length(h_lm68hv_max$maxl$x)
  
  if(leMax > 0 ){
    
    if ( leMax >= 3 ) {
      
      cutpointIndex = h_lm68hv_max$mins$indx[2]
      cutpoint = h_lm68hv_max$mins$x[2]
      if(PLOTS){
        
        text( cutpoint - 0.1,h_lm68hv_max$maxl$y[2],"Tail")
      }
      
      #no tail 
    }else {
      
      # limit the cut by 1
      cutpoint = h_lm68hv_max$mins$x[1]
      cutpointIndex = h_lm68hv_max$mins$indx[1]
      if( PLOTS){
        
        text( h_lm68hv_max$maxl$x[1] -0.2 ,h_lm68hv_max$maxl$y[1],"CD4")
      }
      
    }
    
    if( PLOTS){
      abline(v = cutpoint,col = 1,lwd = 2,lty = 2 )
      dev.off()
    }
  }else{
    
    cutpoint = min(vps34)
    if( PLOTS){
      
      abline(v = cutpoint,col = 1,lwd = 2,lty = 2 )
      
      dev.off()
    }
    
  }  
  
  
  out =  vps34  >  cutpoint
  out
  
}

SelType40r8f <-function(){
  
  bUseRatio = FALSE
  ratio3vd4 = filedatalCD3_sel_tail$Area3 - filedatalCD3_sel_tail$Area4
  h34 = hist(ratio3vd4,100,plot = F)
  h_34max = maxiNu(h34$mids,h34$counts,PLOTS,plotpath,"Detect h34 by ratio ch3-ch4",TRUE,50,10)
  if( length(h_34max$maxl$indx) == 2 ){
    
    
    middelValue = mean(h_34max$mins$x)
    ratio3vd4 = filedatalCD3$Area3 - filedatalCD3$Area4
    out = ratio3vd4 < middelValue
    
    if(PLOTS){
      png(filename= paste0(plotpath,"CD3_Selection CD4_CD8 USE Ratio",".png"))
      with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = 0.2,col = colbr[1+!out[samp]],main = paste0(CarName," CD3_Selection CD4_CD8 USE Ratio")))
      dev.off()  
    }
    
    bUseRatio = TRUE
  }
  
  if(!bUseRatio ) {
    
    x2 = lm34$coefficients[2]
    vp = c(1,x2)/sqrt(1 + x2^2)
    vv  = c(1,-1/x2)/sqrt(1 + (1/x2)^2)
    vvs34woTail = with(filedatalCD3_sel_tail,cbind(Area4,Area3) %*% vv)
    
    if(PLOTS){
      #change new max 210617DD
      #png(filename= paste0(plotpath,"vvs34woTail",".png"))
      #h_vvs34woTail  = hist(vvs34woTail,100)
      h_vvs34woTail  = hist(vvs34woTail,150,plot = F)
      
    }else{
      
      h_vvs34woTail  = hist(vvs34woTail,100,plot = F)
      
    } 
    
    fr =  filter1(bf,1,h_vvs34woTail$mids,h_vvs34woTail$counts)
    h_vvs34woTail$mids = fr$x
    h_vvs34woTail$counts = fr$y
    
    h_vvs34woTailvvs_max = maxiNu(h_vvs34woTail$mids,h_vvs34woTail$counts,PLOTS,plotpath,"h_vvs34woTail Detect CD4 CD8",FALSE,10,6)
    
    minbtween = 0 #Version4 Continue 
    #find the minmum between max ( asume 2 )
    #change new max 250617DD
    # if( length(h_vvs34woTailvvs_max$maxl$indx) >= 2){
    #   
    #   minbtween = h_vvs34woTailvvs_max$mins$x[2]
    # }
    
    lenMax = length(h_vvs34woTailvvs_max$maxl$indx)
    
    if( lenMax > 2){
      
      #Algcorrection 
      #bounary4 = h_vvs34woTailvvs_max$mins$x[lenMax]
      bounary4 = h_vvs34woTailvvs_max$mins$x[lenMax - 1]
      bounary8 = h_vvs34woTailvvs_max$mins$x[2]
      
      if( PLOTS ){
        
        abline(v = bounary4, col = 4, lwd = 2 ,lty = 2)
        abline(v = bounary8, col = 4, lwd = 2 ,lty = 2)
      }
      
      sel4  = vvs34woTail > bounary4
      sel8  = vvs34woTail < bounary8
      
      me4  = mean(vvs34woTail[sel4])
      sd4  = sd(vvs34woTail[sel4])
      me8  = mean(vvs34woTail[sel8])
      sd8  = sd(vvs34woTail[sel8])
      
      sdall = sd8+sd4
      minbtween = me8 + (me4  - me8)*sd4/sdall   
    }
    
    if( lenMax == 2){
      
      minbtween = h_vvs34woTailvvs_max$mins$x[2]
      
      sel4  = vvs34woTail > minbtween
      sel8  = vvs34woTail < minbtween
      
    }
    
    #in thecase only one is detect
    if( lenMax == 1){
      
      ra = range(vvs34woTail)
      
      if(h_vvs34woTailvvs_max$maxl$x > ( ra[1] + 0.6 * ( ra[2] - ra[1]) ) ){
        
        minbtween = findbycut(h_vvs34woTailvvs_max$maxl$indx[1],h_vvs34woTail,5e-2,1)
        minbtween = h_vvs34woTail$mids[minbtween]
        
      }else if(h_vvs34woTailvvs_max$maxl$x < ( ra[1] + 0.4 * ( ra[2] - ra[1]) ) ){
        
        minbtween = h_vvs34woTail$mids[findbycut(h_vvs34woTailvvs_max$maxl$indx[1],h_vvs34woTail,5e-2,0)] 
      }
      else{
        
        minbtween = h_vvs34woTailvvs_max$mins$x[2]
      }
      
    }
    
    if( lenMax == 0  ) {
      
      med3  = median( filedatalCD3_sel_tail$Area3 )
      med4  = median( filedatalCD3_sel_tail$Area4 )
      
      diffmedian = med3 - med4
      
      # 4 is bigger than 3
      if ( diffmedian < 0 ){
        
        #All CD8 
        minbtween =  min(vvs34woTail)
        
      }else{
        
        #All CD8 
        minbtween =  max(vvs34woTail)
        
      }
      
    }
    
    if(PLOTS){
      
      abline(v = minbtween, col = 3, lwd = 2 ,lty = 2)
      dev.off()
      
    }
    
    #Create CD4 CD8 selection 
    vvs34 = with(filedatalCD3,cbind(Area4,Area3) %*% vv)
    #create between 4 and 8 selection variable
    out = vvs34 > minbtween
    
    le = dim(filedatalCD3)[1]
    samp  = sample(1:le,min(le,1e4))
    
    if(PLOTS){
      
      png(filename= paste0(plotpath,"CD3_Selection CD4_CD8",".png"))
      with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = 0.2,col = colbr[1+!out[samp]],main = paste0(CarName," First Cut for CD4 and CD8")))
      dev.off()  
    }
    
  }
  
  
  out
  
}

sel_doubleNegtovef <- function(){
  
  if(!bType8) 
  {
    
    
    le = dim(filedatalCD3only4)[1]
    samp  = sample(1:le,min(le,1e4))
    da1 = rbind(filedatalCD3only4[,c("Area1","Area3")],data.frame(Area1 =rep(1,1e4),Area3 = rep(1,1e4)))
    lmonly4 = with(da1, lm(Area3 ~ Area1))
    
    
    x2 = lmonly4$coefficients[2]
    vp = c(1,x2)/sqrt(1 + x2^2)
    vv  = c(1,-1/x2)/sqrt(1 + (1/x2)^2)
    
    vvslmonly4 = with(filedatalCD3only4,cbind(Area1,Area3) %*% vv)
    if(PLOTS){
      h_vvslmonly4  = hist(vvslmonly4,200,plot = F)
    }else{
      h_vvslmonly4  = hist(vvslmonly4,200,plot = F)
    }
    
    fr =  filter1(bf4,1,h_vvslmonly4$mids,h_vvslmonly4$counts)
    h_vvslmonly4$mids = fr$x
    h_vvslmonly4$counts = fr$y
    
    h_vvslmonly4_max = maxiNu(h_vvslmonly4$mids,h_vvslmonly4$counts,PLOTS,plotpath,"h_vvslmonly4 Detection DN",FALSE,10,10)
    
    if ( length(h_vvslmonly4_max$maxl) > 0 )
    {
      
      
      if (length(h_vvslmonly4_max$maxl$indx) > 2 ) {
        cutIndex = length(h_vvslmonly4_max$mins$indx) - 1
        doubleNegtoveBoundary = h_vvslmonly4_max$mins$x[cutIndex]
        if(PLOTS){
          text( h_vvslmonly4_max$maxl$x[cutIndex]- 0.1,h_vvslmonly4_max$maxl$y[cutIndex],"DN")
        }  
        
      }else if  (length(h_vvslmonly4_max$maxl$indx) == 2 ) {
        
        cutIndex = 2
        doubleNegtoveBoundary = h_vvslmonly4_max$mins$x[cutIndex]
        
        if(PLOTS == TRUE){
          
          if(length(h_vvslmonly4_max$maxl$indx) == 2 ){
            
            text( h_vvslmonly4_max$maxl$x[cutIndex] - 0.05,h_vvslmonly4_max$maxl$y[cutIndex],"DN")  
          }
        }
      }else{
        
        doubleNegtoveBoundary  = h_vvslmonly4_max$mins$x[2]
        
      }
      
      
    }
    else
    {
      doubleNegtoveBoundary = max(vvslmonly4)
    }
    
    if(PLOTS){
      abline(v = doubleNegtoveBoundary, col = 2 , lwd =2 , lty = 2 )
      dev.off()
    }
    
    sel_doubleNegtoveBoundary  = vvslmonly4 > doubleNegtoveBoundary
    
    if(PLOTS){
      
      png(filename= paste0(plotpath,"CD3only4 DN Detection",".png"))#plot24
      le = dim(filedatalCD3only4)[1]
      samp  = sample(1:le,min(le,1e4))
      with(filedatalCD3only4[samp,], plot(Area1,Area3,
                                          pch = 19, cex = 0.2,col = colbr[2- !sel_doubleNegtoveBoundary[samp]],
                                          main = paste0( CarName , " Select Double Negative ",CarName)))
      legend("topleft",c("CD4","DN"),col = c(1,2),pch=19 )
      dev.off()
    }
    
    
    SumdoubleNegtoveBoundary  = sum(sel_doubleNegtoveBoundary) 
    
    vvslmonly4forall = with(filedatalCD3,cbind(Area1,Area3) %*% vv)
    out  = vvslmonly4forall > doubleNegtoveBoundary
    
    
    
  }
  else
  {
    
    out = rep( FALSE, dim(filedatalCD3)[1] )
    
  }
  
  out
}

FindDoublepositive <- function(){
  
  selpreDP  = !sel_doubleNegtoveBoundaryIndForAll
  
  if(PLOTS){
    png(filename= paste0(plotpath,"FinalArea4Area3",".png"))
    le = dim(filedatalCD3)[1]
    samp  = sample(1:le,min(le,1e4))
    
    with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = 0.2,col = colbr[1+(selpreDP + 2*selType40r8 )[samp]],main = paste0(CarName," Pre Double Positive ") ))
    legend("topleft",c("CD4DN","CD4","LoWCD8","CD8"),col = c(1,2,3,4),pch=19 )
    dev.off()
  }
  
  if( !bOnlyOneType48 )
  {
    
    
    preDP =  filedatalCD3[selpreDP,c("Area4","Area3")]
    preDPselmin = selType40r8[selpreDP]
    
    g1 =  preDP[preDPselmin,]
    g2 =  preDP[!preDPselmin,]
    
    g1m = colMeans(g1,na.rm = T)
    g2m = colMeans(g2,na.rm = T)
    
    us = cbind(g1m,g2m)
    fl = flda(us[,1],us[,2],g1,g2)
    
    if ( sum(is.na(fl[,1])) == 0 ){
      
      fldaP = fldaProjection(preDP,us[,1],us[,2],fl)
      
      
      if(PLOTS){
        png(filename= paste0(plotpath,"fldaPForDoublePositive",".png"))
        h_fldaP  = hist(fldaP,200,main = paste0(CarName," LDA for finding the Double Positive points"),plot = T )
      }else{
        h_fldaP  = hist(fldaP,200,plot = F )
      }  
      d1 = dim(g1)[1]
      d2 = dim(g2)[1]
      p1 = d1/(d1+d2)
      p2 = d2/(d1+d2)
      vp = log(p2/p1)
      
      if(PLOTS){
        abline(v = vp, col = 3, lwd = 2 , lty = 2)
        
      }
      
      ind = min(which(h_fldaP$mids > vp )):(length( h_fldaP$mids))
      cutvalueUp =  min(max(h_fldaP$counts[ind])* 4e-2 + min(h_fldaP$counts[ind]),30)
      upind  = min(which(h_fldaP$counts[ind] > cutvalueUp ) ) 
      upind = ind[1] + upind
      
      if(PLOTS){
        abline(v = h_fldaP$mids[upind],col =2 , lwd = 2 ,lty = 2 )
      }
      
      ind = max(which(h_fldaP$mids < vp )) : 1
      cutvalueUDown =  min(max(h_fldaP$counts[ind])* 4e-2 + min(h_fldaP$counts[ind]),30)
      downind  = min(which(h_fldaP$counts[ind] > cutvalueUDown ) ) 
      downind = ind[1] - downind
      
      if(PLOTS){
        abline(v = h_fldaP$mids[downind],col =2 , lwd = 2 ,lty = 2 )
        dev.off()
      }
      
      filedatalCD3forProojection =  filedatalCD3[,c("Area4","Area3")]
      
      
      fldapreDPCD3 = fldaProjection(filedatalCD3forProojection,us[,1],us[,2],fl)
      
      out  = !((h_fldaP$mids[downind] < fldapreDPCD3 ) & (h_fldaP$mids[upind] > fldapreDPCD3 ))
      out = selpreDP & !out
      
      if(PLOTS){
        
        png(filename= paste0(plotpath,"DoublePositive_Area4_Area3",".png"))
        le = dim(filedatalCD3)[1]
        samp  = sample(1:le,min(le,2e4))
        with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = c(0.2,0.6)[1  + out[samp]],col = colbr[1  + out[samp]  ],main = paste0(CarName," Select Double Positive") ))
        dev.off()
      }
    }else{
      out = rep(FALSE,dim(filedatalCD3)[1])
    }
    
  }else{
    out = rep(FALSE,dim(filedatalCD3)[1])
  }
  
  out
  
}  

#--------------------------------  

#  ini Param -----------------------
  
CD45num = 0
CD45Livenum = 0 
CD45deadnum = 0
CD3num = 0 
CD3Livenum = 0 
CD3deadnum = 0 
selcleanCD4_NoDPnum = 0
selcleanCD8_NoDPnum = 0
Cd45liveRatio = 0 
Cd3liveRatio = 0
Cd3Cd45Ratio = 0  

#wrkingFilepath = fl#filepath1
#Load file-----------------------

  filedata = read.csv(wrkingFilepath,header = T)
  
  di = dim(filedata)
  
  #Extract Cartrige num  
  runInformation  = getCartrigeNameAdnDevice(wrkingFilepath)
  
  CarName = paste0("C",runInformation$Cartnum)
  
  Cartnum = runInformation$Cartnum
  
  plotpath = paste0("C:/Project/LeukoDx/LudaFacsValidation/Debug1207/",Cartnum,"/")
  
  if(PLOTS){  
    if(!dir.exists(plotpath)){
      dir.create(plotpath)
    }
  }
  
  #corrupt file Error
  if( !( ( di[1] > 1e3 ) & ( di[2] > 28 ) ) ){
    #Error Corrupt input file
    QCResults$ErrorString  =  "Corrupt input file"
    QCResults$ErrorNum = 99
    return  ( errorH(UseValidationMode,QCResults$ErrorString, plotpath,QCResults ) )
    
  }
  
#Remove low width--------------
  widthb = filedata$Width > 4
  
  if(PLOTS) {
    png(filename= paste0(plotpath,"Width",".png"))
    hist(filedata$Width,500,xlim = c(0,30),xlab = "Width", main = as.character(CarName))
    abline(v = 4,col ="red",lwd = 2 )
    dev.off()
  }
  
  filedata = filedata[widthb,]
#------------------------------  
  
  filedatals =  data.frame(apply(filedata[,ty[1:9]],2,log10))
  filedatals$FCS <- filedata$Peak9
  filedatals$sumCh =  apply(filedatals[,1:8],1,sum)
  
#Sumch Detection--------------
  #sel_sumCh = RemoveLowSumCH()
  
  sel_sumCh = geG(filedatals,"sumCh","FCS",TRUE) == 2

  Num_sel_sumCh  = sum( sel_sumCh )
  
  if(  Num_sel_sumCh < sizeLimits[sizeLimits$Type == "NumberOf45","Value"] )  { 
    #Return Error 1 low number of 45
    QCResults$ErrorNum = 1
    QCResults$ErrorString = paste0("ErrorNum =",QCResults$ErrorNum ,"  QCNumberOf45 ( Num_sel_sumCh ) =",Num_sel_sumCh )
    return(errorH(UseValidationMode,QCResults$ErrorString , plotpath,QCResults ) )
  }
  
  filedatal = filedatals[sel_sumCh,]
  
  if(FALSE){
    
    dim(filedatal)
    
    plot(d$y)
    d= density(filedatal$Area6)
    
    plot(d)
    
    h = with(filedatal,density(Area6))
    ro = HighDensityGroup(h$y,100)
    
    plot(ro[,2])
    
    HighDensityGroup(h$count,10)
    
    with(filedatal,plot(Area6,Area1))
    
    k <- with(filedatal,MASS::kde2d(Area6,Area1, n=50))
    h = hist(k$z,plot = F)
    ind = which(max(h$counts) == h$counts)
    q = 2*h$mids[ind]# + (h$mids[2]  + h$mids[1])/2
    #q = quantile(k$z,0.99)
    
    k$z[k$z < q ] = 0
    k$z[k$z != 0 ] = 3
    kz = GetLabel(k$z)
    le = length(k$x)
    za = matrix(kz,1,)
        
    df = data.frame( x = rep(k$x,each  = le),y = rep(k$y,le), za = t(matrix(kz,1,)))
      
    head(df)  
    
    fn <- function(l){c(median(l$x),median(l$y))}    
        
        sp = split(df,list(za))
        
    rs = sapply(sp,fn)[,-1]
              
      
      
    
    
    
str(za)    
za1 = data.frame(za = t(za))
unique(za1)



    image(k$x,k$y,kz)
    hist(kz)
  }
  
  #Orignal Algorithm 100717DD---------------------
  #sel45 =  Sel45f()
  sel45   = geG(filedatal,"Area1","Peak9",FALSE) == 2
  
  filedatal_sel45 = filedatal[sel45,]
  #sel6   = Sel6f()
  sel6   = geG(filedatal_sel45,"Area6","Area1",FALSE) == 2
  
  #selpreDead   = geG(filedatal_sel45,"Area6","Area8",FALSE) == 2
  #selDead   = geG(filedatal_sel45[selpreDead,],"Area6","Area8",FALSE) == 2
  
#Hirsh and Luda Code    
if(FALSE){
  
  #HirshAlgorithm 070517DD -----------------
  
  FCS_Area1 <- log10(with(filedatal,(10^FCS)*(10^Area1)))
  h_FCS_Area1  = hist((FCS_Area1),100,plot = F)
  fr = filter1(bf,2,h_FCS_Area1$mids, h_FCS_Area1$counts)
  h_FCS_Area1$mids = fr$x
  h_FCS_Area1$counts = fr$y
  fmax_h1 = maxiNu(h_FCS_Area1$mids,h_FCS_Area1$counts,PLOTS,plotpath,"h1",FALSE,30,5)
  cutmax = fmax_h1$maxl$x[length(fmax_h1$maxl$x)]
  if(PLOTS){
    xcut = cutmax - sd(FCS_Area1[FCS_Area1 > cutmax])
    abline(v = xcut)
    dev.off()
  }
 
  
  sel45Hirsh = FCS_Area1 > xcut
  h6 = with(filedatal[sel45Hirsh,],hist(Area6,100,plot = F))
  fr = filter1(bf,5,h6$mids, h6$counts)
  h6$mids = fr$x
  h6$counts = fr$y
  fmax_h6 = maxiNu(h6$mids,h6$counts,PLOTS,plotpath,"h6",FALSE,30,10)
  
  if(length(fmax_h6$maxl) > 0 ){
    
    if( length(fmax_h6$maxl$x) > 1 ){
  
      cut6 = fmax_h6$min$x[length(fmax_h6$min$x) - 1]
    
    }else{
      
      cut6 = h6$mids[SimpleLeftCut(fmax_h6$maxl$indx[1],h6,0.1)]
      
      
    }
  }else{
    
    cut6 = 0
    
  }
      
  if(PLOTS){
    
    abline(v = cut6)
    
  }
  
  sel6 = with(filedatal,Area6 > cut6)
  h1 = with(filedatal[sel45Hirsh & sel6, ], hist(Area1,100,plot = F) )
  fmax_h1 = maxiNu(h1$mids,h1$counts,PLOTS,plotpath,"h1",FALSE,30,5)
  
  if(length(fmax_h1$maxl) > 0 ){
    cut1 = h1$mids[SimpleLeftCut(fmax_h1$maxl$indx[1],h1,0.1)]
    with(filedatal[sel45Hirsh,],plot(Area6,Area1,pch  =19,cex =0.2))
  }else{
    
    cut1  = 0
    
  }
  
  
  
  #results
  
  
  #cd45pre = sum(sel45Hirsh)/dim(filedatal)[1]
  
  se  = seq(2,-2,-0.01)
  bu = NULL
  for(s in se){
    sel45Hirsh = FCS_Area1 > xcut+ s
    #with(filedatal[sel45Hirsh,],points(Area6,Area1,pch  =19,cex = 0.1, col = 2, ylim  = c(1,4) ))
    A = sum(filedatal[sel45Hirsh,"Area1"] < cut1 )
    B = sum( filedatal[sel45Hirsh,"Area1"] > cut1 )
    bu = c( bu , A/(A+B) ) 
  }
  
  
  xout  = (xcut+se)[min(which(0.02*(max(bu,na.rm = T) - min(bu,na.rm = T)) + min(bu,na.rm = T) <= bu))]
  
  if(PLOTS){
    png(filename= paste0(plotpath,"Check FCS_Area1",".png"))
    plot(xcut+se,bu)
    abline(v  =  xcut)
    dev.off()
  }
  
  
  sel45Hirsh = FCS_Area1 > xout
  cd3_45pre = sum(filedatal[sel45Hirsh,]$Area6 >  cut6)/dim(filedatal[sel45Hirsh,])[1]
  
  if(PLOTS){
    png(filename= paste0(plotpath,"Over shoot cut example ",".png"))
    with(filedatal[FCS_Area1 > (xout-0.5),],plot(Area6,Area1,pch  =19,cex =0.2))
    with(filedatal[sel45Hirsh,],points(Area6,Area1,pch  =19,cex =0.2,col = 2))
    abline(v = cut6)
    abline(h = cut1)
    dev.off()
  }
  
  if(PLOTS) {
    
    png(filename= paste0(plotpath,"FCS vs Area1 Hirsh Algorithm ",".png"))
    le = dim(filedatal)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal[samp,],plot(Area1,10^FCS , pch = 19, cex= 0.2, 
                               xlab = "Area1",ylab = "FCS",
                               col = colbr[1+ sel45Hirsh[samp]],
                               #xlim = c(0,20),ylim = c(1,5),
                               main = paste0("FCS vs Area1 Angle = ", slop ,"  ", CarName)))
    #abline(c(0,3/20))
    dev.off()
  }
  
  sel6Hirsh = filedatal$Area6 > cut6
  
  
  if(PLOTS) {
    
    png(filename= paste0(plotpath,"CD3 selection ",".png"))
    le = dim(filedatal)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal[sel45Hirsh,],plot(Area6,Area1 , pch = 19, cex= 0.2, 
                               xlab = "Area6",ylab = "Area1",
                               col = colbr[1+ sel6Hirsh[sel45Hirsh]],
                               #xlim = c(0,20),ylim = c(1,5),
                               main = paste0("CD3 selection")))
    #abline(c(0,3/20))
    dev.off()
  }
  

  #Luda Algorithm --------------------------------

  df = with(filedatal,data.frame(x = Area1,y = 10^FCS ))
  re = range(df$x)
  seq1 = seq(re[1],re[2],0.01)
  df$interval <- findInterval(df$x,seq1)  
  countinterval = aggregate(df$y,list(df$interval),length)
  inmedain = aggregate(df$y,list(df$interval),median)
  dfo = data.frame(countinterval[,1:2],inmedain[,2],seq1[countinterval[,1]])
  dfo = dfo[dfo[,2] > 30,]
  
  q70 = quantile(dfo[,4],0.70)
  sel75 = dfo[,4] > q70
  x = dfo[sel75,4]+0.05
  y = dfo[sel75,3]
  m1 = lm(y~x)
  #q20 = quantile(dfo[,4],0.20)
  #sel20 = dfo[,4] < q20
  #x = dfo[sel25,4]+0.05
  #y = dfo[sel25,3]
  #m2 = lm(y~x)
  
  #if( m2$coefficients[1]  < -500 | m2$coefficients[1]  > 1000 ){
    
    q20 = quantile(dfo[,4],0.20)
    m2$coefficients[1] = min(median(y),500)
    m2$coefficients[2] = 0
  #}
  
  if(PLOTS) {
    png(filename= paste0(plotpath,"Transition location.png")) 
    plot(dfo[,4]+0.05,dfo[,3],pch=2 )
    abline(m1)
    abline(m2)
    
  }
  
  xre = (m2$coefficients[1] - m1$coefficients[1] ) / (m1$coefficients[2] - m2$coefficients[2])
  
  if(PLOTS) {
    
    abline(v = xre,col = 2,lwd = 2, lty = 2) 
    dev.off()
  }  
  
  if( ( xre < re[1] ) | ( xre > re[2] ) )  {
  
    print( "xre Out of range")
    
  }
  
  #Direction
  slop =  -10^7/(2*m1$coefficients[2])
  b  =   -xre * slop 
  
  #xx = seq(1,3,0.01)
  #points(xx,xx*slop  + b,"l",lwd = 2)
  sel45Luda  = ( 10^filedatal$FCS >=  ( slop*filedatal$Area1 + b ) )  

  h6 = with(filedatal[sel45Luda,],hist(Area6,100,plot = F))
  fr = filter1(bf,5,h6$mids, h6$counts)
  h6$mids = fr$x
  h6$counts = fr$y
  fmax_h6 = maxiNu(h6$mids,h6$counts,PLOTS,plotpath,"h6",FALSE,30,5)
  
  if(length(fmax_h6$mins) > 0 ){
    
    if(length(fmax_h6$mins$x) > 2 ){
      cut6 = fmax_h6$mins$x[length(fmax_h6$mins$x)-1]
    }else{
      
      cut6 = quantile(filedatal[sel45Luda,]$Area6,2e-2)
    }
    
    sel6Luda = ( filedatal$Area6 > cut6 )
    cut6_1 = with(filedatal[sel6Luda,],quantile(Area1,0.05))
  }
  
  if(PLOTS) {
    
    abline(v = cut6)
    dev.off()
  }  

  seq1 = seq(-1.0,1.0,0.1)
  bu = NULL
  xu = NULL
  for(s in seq1 ){
  
    #Direction
    slop =  -10^7/(2*m1$coefficients[2])
    b  =   -(xre+s) * slop 
    
    #xx = seq(1,3,0.01)
    #points(xx,xx*slop  + b,"l",lwd = 2)
    sel45loop  = ( 10^filedatal$FCS >=   ( slop*filedatal$Area1 + b ) )
    bu = c(bu,with(filedatal[!sel6Luda & sel45loop,],length(Area1[Area1 < cut6_1])/length(Area1[Area1 > cut6_1])))
    xu = c(xu ,xre+s)
  }
  
  
  if(sum( bu > 0.05 ) > 0 ){
    xre = xu[min(which(bu < ( min(bu) + 0.05 *(max(bu) - min(bu)))))]
  }
  slop =  -10^7/(2*m1$coefficients[2])
  b  =   - xre * slop 
  #sel45  = ( 10^filedatal$FCS >=   ( slop*filedatal$Area1 + b ) )

  if(PLOTS) {
    
    png(filename= paste0(plotpath,"Transition location.png")) 
    plot(dfo[,4]+0.05,dfo[,3],pch=2 )
    abline(m1)
    abline(m2)
    abline(v = xre,col = 2, lwd = 2 ,lty = 2 )
    dev.off()
    
  }
    
    sel45Luda  = ( 10^filedatal$FCS >=   ( slop*filedatal$Area1 + b ) )
  
    if(PLOTS) {
      
      png(filename= paste0(plotpath,"Area1 vs Area6",".png"))

      slop =  -10^7/(2*m1$coefficients[2])
      b  =   -( xre - 0.5) * slop 
      sel45s  = ( 10^filedatal$FCS >=   ( slop*filedatal$Area1 + b ) )
      
      with(filedatal[sel45s,],plot(Area6,Area1,pch =19 ,col = 2,cex =0.2,xlim = c(0,3),ylim = c(1,4)))
      abline(v = cut6,col = 2,lwd =2 ,lty  = 2 )
      abline(h = cut6_1,col = 2,lwd =2 ,lty  = 2 )
      
      
      with(filedatal[sel45Luda,],points(Area6,Area1,pch =19 ,col = 1 ,cex =0.3,xlim = c(0,3),ylim = c(1,4)))
      abline(v = cut6,col = 2,lwd =2 ,lty  = 2 )
      abline(h = cut6_1,col = 1,lwd =2 ,lty  = 2 )
      
      dev.off()
    }
  

  if(PLOTS) {
    
    png(filename= paste0(plotpath,"FCS vs Area1",".png"))
    le = dim(filedatal)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal[samp,],plot(Area1,10^FCS , pch = 19, cex= 0.2, 
                               xlab = "Area1",ylab = "FCS",
                               col = colbr[1+ sel45Luda[samp]],
                               #xlim = c(0,20),ylim = c(1,5),
                               main = paste0("FCS vs Area1 Angle = ", slop ,"  ", CarName)))
    #abline(c(0,3/20))
    dev.off()
  }
  
  #------------------------------------
  
    
  if(FALSE){    
    if( LUDA ){
      sel45 = sel45Luda
      sel6 =  sel6Luda
    }else{
      sel45 = sel45Hirsh
      sel6 =  sel6Hirsh
    }
    
  }
}  
  
  filedatal45pl = filedatal[sel45,]
  selCD3 =  sel6 & sel45
  #--------------------
  #Old CODE 280617
  if(FALSE){
  #Change FCS to Area7 210617 delete in the end 
    
  if(PLOTS) {
    png(filename= paste0(plotpath,"sumCh_FCS",".png"))
    le = dim(filedatal)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal[samp,],plot(sumCh,Area7 , pch = 19, cex= 0.2, 
                               xlab = "sumCh",ylab = "Area7",
                               col = colbr[1+ sel_sumCh[samp]],
                               xlim = c(0,20),ylim = c(1,5),
                               main = paste0("Area7 vs Area1 Angle = ", round(x2,2) ,"  ", CarName)))
    #abline(c(0,3/20))
    dev.off()
  }
  
  
  #Original
  # if(PLOTS) {
  #   png(filename= paste0(plotpath,"sumCh_FCS",".png"))
  #   le = dim(filedatal)[1]
  #   samp  = sample(1:le,min(le,2e4))
  #   with(filedatal[samp,],plot(sumCh,FCS , pch = 19, cex= 0.2, 
  #                            xlab = "sumCh",ylab = "FCS",
  #                            col = colbr[1+ sel_sumCh[samp]],
  #                            xlim = c(0,20),ylim = c(1,5),
  #                            main = paste0("FCS vs Area1 Angle = ", round(x2,2) ,"  ", CarName)))
  #   #abline(c(0,3/20))
  #   dev.off()
  # }
  # 
  # 
  
  filedatal = filedatal[sel_sumCh,]
  }
  
  
  #Old CODE 280617---------------------
  if(FALSE){
  #detect 45 High Area1 -----------------------
  if( PLOTS) {
    h_Area1 = with(filedatal,hist(Area1,200,plot = F))
  }else{
    h_Area1 = with(filedatal,hist(Area1,200,plot = F))
  }
  
  
  fr = filter1(bf,5,h_Area1$mids, h_Area1$counts)
  h_Area1$mids = fr$x
  h_Area1$counts = fr$y
  #Ver2.1
  #h_Area1_max  = findmax1(h_Area1$mids, h_Area1$counts,20,5,5,20,main = CarName)
  
  #change new max 210617DD
  #mindiff = min(max(h_Area1$counts)/10,100)
  #h_Area1_max = maxdetcetion(h_Area1$counts,0.85,2*mindiff,1.2,mindiff,PLOTS,h_Area1$mids,main = CarName)
  h_Area1_max = maxiNu(h_Area1$mids,h_Area1$counts,PLOTS,plotpath,"h_Area1",FALSE,30,10)
  
  
  le = length(h_Area1_max$maxl$x)
  
  if( le == 2 ) {
    
    ind = h_Area1_max$maxl$indx[1]:h_Area1_max$maxl$indx[2]
    indmin = which(min(h_Area1$counts[ind]) == h_Area1$counts[ind]) + h_Area1_max$maxl$indx[1]
    cutLowArea1 = h_Area1$mids[indmin]
    
    se = seq(cutLowArea1-0.2,cutLowArea1+0.2,0.01)
    bu = NULL
    for(l in se){  
      selArea1 = filedatal$Area1 > l
      filedatal45pl = filedatal[selArea1,]
      bu  = c(bu,dim(filedatal45pl)[1])
    }
    
    diffbu = diff(bu)
    
    if(PLOTS) {
      
      plot(se[2:length(se)],diffbu,'b',pch = 19,cex = 0.5)
    }
    
    
    
    #diffbu = average2(diffbu)
    diffbu = filter1(bf,3,1:length(diffbu), diffbu)$y
    
    if(PLOTS) {
      points(se[3:(length(se)-2)],diffbu,'l',pch = 19,col = 2)
    }
    
    maxlocation = which(max(diffbu) == diffbu)
    
    if(PLOTS) {
      points(se[min(maxlocation+2,length(se))],diffbu[min(maxlocation,length(diffbu))],col =3, lwd = 2,cex =3 )
    }
    
    cutLowArea1 =  se[ min(maxlocation+2,length(se)) ]
    
  }else if ( le == 1 ){
    
    indmin = which.max(cumsum(h_Area1$counts)/sum(h_Area1$counts) > 0.02)
    cutLowArea1 = h_Area1$mids[ indmin ]
    
  }
  
  
  #Debug 280617DD
  
  #selArea1 = filedatal$Area1 > cutLowArea1
  selArea1 = filedatal$Area1 > 2.5
  
  #
  
  if(PLOTS) {
    
    abline(v = cutLowArea1,col = 2 ,lty = 2, lwd =2)
    dev.off()
  }
  
  filedatal45pl = filedatal[selArea1,]
  
  if(PLOTS) {
    
    png(filename= paste0(plotpath,"Area1",".png"))
    le = dim(filedatal)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal[samp,],plot(Area1,FCS,pch =19 ,cex = 0.2 ,col = colbr[1 + selArea1[samp]] ))
    dev.off()
  }
  
  }

#Find Dead---------------
seldeadAll = SelDeadf()

#QC45! 
if(bAddmedianInformation){
  
  NumOfCD3Live = sum((!seldeadAll) & selCD3)
  
  addata = data.frame( testNum = Cartnum,Type = "CD3" , size = dim(filedatal[(!seldeadAll) & selCD3,])[1], t(apply(filedatal[(!seldeadAll) & selCD3,c(paste0("Area",1:8),"FCS")],2,median)) )
  medianDat = rbind( medianDat , addata )
  addata = data.frame( testNum = Cartnum,Type = "NoCD3" , size = dim(filedatal[(!seldeadAll) & (!selCD3),])[1] , t(apply(filedatal[(!seldeadAll) & (!selCD3),c(paste0("Area",1:8),"FCS")],2,median)) )
  medianDat = rbind( medianDat , addata )
  
  QCResults = QcCheck("QCNumberOf45",QCResults) 
  
  if(QCResults$ErrorNum != 0 ){
    
    if ( USE_ACCEPTANCE_CRITERIA){
      return  ( errorH(UseValidationMode,QCResults$ErrorString, plotpath,QCResults ) )
    }
    
  }
  
  addata = data.frame( testNum = Cartnum,Type = "CD45",size = dim(filedatal[sel45,])[1] , t(apply(filedatal[sel45,c(paste0("Area",1:8),"FCS")],2,median)) )
  medianDat = rbind( medianDat , addata )
  
  addata = data.frame( testNum = Cartnum,Type = "Junk",size = dim(filedatal[!sel45,])[1] , t(apply(filedatal[!sel45,c(paste0("Area",1:8),"FCS")],2,median)))
  medianDat = rbind( medianDat , addata )
}

#CD48----------------------------------
if(FALSE){

  #Find CD3
  if(PLOTS){
    #change new max 210617DD
    #png(filename= paste0(plotpath,"FD45pl_Area6",".png"))
    h_hist6 = hist(filedatal45pl$Area6,200,plot = F)
    
  }else{
    
    h_hist6 = hist(filedatal45pl$Area6,200,plot = F)
  }
  
  fr =  filter1(bf4,1,h_hist6$mids,h_hist6$counts)
  h_hist6$mids = fr$x
  h_hist6$counts = fr$y
  
  #change new max 210617DD
  #mindiff = min(max(h_hist6$counts)/10,100)
  #h_lm68hp_max = maxdetcetion(h_hist6$counts,0.85,2*mindiff,1.2,mindiff,PLOTS,h_hist6$mids,main = paste0(CarName," Find CD3 Hist of the CH6 Take the high heap"),xlab = "Area6 Index")
  h_lm68hp_max = maxiNu(h_hist6$mids,h_hist6$counts,PLOTS,plotpath,"h_hist6_CD45_Area6 CD3 Detection",FALSE,30,10)
  h_hist6_max_right = h_lm68hp_max$maxl$indx[length(h_lm68hp_max$maxl$indx)]
  
    
  if( length(h_lm68hp_max$maxl$indx) == 1 ) {
    
    #Change 210617DD
    #cut0.02ind = findbycut(h_hist6_max_right,h_hist6,1.5e-2,1)  
    cut0.02ind  = SimpleLeftCut(h_hist6_max_right,h_hist6,1.5e-2)
    CD3selectionBoundary = h_hist6$mid[cut0.02ind]
    
  }else{# new Code Version 3
    
    minmiddle = h_lm68hp_max$mins$y[2]
    cut0.02ind  = h_lm68hp_max$mins$indx[2]
    CD3selectionBoundary = h_lm68hp_max$mins$x[2]
  
  }
  
  #--
  
  if(PLOTS){
    abline(v = CD3selectionBoundary, col = 3, lwd = 2,lty = 2 )
    dev.off()
    
    png(filename= paste0(plotpath,"CD3 Detection",".png"))
    le = dim(filedatal45pl)[1]
    samp  = sample(1:le,min(le,2e4))
    
    plot(filedatal45pl$Area6[samp],filedatal45pl$Area1[samp] , pch = 19, cex= 0.2, 
         xlab = "Area6",ylab = "Area1",
         xlim = c( 0,quantile(filedatal45pl$Area6[samp],0.999)) , ylim = c( 0.1,max(filedatal45pl$Area1[samp]) ) ,
         main = paste0("Area1 vs Area6 ", CarName," Selected CD3"),
         col = colbr[1+( filedatal45pl$Area6[samp] > CD3selectionBoundary ) ] )
    
    abline(v = CD3selectionBoundary, col = 3, lwd = 2 ,lty = 2)
    legend("topleft",c("~CD3","CD3"),col = c(1,2),pch=19 )
    dev.off()
  }
  
  selCD3 = filedatal45pl$Area6 > CD3selectionBoundary
}

#QC3!
if(bAddmedianInformation){
  
  NumOfCD3LiveNegative = sum((!seldeadAll) & (!selCD3))
  
  QCResults = QcCheck("QC3",QCResults) 
  
  if(QCResults$ErrorNum != 0 ){
    
    if ( USE_ACCEPTANCE_CRITERIA){
      
      return(errorH(UseValidationMode,QCResults$ErrorString , plotpath,QCResults))
    }
  }
}

selCD3Live =  (!seldeadAll) & selCD3
selCD3Dead =   seldeadAll & selCD3

#QCDraq7!
if(bAddmedianInformation){
  
  addata = data.frame( testNum = Cartnum,Type = "CD3Live" , size = dim(filedatal[selCD3Live,])[1], t(apply(filedatal[selCD3Live,c(paste0("Area",1:8),"FCS")],2,median)) )
  medianDat = rbind( medianDat , addata )
  addata = data.frame( testNum = Cartnum,Type = "CD3Dead" , size = dim(filedatal[selCD3Dead,])[1], t(apply(filedatal[selCD3Dead,c(paste0("Area",1:8),"FCS")],2,median)) )
  medianDat = rbind( medianDat , addata )
  
  NumCD3Dead = sum( selCD3Dead )
  
  QCResults = QcCheck("QCDraq7",QCResults)
  if(QCResults$ErrorNum != 0 ){
    
    if ( USE_ACCEPTANCE_CRITERIA){
      return(errorH(UseValidationMode,QCResults$ErrorString, plotpath,QCResults ) ) 
    }
  }
}  

CD3Livenum = sum( selCD3Live )
CD3deadnum = sum( selCD3Dead )
CD3num  = sum(selCD3)
CD45num = sum(sel45)
CD45deadnum = sum(seldeadAll & sel45 )
CD45Livenum  = CD45num   - CD45deadnum 
Cd3liveRatio =  CD3Livenum/CD3num
Cd45liveRatio =  CD45Livenum/CD45num
Cd3Cd45Ratio = CD3Livenum/CD45Livenum

#Area3 vs Area4 select the high Area6
filedatalCD3 = filedatal[selCD3Live,]

if(PLOTS){
  
  png(filename= paste0(plotpath,"CD3 Direction",".png"))
  le = dim(filedatalCD3)[1]
  samp  = sample(1:le,min(le,2e4))
  with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = 0.2,main = paste0(CarName," CD4 CD8 Create speration direction" )))
}

da1 = rbind(filedatal45pl[selCD3,c("Area4","Area3")],data.frame(Area4 =rep(1,1e4),Area3 = rep(1,1e4)))
lm34 = with(da1, lm(Area3 ~ Area4))

#Create parallel for to cut tail for better differ cd4 and cd8 ----------------------
selNonTailCD48 = SelNonTailCD48f()
filedatalCD3_sel_tail = filedatalCD3[selNonTailCD48,]

#Detect 4 and 8-------------------------------------------------------------
# Use filedatalCD3_sel_tail to find 4 and by dividing
selType40r8 = SelType40r8f()
 
#Clean ONLY CD4 from Double Negtive --------------------------------------
filedatalCD3only4 = filedatalCD3[!selType40r8,]
bType8 = all( selType40r8 )
bType4 = all( !selType40r8 )
bOnlyOneType48 = bType8 | bType4

sel_doubleNegtoveBoundaryIndForAll = sel_doubleNegtovef()

#CD4 without Double Negtive
selcleanCD4 = (!sel_doubleNegtoveBoundaryIndForAll  & (!selType40r8))

#Find the Double positive ------------------------------------
sel_DP = FindDoublepositive()

selcleanCD4_NoDP =  ( selcleanCD4 ) & ( !sel_DP )
selcleanCD8_NoDP  = selType40r8 & (!sel_DP )

selcleanCD4_NoDPnum = sum(selcleanCD4_NoDP)
selcleanCD8_NoDPnum = sum(selcleanCD8_NoDP)


#QC4!
if(bAddmedianInformation){
  
  NumCD4 = sum(selcleanCD4_NoDP)
  NumCD8 = sum(selcleanCD8_NoDP)
  
  
  addata = data.frame( testNum = Cartnum,Type = "CD4" , size = dim(filedatalCD3[selcleanCD4_NoDP,])[1], t( apply(filedatalCD3[selcleanCD4_NoDP,c(paste0("Area",1:8),"FCS")],2,median)) )
  medianDat = rbind( medianDat , addata )
  addata = data.frame( testNum = Cartnum,Type = "CD8" , size = dim(filedatalCD3[selcleanCD8_NoDP,])[1] ,  t( apply(filedatalCD3[selcleanCD8_NoDP,c(paste0("Area",1:8),"FCS")],2,median)) )
  medianDat = rbind( medianDat , addata )
  
  QCResults = QcCheck("QC4",QCResults) 
  
  if(QCResults$ErrorNum != 0 ){
    if ( USE_ACCEPTANCE_CRITERIA){
      
      return(errorH(UseValidationMode,QCResults$ErrorString , plotpath,QCResults ) )  
    }
  }
  
  
  #QC8!
  QCResults = QcCheck("QC8",QCResults) #c("QCNumberOf45","QC45","QC3","QCDraq7","QC4","QC8")
  if(QCResults$ErrorNum != 0 ){
    if ( USE_ACCEPTANCE_CRITERIA ){
      return(errorH( UseValidationMode, QCResults$ErrorString , plotpath,QCResults ) )
    }
  }
  
}

if(bAddmedianInformation){
  
  if(!dir.exists(plotpath)){
    dir.create(plotpath)
  }
  medianDatpaath = paste0(plotpath,"medianDat.csv")
  write.csv(medianDat,medianDatpaath,row.names = F)
}

if(PLOTS){graphics.off()}

#result
if(UseValidationMode == TRUE)
{
  
  Values = c(
    runInformation$DeviceNum,
    Cartnum,
    CD45num,
    CD45Livenum, 
    CD45deadnum,
    CD3num,
    CD3Livenum, 
    CD3deadnum,
    selcleanCD4_NoDPnum,
    selcleanCD8_NoDPnum,
    100*round(Cd45liveRatio,3),
    100*round(Cd3liveRatio,3),
    100*round(Cd3Cd45Ratio,3),
    #100*round(CD4CD8Ratio,2),
    #100*round((selcleanCD4_NoDPnum - SumdoubleNegtoveBoundaryInd)/CD3Livenum,2),
    100*round(selcleanCD4_NoDPnum/CD3Livenum,3),
    100*round(selcleanCD8_NoDPnum/CD3Livenum,3),
    round(selcleanCD4_NoDPnum/selcleanCD8_NoDPnum,2),
    0,#Error
    0)#"Location Warning"
  
  
  Values  = as.character( Values )
  if( length(QCResults$QCarra) > 0 ){
    WarningStringOut = paste0( QCResults$QCarra , collapse = ";")
    Values[length(Values)] = WarningStringOut
  }
  
  
  Variable = c(
    "Device",
    "Cartrige",
    "CD45_total",
    "CD45_live",
    "CD45_dead",
    "CD3_total",
    "CD3_live",
    "CD3_dead",
    "CD4_live",
    "CD8_live",
    "Cd45liveRatio",
    "Cd3liveRatio",
    "Cd3Cd45Ratio",
    #"CD4CD8Ratio",
    "CD4CD3Ratio",
    "CD8CD3Ratio",
    "CD4CD8Ratio",
    "ErrorNum",
    "WarningNums")
  
  df = data.frame(Values,Variable,stringsAsFactors = F)
  df
  
}
else
{  
  
  Values = c(
    #runInformation$DeviceNum,
    Cartnum,
    CD45num,
    CD45Livenum, 
    CD45deadnum,
    CD3num,
    CD3Livenum, 
    CD3deadnum,
    selcleanCD4_NoDPnum,
    selcleanCD8_NoDPnum,
    100*round(Cd45liveRatio,3),
    100*round(Cd3liveRatio,3),
    100*round(Cd3Cd45Ratio,3),
    #100*round(CD4CD8Ratio,2),
    #100*round((selcleanCD4_NoDPnum - SumdoubleNegtoveBoundaryInd)/CD3Livenum,2),
    100*round(selcleanCD4_NoDPnum/CD3Livenum,3),
    100*round(selcleanCD8_NoDPnum/CD3Livenum,3),
    round(selcleanCD4_NoDPnum/selcleanCD8_NoDPnum,2),
    0,#NO Error
    0)#Warning Location
  
  
  Values  = as.character( Values )
  if( length(QCResults$QCarra) > 0 ){
    WarningStringOut = paste0( QCResults$QCarra , collapse = ";")
    Values[length(Values)] = WarningStringOut
  }
  
  
  
  Variable = c(
    #"Device",
    "Cartrige",
    "CD45_total",
    "CD45_live",
    "CD45_dead",
    "CD3_total",
    "CD3_live",
    "CD3_dead",
    "CD4_live",
    "CD8_live",
    "Cd45liveRatio",
    "Cd3liveRatio",
    "Cd3Cd45Ratio",
    #"CD4CD8Ratio",
    "CD4CD3Ratio",
    "CD8CD3Ratio",
    "CD4CD8Ratio",
    "ErrorNum",
    "WarningNums")
  
  
  df = data.frame(Values,Variable,stringsAsFactors = F)
  #write.csv(df,"C:/Project/LeukoDx/TCL/TCL_Results/luda/re288.csv")
  
  #write.csv(df,"c:/Temp/Results283.csv")
  
  paste0(df[,c("Values")],collapse = ",")
}

}

#runEnv ----- 
#main(files[1])
#main(args[1])
