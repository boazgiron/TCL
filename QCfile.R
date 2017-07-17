#Ver4.1
require(signal)

#Constants-------------------------------------------
col  = topo.colors(255)
colbr = c("black","red","green","blue")
ty  = c(paste0("Area",1:8),"Peak9")
bf <- butter(2, 0.2, type="low")
args<-commandArgs(TRUE)
PLOTS = TRUE
UseValidationMode = TRUE
bAddmedianInformation = TRUE
USE_ACCEPTANCE_CRITERIA  = FALSE
#FUNCTIONS-------------------------------------------

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
  vcut = vmax * alfa
  if (type1 == 0 ) {
    
    mincutind = min( which( h$counts[1:maxrv] > vcut ) ) 
    out = 2 * maxrv - mincutind
    
  }else {
    
    wv = which( h$counts[maxrv:length(h$counts)] < vcut )
    if(length(wv) > 0 ){
      
      mincutind = min( wv )
      
    }else{
      
      mincutind = length( maxrv:length(h$counts) )
      
    }
    
    
    out = maxrv - mincutind
    
    
  }
  
  out
  #h_vpArea1FCS$mids[maxrv + mincutind]
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

#source("tools.R")
filter1 <- function(flitercoff,delay,x,y){
  out  = list(x = x[(delay+1):length(x)]  ,y = filter( flitercoff, y )[(delay+1):length(x)])
  out 
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
  
  
  
  #browser()
  
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

#Error Handle function ----------------------------------------------------------
checkType  = c("QCNumberOf45","QC45","QC3","QCDraq7","QC4","QC8")
RatioMaxlimits = c(1e8,1e8,1e8,1e8,1e8,1e8)
RatioMinlimits = c(1.0,1.6,1.8,1.5,1.0,4.0)
sizelimits = c(1000,300,300,300,300,300)
Ratiolimits =  data.frame(checkType,RatioMinlimits,RatioMaxlimits)

sizeLimits = data.frame(Type = c("NumberOf45","NumOfCD3Live","NumOfCD3LiveNegative","NumCD3Dead","NumCD4","NumCD8"),Value = c(1000,300,300,300,300,300) ,stringsAsFactors = F)

#190617 NewCorrection to check

#Main function ------------------------------------------------------------
runAlgo_shortData <- function(wrkingFilepath){

#InFunction  
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

RemoveLowSumCH <- function(filedatal){
  
  
  
  qi =  quantile(filedatal$sumCh,c(0.1,0.9))
  x1 = median(filedatal$sumCh[filedatal$sumCh < qi[1]])
  x2 = median(filedatal$sumCh[filedatal$sumCh > qi[2]])
  y1 = median(filedatal$FCS[filedatal$sumCh < qi[1]])
  y2 = median(filedatal$FCS[filedatal$sumCh > qi[2]])
  x2 = 1/(4*(y2 - y1)/(x2 - x1))
  x2 = min(10,x2)
  
  #x2 = 10  # 3/20
  vp = c(1,x2)/sqrt(1 + x2^2 )
  vv  = c(1,-1/x2)/sqrt(1 +  ( 1/x2 )^2 )
  
  sumCh = cbind( filedatal$sumCh,filedatal$FCS ) %*% (vp)
  
  
  if(PLOTS){
    png(filename= paste0(plotpath,"sumCh",".png"))
    hsumch  =  hist(sumCh,200,xlab = "Sum Ch1-Ch8",main = paste0("Sum Ch1-Ch8 " ,CarName ),plot = T )
  }else{
    hsumch  =  hist(sumCh,200,plot = F )
  }
  
  
  fr = filter1(bf,5,hsumch$mids, hsumch$counts)
  hsumch$mids = fr$x
  hsumch$counts = fr$y
  
  mindiff = min(max(hsumch$counts)/10,100)
  fmax_hsumch = maxdetcetion(hsumch$counts,0.85,2*mindiff,1.2,mindiff,PLOTS,hsumch$mids,main = CarName)
  
  if ( length(fmax_hsumch$maxl$x)  %in% c(2,3) ) {
    
    
    bondaryv = fmax_hsumch$mins$x[2]
    
  }else{
    
    cutv = fmax_hsumch$maxl$y[1]* 0.05
    bondaryv  = hsumch$mids[fmax_hsumch$maxl$indx[1] + min(which(hsumch$counts[fmax_hsumch$maxl$indx[1]:length(hsumch$counts)] < cutv)) - 1]
    
  }
  
  if(PLOTS) {
    
    abline(v = bondaryv , col = 2 ,lty = 2 ,lwd =2)
    dev.off()
  }
  
  sel_sumCh = sumCh > bondaryv
  
  # if(  sum(sel_sumCh) < sizeLimits[sizeLimits$Type == "NumberOf45","Value"] )  { 
  #   #Return Error 1 low number of 45
  #   QCResults$ErrorNum = 1
  #   QCResults$ErrorString = paste0("ErrorNum =",QCResults$ErrorNum ,"  QCNumberOf45 (sel_sumCh ) =",sel_sumCh )
  #   return(errorH(UseValidationMode,QCResults$ErrorString , plotpath,QCResults ) )
  # }
  
  
  if(PLOTS) {
    png(filename= paste0(plotpath,"sumCh_FCS",".png"))
    le = dim(filedatal)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal[samp,],plot(sumCh,FCS , pch = 19, cex= 0.2, 
                               xlab = "sumCh",ylab = "FCS",
                               col = colbr[1+ sel_sumCh[samp]],
                               xlim = c(0,20),ylim = c(1,5),
                               main = paste0("FCS vs Area1 Angle = ", round(x2,2) ,"  ", CarName)))
    #abline(c(0,3/20))
    dev.off()
  }
  
  
  filedatal = filedatal[sel_sumCh, ]
  
  Num_sel_sumCh = sum(sel_sumCh) 
  out  = list( filedatal = filedatal, Num_sel_sumCh  = Num_sel_sumCh )
  out 
}

Select45ByHighArea1 <- function(){
  
  if(PLOTS) {
    
    png(filename= paste0(plotpath,"Area1_FCS",".png"))
    le = dim(filedatal)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal[samp,],plot(Area1,FCS , pch = 19, cex= 0.2, 
                               xlab = "Area1",ylab = "FCS",
                               #xlim = c( 0,quantile(filedatal$Area1,0.999)) , ylim = c( 0.1,max(filedatal$FCS) ) ,
                               main = paste0("FCS vs Area1 ", CarName)))
    dev.off()
  }
  
  #detect 45 High Area1
  if( PLOTS) {
    png(filename= paste0(plotpath,"Area1",".png"))
    h_Area1 = with(filedatal,hist(Area1,200,plot = T))
  }
  else{
    h_Area1 = with(filedatal,hist(Area1,200,plot = F))
  }
  
  fr = filter1( bf,5,h_Area1$mids, h_Area1$counts )
  h_Area1$mids =   fr$x
  h_Area1$counts = fr$y
  
  mindiff = min(max(h_Area1$counts)/10,100)
  h_Area1_max = maxdetcetion(h_Area1$counts,0.85,2*mindiff,1.2,mindiff,PLOTS,h_Area1$mids,main = CarName)
  
  
  le = length(h_Area1_max$maxl$x)
  
  if( le == 2 ) {
    
    ind = h_Area1_max$maxl$indx[1]:h_Area1_max$maxl$indx[2]
    indmin = which(min(h_Area1$counts[ind]) == h_Area1$counts[ind]) + h_Area1_max$maxl$indx[1]
    
    if(length(indmin) > 1){
      
      indmin  = indmin[1]
    }
    
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
    
    diffbu = average2(diffbu)
    if(PLOTS) {
      points(se[3:length(se)],diffbu,'l',pch = 19,col = 2)
    }
    
    maxlocation = which(max(diffbu) == diffbu)
    
    if(length(maxlocation) > 1){
      
      maxlocation  = maxlocation[1]
    }
    
    if(PLOTS) {
      points(se[maxlocation+1],diffbu[maxlocation],col =3, lwd = 2,cex =3 )
    }
    
    cutLowArea1 =  se[ maxlocation + 1 ]
    
  }else if ( le == 1 ){
    
    indmin = which.max(cumsum(h_Area1$counts)/sum(h_Area1$counts) > 0.02)
    
    #190617 NewCorrection to check
    if(length(indmin) > 1){
      
      indmin  = indmin[1]
    }
    cutLowArea1 = h_Area1$mids[ indmin ]
    
  }
  
  selArea1 = filedatal$Area1 > cutLowArea1
  
  if(PLOTS) {
    
    #h_Area1 = with(filedatal,hist(Area1,100),plot = T)
    abline(v = cutLowArea1,col = 2 ,lty = 2, lwd =2)
    dev.off()
  }
  
  #filedatal45pl = filedatal[selArea1,]
  out = filedatal[selArea1,]
  
  if(PLOTS) {
    png(filename= paste0(plotpath,"Area1",".png"))
    le = dim(filedatal)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal[samp,],plot(Area1,FCS,pch =19 ,cex = 0.2 ,col = colbr[1 + selArea1[samp]] ))
    dev.off()
  }
  
  
  out = list(filedatal45pl = filedatal[selArea1,],selArea1 = selArea1 )
  
  out
  
}

FindPreDead <-function(){    
  lm68 = with(filedatal45pl[filedatal45pl[,"Area6"] > 1,],lm(Area8 ~ Area6 ))
  x2 = lm68$coefficients[2]
  vp = c(1,x2)/sqrt(1 + x2^2 )
  vv  = c(1,-1/x2)/sqrt(1 +  ( 1/x2 )^2 )
  
  # Create parallel for cutting low dots
  vpArea6Area8 =  cbind( filedatal45pl$Area6,filedatal45pl$Area8 ) %*% (vp)
  
  
  if(PLOTS){
    png(filename= paste0(plotpath,"vpArea6Area8",".png"))
    h_vpArea6Area8 = hist(vpArea6Area8,100,plot = T)
  }else{
    h_vpArea6Area8 = hist(vpArea6Area8,100,plot = F)
  }
  
  fr = filter1(bf,5,h_vpArea6Area8$mids, h_vpArea6Area8$counts)
  h_vpArea6Area8$mids = fr$x
  h_vpArea6Area8$counts = fr$y
  #h_vpArea6Area8_max  = findmax1(h_vpArea6Area8$mids, h_vpArea6Area8$counts,20,5,5,20,main = CarName)
  mindiff = min(max(h_vpArea6Area8$counts)/10,100)
  h_vpArea6Area8_max = maxdetcetion(h_vpArea6Area8$counts,0.85,2*mindiff,1.2,mindiff,PLOTS,h_vpArea6Area8$mids,main = CarName)
  
  
  le = length(h_vpArea6Area8_max$maxl$indx)
  if(le >= 2){
    
    h_lm68hpmiddsleMin = h_vpArea6Area8_max$mins$x[2]
    
  }else if ( le == 1 ) {
    
    hmaxs  = findbycut(h_vpArea6Area8_max$maxl$indx[1],h_vpArea6Area8,2e-2,1)
    h_lm68hpmiddsleMin =  h_vpArea6Area8$mids[hmaxs]
  }else{
    
    print(paste0( "le = length(h_vpArea6Area8_max$maxl$indx )le =",le ))
    #Problem
  }
  
  if(PLOTS){
    abline(v = h_lm68hpmiddsleMin,col  = 3,lwd = 2,lty =2)
    dev.off()
  }
  
  #Remove the spread dots
  selpredead  = vpArea6Area8 > h_lm68hpmiddsleMin
  
  
  if(PLOTS){
    png(filename= paste0(plotpath,"SelPreDead",".png"))#plot11
    le = dim(filedatal45pl)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatal45pl[samp,],plot(Area6,Area8 , pch = 19, cex= 0.2, 
                                   xlab = "Area6",ylab = "Area8",
                                   #xlim = c( 0,quantile(filedatal45pl$Area6,0.999)) , ylim = c( 0.1,max(filedatal45pl$Area8) ) ,
                                   main = paste0("Area8 vs Area6 ", CarName),
                                   col = colbr[1+selpredead[samp,]]))
    
    abline(lm68,col = "blue",lwd = 2)
    
  }
  
  #Cut tail select up Area6 and Are8 
  filedatal45pl_predead = filedatal45pl[selpredead,]
  
  out = list(filedatal45pl_predead = filedatal45pl[selpredead,],selpredead = selpredead) 
  
  out
}

FindDead <- function(){
  
  #FIRST iteration
  le = dim(filedatal45pl_predead)[1]
  samp  = sample(1:le,min(le,4e4))
  Area6 = c(rep(-2,1e4),filedatal45pl_predead[samp,"Area6"])
  Area8 = c(rep(-2,1e4),filedatal45pl_predead[samp,"Area8"])
  lm68s = lm(Area8 ~ Area6 )
  
  if(PLOTS){
    abline(lm68s,col = 2,lwd = 2,lty = 2)
    legend("topleft",c("first","second"),col = c("blue","red"),pch = 19 ) #SecondVer
    dev.off() #SecondVer
  }
  
  x2 = lm68s$coefficients[2]
  vp = c(1,x2)/sqrt(1 + x2^2)
  vv  = c(1,-1/x2)/sqrt(1 + (1/x2)^2)
  vvs = cbind( filedatal45pl_predead$Area6,filedatal45pl_predead$Area8) %*% vv
  
  
  if(PLOTS){
    
    png(filename= paste0(plotpath,"h_lm68hv",".png"))
    h_lm68hv = hist(vvs,200,plot = T )
    
  }else{
    
    h_lm68hv = hist(vvs,200,plot = F )
  }  
  
  fr =  filter1(bf,1,h_lm68hv$mids,h_lm68hv$counts)
  h_lm68hv$mids = fr$x
  h_lm68hv$counts = fr$y
  
  mindiff = min(max(h_lm68hv$counts)/10,100)
  h_lm68hv_max = maxdetcetion(h_lm68hv$counts,0.85,2*mindiff,1.2,mindiff,PLOTS,h_lm68hv$mids,main = CarName)
  h_lm68hv_max_right = h_lm68hv_max$maxl$indx[length(h_lm68hv_max$maxl$indx)]
  cut0.02ind = findbycut(h_lm68hv_max_right,h_lm68hv,2e-2,1)
  
  seldeadv = h_lm68hv$mids[cut0.02ind]
  
  if(PLOTS){
    
    abline(v = seldeadv ,col = 2,lwd = 2,lty = 2)
    dev.off()
  }
  
  seldead = vvs < seldeadv
  vvsAll = cbind( filedatal45pl$Area6,filedatal45pl$Area8) %*% vv
  
  seldeadAll = ( vvsAll < seldeadv ) & selpredead
  
  out = list( seldeadAll = seldeadAll, seldead = seldead )
  out
}

FindCD3 <-function(){
  
  if(PLOTS){
    
    png(filename= paste0(plotpath,"FD45pl_Area6",".png"))
    h_hist6 = hist(filedatal45pl$Area6,200,plot = T)
    
  }
  else{
    
    h_hist6 = hist(filedatal45pl$Area6,200,plot = F)
  }
  
  fr =  filter1(bf,1,h_hist6$mids,h_hist6$counts)
  h_hist6$mids = fr$x
  h_hist6$counts = fr$y
  
  mindiff = min(max(h_hist6$counts)/10,100)
  h_lm68hp_max = maxdetcetion(h_hist6$counts,0.85,2*mindiff,1.2,mindiff,PLOTS,h_hist6$mids,main = paste0(CarName," Find CD3 Hist of the CH6 Take the high heap"),xlab = "Area6 Index")
  h_hist6_max_right = h_lm68hp_max$maxl$indx[length(h_lm68hp_max$maxl$indx)]
  
  
  if( length(h_lm68hp_max$maxl$indx) == 1 ) {
    cut0.02ind = findbycut(h_hist6_max_right,h_hist6,1.5e-2,1)
    CD3selectionBoundary = h_hist6$mid[cut0.02ind]
  }
  else{# new Code Version 3
    
    minmiddle = h_lm68hp_max$mins$y[2]
    cut0.02ind  = h_lm68hp_max$mins$indx[2]
    CD3selectionBoundary = h_lm68hp_max$mins$x[2]
  }
  
  if(PLOTS){
    abline(v = CD3selectionBoundary, col = 3, lwd = 2,lty = 2 )
    dev.off()
    
    png(filename= paste0(plotpath,"Area6_Area1",".png"))#plot16
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
  selCD3countes  = sum(selCD3)
  
  selCD3
}

selCD48 <- function(){
  
  #Create parallel for to cut tail for better differ cd4 and cd8
  
  da1 = rbind(filedatal45pl[selCD3,c("Area4","Area3")],data.frame(Area4 =rep(1,1e4),Area3 = rep(1,1e4)))
  lm34 = with(da1, lm(Area3 ~ Area4))
  
  if(PLOTS){  
    
    abline(lm34,col = 3,lwd = 2)
    dev.off()
    
  }
  
  x2 = lm34$coefficients[2]
  vp = c(1,x2)/sqrt(1 + x2^2)
  vv  = c(1,-1/x2)/sqrt(1 + (1/x2)^2)
  vps34 = with(filedatalCD3,cbind(Area4,Area3) %*% vp)
  
  if(PLOTS){
    png(filename= paste0(plotpath,"h_vps34",".png"))#plot20
    h_vps34  = hist(vps34,200,plot = T)
  }
  else{
    
    h_vps34  = hist(vps34,200,plot = F)
  }
  
  fr =  filter1(bf,1,h_vps34$mids,h_vps34$counts)
  h_vps34$mids = fr$x
  h_vps34$counts = fr$y
  
  #Ver2.1
  #h_lm68hv_max  = findmax1(h_vps34$mids, h_vps34$counts,21,5,5,20,main = paste0(CarName," CH3 vs CH4 Parallel Direction Histogram ( Cut low tail )"))
  mindiff = min(max(h_vps34$counts)/10,100)
  h_lm68hv_max = maxdetcetion(h_vps34$counts,0.85,2*mindiff,1.2,mindiff,PLOTS,h_vps34$mids,main = paste0(CarName," CH3 vs CH4 Parallel Direction Histogram ( Cut low tail )"))
  leMax = length(h_lm68hv_max$maxl$indx)
  
  if ( leMax >= 3 ) 
  {
    
    cutpointIndex = h_lm68hv_max$mins$indx[2]
    cutpoint = h_lm68hv_max$mins$x[2]
    text( cutpoint - 0.1,h_lm68hv_max$maxl$y[2],"Tail")
    
    #no tail 
  }
  else 
  {
    
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
  
  
  selNonTailCD48 =  vps34  >  cutpoint
  filedatalCD3_sel_tail = filedatalCD3[selNonTailCD48,]
  
  
  x2 = lm34$coefficients[2]
  vp = c(1,x2)/sqrt(1 + x2^2)
  vv  = c(1,-1/x2)/sqrt(1 + (1/x2)^2)
  vvs34woTail = with(filedatalCD3_sel_tail,cbind(Area4,Area3) %*% vv)
  if(PLOTS){
    png(filename= paste0(plotpath,"vvs34woTail",".png"))
    h_vvs34woTail  = hist(vvs34woTail,100)
  }else{
    h_vvs34woTail  = hist(vvs34woTail,100,plot = F)
    
  } 
  
  fr =  filter1(bf,1,h_vvs34woTail$mids,h_vvs34woTail$counts)
  h_vvs34woTail$mids = fr$x
  h_vvs34woTail$counts = fr$y
  
  mindiff = min(max(h_vvs34woTail$counts)/10,100)
  h_vvs34woTailvvs_max = maxdetcetion(h_vvs34woTail$counts,0.85,2*mindiff,1.2,mindiff,PLOTS,h_vvs34woTail$mids,main = paste0(CarName," CH4 vs CH3 vertical histogram"))
  
  minbtween = 0 #Version4 Continue 
  #find the minmum between max ( asume 2 )
  if( length(h_vvs34woTailvvs_max$maxl$indx) >= 2){ 
    minbtween = h_vvs34woTailvvs_max$mins$x[2]
  }
  
  if(PLOTS){
    abline(v = minbtween, col = 3, lwd = 2 ,lty = 2)
    dev.off()
  }
  
  #Create CD4 CD8 selection 
  vvs34 = with(filedatalCD3,cbind(Area4,Area3) %*% vv)
  #create between 4 and 8 selection variable
  selType40r8 = vvs34 > minbtween
  
  selType40r8
}

CleanCD4DoubleNegtive <- function(){
  
  le = dim(filedatalCD3only4)[1]
  samp  = sample(1:le,min(le,1e4))
  da1 = rbind(filedatalCD3only4[,c("Area1","Area3")],data.frame(Area1 =rep(1,1e4),Area3 = rep(1,1e4)))
  lmonly4 = with(da1, lm(Area3 ~ Area1))
  
  x2 = lmonly4$coefficients[2]
  vp = c(1,x2)/sqrt(1 + x2^2)
  vv  = c(1,-1/x2)/sqrt(1 + (1/x2)^2)
  
  si  = dim(filedatalCD3only4)[1]
  
  # if(si < 1e3 ){
  #   return(errorH(UseValidationMode,paste0("filedatalCD3only4 =", si),plotpath))
  # }
  
  vvslmonly4 = with(filedatalCD3only4,cbind(Area1,Area3) %*% vv)
  if(PLOTS){
    png(filename= paste0(plotpath,"vvslmonly4",".png"))#plot23
    h_vvslmonly4  = hist(vvslmonly4,200,plot = T)
  }
  else{
    h_vvslmonly4  = hist(vvslmonly4,200,plot = F)
  }
  
  fr =  filter1(bf,1,h_vvslmonly4$mids,h_vvslmonly4$counts)
  h_vvslmonly4$mids = fr$x
  h_vvslmonly4$counts = fr$y
  
  mindiff = min(max(h_vvslmonly4$counts)/10,100)
  h_vvslmonly4_max = maxdetcetion(h_vvslmonly4$counts,0.85,2*mindiff,1.2,mindiff,PLOTS,h_vvslmonly4$mids,main = paste0(CarName," Sperate direction CH3 Vs CH1 histogram \n Clean CD4 from Double Negtive" ))
  
  if(length(h_vvslmonly4_max$maxl$indx) > 2 ) 
  {
    cutIndex = length(h_vvslmonly4_max$mins$indx) - 1
    doubleNegtoveBoundary = h_vvslmonly4_max$mins$x[cutIndex]
    text( h_vvslmonly4_max$maxl$x[cutIndex]- 0.1,h_vvslmonly4_max$maxl$y[cutIndex],"DN")
    
  }
  else if(length(h_vvslmonly4_max$maxl$indx) == 2 )
  {
    cutIndex = 2
    doubleNegtoveBoundary = h_vvslmonly4_max$mins$x[cutIndex]
    if(PLOTS == TRUE){
      if(length(h_vvslmonly4_max$maxl$indx) == 2 ){
        text( h_vvslmonly4_max$maxl$x[cutIndex] - 0.05,h_vvslmonly4_max$maxl$y[cutIndex],"DN")  
      }
    }
  }
  else
  {
    
    doubleNegtoveBoundary  = h_vvslmonly4_max$mins$x[2]
    
  }
  
  if(PLOTS){
    abline(v = doubleNegtoveBoundary, col = 2 , lwd =2 , lty = 2 )
    dev.off()
  }
  
  sel_doubleNegtoveBoundary  = vvslmonly4 > doubleNegtoveBoundary
  
  
  if(PLOTS){
    png(filename= paste0(plotpath,"CD3only4_Area1_Area3",".png"))#plot24
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
  sel_doubleNegtoveBoundaryIndForAll  = vvslmonly4forall > doubleNegtoveBoundary
  
  sel_doubleNegtoveBoundaryIndForAll
  
}

FindDoublepositive <-function(){
  
  selpreDP  = !sel_doubleNegtoveBoundaryIndForAll
  
  if(PLOTS){
    png(filename= paste0(plotpath,"FinalArea4Area3",".png"))
    le = dim(filedatalCD3)[1]
    samp  = sample(1:le,min(le,1e4))
    with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = 0.2,col = colbr[1+(selpreDP + 2*selType40r8 )[samp]],main = paste0(CarName," Pre Double Positive ") ))
    legend("topleft",c("CD4DN","CD4","LoWCD8","CD8"),col = c(1,2,3,4),pch=19 )
    dev.off()
  }
  
  preDP =  filedatalCD3[selpreDP,c("Area4","Area3")]
  preDPselmin = selType40r8[selpreDP]
  
  g1 =  preDP[preDPselmin,]
  g2 =  preDP[!preDPselmin,]
  
  g1m = colMeans(g1,na.rm = T)
  g2m = colMeans(g2,na.rm = T)
  
  us = cbind(g1m,g2m)
  
  fl = flda(us[,1],us[,2],g1,g2)
  fldaP = fldaProjection(preDP,us[,1],us[,2],fl)
  
  if(PLOTS){
    png(filename= paste0(plotpath,"fldaPForDoublePositive",".png"))
    h_fldaP  = hist(fldaP,200,main = paste0(CarName," LDA for finding the Double Positive points") )
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
  
  sel_DP  = !((h_fldaP$mids[downind] < fldapreDPCD3 ) & (h_fldaP$mids[upind] > fldapreDPCD3 ))
  sel_DP = selpreDP & !sel_DP
  
  
  if(PLOTS){
    png(filename= paste0(plotpath,"DoublePositive_Area4_Area3",".png"))#plot28
    le = dim(filedatalCD3)[1]
    samp  = sample(1:le,min(le,2e4))
    with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = c(0.2,0.6)[1  + sel_DP[samp]],col = colbr[1  + sel_DP[samp]  ],main = paste0(CarName," Select Double Positive") ))
    dev.off()
  }
  
  sel_DP
}

#ini Param-----------------------------------------------------------------
  
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

#Constant--------------------------------------------------------------------
bAddmedianInformation = TRUE
medianDat = NULL
QCResults = list(ErrorNum = 0,QCarra = NULL, ErrorString = NULL )
PLOTS = TRUE
UseValidationMode = TRUE
bAddmedianInformation = TRUE
#Load file ------------------------------------------------------------------
  
  #wrkingFilepath = ValidationFiles_[grep("257",ValidationFiles_)]

  filedata = read.csv(wrkingFilepath,header = T)
  di = dim(filedata)
  
  #Check size of file
  #sizeLimits[sizeLimits$Type == "NumCD4","Value"]
  
  if( ( di[1] < sizeLimits[sizeLimits$Type == "NumberOf45","Value"] ) | ( di[2] < 20 ) ) 
  { 
    #Return Error 1 low number of 45
    QCResults$ErrorNum = 1
    QCResults$ErrorString = paste0("ErrorNum = ",QCResults$ErrorNum ,"  QCNumberOf45(preAnaysis) = " ,  di[1] )
    return(errorH(UseValidationMode, QCResults$ErrorString , plotpath,QCResults ) )
  }
    
  runInformation  = getCartrigeNameAdnDevice(wrkingFilepath)
  
  CarName = paste0("C",runInformation$Cartnum)
  
  Cartnum = runInformation$Cartnum
  
  if(PLOTS){
    plotpath = paste0("C:/Project/LeukoDx/LudaFacsValidation/ValidationVer3/",Cartnum,"/")
    if(!dir.exists(plotpath)){
      dir.create(plotpath)
    }
  }
  
  #Prepare result file ONLY IN DEBUG MODE
  if(bAddmedianInformation){
    medianDatpaath = paste0("C:/Project/LeukoDx/LudaFacsValidation/ValidationVer3/",Cartnum,"/medianDat.csv")
  }
  
#Remove low width ------------------------------------------------------
  widthb = filedata$Width > 4
  
  if(PLOTS) {
    png(filename= paste0(plotpath,"Width",".png"))
    hist(filedata$Width,500,xlim = c(0,30),xlab = "Width", main = as.character(CarName))
    abline(v = 4,col ="red",lwd = 2 )
    dev.off()
  }
  
  filedata = filedata[widthb,]
  #Create log10 Variables-------------------------------------------------
  filedatal =  data.frame(apply(filedata[,ty[1:9]],2,log10))
  filedatal$FCS <- log10(filedata$Peak9) 
  filedatal$sumCh =  apply(filedatal[,1:8],1,sum)

  
#Remove low sumCH events -------------------------------------------------
  results = RemoveLowSumCH(filedatal)
  
  filedatal = results$filedatal
  Num_sel_sumCh  = results$Num_sel_sumCh
  
  if(  results$Num_sel_sumCh < sizeLimits[sizeLimits$Type == "NumberOf45","Value"] )  { 
    #Return Error 1 low number of 45
    QCResults$ErrorNum = 1
    QCResults$ErrorString = paste0("ErrorNum =",QCResults$ErrorNum ,"  QCNumberOf45 ( Num_sel_sumCh ) =",Num_sel_sumCh )
    return(errorH(UseValidationMode,QCResults$ErrorString , plotpath,QCResults ) )
  }
  
#Select 45 by Area1-------------------------------------------------------
  results = Select45ByHighArea1()
  filedatal45pl = results$filedatal45pl
  selArea1  = results$selArea1

#QCNumberOf45! 
if(bAddmedianInformation){
  
  QCResults = QcCheck("QCNumberOf45",QCResults) 
  
  if(QCResults$ErrorNum != 0 ){
    
    if ( USE_ACCEPTANCE_CRITERIA){
      return  ( errorH(UseValidationMode,QCResults$ErrorString, plotpath,QCResults ) )
    }
    
  }
  
  addata = data.frame( testNum = Cartnum,Type = "CD45",size = dim(filedatal45pl)[1] , t(apply(filedatal45pl[,c(paste0("Area",1:8),"FCS")],2,median)) )
  medianDat = rbind( medianDat , addata )
  
  addata = data.frame( testNum = Cartnum,Type = "Junk",size = dim(filedatal45pl[!selArea1,])[1] , t(apply(filedatal[!selArea1,c(paste0("Area",1:8),"FCS")],2,median)))
  medianDat = rbind( medianDat , addata )
}

#Find Dead-----------------------------------------------------------------

  result = FindPreDead() 
  
  filedatal45pl_predead   = result$filedatal45pl_predead

  selpredead = result$selpredead

  results = FindDead()
  
  seldeadAll = results$seldeadAll
  seldead = results$seldead
  
  if(bAddmedianInformation){
  
    addata = data.frame( testNum = Cartnum,Type = "CD45Live" ,size = dim(filedatal45pl[!seldeadAll,])[1],  t(apply(filedatal45pl[!seldeadAll,c(paste0("Area",1:8),"FCS")],2,median)) )
    medianDat = rbind( medianDat , addata )
    
    
    addata = data.frame( testNum = Cartnum,Type = "CD45Dead" , size = dim(filedatal45pl[seldeadAll,])[1], t(apply(filedatal45pl[seldeadAll,c(paste0("Area",1:8),"FCS")],2,median)) )
    medianDat = rbind( medianDat , addata )
  }

  DeadNum  = sum(seldeadAll)

  if(PLOTS){
    
    png(filename= paste0(plotpath,"Dead",".png"))#plot14
    le = dim(filedatal45pl_predead)[1]
    samp  = sample(1:le,min(le,2e4))
    plot(filedatal45pl_predead$Area6[samp],filedatal45pl_predead$Area8[samp] , pch = 19, cex= 0.2, 
         xlab = "Area6",ylab = "Area8",
         xlim = c( 0,quantile(filedatal45pl_predead$Area6[samp],0.999,na.rm = T)) , ylim = c( 0.1,max(filedatal45pl_predead$Area8[samp],na.rm = T) ) ,
         main = paste0("Area8 vs Area6 ", CarName),
         col = colbr[1+seldead[samp]])
    legend("topleft",c("live","Dead"),col = c(1,2),pch=19 )
    #abline(lm68live,col  = 5, lwd = 2 )
    dev.off()
  
  }

#Find CD3-----------------------------------------------------------------
  
selCD3 = FindCD3()
  
#QC45!
  if(bAddmedianInformation){
  
    #selCD3live = filedatal45pl[!seldeadAll,]$Area6 > CD3selectionBoundary
    NumOfCD3Live = sum((!seldeadAll) & selCD3)
    
    addata = data.frame( testNum = Cartnum,Type = "CD3" , size = dim(filedatal45pl[(!seldeadAll) & selCD3,])[1], t(apply(filedatal45pl[(!seldeadAll) & selCD3,c(paste0("Area",1:8),"FCS")],2,median)) )
    medianDat = rbind( medianDat , addata )
    addata = data.frame( testNum = Cartnum,Type = "NoCD3" , size = dim(filedatal45pl[(!seldeadAll) & (!selCD3),])[1] , t(apply(filedatal45pl[(!seldeadAll) & (!selCD3),c(paste0("Area",1:8),"FCS")],2,median)) )
    medianDat = rbind( medianDat , addata )
    
    QCResults = QcCheck("QC45",QCResults)
    
    if ( QCResults$ErrorNum != 0 ) {
      
      if ( USE_ACCEPTANCE_CRITERIA){
        return( errorH( UseValidationMode , QCResults$ErrorString , plotpath , QCResults ) )
      }
    }
    #QC3!  
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
    
    addata = data.frame( testNum = Cartnum,Type = "CD3Live" , size = dim(filedatal45pl[selCD3Live,])[1], t(apply(filedatal45pl[selCD3Live,c(paste0("Area",1:8),"FCS")],2,median)) )
    medianDat = rbind( medianDat , addata )
    addata = data.frame( testNum = Cartnum,Type = "CD3Dead" , size = dim(filedatal45pl[selCD3Dead,])[1], t(apply(filedatal45pl[selCD3Dead,c(paste0("Area",1:8),"FCS")],2,median)) )
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
  CD45num = length(seldeadAll)
  CD45deadnum = sum(seldeadAll)
  CD45Livenum  = CD45num   - CD45deadnum 
  Cd3liveRatio =  CD3Livenum/CD3num
  Cd45liveRatio =  CD45Livenum/CD45num
  Cd3Cd45Ratio = CD3Livenum/CD45Livenum


#Area3 vs Area4 select the high Area6
filedatalCD3 = filedatal45pl[selCD3Live,]

if(PLOTS){
  png(filename= paste0(plotpath,"FdCD3",".png"))
  le = dim(filedatalCD3)[1]
  samp  = sample(1:le,min(le,2e4))
  with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = 0.2,main = paste0(CarName," CD4 CD8 Create speration direction" )))
}

#Find CD48--------------------------------------------------------------

selType40r8  = selCD48()

le = dim(filedatalCD3)[1]
samp  = sample(1:le,min(le,1e4))

if(PLOTS){
  png(filename= paste0(plotpath,"FdCD3_Area4Area3",".png"))
  with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = 0.2,col = colbr[1+!selType40r8[samp]],main = paste0(CarName," First Cut for CD4 and CD8")))
  dev.off()  
}

#Clean ONLY CD4 from Double Negtive------------------------------------
filedatalCD3only4 = filedatalCD3[!selType40r8,]
sel_doubleNegtoveBoundaryIndForAll  = CleanCD4DoubleNegtive()

#CD4 without Double Negtive
selcleanCD4 = (!sel_doubleNegtoveBoundaryIndForAll  & (!selType40r8))

#Find the Double positive---------------------------------------------
sel_DP = FindDoublepositive()

#Summerize------------------------------------------------------------
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

# rb = as.data.frame(rbind(Area1DArea2,DeadCD3Area7dArea8,Area4DArea8,CD4Area2dArea4,CD8Area4dArea5))
# rb$VariableName <-  rownames(rb)
# rownames(rb) <- NULL
# colnames(rb)  = c("Value","NotInRange","VariableName")
# rb = rb[, c("VariableName","Value","NotInRange")]
# if(PLOTS){
#   write.csv(rb,paste0(plotpath,"Results.txt"),row.names = F)
# }

if(bAddmedianInformation){
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

#Run script---
#runAlgo_shortData(args[1])
#----
