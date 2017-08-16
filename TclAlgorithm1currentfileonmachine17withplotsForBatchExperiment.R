#Ver7
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
col  = topo.colors(255)
colbr = c("black","red","green","blue")
ty  = c(paste0("Area",1:8),"Peak9")
bf <- butter(2, 0.2, type="low")
bf4 <- butter(4, 0.2, type="low")
args<-commandArgs(TRUE)
PLOTS = TRUE
PLOTS2 = TRUE
UseValidationMode = TRUE
WRITE_RESUTLS = FALSE
LUDA = FALSE


#FUNCTIONS------------------

DEV.OFF <- function(){
  
  if(TRUE){
    
    result = tryCatch({
      dev.off()
      
    }, warning = function(w) {
      
    }, error = function(e) {
      FALSE
    }, finally = {
      
    })
    
    
  }
  
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
  #ChngefromOriginal
  #filename = strsplit(filePath,"/")[[1]]
  filename = basename(filePath)
  #filepart = strsplit(filePath,"_")[[1]]
  filepart = strsplit(filename,"_")[[1]]
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

errorH <- function(UseValidationMode1,ErrorString,plotpath1){
  
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
      "CD4CD8Ratio")
    
    Values  = rep(0,length(Variable))
    
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
      "CD4CD8Ratio") 
    
    Values  = rep(0,length(Variable))
    
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
      
      DEV.OFF()
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

#main function
runAlgo_shortData <- function(wrkingFilepath , datadir,UseCheck ){
  
  #ini Param-----------------------
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
  
  #wrkingFilepath = fll#filepath1
  #Load file-----------------------
  
  #wrkingFilepath = fll
  #datadir = dataDirname
  #UseCheck = TRUE
  
  filedata = read.csv(wrkingFilepath,header = T)
  
  di = dim(filedata)
  
  if( ( di[1] > 1e3 ) & ( di[2] > 28 ) ) {
    
    runInformation  = getCartrigeNameAdnDevice(wrkingFilepath)
    
    CarName = paste0("C",runInformation$Cartnum)
    
    Cartnum = runInformation$Cartnum
    
    plotpath = paste0(datadir,"/",Cartnum,"/")
    
    if(PLOTS){  
      if(!dir.exists(plotpath)){
        dir.create(plotpath)
      }
    }
    
    #Remove low width
    widthb = filedata$Width > 4
    
    if(PLOTS) {
      png(filename= paste0(plotpath,"Width",".png"))
      hist(filedata$Width,500,xlim = c(0,30),xlab = "Width", main = as.character(CarName))
      abline(v = 4,col ="red",lwd = 2 )
      DEV.OFF()
    }
    
    filedata = filedata[widthb,]
    
    filedatals =  data.frame(apply(filedata[,ty[1:9]],2,log10))
    filedatals$FCS <- filedata$Peak9
  
    #DD
     
    
    for( i in 1:dim(filedatals)[2]){
      filedatals[is.infinite(filedatals[,i]),i] = NA
    }
    
    sel4 = complete.cases(filedatals)
    filedatals = filedatals[sel4,]
    

    filedatals$sumCh =  apply(filedatals[,1:8],1,sum)
    
    #Sumch Detection--------------------------------------------
    
    SumCH_MinMAX = 10
    sumCh = filedatals$sumCh 
    hsumch  = hist(sumCh,200,plot = F )
    fr = filter1(bf,5,hsumch$mids, hsumch$counts)
    hsumch$mids = fr$x
    hsumch$counts = fr$y
    fmax_hsumch = maxiNu(hsumch$mids,hsumch$counts,PLOTS,plotpath,"sumCh",FALSE,30,10)
    if(length(fmax_hsumch$maxl) ==  0 ){
      DEV.OFF()
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
    
    #--------------------
    
    if(PLOTS) {
      
      abline(v = bondaryv1 , col = 2 ,lty = 2 ,lwd =2)
      DEV.OFF()
    }
    
    sel_sumCh = sumCh > bondaryv1
    filedatal = filedatals[sel_sumCh,]
    
    if(sum(sel_sumCh ) < 1e3 ){
      
      return(errorH(UseValidationMode,paste0("sel_sumCh =", sum(sel_sumCh )),plotpath))
    }
    
    if(PLOTS) {
      
      png(filename= paste0(plotpath,"SumCh_Area7 boundary bondaryv1",".png"))
      le = dim(filedatals)[1]
      samp  = sample(1:le,min(le,2e4))
      with(filedatals[samp,],plot(sumCh,Area7 , pch = 19, cex= 0.2, 
                                  xlab = "sum FL",ylab = "Area7",
                                  col = colbr[1+ sel_sumCh[samp]],
                                  ylim = c(0,5),
                                  main = paste0("SumCh_Area7 boundary bondaryv1 ", CarName)))
      DEV.OFF()
    }
    
    #Orignal Algorithm 100717DD---------------------
    
    
    qi =  quantile(filedatal$Area1,c(0.05,0.95))
    x1 = median( filedatal$Area1[filedatal$Area1 < qi[1]] )
    x2 = median( filedatal$Area1[filedatal$Area1 > qi[2]] )
    y1 = median( filedatal$FCS[filedatal$Area1 < qi[1]] )
    y2 = median( filedatal$FCS[filedatal$Area1 > qi[2]] )
    x2 = 1/(4*(y2 - y1)/(x2 - x1))
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
      points(c(x1,x2),c(y1,y2),pch  = 19,cex = 2,col = 5) 
      DEV.OFF()
    }
    
    
    m1 = mean(filedatal$Area1[cl$cluster == 1])
    m2 = mean(filedatal$Area1[cl$cluster == 2])
    if(m1 > m2){
      
      sel45 = cl$cluster == 1
      
    }else{
      
      sel45 = cl$cluster == 2
    }
    
    #DD
    fsum <- function(x){sum(x,na.rm = T)}
    
    if( UseCheck){
    
      filedatals$sumCh =  apply(filedatals[,2:8],1,fsum )
      ou = checkfs(plotpath,CartNum)
      sel45 = ou == 2
      filedatal = filedatals
      filedatal_sel45 = filedatal[sel45,]
      
    }else{
      
    
      filedatal_sel45 = filedatal[sel45,]
    }
    #DDf
  
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
      DEV.OFF()
    }

    
    
    sel6 = filedatal$Area6  > bounadry_h_Area1_fcs
    
    if(PLOTS) {
      
      png(filename= paste0(plotpath,"cd3 Area6_Area1 ",".png"))
      sel6_ = sel6[sel45]  
      le = dim(filedatal_sel45)[1]
      samp  = sample(1:le,min(le,2e4))
      with(filedatal_sel45[samp,],plot(Area6,Area1 , pch = 19, cex= 0.2, 
                                       xlab = "Area6",ylab = "Area1",
                                       col = colbr[1 + sel6_[samp]] ,
                                       #ylim = c(0,5) ,
                                       main = paste0("CD3 Area6_Area1 ", CarName)) )
      DEV.OFF()
    } 
    
    filedatal45pl = filedatal[sel45,]
    selCD3 =  sel6 & sel45
    
    #Old CODE 280617---------------------
    
    dimfiledatal45pl  = dim(filedatal45pl)[1]
    
    if(dimfiledatal45pl < 1e3 ){
      
      return(errorH(UseValidationMode,paste0("CD45 countes =", dimfiledatal45pl),plotpath))
    }    
    
    
    
    
    Area1DArea2b = 0
    Area1DArea2  = round(median(filedatal[sel45,]$Area1)/median(filedatal[sel45,]$Area2),2)
    
    #DD
    Area1DArea2b = tryCatch( { ifelse(Area1DArea2 < 1.22,1,0)
    }, warning = function(w) {
      0#warning-handler-code
    }, error = function(e) {
      0#error-handler-code
    }, finally = {
      0#cleanup-code
    })
    
    # if(Area1DArea2 < 1.22 ) {
    #   Area1DArea2b = 1
    #   #return(errorH(UseValidationMode,paste0("CD45Area1/Area2 =", Area1DArea2),plotpath))
    # }
    # 
    
    
    Area1DArea2 = c(Area1DArea2,Area1DArea2b)
    
    
    #Find Dead---------------
    
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
      DEV.OFF()
    }
    
    vpArea6Area8 =  cbind( filedatal$Area6,filedatal$Area8 ) %*% (vp)
    selpredead  = vpArea6Area8 > h_lm68hpmiddsleMin
    
    if(PLOTS){
      
      png(filename= paste0(plotpath,"SelPreDead",".png"))#plot11
      le = dim(filedatal45pl)[1]
      samp  = sample(1:le,min(le,2e4))
      
      with(filedatal[samp,],plot(Area6,Area8 , pch = 19, cex= 0.2, 
                                 xlab = "Area6",ylab = "Area8",
                                 main = paste0("SelPreDead Cut tail ", CarName),
                                 col = colbr[1+selpredead[samp,]]))
      
      abline(lm68,col = "blue",lwd = 2)
      DEV.OFF()
    }
    
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
      
      DEV.OFF()
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
      DEV.OFF()
    }
    
    
    vvs = cbind( filedatal$Area6,filedatal$Area8) %*% vv
    seldead = vvs < seldeadv
    
    vvsAll = cbind( filedatal$Area6,filedatal$Area8) %*% vv
    seldeadAll = ( vvsAll < seldeadv ) & selpredead
    
    DeadNum  = sum(seldeadAll & sel45 )
    
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
      DEV.OFF()
      
    }
    
    #CD48-------------------------------------------------------------
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
        DEV.OFF()
        
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
        DEV.OFF()
      }
      
      selCD3 = filedatal45pl$Area6 > CD3selectionBoundary
    }
    
    selCD3countes  = sum(selCD3)
    selCD3Live =  (!seldeadAll) & selCD3
    selCD3Dead =   seldeadAll & selCD3
    
    DeadCD3Area7dArea8b = 0
    DeadCD3Area7dArea8 = 0 
    if( sum(selCD3Dead) > 3 ){
      DeadCD3Area7dArea8  = round(median(filedatal[selCD3Dead,]$Area7) / median( filedatal[selCD3Dead,]$Area8 ),2)
      if( ( DeadCD3Area7dArea8 < 1.2 ) | ( DeadCD3Area7dArea8 > 2.1 ) )  {
        DeadCD3Area7dArea8b = 1
        #return(errorH(UseValidationMode,paste0("CD45Area1/Area2 =", deadCD3Area7dArea8),plotpath))
      }
    }else{
      
      DeadCD3Area7dArea8b = 1
      
    }  
    
    DeadCD3Area7dArea8 = c(DeadCD3Area7dArea8,DeadCD3Area7dArea8b)
    
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
    
    if(PLOTS){  
      
      abline(lm34,col = 3,lwd = 2)
      DEV.OFF()
      
    }
    
    
    
    #Create parallel for to cut tail for better differ cd4 and cd8 ----------------------
    
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
      
      DEV.OFF()
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
        DEV.OFF()
      }
    }else{
      
      cutpoint = min(vps34)
      if( PLOTS){
        
        abline(v = cutpoint,col = 1,lwd = 2,lty = 2 )
        
        DEV.OFF()
      }
      
    }  
    
    
    selNonTailCD48 =  vps34  >  cutpoint
    filedatalCD3_sel_tail = filedatalCD3[selNonTailCD48,]
    
    #Addtion120717----     
    #   }
    #   
    #   out = "x1"
    #   out
    # }
    # 
    # if(FALSE){
    #   if(FALSE){
    #     
    #------------------
    
    
    #Detect 4 and 8-------------------------------------------------------------
    # Use filedatalCD3_sel_tail to find 4 and by dividing
    bUseRatio = FALSE
    ratio3vd4 = filedatalCD3_sel_tail$Area3 - filedatalCD3_sel_tail$Area4
    h34 = hist(ratio3vd4,100,plot = F)
    h_34max = maxiNu(h34$mids,h34$counts,PLOTS,plotpath,"Detect h34 by ratio ch3-ch4",TRUE,50,10)
    if( length(h_34max$maxl$indx) == 2 ){
      
      
      middelValue = mean(h_34max$mins$x)
      ratio3vd4 = filedatalCD3$Area3 - filedatalCD3$Area4
      selType40r8 = ratio3vd4 < middelValue
      
      if(PLOTS){
        png(filename= paste0(plotpath,"CD3_Selection CD4_CD8 USE Ratio",".png"))
        with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = 0.2,col = colbr[1+!selType40r8[samp]],main = paste0(CarName," CD3_Selection CD4_CD8 USE Ratio")))
        DEV.OFF()  
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
        
        bounary4 = h_vvs34woTailvvs_max$mins$x[lenMax]
        bounary8 = h_vvs34woTailvvs_max$mins$x[2]
        #ADD IF 180717
        if (PLOTS){
          
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
        DEV.OFF()
        
      }
      
      #Create CD4 CD8 selection 
      vvs34 = with(filedatalCD3,cbind(Area4,Area3) %*% vv)
      #create between 4 and 8 selection variable
      selType40r8 = vvs34 > minbtween
      
      le = dim(filedatalCD3)[1]
      samp  = sample(1:le,min(le,1e4))
      
      if(PLOTS){
        
        png(filename= paste0(plotpath,"CD3_Selection CD4_CD8",".png"))
        with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = 0.2,col = colbr[1+!selType40r8[samp]],main = paste0(CarName," First Cut for CD4 and CD8")))
        DEV.OFF()  
      }
      
    }
    #Clean ONLY CD4 from Double Negtive --------------------------------------
    filedatalCD3only4 = filedatalCD3[!selType40r8,]
    bType8 = all( selType40r8 )
    bType4 = all( !selType40r8 )
    bOnlyOneType48 = bType8 | bType4
    
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
        DEV.OFF()
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
        DEV.OFF()
      }
      
      
      SumdoubleNegtoveBoundary  = sum(sel_doubleNegtoveBoundary) 
      
      vvslmonly4forall = with(filedatalCD3,cbind(Area1,Area3) %*% vv)
      sel_doubleNegtoveBoundaryIndForAll  = vvslmonly4forall > doubleNegtoveBoundary
      
      
      
    }
    else
    {
      
      sel_doubleNegtoveBoundaryIndForAll = rep( FALSE, dim(filedatalCD3)[1] )
      
    }
    
    #CD4 without Double Negtive
    selcleanCD4 = (!sel_doubleNegtoveBoundaryIndForAll  & (!selType40r8))
    
    #Find the Double positive ------------------------------------
    selpreDP  = !sel_doubleNegtoveBoundaryIndForAll
    
    if(PLOTS){
      png(filename= paste0(plotpath,"FinalArea4Area3",".png"))
      le = dim(filedatalCD3)[1]
      samp  = sample(1:le,min(le,1e4))
      
      with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = 0.2,col = colbr[1+(selpreDP + 2*selType40r8 )[samp]],main = paste0(CarName," Pre Double Positive ") ))
      legend("topleft",c("CD4DN","CD4","LoWCD8","CD8"),col = c(1,2,3,4),pch=19 )
      DEV.OFF()
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
          DEV.OFF()
        }
        
        filedatalCD3forProojection =  filedatalCD3[,c("Area4","Area3")]
        
        
        fldapreDPCD3 = fldaProjection(filedatalCD3forProojection,us[,1],us[,2],fl)
        
        sel_DP  = !((h_fldaP$mids[downind] < fldapreDPCD3 ) & (h_fldaP$mids[upind] > fldapreDPCD3 ))
        sel_DP = selpreDP & !sel_DP
        
        if(PLOTS){
          
          png(filename= paste0(plotpath,"DoublePositive_Area4_Area3",".png"))
          le = dim(filedatalCD3)[1]
          samp  = sample(1:le,min(le,2e4))
          with(filedatalCD3[samp,], plot(Area4,Area3,pch = 19, cex = c(0.2,0.6)[1  + sel_DP[samp]],col = colbr[1  + sel_DP[samp]  ],main = paste0(CarName," Select Double Positive") ))
          DEV.OFF()
        }
      }else{
        sel_DP = rep(FALSE,dim(filedatalCD3)[1])
      }
      
    }else{
      sel_DP = rep(FALSE,dim(filedatalCD3)[1])
    }
    
    sum(( selcleanCD4))
    selcleanCD4_NoDP =  ( selcleanCD4 ) & ( !sel_DP )
    selcleanCD8_NoDP  = selType40r8 & (!sel_DP )
    
    selcleanCD4_NoDPnum = sum(selcleanCD4_NoDP)
    selcleanCD8_NoDPnum = sum(selcleanCD8_NoDP)
    
    Area4DArea8b = 0
    Area4DArea8  = round(( median(filedatalCD3[selcleanCD4_NoDP,]$Area4) + median( filedatalCD3[selcleanCD4_NoDP,]$Area5 ) )/( median(filedatalCD3[selcleanCD8_NoDPnum,]$Area4) + median(filedatalCD3[selcleanCD8_NoDP,]$Area5) ),2)
    if ( (Area4DArea8 < 0.0674 )| (Area4DArea8 < 0.4385 ) ) {
      Area4DArea8b = 1
      #return(errorH(UseValidationMode,paste0("CD45Area1/Area2 =", Area4DArea8),plotpath))
    }
    Area4DArea8 = c(Area4DArea8,Area4DArea8b)
    
    CD4Area2dArea4b = 0
    CD4Area2dArea4  = round(median( filedatalCD3[selcleanCD4_NoDP,]$Area2 ) / median( filedatalCD3[selcleanCD4_NoDP,]$Area4 ),2)
    if( (CD4Area2dArea4 < 1.323 ) | (CD4Area2dArea4 >  5.3 ) ) {
      CD4Area2dArea4b = 1
      #return(errorH(UseValidationMode,paste0("CD45Area1/Area2 =", CD4Area2dArea4),plotpath))
    }
    CD4Area2dArea4 = c(CD4Area2dArea4,CD4Area2dArea4b)
    
    CD8Area4dArea5b = 0
    CD8Area4dArea5  = round(median(filedatalCD3[selcleanCD8_NoDP,]$Area4) / median( filedatalCD3[selcleanCD8_NoDP,]$Area5 ),2)
    if( ( CD8Area4dArea5 < 1.5 ) | ( CD8Area4dArea5 > 3.71 ) )  {
      CD8Area4dArea5b = 1
      #return(errorH(UseValidationMode,paste0("CD45Area1/Area2 =", CD8Area4dArea5),plotpath))
    }
    
    CD8Area4dArea5 = c(CD8Area4dArea5,CD8Area4dArea5b)
    
    rb = as.data.frame(rbind(Area1DArea2,DeadCD3Area7dArea8,Area4DArea8,CD4Area2dArea4,CD8Area4dArea5))
    rb$VariableName <-  rownames(rb)
    rownames(rb) <- NULL
    colnames(rb)  = c("Value","NotInRange","VariableName")
    rb = rb[, c("VariableName","Value","NotInRange")]
    
    if(PLOTS){
      
      write.csv(rb,paste0(plotpath,"Results.txt"),row.names = F)
    }
    
    
    if(PLOTS){graphics.off()}
    
    if(UseValidationMode == TRUE){
      
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
        round(selcleanCD4_NoDPnum/selcleanCD8_NoDPnum,2)) 
      
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
        "CD4CD8Ratio")
      
      df = data.frame(Values,Variable)
      
      if(WRITE_RESUTLS){
        
        write.csv(df,paste0(plotpath,"ReadMe.txt"))
        
      }
      
      df
      
    }else{  
      
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
        round(selcleanCD4_NoDPnum/selcleanCD8_NoDPnum,2)) 
      
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
        "CD4CD8Ratio")
      
      
      df = data.frame(Values,Variable)
      
      
      if(WRITE_RESUTLS){
        
        write.csv(df,paste0(plotpath,"ReadMe.txt"))
        
      }
      
      paste0(df[,c("Values")],collapse = ",")
    }
    
    
    #GENRAL   
  }else{
    
    if(UseValidationMode == TRUE) {
      
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
        "CD4CD8Ratio")
      
      Values  = rep(0,length(Variable))
      
      df = data.frame(Values,Variable)
      
      if(WRITE_RESUTLS){
        
        write.csv(df,paste0(plotpath,"ReadMe.txt"))
        
      }
      
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
        "CD4CD8Ratio") 
      
      Values  = rep(0,length(Variable))
      
      df = data.frame(Values,Variable)
      
      if(WRITE_RESUTLS){
        
        write.csv(df,paste0(plotpath,"ReadMe.txt"))
        
      }
      
      paste0(df[,c("Values")],collapse = ",")
      
    }
    
  }
}

dirpath = choose.dir(caption = "Select floder wit the list of ..._events.csv files")
timeStemp = gsub("-","_",gsub(":","_",Sys.time()))
fileList = dir(dirpath,full.names = T)
dataDirname  = paste0(dirpath,"\\data", timeStemp )
dir.create(dataDirname)
fileList = fileList[!(file.info(fileList))$isdir]

fll = fileList[2]

br = NULL
runfile = NULL
i = 0
for(fll in fileList){
  i = i+1
  #print(i)
  #print(basename(fl))
  #undebug(runAlgo_shortData)
  #re=runAlgo_shortData(fll,dataDirname,TRUE)
  result = TRUE
  result = tryCatch({
    re=runAlgo_shortData(fll,dataDirname,TRUE)
    TRUE
  }, warning = function(w) {
    TRUE
  }, error = function(e) {
    print("Error")
    FALSE
    
  }, finally = {
    print("final")
  })
  
  if( result ){
    
    br = rbind(br,re[,1])
    runfile = c(runfile,basename(fll))
  }
}

colnames(br) = re[,2]
br = data.frame(br,FILE.NAME = runfile)
write.csv(br,paste0(dataDirname,"\\Results",timeStemp,".csv"))

#runEnv ----- 
#runAlgo_shortData(args[1])
#library(readr)
#filelist <- read.csv("C:/Temp/AronRuns/filelist.csv",header  = T )

# br = NULL
# runfile = NULL
# 
# 
# 
# re = runAlgo_shortData(file.choose())
# 
# 
# i = 0 
# for(fl in filelist[1:7,1]){
#   i = i+1
#   print(i) 
#   d = dir(as.character(fl),full.names = T)
#   l = grep(".events.csv",d)
#   if(length(l) > 0 ){
#     print( basename( d[l] ) )
#     re = runAlgo_shortData( d[l] )
#     br = rbind(br,re[,1] )
#     runfile = c(runfile,basename(d[l]) )
#   }
# }
# 
# colnames(br) = re[,2]
# br = br1
# br <- data.frame(br)
# br = within( br,FILE.NAME <- (runfile) )
# 
# write.csv(br,"c:/Temp/AronRuns/Reuslts.csv")
# 
# re = runAlgo_shortData(file.choose())
# 
# View(br)
# 
# #re = runAlgo_shortData(args[1])
# #write.csv(re,"C:/Accellix/CheckResult.csv")