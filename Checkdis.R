library("e1071")
createdfz <-function(h){
  zr =NULL
  for(i in 1:length(h$x)){
    for(j in 1:length(h$y)){
      zr = rbind(zr , c(i,j,h$x[i],h$y[j],h$z[i,j]))    
    }
  }
  
  zr
}
fsum <- function(x){sum(x,na.rm = T)}
is.notinfinite.any <- function(x){!any( is.infinite(x) ) }

grouping <-function(clnames,plotpathg,lastdev){
  
  s = sample(nrow(filedatals),min(1e4,nrow(filedatals)))
  
  i = 0
  cont =  TRUE
  h2 = NULL
  n = 25
  while( cont & (i < n ) ) {
    
    #plot(filedatals[s,clnames],pch =19,cex= 0.3)
    
    hlast = h2
    h2 <- kde2d(filedatals[s,clnames[1]], filedatals[s,clnames[2]], n = 30,lims = c(0, max(filedatals[s,clnames[1]]) ,c(1e2,quantile(filedatals[s,clnames[2]],0.99))))
    h2$z = 2*(n - i) *(h2$z- min(h2$z))/diff(range(h2$z))
    
    h2$z[h2$z < 1] = 0
    h2$z[h2$z > 1] = 1
    h2$z = GetLabel(h2$z)
    unique(as.vector(h2$z))
    
    dfh = as.data.frame(createdfz(h2))
    colnames(dfh) <- c("indx","indy","x","y","z")
    
    ag = aggregate(dfh$x,list(dfh$z),mean)
    agl = aggregate(dfh$x,list(dfh$z),length)
    agl = agl[!(agl$Group.1 == 0),]
    agl = agl[order(agl$x,decreasing = T),]
    agl = agl[agl$x > 10,]
    
    cont  = length(agl$x) < 2
    i= i + 1
    
  }
  
  #getMaximumGroup
  ta = table(as.vector(h2$z))  
  ta = ta[!names(ta) =="0"]
  wh = which(max(ta) == ta)
  maxgroupnumber = as.numeric(names(ta)[wh])
  
  if(PLOTS) {
    if(is.null(hlast)){
      hlast  = h2
    }
    png(filename= paste0(plotpathg,"LastChekcdis","_CD45_", clnames[2],"_vs_",clnames[1],".png"))
    image(hlast,breaks = 1:(n+1),col = rainbow(n) )
    dev.off()
    
    png(filename= paste0(plotpathg,"Chekcdis","_CD45_", clnames[2],"_vs_",clnames[1],".png"))
    image(h2,breaks = 1:(n+1),col = rainbow(n) )
    points(ag[2:3,"x"],ag[2:3,"y"])
    dev.off()
  }
  
  
  dfh = dfh[dfh$z != 0,]
  dfh$z = factor(ifelse( dfh$z ==  maxgroupnumber,1,2))
  dfh = dfh[,c("x","y","z")]
  
  svm_model <- svm(z ~ ., data=dfh)
  summary(svm_model)
  df = data.frame(x = filedatals[s,clnames[1]],y= filedatals[s,clnames[2]])
  pre = predict(svm_model,df)
  
  if(PLOTS){
    
    png(filename= paste0(plotpathg,"CD45_", clnames[2],"_vs_",clnames[1],"_SVM.png"))
    with(df,plot(x,y, col = unlist(pre),pch = 19,cex = 0.3))
    if(lastdev){
      dev.off()
    }
  }
  
  
  df = data.frame(x = filedatals[,clnames[1]],y= filedatals[,clnames[2]])
  pre = predict(svm_model,df)
}

checkdis <-function(wrkingFilepath,datadir){
  
  
  #wrkingFilepath = fll
  filedata = read.csv(wrkingFilepath,header = T)
  
  runInformation  = getCartrigeNameAdnDevice(wrkingFilepath)
  
  CarName = paste0("C",runInformation$Cartnum)
  
  Cartnum = runInformation$Cartnum
  
  plotpathin = paste0(datadir,"/",Cartnum,"/")
  
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

  
  filedatals= filedatals[unlist(apply(filedatals,1,is.notinfinite.any)),]
  filedatals= filedatals[complete.cases(filedatals),]
  
  
  filedatals$sumCh =  apply(filedatals[,1:8],1,fsum )
  
  #filedata[which(apply(filedata,1,is.infinite.any)),]
  
  clnames = c("sumCh","FCS")
  m1  = grouping(clnames,plotpathin,FALSE)
  
  #df = with(filedatals,data.frame(x=sumCh,y=FCS))
  #pr1 = predict(m1,df) 
  me = aggregate(filedatals$sumCh,list(m1),mean)
  AverageLowGroup = min(me$x) 
  abline(v = AverageLowGroup ,col = 3 ,lwd =2,lty =2)
  dev.off()
  filedatals  = filedatals[filedatals$sumCh >  AverageLowGroup,]
  
  clnames = c("Area1","FCS")
  m2 = grouping(clnames,plotpathin,TRUE)
  #df = with(filedatals,data.frame(x=Area1,y=FCS))
  #pr2 = predict(m2,df) 
  ##DD210817xxx

  #hist(filedatals$Area1)
  #hist(filedatals$Peak9)

}


dirpath = choose.dir(caption = "Select floder wit the list of ..._events.csv files")
timeStemp = gsub("-","_",gsub(":","_",Sys.time()))
fileList = dir(dirpath,full.names = T)
dataDirname  = paste0(dirpath,"\\data", timeStemp )
dir.create(dataDirname)
fileList = fileList[!(file.info(fileList))$isdir]

for(fll in fileList[1:1]){
  
  #re= checkdis(fll,dataDirname)
  checkdis(fll,dataDirname)
  #br = rbind(br,re[,1])
  #runfile = c(runfile,basename(fll))
  
}



