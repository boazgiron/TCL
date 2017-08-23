library("grDevices", lib.loc="C:/Program Files/R/R-3.2.1/library")
#Tools
isclose <- function(rein){
  
  x1 = c(1,length(rein$x))
  de = sqrt(diff(rein$x[x1])^2 + diff(rein$y[x1])^2)
  out = de == 0
  out
}

closeOpen <- function(rein){
  rein = re
  rout = rein
  if(!isclose(rein)){
    
    rout$x = c(rout$x,rout$x[1]) 
    rout$y = c(rout$y,rout$y[1]) 
  
  }
  
  rout

}

contourLength <- function(rein){
  
  c(sum(sqrt(diff(rein$x)^2 + diff(rein$y/100)^2)),rein$level)
}

contourLength1 <- function(rein){
  
  c(sqrt(diff(range(rein$x))^2 + diff(range(rein$y)/1000)^2),rein$level )
}



Getlevles <- function(rein){
  rein$level
}

#function-----
conture<-function(plotpathin,h){
  
  #plotpathg = plotpathin
  zseq = quantile(z,seq(0.7,0.989,0.01))
  #h2$z = sigmoid(h2$z/quantile(h2$z,0.8))
  #z= h2$z
  nlevels = 1
  le = NULL
  
  png(filename = paste0(plotpathin,"image","_CD45_", "Area1","_vs_","FCS",".png"))
  image(h)
  dev.off()

  png(filename = paste0(plotpathin,"LastChecounter","_CD45_", "Area1","_vs_","FCS",".png"))
  first = T
  j=0
  for( z1 in zseq ){
  
    j= j+1
    #z1 = zseq[i]
    re = contourLines(x = min(h$x) + diff(range(h$x))*seq(0, 1, length.out = nrow(h$z)),
                      y = min(h$y) + diff(range(h$y))*seq(0, 1, length.out = ncol(h$z)),
                      h$z, nlevels = nlevels,
                      #levels = pretty(range(z, na.rm = TRUE), nlevels),
                      levels = z1)
    relengths = t(sapply(re,contourLength1))[,1]
    
    le = c(le ,length(re[relengths > 1]))
    #table(gl)
    
    if(first){
      
      plot(re[[1]]$x,re[[1]]$y,"l",xlim = c(-1,3.5),ylim = c(-1,3e3))#,col = cl$cluster[1])
      first = FALSE
      
    }else{
      
      if(length(re) > 0){
        if(j%%3 == 0){
          for( i in 1:length(re)){
          
            points(re[[i]]$x,re[[i]]$y,"l",xlim = c(-1,3.5),ylim = c(-1,3e3))#,col = cl$cluster[1])
          }
        }
      }
    }
    
  }
  
  z1 = NULL
  z1 = zseq[min(which(max(le) == le))]
  
  # wh = which(diff(diff(ifelse(le == 2,1,0)) ) == -1 ) + 1
  # if( length(wh) > 0 ){
  #   
  #   z1 = zseq[min(wh)]  
  #   
  # }else{
  #   
  #   z1 = zseq[min(which(le == 2))]
  #   
  # }
  # 
  # if(length(z1) > 0 ){
  #   
  #   wh = which(diff(diff(ifelse(le >= 2,1,0)) ) == -1 ) + 1
  #   if( length(wh) > 0 ){
  #     
  #     z1 = zseq[min(wh)]  
  #     
  #   }else{
  #     
  #     z1 = zseq[min(which(le == 3))]
  #     
  #   } 
  #   
  # }
  
  
  #diff(which(le == 2)) == 0
  
  #z1 = zseq[min(which(le == 2))]
  
  if(length(z1) > 0 ){
    
    re = contourLines(x = min(h$x) + diff(range(h$x))*seq(0, 1, length.out = nrow(z)),
                    y = min(h$y) + diff(range(h$y))*seq(0, 1, length.out = ncol(z)),
                    h$z, nlevels = nlevels,
                    #levels = pretty(range(z, na.rm = TRUE), nlevels),
                    levels = z1)
    if(length(re)>0){  
      for( i in 1:length(re)){
      
        points(re[[i]]$x,re[[i]]$y,"l",xlim = c(-1,3.5),ylim = c(-1,3e3),col = 2,lwd = 2)#,col = cl$cluster[1])
      }
    }
    
  }  
  dev.off()
  
  out = z1
  
  out
}





CheckConture<-function(wrkingFilepath,datadir){
    
    #datadir = dataDirname
    #wrkingFilepath = fll
    filedata = read.csv(wrkingFilepath,header = T)
    
    runInformation  = getCartrigeNameAdnDevice(wrkingFilepath)
    
    CarName = paste0("C",runInformation$Cartnum)
    
    Cartnum = runInformation$Cartnum
    
    plotpathin = paste0(datadir,"/",Cartnum,"/")
    
    if(PLOTS){  
      if(!dir.exists(plotpathin)){
        dir.create(plotpathin)
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
    
    
    print(dim(filedatals))
    
    s = sample(nrow(filedatals),min(1e4,nrow(filedatals)))
    
    hin <- kde2d(filedatals[s,"Area1"], filedatals[s,"FCS"], n = 30,
                lims = c(-1, max(filedatals[s,"Area1"]) ,
                         c(-1,quantile(filedatals[s,"FCS"],0.99))))
    
    conture(plotpathin,hin)
    
    rm(filedatals)
}



#---------------------------------

gl = sapply(re,Getlevles)
minindex  = min(length(re) - which(table(gl) == 1))

glb = (table(gl) == 1)
glb[9]

table(gl)


plot(t(sapply(re,contourLength)),col = factor(table(gl)))




s = sample(nrow(filedatals),min(1e4,nrow(filedatals)))
h2 <- kde2d(filedatals[s,"Area1"], filedatals[s,"FCS"], n = 30,
                       lims = c(-1, max(filedatals[s,"Area1"]) ,
                                c(-1,quantile(filedatals[s,"FCS"],0.99))))
#h2$z = 40*(h2$z- min(h2$z))/diff(range(h2$z))
#h2$z = sigmoid(h2$z/5)

image(h2)

h2$z[h2$z < 1] = 0
h2$z[h2$z > 1] = 1
z = h2$z
h2$z = 40*sigmoid(z/5)
 
nlevels = 3


z1 = quantile(z,seq(0.6,0.989,0.05))


#z1 = quantile(z,seq(0.95))
#z1  = max(z)-seq(0.1,0.5,0.05)
nlevels = length(z1)

re = contourLines(x = min(h2$x) + diff(range(h2$x))*seq(0, 1, length.out = nrow(z)),
             y = min(h2$y) + diff(range(h2$y))*seq(0, 1, length.out = ncol(z)),
             z, nlevels = nlevels,
             #levels = pretty(range(z, na.rm = TRUE), nlevels),
             levels = z1)

re = lapply(re,closeOpen)

points(re[[13]]$x,re[[13]]$y,"l",col = ifelse(isclose(re[[minindex-1]]),2,1))

re = re[which(sapply(re,contourLength)[1,] < 6100)]

#plot(re[[1]]$x,re[[1]]$y,"l",xlim = c(-1,3.5),ylim = c(-1,3e3),col = ifelse(isclose(re[[1]]),2,1))
plot(re[[1]]$x,re[[1]]$y,"l",xlim = c(-1,3.5),ylim = c(-1,3e3))#,col = cl$cluster[1])
mea = NULL

for(i in 2:length(re)){
  #points(re[[i]]$x,re[[i]]$y,"l",col = ifelse(isclose(re[[i]]),2,1))
  points(re[[i]]$x,re[[i]]$y,"l",col = cl$cluster[i])
  mes = c( mean(re[[i]]$x),mean(re[[i]]$y))
  mea = rbind(mea ,mes)
  #points(mes[1],mes[2],pch = 19,col = 2 + ifelse(isclose(re[[i]]),2,1))
  points(mes[1],mes[2],pch = 19,col = cl$cluster[i])
}

cl = kmeans(mea,4)
plot(mea,col = cl$cluster)


h2$z = hist(as.vector(sigmoid(h2$z/1e-3)),100)
h2$z = sigmoid(h2$z)
h2$z = hist(as.vector(h2$z),100)

bu =NULL
for( i in 1:length(re)){
  da = sum(sqrt(diff(re[[i]]$x)^2 + diff(re[[i]]$y)^2))
  x1 = c(1,length(re[[i]]$x))
  de = sqrt(diff(re[[i]]$x[x1])^2 + diff(re[[i]]$y[x1])^2)
  bu  = c(bu,de/da)
  
}


sapply(re,isclose)

i = 2
points(re[[i]]$x,re[[i]]$y,"l")
points(mean(re[[i]]$x),mean(re[[i]]$y),pch = 2,cex = 2,col =5)

plot(diff(re[[2]]$x))

plot(re[[2]]$x)

table(unlist(lapply(re,function(x){x$level})))

aggregate(re,list(re$level),length)

dirpath = choose.dir(caption = "Select floder wit the list of ..._events.csv files")
timeStemp = gsub("-","_",gsub(":","_",Sys.time()))
fileList = dir(dirpath,full.names = T)
dataDirname  = paste0(dirpath,"\\data", timeStemp )
dir.create(dataDirname)
fileList = fileList[!(file.info(fileList))$isdir]

fll = fileList[1]

for(fll in fileList[10:20]){
  print(fll)
  #re= checkdis(fll,dataDirname)
  CheckConture(fll,dataDirname)
  
  #br = rbind(br,re[,1])
  #runfile = c(runfile,basename(fll))
  
}



