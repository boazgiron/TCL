args<-commandArgs(TRUE)

s = sample(nrow(filedatals),1e3)

# plot_ly(filedatals[s,], x = ~Area1, y = ~Area2, z = ~Area3)#, color = 1,colors  = c('#BF382A'))
# chart_link = plotly_POST(p, filename="scatter3d/colorscale")
# chart_link

# with(filedatals[s,],scatterplot3d(Area1, Area2, Area3, highlight.3d=TRUE, col.axis="blue",
#               col.grid="lightblue", main="scatterplot3d - 1", pch=20))

dff = cbind(x1=scale(filedatals$Area1),
            y1=scale(filedatals$FCS),
            z1=scale(filedatals$sumCh)) 

x1=scale(dff[s,1])
y1=scale(dff[s,2])
z1=scale(dff[s,3])     




col=rainbow(1000)
spo  = cbind(x1,y1,z1)
sco <- specc(spo,kernel = "laplacedot",kpar = list(sigma= 0.4), centers = 2 )
plot3d(x1, y1, z1,xlab = "Area1", ylab  = "FCS", zlab = "sumCh",col = factor(sco@.Data)) 


df = data.frame(spo,cl = factor(sco@.Data))
#with(df,plot(X1,X2,col = cl,pch = 19 ,cex = 0.3,ylim = c(-1,7)))
svm_model <- svm(cl ~ ., data = df)# ,probability = TRUE)
str(svm_model)


pre = predict(svm_model,dff)#,probability = TRUE)
dff = as.data.frame(dff)




for( i in 1:3){
  
  filedatals[is.infinite(filedatals[,i]),i] = NA
}

sel4 = complete.cases(filedatals)
filedatals = filedatals[sel4,]

#filedatals$Area11 <- filedatals$sumCh *1/10

selv = which(filedatals$sumCh > 1)
filedatals1 = filedatals[selv,]
out = geG(filedatals1,c("Area1", "FCS" ,"sumCh"),T,plotpath)

out1  = rep(1,dim(filedatals)[1])

for(i in 1:length(selv)){
  
  out1[selv[i]]  <- out[i]
}


with(filedatals,plot3d(Area1, FCS , sumCh,xlab = "Area1", ylab  = "FCS", zlab = "sumCh",col = factor(out1)) )

dff = as.data.frame(dff)

  geG(filedatals,c("Area1", "FCS" ,"sumCh"),T)  
  

par(mfrow = c(2,1))
with( dff, plot(z1,y1,pch = 19, cex =0.2,col  = factor(pre) ) )
with( dff, plot(x1,y1,pch = 19, cex =0.2,col  = factor(pre) ) )
#with( filedatals, plot(sumCH,FCS,pch = 19, cex =0.2,col  = factor(sco@.Data)))
#with( filedatals, plot(Area1,FCS,pch = 19, cex =0.2,col  = factor(sco@.Data)))


rdirname = dir("L:/Scientific/Raw Data/Reader_Data/TCL Validation/accellix17/",full.names = T)

f <- function(x){getCartrigeNameAdnDevice(x)$Cartnum}
filePath =file.choose()
fcsf = read.csv(filePath,header = T)
colnames(fcsf)
#fcsfcn = sapply(basename(as.character(fcsf$FILE.NAME)),f)
fcsfcn =  sapply(as.character(fcsf$FILE.NAME),f)

names(fcsfcn) = NULL
ludaResult = data.frame(fileName = fcsf$FILE.NAME, Type = fcsf$Type, FILE.NAME = fcsfcn,facsbu45 = fcsf$X.CD45.TOT)



#ludaResult_cn = sapply(as.character(basename(ludaResult$FILE.NAME)),f)
list_cn = sapply(basename(rdirname),f)
names(list_cn) = NULL
rdirname = rdirname[ list_cn %in%  ludaResult$FILE.NAME ]

fl = rdirname[grep("3300",basename(rdirname))]

fl = rdirname[1]
ou = checkf(rdirname[34])
sel45 = ou == 2
sum(sel45)

#bu45 = NULL
#bufl = NULL
for(fl in rdirname[101:112] ){
  
   ou = checkf(fl)
   sel45 = ou == 2
   sum45 = sum(sel45)
   bu45 = c(bu45,sum45)
   print(paste0(getCartrigeNameAdnDevice(basename(fl))$Cartnum," - ", sum45)) 
   bufl = c(bufl,fl)
}

keepBu45  = bu45
keepbufl  = bufl
#le = length(bu45)
#bu45 = bu45[-le]
#bufl = bufl[-le]
cn = sapply(basename(bufl),f)
names(cn) = NULL
cn





df = data.frame(FILE.NAME = cn ,auto   =  bu45 )
#cn = sapply(basename(keepbufl),f)
#names(cn) = NULL
#df = data.frame(FILE.NAME = cn ,auto   =  keepBu45 )



me = merge(ludaResult,df)
me = data.frame(me)

resultWitFilter.path = file.choose()
resultWitFilter.data = read.csv(resultWitFilter.path,header = T ,stringsAsFactors = T)

View(resultWitFilter.data)

write.csv(me,"C:/Project/LeukoDx/LudaFacsValidation/Debug1408/ResultwithShortEvents.csv")

View(me)

results2 = data.frame(fullname  = basename(bufl),FILE.NAME = cn,bu45)
write.csv(results2,"C:/Project/LeukoDx/LudaFacsValidation/Debug0708/Results080917_1129sumCh_ch1_.csv")
re2 = results2[ list_cn %in% ludaResult$FILE.NAME, ]

result_1 = read.csv(file.choose(),header = T,stringsAsFactors = F)

re2 = results2[ list_cn %in% ludaResult$FILE.NAME, ]
me2 = merge(ludaResult, re2 ,by = "FILE.NAME")
dim(me2)
re2r  = round(100*(me2$bu45 - me2$facsbu45 )/me2$facsbu45,2)
hist(re2r)

length(rdirname)

me2[me2[re2r < -20,]$FILE.NAME %in% ludaResult$FILE.NAME,]$FILE.NAME

#ludaResult[grep("C427",as.character(fcsf$FILE.NAME)),]



length(me2[re2r < -20,]$FILE.NAME)


me = merge(ludaResult, re2 ,by = "FILE.NAME")
re1 = result_1[ list_cn %in% ludaResult$FILE.NAME, ]
me1 = merge(ludaResult, re1 ,by = "FILE.NAME")

me = merge(ludaResult, re2 ,by = "FILE.NAME")


rer =  round(100*(me$bu45 - me$facsbu45 )/me$facsbu45,2)
re1r  = round(100*(me1$bu45 - me1$facsbu45 )/me$facsbu45,2)


library("kernlab", lib.loc="~/R/win-library/3.2")
library("e1071", lib.loc="~/R/win-library/3.2")
library("rgl", lib.loc="~/R/win-library/3.2")
library("MASS",lib.loc="~/R/win-library/3.2")

checkf <- function(wrkingFilepath){
  
  #wrkingFilepath = rdirname[32]
  
  filedata = read.csv(wrkingFilepath,header = T)
  di = dim(filedata)
  
  #Extract Cartrige num  
  runInformation  = getCartrigeNameAdnDevice(wrkingFilepath)
  
  CarName = paste0("C",runInformation$Cartnum)
  
  Cartnum = runInformation$Cartnum
  
  plotpath = paste0("C:/Project/LeukoDx/LudaFacsValidation/Debug1608/",Cartnum,"/")
  
  if( PLOTS ){  
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
  
  filedatals = filedata[widthb,]
  
  filedatals =  data.frame(apply(filedata[,ty[1:9]],2,log10))
  filedatals$FCS <- filedata$Peak9
  fsum <- function(x){sum(x,na.rm = T)} 
  
  
  #filedata = filedata[widthb,]
  #------------------------------  
  
  for( i in 1:dim(filedatals)[2]){
    
    filedatals[is.infinite(filedatals[,i]),i] = NA
  }
  
  sel4 = complete.cases(filedatals)
  filedatals = filedatals[sel4,]
  
  filedatals$sumCh =  apply(filedatals[,2:8],1,fsum )
  
  if( dim(filedatals)[1] > 100 ){
    
    #out = geG(filedatals,c("Area1", "Area7" ,"sumCh"),T,plotpath)
    selv = which(filedatals$sumCh > 1 & filedatals$Area1 > 1)
    
    filedatals1 = filedatals[selv,]
    CNAME_in = c("Area1","FCS","sumCh", "Area7" )
    CNAME_in = c("Area1","FCS","sumCh", "Area7" )
    out = geG(filedatals1,CNAME_in,T,plotpath,c(1,1,1,1))
    out1  = rep(1,dim(filedatals)[1])
    
    for(i in 1:length(selv)){
      
      out1[selv[i]]  <- out[i]
    }
    
    png(filename= paste0(plotpath,CNAME_in[1],"_",CNAME_in[2],"_SVM.png"))
    plot(filedatals[,c(CNAME_in[1],CNAME_in[2])],col = out1,pch = 19 ,cex = 0.2,xlab = CNAME_in[1],ylab = CNAME_in[2])
    #points(c(g1x,g2x),c(g1y,g2y),col = 3, cex = 2,pch = 19)
    dev.off()
    
    png(filename= paste0(plotpath,CNAME_in[3],"_",CNAME_in[2],"_SVM.png"))
    plot(filedatals[,c( CNAME_in[3] , CNAME_in[2] )],col = out1,pch = 19 ,cex = 0.2,xlab = CNAME_in[3],ylab = CNAME_in[2])
    plot(filedatals[,c( CNAME_in[3] , CNAME_in[2] )],col = out1,pch = 19 ,cex = 0.2,xlab = CNAME_in[3],ylab = CNAME_in[2])
    dev.off()
    
  }else{
    
    print(paste0("No points !!!",dim(filedatals)[1]))
  
  }
  
  out1
  
}


checkfs <- function(plotpath, CartNum ){
  
  if( dim(filedatals)[1] > 100 ){
    
    #out = geG(filedatals,c("Area1", "Area7" ,"sumCh"),T,plotpath)
    selv = which(filedatals$sumCh > 1 & filedatals$Area1 > 1)
    
    filedatals1 = filedatals[selv,]
    CNAME_in = c("Area1","FCS","sumCh", "Area7" )
    CNAME_in = c("Area1","FCS","sumCh", "Area7" )
    out = geG(filedatals1,CNAME_in,T,plotpath,c(1,1,1,1))
    out1  = rep(1,dim(filedatals)[1])
    
    for(i in 1:length(selv)){
      
      out1[selv[i]]  <- out[i]
    }
    
    png(filename= paste0(plotpath,CNAME_in[1],"_",CNAME_in[2],"_SVM.png"))
    plot(filedatals[,c(CNAME_in[1],CNAME_in[2])],col = out1,pch = 19 ,cex = 0.2,xlab = CNAME_in[1],ylab = CNAME_in[2])
    #points(c(g1x,g2x),c(g1y,g2y),col = 3, cex = 2,pch = 19)
    dev.off()
    
    png(filename= paste0(plotpath,CNAME_in[3],"_",CNAME_in[2],"_SVM.png"))
    plot(filedatals[,c( CNAME_in[3] , CNAME_in[2] )],col = out1,pch = 19 ,cex = 0.2,xlab = CNAME_in[3],ylab = CNAME_in[2])
    plot(filedatals[,c( CNAME_in[3] , CNAME_in[2] )],col = out1,pch = 19 ,cex = 0.2,xlab = CNAME_in[3],ylab = CNAME_in[2])
    dev.off()
    
  }else{
    
    print(paste0("No points !!!",dim(filedatals)[1]))
    
  }
  
  out1
  
}



#filedatals$Area11 <- filedatals$Peak9
#Working 140817

# s = sample(nrow(filedatals),2e3)
# f1 <- with(filedatals,kde2d(Area1, FCS, n = 30 ) )# , lims = c(0.5, 6, 40, 100))
# s = scale(f1$z) * 10
# s= s+min(s)
# f1$z <- s
# with(f1,image(x,y,z)) #, zlim = c(0, 0.05))
# persp(f1, phi = 30, theta = 20, d = 5)
# out1 = NULL
# 
#with(filedatals, plot3d(Area1, sumCh, FCS,pch = 19,cex =0.2,col = out1))#, 

#f1 <- with(filedatals,kde2d(Area1, FCS, 10, lims = c(0,5,0,1e4 ) ) )# , lims = c(0.5, 6, 40, 100))
#f1$z = filedatalsf1$z/max(f1$z)
#f1$z[f1$z < 0.1] = 0
#image(f1,col = c(1,2,3,4))
