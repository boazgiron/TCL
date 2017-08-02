## Cluster the spirals data set.
library("kernlab", lib.loc="~/R/win-library/3.2")
library("e1071")
#data(spirals)
#write.csv(filedatals,"C:/Temp/filedatals.csv")
filedatals  = read.csv("C:/Temp/filedatals.csv",header = T,stringsAsFactors = F)
filedatalsSample = filedatals[sample(nrow(filedatals),2e4),]
write.csv(filedatalsSample,"C:/Project/LeukoDx/R/TCL/filedatals.csv")

bu1 = NULL
bu2 = NULL
st = Sys.time()
for(i in 1:10){
 
  filedatalsSample = filedatals[sample(nrow(filedatals),5e2),]
  sp = scale(cbind(filedatalsSample$sumCh,filedatalsSample$FCS))
  sc <- specc(sp,kernel = "laplacedot",kpar = list(sigma= 0.4), centers=2)
  
  bu1  = cbind(bu1,aggregate(sp[,1],list(sc),mean)[,2])
  bu2  = cbind(bu2,aggregate(sp[,2],list(sc),mean)[,2])

}

st = Sys.time()
dt = Sys.time() -st
dt
x = cbind(apply(t(bu1),1,max),apply(t(bu1),1,min))
y = cbind(apply(t(bu2),1,max),apply(t(bu2),1,min))

plot(x[,1],y[,1],xlim = c(-1,2.5),ylim = c(-1,2.5))
points(x[,2],y[,2])

plot(sp)

sc <- specc(sp, centers=2)

sc
centers(sc)
size(sc)
withinss(sc)

plot(sp, col=sc,pch = 19)





filedatalsSample = filedatals[sample(nrow(filedatals),5e2),]
#sp = scale(cbind(filedatalsSample$sumCh,filedatalsSample$Area7))
sp = (cbind(filedatalsSample$sumCh,filedatalsSample$Area7))
sc <- specc(sp,kernel = "laplacedot",kpar = list(sigma= 1), centers=2)
df = data.frame(sp,cl = factor(sc@.Data))
colnames(df)[1:2] <- c("Area7","sumCh") 

svm_model <- svm(cl ~ ., data=df)
summary(svm_model)
plot(svm_model,data=df)


filedatalsSample = filedatals[sample(nrow(filedatals),1e2),]
sp = scale(cbind(filedatalsSample$sumCh,filedatalsSample$FCS))
sc <- specc(sp,kernel = "laplacedot",kpar = list(sigma= 0.4), centers=2)
df = data.frame(sp,cl = factor(sc@.Data))


filedatalsSample = filedatals[sample(nrow(filedatals),1e3),]
sp = scale(cbind(filedatalsSample$sumCh,filedatalsSample$FCS))
df = data.frame(sp)


filedatalsSample = filedatals[sample(nrow(filedatals),5e2),]
#sp = scale(cbind(filedatalsSample$sumCh,filedatalsSample$Area7))
sp = scale(cbind(filedatalsSample$sumCh,filedatalsSample$FCS))
sc <- specc(sp,kernel = "laplacedot",kpar = list(sigma= 0.4), centers=2)
df = data.frame(sp,cl = factor(sc@.Data))

svm_model <- svm(cl ~ ., data=df)
summary(svm_model)
plot(svm_model,data=df)

filedatalsSample = filedatals[sample(nrow(filedatals),1e4),]
sp = scale(cbind(filedatalsSample$sumCh,filedatalsSample$FCS ))
df = data.frame( sp )
cl = predict(svm_model,df)
plot(sp,col = cl,pch =19)


filedatalsSample = filedatals[sample(nrow(filedatals),5e2),]
sp = scale(cbind(filedatalsSample$sumCh,filedatalsSample$FCS))
sc <- specc(sp,kernel = "laplacedot",kpar = list(sigma= 0.4), centers = 2 )
df = data.frame(sp,cl = factor(sc@.Data))

svm_model <- svm(cl ~ ., data=df)
summary(svm_model)
plot(svm_model,data=df)


filedatalsSample = filedatals[sample(nrow(filedatals),2e4),]
filedatals$sumCh[!is.finite(filedatals$sumCh)] <- NA
sp = scale(cbind(filedatalsSample$sumCh,filedatalsSample$FCS))
df = data.frame(sp)
cl = predict(svm_model,df)

plot(filedatalsSample$sumCh,filedatalsSample$FCS,col = cl,pch = 19 ,cex =0.2)



getsumCH <-function(){
  
  
  


filedatals$sumCh[!is.finite(filedatals$sumCh)] <- NA
filedatalsSample = filedatals[sample(nrow(filedatals),2e4),]
sp1 = (cbind(filedatals$sumCh,filedatals$FCS))
dfa = data.frame(sp1)

  
cla = NULL
st = Sys.time()
for(i in 1:5){
  
  filedatalsSample = filedatals[sample(nrow(filedatals),1e3),]
  df = cbind(filedatalsSample$sumCh,filedatalsSample$FCS)
  mes = colMeans(df)
  mds = apply(df,2,sd)
  sp = scale(df)
  
  dff = subformcol(dfa,mes)
  dff = divformcol(dff,mds)
  sc <- specc(sp,kernel = "laplacedot",kpar = list(sigma= 0.4), centers = 2 )
  df = data.frame(sp,cl = factor(sc@.Data))
  
  
  svm_model <- svm(cl ~ ., data = df)
  plot(svm_model,data=df)
  
  pre = predict(svm_model,df)
  me = with(df,aggregate(X1,list(pre),mean))[,2]
  men = which(max(me) == me)
  if( men ==  2){
    pre = predict(svm_model,dff)  
    pre = ifelse(pre == 1,2,1)  
  }
  cla = rbind(cla,pre)
  
}

head(dff)

dt = Sys.time() - st
dt

cl  = ifelse(sapply(cla,mean) > 1.5,2,1)

plot(dfa$X1,dfa$X2,col = cla[2,],pch = 19 ,cex = 0.3,ylim = c(-1,7))

plot(dff$X1,dff$X2,col = pre,pch = 19 ,cex = 0.3 )

plot(df[,1:2],col = df$cl)

#020817DD---------------------------------------

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



## Cluster the spirals data set.
library("kernlab", lib.loc="~/R/win-library/3.2")
library("e1071")
#data(spirals)
#write.csv(filedatals,"C:/Temp/filedatals.csv")
filedatals  = read.csv("C:/Temp/filedatals.csv",header = T,stringsAsFactors = F)

filedatals$sumCh[!is.finite(filedatals$sumCh)] <- NA
filedatalsSample = filedatals[sample(nrow(filedatals),5e4),]


#sp1 = (cbind(filedatals$sumCh,filedatals$FCS))
#dfa = data.frame(sp1)


cla = NULL
st = Sys.time()

geG <- function(){


  filedatalsSample = filedatals[sample(nrow(filedatals),1e3),]
  df = cbind(filedatalsSample$sumCh,filedatalsSample$FCS)
  mes = colMeans(df,na.rm = T)
  mds = apply(df,2,sd,na.rm = T)
  sp = scale(df)

  sc <- specc(sp,kernel = "laplacedot",kpar = list(sigma= 0.4), centers = 2 )
  df = data.frame(sp,cl = factor(sc@.Data))
  #with(df,plot(X1,X2,col = cl,pch = 19 ,cex = 0.3,ylim = c(-1,7)))
  svm_model <- svm(cl ~ ., data = df)# ,probability = TRUE)
  sp1 = (cbind(filedatals$sumCh,filedatals$FCS))
  dfa = data.frame(sp1)
  dff = subformcol(dfa,mes)
  dff = divformcol(dff,mds)
  dff = dff[complete.cases(dff),]
  pre = predict(svm_model,dff)#,probability = TRUE)
  
  g1xx = mean(dff[pre == 1,1])
  g1yy = mean(dff[pre == 1,2])
  
  g2xx = mean(dff[pre == 2,1])
  g2yy = mean(dff[pre == 2,2])
  
  
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
  
  if(g1xx > g2xx ){
    pre = ifelse(pre == 1,2,1)
  }
  
  if(PLOTS){
    
    png(filename= paste0(plotpath,sumCH_SVM,".png"))
    plot(dff,col = pre,pch = 19 ,cex = 0.3,ylim = c(-1,10))
    dev.off()
  }
      
  pre
  
}  
#dff[which(diff(as.numeric(names(pre)))>1)+1,]

str(pre)

dt = Sys.time() - st
dt

#colnames(dfa)[1:2] <- c("sumCH","FCS")



# pre = predict(svm_model,df)
# me = with(df,aggregate(X1,list(pre),mean))[,2]
# men = which(max(me) == me)
# if( men ==  2){
#   pre = predict(svm_model,dff)  
#   pre = ifelse(pre == 1,2,1)  
# }
# cla = rbind(cla,pre)
  

cl  = ifelse(sapply(cla,mean) > 1.5,2,1)

plot(dfa$X1,dfa$X2,col = cla[2,],pch = 19 ,cex = 0.3,ylim = c(-1,7))

plot(dff$X1,dff$X2,col = pre,pch = 19 ,cex = 0.3 )


