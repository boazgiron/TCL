## Cluster the spirals data set.
library("kernlab", lib.loc="~/R/win-library/3.2")
library("e1071")
#data(spirals)
#write.csv(filedatals,"C:/Temp/filedatals.csv")
filedatals  = read.csv("C:/Temp/filedatals.csv",header = T,stringsAsFactors = F)

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
sc <- specc(sp,kernel = "laplacedot",kpar = list(sigma= 0.4), centers=2)
df = data.frame(sp,cl = factor(sc@.Data))

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
  sp = (cbind(filedatals$sumCh,filedatals$FCS))
  dfa = data.frame(sp)

  
subformcol <-function(df,v){
  
  for(i in 1:dim(df)[2]){
   
    df[,i] <- df[,i] + v[i] 
  }
  
  df
  
}

divformcol <-function(df,v){
  
  for(i in 1:dim(df)[2]){
    
    df[,i] <- df[,i]/v[i] 
  }
  
  df
  
}
  
cla = NULL
st = Sys.time()
for(i in 1:5){
  
  filedatalsSample = filedatals[sample(nrow(filedatals),1e3),]
  df = cbind(filedatalsSample$sumCh,filedatalsSample$FCS)
  sp = scale(df)
  mes = colMeans(df)
  mds = apply(df,2,sd)
  
  dff = subformcol(dfa,mes)
  dff = subformcol(dff,mds)
  sc <- specc(sp,kernel = "laplacedot",kpar = list(sigma= 0.4), centers = 2 )
  df = data.frame(sp,cl = factor(sc@.Data))
  svm_model <- svm(cl ~ ., data = df)
  df = data.frame(sp)
  pre = predict(svm_model,df)
  me = with(df,aggregate(X1,list(pre),mean))[,2]
  men = which(max(me) == me)
  if( men ==  2){
    pre = predict(svm_model,dff)  
    pre = ifelse(pre == 1,2,1)  
  }
  cla = rbind(cla,pre)
  
}

dt = Sys.time() - st
dt

cl  = ifelse(sapply(cla,mean) > 1.5,2,1)

plot(dfa$X1,dfa$X2,col = cla[2,],pch = 19 ,cex = 0.3,ylim = c(-1,7))

plot(dff$X1,dff$X2,col = pre,pch = 19 ,cex = 0.3 )
