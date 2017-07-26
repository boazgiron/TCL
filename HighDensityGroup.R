x= c(rep(0,5),rep(5,3),rep(0,2),rep(17,10),rep(0,10),rep(15,4),rep(10,2) ,rep(12,4) )

HighDensityGroup <- function(x,n){
  
  #x = d$y
  #n =100
  st = max(x)
  end = min(x)
  
  if( n > 0 ){
    
    stp = ( min(x) - max(x) )/n
    
  }else{
    
    stp = n
  }
  
  se = seq(st,end,stp)
  
  r = NULL
  for(s in se ){
    
    bb = (x >= s)
    z = 1
    i = 0
    r1 = NULL
    
    for(b in bb){
      if( b ){
        if( z == 0  ){
          i =i+1
          z = 1
        }
      }else{
        
        z = 0
      
      }
      
      r1 = c(r1,i*z)
    }
    
    r = rbind( r, c(s,i,r1) )
  }
  
  r
}


ClustringBylevel <- function(x,s){
  
  
  bb = (x >= s)
  z = 1
  i = 0
  r =NULL
  
  for(b in bb){
    
    if( b ){
      if( z == 0 ){
        i =i+1
        z = 1
      }
    }else{
      
      z = 0
      
    }
    
    r = c(r,i*z)
  }
  
  r
}  



upfromMin <- function(x,vl,n ,stNumOfGroup){
  
  # vl = vlm
  # n  = 20
  # stNumOfGroup = 2
  # x = df$y

  end = max(x)
  st = vl
  
  if( n > 0 ){
    
    stp = ( max(x) - st )/n
    
  }else{
    
    stp = n
  }
  
  se = seq(st,end,stp)
  
  j = 1
  le = length(se)
  i = stNumOfGroup
  
  while( (j <=  le ) & ( i == stNumOfGroup  ) ){
    
    s  = se[j]
    j=j+1
    bb = (x >= s)
    z = 1
    i = 0
    r1 = NULL
    
    for(b in bb){
      if( b ){
        if( z == 0  ){
          i =i+1
          z = 1
        }
      }else{
        
        z = 0
        
      }
      
    }
  
  }    
  
  se[max(j-2,1)]

}


maxlocation <- function(fs){
  
  fs$x[which( max(fs$y) == fs$y )]
  
}


#Check

df = data.frame(x = h_vpArea6Area8$mids,y  = h_vpArea6Area8$counts)
df = read.csv("C:/Project/LeukoDx/R/TCL/h_vpArea6Area8_Example.csv")



#write.csv(df,"C:/Project/LeukoDx/R/TCL/h_vpArea6Area8_Example.csv")

plot(df$x,df$y,"l")
ro = HighDensityGroup(df$y,100)

rol  = ro[,3:dim(ro)[2]]

for(i in 1:dim(rol)[1]){
  
  getmaxlocation(rol)
  
}

sapply(t(rol[1:2,]),getmaxlocation)

cl = rol[1,]

bu = NULL
for(i in 1:dim(rol)[1]){#1:3){#
  
  getmaxlocation(rol[1,])
  
  bu = rbind(bu , getmaxlocation(rol[i,]) )

}



plot(bu$h.Group.1,bu$h.x)

getmaxlocation <- function(cl){

  #cl = rol[1,]
  hs = aggregate(df$y,list(cl),max )
  dfs = split(df,list(cl))
  yy = unlist(lapply(dfs,max))
  xx = unlist(lapply(dfs,maxlocation))
  ddf = data.frame(h  = hs,x = xx)
}

plot(df$x,df$y,"l")
vl = min(ro[ro[,2] == 2,1])
points(df$x,max(df$y)*ClustringBylevel(df$y,vl)/2,pch = 19,type = "l",col =2)
vlu =  upfromMin(df$y,vl,20,2)
points(df$x,max(df$y)*ClustringBylevel(df$y,vlu)/2.1,pch = 19,type = "l",col =2)

#plot(ClustringBylevel(df$y,vl))

vlm = max(ro[ro[,2] == 2,1])
#plot(ClustringBylevel(df$y,vl3))


points(df$x,max(df$y)*ClustringBylevel(df$y,(vl + vlu)/2)/2,pch = 19,type = "l", col = 3)
vl = min(ro[ro[,2] == 2,1])
cl = ClustringBylevel(df$y,vl)
ar = aggregate(df$y,list(cl),max)

dfs = split(df,cl)

sapply(dfs,maxlocation)

ar = aggregate(df,list(cl),maxlocation)

abline(h = ar$x[2:3])



paste0( "C:/YossiResults/CheckResult",gsub("-","",gsub(":","",Sys.time())),".csv")

plot(ro)


lo = ro[ind,3:dim(ro)[2]]
tlo  = t(lo)

gy <- function(y){

sapply(tlo,aggregate(x,))

View(ro)
  
plot(ro[,1],ro[,2],"b")

ro[diff(ro[,2]) != 0,1]

ind = which(diff(ro[,2]) < 0 )

le1 = ro[ind ,2]
ro[ind,3:dim(ro)[2]]


ClustringBylevel(x,le1[1])

df = data.frame(x,cl = ClustringBylevel(x,le1[1])) 






with(df,aggregate( x, list(cl), length ))  
with(df,aggregate(x, list(cl), mean ))  
with(df,aggregate(x, list(cl), sd ))  
with(df,aggregate(1:length(x), list(cl), range ))  
  
  
with(df,aggregate( x, list(cl), length ))  

  
  
