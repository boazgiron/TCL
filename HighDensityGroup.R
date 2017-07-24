x= c(rep(0,10),rep(17,10),rep(0,10),rep(15,10) )

HighDensityGroup <- function(x,n){
  
  #x = d$y
  #step1 = 1e-4
  st = max(x)
  end = min(x)
  st = ( min(x) - max(x) )/n
  se = seq(st,end,st)
  r = NULL
  for(s in se ){
    
    bb = (x >= s)
    z = T
    i = 0
    
    for(b in bb){
      if( b ){
        if( !z ){
          i =i+1
          z = T
        }
      }else{
        
        z = F
      
      }
    }
    
    r = rbind( r, c(s,i ) )
  }
  
  r
}

