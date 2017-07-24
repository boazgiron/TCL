args<-commandArgs(TRUE)

getheadname <- function(x){(strsplit(basename(x),"_EDv")[[1]])[1]}


runondatabase(comline){
  
  files = dir(comline[1])
  
  fieldBur = NULL
  FILENAME = NULL
  for (fd in fielddir){
    
    fi = dir(fd,full.names = T)
    fi = fi[grep(".events.csv",fi)]
    re  = main(fi)
    FILENAME = c(FILENAME,basename(fi))
    fieldBur  = rbind(fieldBur, re[,1] )
  }
  
  colnames(fieldBu) =  re[,2]
  
  fieldBur <- data.frame(fieldBur)
  
  fieldBur$FILENAME <- FILENAME
  
  fieldBur <- within( fieldBur, FILE.NAME.ID <-unlist( lapply(FILENAME, getheadname ) ) ) 
  write.csv(fieldBur,"C:/Project/LeukoDx/LudaFacsValidation/fourCasesForValidation.csv")
  
  AcceptFiles = read.csv("L:\\Personal\\Boaz\\TCL\\accellix 17 MANUAL DATA FROM 201609_Accept.csv",header = T ,stringsAsFactors = F)
  AcceptFiles$FILE.NAME.ID <- unlist(lapply(AcceptFiles$FILE.NAME,getheadname))
  AcceptFiles  = AcceptFiles$FILE.NAME.ID[ AcceptFiles$FILE.NAME.ID %in% fieldBur$FILE.NAME.ID, ]
  me = merge(fieldBur,AcceptFiles,by  = "FILE.NAME.ID")
  
}


runondatabase(args[1])


