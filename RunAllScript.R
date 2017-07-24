#parameter
#TCL Script location
#cmd /c Rscript.exe "C:\Temp\example.R" 
args<-commandArgs(TRUE)
batchLocation = args[1]
library("tcltk", lib.loc="C:/Program Files/R/R-3.2.1/library")
getheadname <- function(x){(strsplit(basename(x),"_EDv")[[1]])[1]}

#FileTocompare = "L:\\Personal\\Boaz\\TCL\\accellix 17 MANUAL DATA FROM 201609_Accept.csv"
ResultLocation = paste0(batchLocation,"/Results.csv")

plotdir = paste0(,"/plots")

if(!dir.exists(plotdir)){
  
  dir.create(plotdir)
  
}
ResultLocation = "C:/Project/LeukoDx/LudaFacsValidation/fourCasesForValidation.csv"
#ResultPlotsLocation = "C:/Project/LeukoDx/LudaFacsValidation/Debug240717"


dirpath  = tk_choose.dir(default = "", caption = "Select CSV files directory")
files = dir( dirpath )
#ResultLocationfile  = tk_choose.files(default = "", caption = "Select Compare file")
#ResultPlotsLocation  = tk_choose.dir(default = "", caption = "Select Plots Directory")

fieldBur = NULL
FILENAME = NULL

for (fd in fielddir){
  
  fi = dir(fd,full.names = T)
  fi = fi[grep(".events.csv",fi)]
  re  = main(fi,plotdir)
  FILENAME = c(FILENAME,basename(fi))
  fieldBur  = rbind(fieldBur, re[,1] )
}

colnames(fieldBu) =  re[,2]

fieldBur <- data.frame(fieldBur)

fieldBur$FILENAME <- FILENAME

fieldBur <- within( fieldBur, FILE.NAME.ID <-unlist( lapply(FILENAME, getheadname ) ) ) 
write.csv(fieldBur,ResultLocation)


if( !is.null(FileTocompare) ){
  
  AcceptFiles = read.csv(FileTocompare,header = T ,stringsAsFactors = F)
  AcceptFiles$FILE.NAME.ID <- unlist(lapply(AcceptFiles$FILE.NAME,getheadname))
  AcceptFiles  = AcceptFiles$FILE.NAME.ID[ AcceptFiles$FILE.NAME.ID %in% fieldBur$FILE.NAME.ID, ]
  me = merge(fieldBur,AcceptFiles,by  = "FILE.NAME.ID")
  
}  






