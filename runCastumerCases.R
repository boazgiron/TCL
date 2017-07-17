fielddir = c("L:/Scientific/Raw Data/Reader_Data/Accellix 19/TCL/20170710/TCL CID-320",
             "L:/Scientific/Raw Data/Reader_Data/Accellix 19/TCL/20170710/TCL CID-321",
             "L:/Scientific/Raw Data/Reader_Data/Accellix 19/TCL/20170711/TCL CID-318",
             "L:/Scientific/Raw Data/Reader_Data/Accellix 19/TCL/20170711/TCL CID-319")
fieldBu1 = NULL
FILENAME = NULL
for (fd in fielddir){
  
  fi = dir(fd,full.names = T)
  fi = fi[grep(".events.csv",fi)]
  re  = runAlgo_shortData(fi)
  FILENAME = c(FILENAME,basename(fi))
  fieldBu1  = rbind(fieldBu1, re[,1] )
}

colnames(fieldBu) =  re[,2]

#fieldBu <- fieldBur

fieldBur <- data.frame(fieldBur)

fieldBur$FILENAME <- FILENAME

write.csv(fieldBur,"C:/Project/LeukoDx/LudaFacsValidation/fourCasesForValidation.csv")


basename(files[3])


AronFilespath = "L:/Personal/Boaz/TCL Data/Aron170717"
files = dir(AronFilespath,full.names = T)
files = files[-3]
files[-3]

fieldBu1 = NULL
FILENAME = NULL
for (fd in files){
  
  #fi = dir(fd,full.names = T)
  #fi = fi[grep(".events.csv",fi)]
  re  = runAlgo_shortData(fd)
  FILENAME = c(FILENAME,basename(fd))
  fieldBu1  = rbind(fieldBu1, re[,1] )
}

colnames(fieldBu1) =  re[,2]
#re  = runAlgo_shortData(files[3])

fieldBu1 <- data.frame(fieldBu1)
fieldBu1$FILENAME <- FILENAME
write.csv(fieldBu1,"C:/Project/LeukoDx/LudaFacsValidation/Aron8cases170717.csv")

re  = runAlgo_shortData(files[1])

