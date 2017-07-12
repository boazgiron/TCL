#Run script---
re  = runAlgo_shortData(runFiles[4])
#filepath  = "C:/working/data/20170619_17-16-19_AX019_TCL_L003_C427_q1543_690_SP_EDv38TCLv16_nofft.events.csv"
#filepath1  = "C:/working/data/20170619_16-16-19_AX019_TCL_L003_C429_q1618_690_SP_EDv38TCLv16_nofft.events.csv"

filepath  = "L:/Personal/Boaz/TCL/WorkingOnBug/data/20170619_17-16-19_AX019_TCL_L003_C427_q1543_690_SP_EDv38TCLv16_nofft.events.csv"
filepath1  = "L:/Personal/Boaz/TCL/WorkingOnBug/data/20170619_16-16-19_AX019_TCL_L003_C429_q1618_690_SP_EDv38TCLv16_nofft.events.csv"


#----


pathd = "L:/Scientific/Raw Data/Reader_Data/TCL Validation/accellix17"
files = dir(pathd,full.names = T)

#keepbr = br
getn <- function(st){
  arst = strsplit(st,"/")[[1]]
  st1 = paste0((strsplit(arst[length(arst)],"[.]")[[1]])[1],".fcs")
  st1
}

sa = 1:length(files)#sample(1:length(files),50)


#sal = (1:length(files))[!((1:length(files)) %in% sa)]

#saa = sample(sa,40)

#C3293,C3377,C233,C3288,C131,C3459,C233
#fl = files[grep("C131",files)]
fl = "L:/Scientific/Raw Data/Reader_Data/TCL Validation/accellix17/20170119_03-06-17_AX017_GT0_L016_C3442_22172_cd3f_45pe_draq_650_EDv38TCLv16_nofft.events.csv"

fl = AcceptedrunFiles[9]
re = runAlgo_shortData(fl)

br = NULL
runfile = NULL
i = 0
for(fl in AcceptedrunFiles){
  i = i+1
  print(i)
  print(basename(fl))
  #undebug(runAlgo_shortData)
  re=runAlgo_shortData(fl)
  br = rbind(br,re[,1])
  runfile = c(runfile,basename(fl))
}

write.csv(br,"C:/Project/LeukoDx/LudaFacsValidation/Results12cases.csv")
colnames(br) = re[,2]
br1 <- br
#colnames(br1)

br1 = data.frame(br1,FILE.NAME = runfile)
br1$FILE.NAME <- as.character(br1$FILE.NAME)
br1 = data.frame(br1, FILE.NAME.ID <- unlist(lapply(br1$FILE.NAME,getheadname) ) )
br1$FILE.NAME.ID <- as.character(br1$FILE.NAME.ID)

merge_br1_AcceptFiles = merge(AcceptFiles,br1,by = "FILE.NAME.ID")
diff = merge_br1_AcceptFiles[,c("Cd45liveRatio","Cd3liveRatio","Cd3Cd45Ratio","CD4CD3Ratio","CD8CD3Ratio")]  - merge_br1_AcceptFiles[,c("Ratio45","Ratio3","Ratio3vs5","Ratio4","Ratio8" ) ]                                        
colnames(diff) <- c("d_45","d_3","d_3vd45","d_4","d_8")
merge_br1_AcceptFiles_  = data.frame( merge_br1_AcceptFiles, diff )

View(merge_br1_AcceptFiles_)



#notInuseCasefile = read.csv(file.choose(),header = T ,stringsAsFactors = F)

getheadname <- function(x){(strsplit(basename(x),"_EDv")[[1]])[1]}
p17 = "L:/Scientific/Raw Data/Reader_Data/TCL Validation/accellix17/"
f17  = dir(p17,full.names = T)
head_f17 = unlist(lapply(basename(f17),getheadname))
rf = c(
  "20160903_04-24-43_AX017_TCL_L002_C152_170_D0_CD4_5PERC_650_EDv38_nofft.events.csv",
  "20160903_04-59-44_AX017_TCL_L002_C154_170_D0_CD8_5PERC_650_EDv38_nofft.events.csv",
  "20160903_05-35-11_AX017_TCL_L002_C156_158_D0_CD4_POS_650_EDv38_nofft.events.csv",
  "20160903_06-07-41_AX017_TCL_L002_C158_170_D0_CD4_POS_650_EDv38_nofft.events.csv",
  "20160903_06-40-53_AX017_TCL_L002_C160_170_D0_CD4_POS_650_EDv38_nofft.events.csv",
  "20160903_07-16-09_AX017_TCL_L002_C162_170_D0_CD4_POS_650_EDv38_nofft.events.csv",
  "20150908_10-39-00_AX017_TCL_L002_C113_215D1_650_EDv38_nofft.events.csv",
  "20150908_13-02-25_AX017_TCL_L002_C193_215D1_LOW CD8_650_EDv38_nofft.events.csv",
  "20151028_15-35-20_AX017_TCL_L002_C134_rc_650_EDv38_nofft.events.csv",
  "20160903_08-43-36_AX017_TCL_L002_C172_158_D0_CD4_NEG_1M_650_EDv38_nofft.events.csv",
  "20150911_14-18-32_AX017_TCL_L002_C237_267_WO_PhenRed_650_EDv38_nofft.events.csv",
  "20160905_12-11-16_AX017_TCL_L002_C119_158_D3_650_EDv38_nofft.events.csv")
head_rf = unlist(lapply(rf,getheadname))
AcceptedrunFiles  = f17[ head_f17 %in%  head_rf ]
length(head_rf)
basename(AcceptedrunFiles)

#AcceptedrunFiles = AcceptedrunFiles[-6]
#prepare file list
if(FALSE){
  
  getheadname <- function(x){(strsplit(basename(x),"_EDv")[[1]])[1]}
  AcceptFiles = read.csv("L:\\Personal\\Boaz\\TCL\\accellix 17 MANUAL DATA FROM 201609_Accept.csv",header = T ,stringsAsFactors = F)
  
  unique(AcceptFiles$Type)
  
  AcceptFilesList =  unique(AcceptFiles$FILE.NAME)
  p17 = "L:/Scientific/Raw Data/Reader_Data/TCL Validation/accellix17/"
  AcceptFilesList <- paste0(p17,AcceptFilesList[grep("AX017",AcceptFilesList)]) 
  #write.csv(AcceptFilesList,"L:/Personal/Boaz/TCL/fileList110.csv")
  f17  = dir(p17,full.names = T)
  
  head_AcceptFilesList = basename(AcceptFilesList)
  head_AcceptFilesList = unlist(lapply(head_AcceptFilesList,getheadname))
  head_f17 = unlist(lapply(basename(f17),getheadname))
  AcceptedrunFiles  = f17[ head_f17 %in%  head_AcceptFilesList ]
  length(head_AcceptFilesList)
}


AcceptFiles$FILE.NAME.ID <- unlist(lapply(AcceptFiles$FILE.NAME,getheadname))
colnames(AcceptFiles)
AcceptFiles$Ratio45 <- with(AcceptFiles,round(100*X.CD45.LIVE/X.CD45.TOT,2))
AcceptFiles$Ratio3 <- with(AcceptFiles,round(100*X.CD3.LIVE/X.CD3.TOTAL,2))
AcceptFiles$Ratio3vs5 <- with(AcceptFiles,round(100*X.CD3.LIVE/X.CD45.LIVE,2))
AcceptFiles$Ratio4 <- with(AcceptFiles,round(100*X.CD4.LIVE/X.CD3.LIVE,2))
AcceptFiles$Ratio8 <- with(AcceptFiles,round(100*X.CD8.LIVE/X.CD3.LIVE,2))
#write.csv(AcceptFiles,"C:/Project/LeukoDx/LudaFacsValidation/result100717.csv")
AcceptFiles = read.csv("C:/Project/LeukoDx/LudaFacsValidation/result100717.csv")
AcceptedrunFiles  = f17[ head_f17 %in% AcceptFiles$FILE.NAME.ID ]

head(AcceptFiles)
colnames(br) = re[,2]
br1 = data.frame(br)
br1$FILE.NAME <- unlist(lapply(runfile,getn))
sel  = results$FILE.NAME %in% br1$FILE.NAME
br2 = merge(br1,results[sel,],by  = "FILE.NAME")
#write.csv(br2,"C:/Project/LeukoDx/LudaFacsValidation/Debug/firstAll.csv")

br2[br2$FILE.NAME %in% verifiedFile, ]

#problem with 235 (17)

colnames(br) = re[,2]
View(br)

results = read.csv("C:/Project/LeukoDx/LudaFacsValidation/Debug/results.csv",header = T ,stringsAsFactors = F)
colnames(results)[dim(results)[2]] = "Cartrige"

br1 = data.frame(br)
sum(results$Cartrige %in% br1$Cartrige)

results$Cartrige[results$Cartrige %in% br1$Cartrige]

results = results[results$Cartrige %in% br1$Cartrige, ]

br2 = merge(br1,results,by  = "Cartrige")
dim(br2)
View(br2)


verifiedFile = c(
  "20150912_09-31-50_AX017_TCL_L002_C227_215D5_WO_PhenRed_650_EDv38_nofft.fcs",
  "20150912_10-12-14_AX017_TCL_L002_C229_215D5_WO_PhenRed_650_EDv38_nofft.fcs",
  "20150912_11-02-38_AX017_TCL_L002_C231_215D5_WO_PhenRed_650_EDv38_nofft.fcs",
  
  "20150908_11-12-55_AX017_TCL_L002_C191_215D1_LOW CD8_650_EDv38_nofft.fcs",
  "20150908_13-02-25_AX017_TCL_L002_C193_215D1_LOW CD8_650_EDv38_nofft.fcs",
  "20150908_13-36-27_AX017_TCL_L002_C195_215D1_LOW CD8_650_EDv38_nofft.fcs",
  "20150908_14-25-35_AX017_TCL_L002_C197_267D1_LOW CD8_650_EDv38_nofft.fcs",
  "20150908_14-58-32_AX017_TCL_L002_C199_267D1_LOW CD8_650_EDv38_nofft.fcs",
  "20150908_15-35-07_AX017_TCL_L002_C248_267D1_LOW CD8_650_EDv38_nofft.fcs",
  
  "20150908_09-24-16_AX017_TCL_L002_C109_215D1_650_EDv38_nofft.fcs",
  "20150908_09-58-20_AX017_TCL_L002_C111_215D1_650_EDv38_nofft.fcs",
  "20150908_10-39-00_AX017_TCL_L002_C113_215D1_650_EDv38_nofft.fcs"
)


verifiedcart = c(227,229,231,191,193,195,197,199,248,109,111,113)

write.csv(br2,"C:/Project/LeukoDx/LudaFacsValidation/Debug/first80.csv")


brv = br1[br1$Cartrige %in% verifiedcart, ]
results1 = results[results$FILE.NAME %in% verifiedFile,] 
write.csv(results1,"C:/Project/LeukoDx/LudaFacsValidation/Debug/result1.csv")
results1$Cartrige <- results1$CartNum
View(results1)

vmer = merge(brv,results1,by = "Cartrige")
View(vmer)
write.csv(vmer,"C:/Project/LeukoDx/LudaFacsValidation/Debug/vmer.csv")
