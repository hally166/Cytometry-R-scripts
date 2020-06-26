#Load the packages
library(flowCore)
library(flowAI)
library(ggcyto)

#How to get help
??flowCore

#Load a single file
myfile <- "C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZU_files/0004.FCS"
fcsfile <- read.FCS(myfile)
fcsfile
names(fcsfile)
exprs(fcsfile)
each_col(fcsfile, median)
keyword(fcsfile)

#Compensation
spillover(fcsfile)
fcsfile_comp <-compensate(fcsfile, spillover(fcsfile)$SPILL)
fcsfile_comp

#Cleaning
fcsfile_comp_clean <- flow_auto_qc(fcsfile_comp)
fcsfile_comp_clean
keyword(fcsfile_comp_clean) <- keyword(fcsfile)
fcsfile_comp_clean
??flowAI

#Transformation
trans <- estimateLogicle(fcsfile_comp_clean, colnames(fcsfile_comp_clean[,3:10]))
fcsfile_comp_clean_trans <- transform(fcsfile_comp_clean, trans)

#Visualise the results
??ggcyto
autoplot(fcsfile_comp_clean)
autoplot(fcsfile_comp_clean_trans)
autoplot(fcsfile_comp_clean_trans, x="PE-Cy7-A", y="PerCP-Cy5-5-A", bins = 256)
autoplot(fcsfile_comp_clean_trans, x="Time", y="FSC-A", bins = 128)
autoplot(transform(fcsfile_comp,trans), x="Time", y="FSC-A", bins = 128)

#In a flowSet
myfiles <- list.files(path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZU_files", pattern=".FCS$")
fs <- flowCore::read.flowSet(myfiles, path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZU_files/")
fs
fs[[1]]
spillover(fs[[1]])
fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)
fs_comp_clean <- flow_auto_qc(fs_comp)
trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[[1]][,3:10]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)

#fsApply
??fsApply
fsApply(fs,each_col,median)