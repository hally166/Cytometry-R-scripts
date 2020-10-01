#Load packages
library(flowCore)
library(flowCut)
library(flowAI)
library(ggcyto)
library(flowWorkspace)
library(openCyto)
setwd("C:\\Users\\chall\\OneDrive\\EMBL\\Test FCSfiles\\Attune\\030220")

#load and clean files usinf flowAI
files <- list.files("C:\\Users\\chall\\OneDrive\\EMBL\\Test FCSfiles\\Attune\\030220", pattern=".fcs", full.names = TRUE,ignore.case = TRUE)
fs <- read.flowSet(files, emptyValue = TRUE) # add this if you are using Attune data: emptyValue = FALSE)
pd <- pData(fs) #this bit is because flowAI currently messed with the metadata in a bad way
resQCfs <- flowAI::flow_auto_qc(fs) 
pData(resQCfs) <- pd

#What has flowAI done
QC_report<-cbind(as.data.frame(fsApply(fs, nrow)),fsApply(resQCfs, nrow)[,1])
names(QC_report)<-c("Pre_Clean","Post_Clean")
QC_report<-cbind(QC_report, 100-round(QC_report$Post_Clean/QC_report$Pre_Clean*100))
names(QC_report)<-c("Pre_Clean","Post_Clean", "% removed")
QC_report

#compensate files
comp <- fsApply(resQCfs, function(x) spillover(x)[[3]], simplify=FALSE) #change the [[3]] dependent on the spillover keyword
fs_comp <- compensate(resQCfs, comp)

#transform files - this function chooses all the paramaters, except FSC SSC and Time, you can add more if you wish
autovect_verbose<- function(ff){
  c<- data.frame(ff@parameters@data)
  d<- grep("FSC|SSC|Time", c$name, invert = TRUE, value = TRUE)
  return(unname(d))
}
x<-autovect_verbose(fs_comp[[1]])
tf<-estimateLogicle(fs_comp[[1]], channels =  x)
fs_trans<-transform(fs_comp,tf)

#analyse data using flowWorkspace
gs<-GatingSet(fs_trans)
thisData<-gs_pop_get_data(gs)
nonDebris<-fsApply(thisData, function(fr)openCyto:::.mindensity(fr,channels = "FSC-A", max=1000000)) #change or remove max to help it find the population
gs_pop_add(gs,nonDebris, parent="root", name="nonDebris")
recompute(gs)

thisData<-gs_pop_get_data(gs)
nonDebris<-fsApply(thisData, function(fr)openCyto:::.singletGate(fr,channels = c("FSC-H","FSC-W")))
gs_pop_add(gs,nonDebris, parent="nonDebris", name="singlets")
recompute(gs)

thisData<-gs_pop_get_data(gs)
nonDebris<-fsApply(thisData, function(fr)openCyto:::.mindensity(fr,channels = "BL1-A"))
gs_pop_add(gs,nonDebris, parent="singlets", name="GFP")
recompute(gs)

f1<-autoplot(gs,x="FSC-A", y="SSC-A", "nonDebris")
f2<-autoplot(gs,x="FSC-H", y="FSC-W", "singlets")
f3<-autoplot(gs,x="BL1-A", y="RL1-A", "GFP")

gs_pop_get_count_fast(gs)
a<-gs_pop_get_count_fast(gs,statistic = c("count", "freq"),subpopulations = "/nonDebris/singlets/GFP")
a<-cbind(a,percentage = round((a$Count/a$ParentCount)*100))
a$Parent<-NULL
a$ParentCount<-NULL

write.csv(a,"shortstats.csv")
write.csv(QC_report,"QC report.csv")
ggsave(plot = f1, width = 15, height = 15, dpi = 300, filename = "root.png")
ggsave(plot = f2, width = 15, height = 15, dpi = 300, filename = "singlets.png")
ggsave(plot = f3, width = 15, height = 15, dpi = 300, filename = "GFP.png")
