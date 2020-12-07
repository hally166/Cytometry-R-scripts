#Load packages
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowAI)
library(gridExtra)

#Load data
myfiles <- list.files(path="C:/Users/hallc/Downloads/FlowRepository_FR-FCM-ZZZV_files", pattern = ".fcs", ignore.case = TRUE)
fs <- read.flowSet(myfiles[1:10], path="C:/Users/hallc/Downloads/FlowRepository_FR-FCM-ZZZV_files", alter.names=TRUE)
matrix<-spillover(fs[[1]])$SPILL
colnames(matrix)<-c("X.FITC.A.", "X.Pacific.Blue.A.", "X.Alexa.680.A.", "X.APC.A.", "X.PE.Cy7.A.", "X.PE.Cy55.A.", "X.PE.Tx.RD.A.", "X.PE.Green.laser.A.")
fs_comp <-compensate(fs,matrix)
fs_comp_clean <- flow_auto_qc(fs_comp)
trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[,c(4,6:12)]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)

#create the empty gating set
auto_gs<-GatingSet(fs_comp_clean_trans)

#cell gate
fs_data<- gs_pop_get_data(auto_gs)
noneDebris_gate<- fsApply(fs_data, function(fr) openCyto:::.flowClust.2d(fr, channels= c("FSC.A","SSC.A")))
gs_pop_add(auto_gs, noneDebris_gate, parent = "root", name="noneDebris_gate")
recompute(auto_gs)
autoplot(auto_gs, x="FSC.A", y="SSC.A", "noneDebris_gate", bins=256)

#Singlet gate
fs_data <- gs_pop_get_data(auto_gs, "noneDebris_gate") #get parent data
singlet_gate <- fsApply(fs_data, function(fr) openCyto:::.singletGate(fr, channels =c("FSC.A", "FSC.H")))
gs_pop_add(auto_gs, singlet_gate, parent = "noneDebris_gate", name = "singlets")
recompute(auto_gs)
autoplot(auto_gs, x = 'FSC.A', y = 'FSC.H', "singlets", bins = 256)

#Quad gate
fs_data <- gs_pop_get_data(auto_gs, "singlets") #get parent data
BGquad_gate <- fsApply(fs_data, function(fr) openCyto:::.quadGate.seq(fr, gFunc="mindensity", min=c(3,3), channels =c('X.FITC.A.', 'X.PE.Tx.RD.A.')))
gs_pop_add(auto_gs, BGquad_gate, parent = "singlets", names = c("CD3+CD4-", "CD3+CD4+", "CD3-CD4+", "CD3-CD4-"))
recompute(auto_gs)
gs_get_pop_paths(auto_gs[[1]])
plot(auto_gs)
autoplot(auto_gs, x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.', gs_get_pop_paths(auto_gs)[4:7], bins = 256)
#Let's take a look at FR-FCM-ZZZV

fs_comp_clean_trans
fsApply(fs_comp_clean_trans, keyword, '$FIL')
colnames(fs_comp_clean_trans)
plot(auto_gs)
autoplot(auto_gs, x="FSC.A", y="SSC.A", "noneDebris_gate", bins=256)
autoplot(auto_gs, x = 'FSC.A', y = 'FSC.H', "singlets", bins = 256)
autoplot(auto_gs, x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.', gs_get_pop_paths(auto_gs)[4:7], bins = 256)

autoplot(auto_gs[[1]], x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.', gs_get_pop_paths(auto_gs)[4:7], bins = 256)

#fix the plot

p<-ggcyto(auto_gs[[1]],aes(x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.'), subset="singlets")
p<- p + geom_hex(bins=256)
p<- p + geom_gate(gs_get_pop_paths(auto_gs)[4:7])
p<- p + geom_stats(gs_get_pop_paths(auto_gs)[4:7],type = c("percent","gate_name"), adjust = c(0.3,0.8))
myPars <- ggcyto_par_set(limits = list(y = c(3,5), x = c(3,5)))
p <- p  + myPars
p <-p+ labs(x="CD4", y="CD3")
p
?ggcyto

#Duplicate the Autoplot - remove the [[1]]

#Faceting

pData(auto_gs)
csvfile<-read.csv("C:/Users/hallc/Downloads/FlowRepository_FR-FCM-ZZZV_files/attachments/Challenge3Metadata.csv")
newpData<-merge(pData(auto_gs),csvfile, by.x =  'name', by.y = 'FCS.File')
rownames(newpData)<- newpData$name
pData(auto_gs)<-newpData
pData(auto_gs)

p2 <- ggcyto(auto_gs, aes(x = `CD8`, y = `CD4`), subset="singlets") 
p2 <- p2 + geom_hex(bins=256)
myPars <- ggcyto_par_set(limits = list(y = c(3,4.5), x = c(2,4)))
p2 <- p2 + myPars
p2 + facet_grid(Sample.Treatment~Sample.Description)

#Backgating

p3 <- ggcyto(auto_gs, aes(x = `CD3`, y = `CD4`), subset = "root")
p3 <- p3 + geom_hex(bins=128)
myPars <- ggcyto_par_set(limits = list(y = c(3,5), x = c(3,5)))
p3 <- p3  + myPars
p3 + geom_overlay(data = "CD3+CD4+",size = 0.01, alpha = 0.05, color = "orange")

#Exporting image

ggsave("plot1.png", p2)
?ggsave

#Playing with grid.arrange

myplot<- function(x){
  p<-ggcyto(auto_gs[[x]],aes(x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.'), subset="singlets")
  p<- p + geom_hex(bins=256)
  p<- p + geom_gate(gs_get_pop_paths(auto_gs)[4:7]) 
  p<- p + geom_stats(gs_get_pop_paths(auto_gs)[4:7],type = c("percent","gate_name"), adjust = c(0.3,0.8))
  myPars <- ggcyto_par_set(limits = list(y = c(3,5), x = c(3,5)))
  p <- p  + myPars
  p <-p+ labs(x="CD4", y="CD3")
  p<-as.ggplot(p)
}

grid.arrange(
  grobs = myplot(1:4),
  widths = c(2, 2),
  layout_matrix = rbind(c(1,2,3,4))
)

grid.arrange(myplot(5:10))





