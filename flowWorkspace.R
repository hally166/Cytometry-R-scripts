library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowAI)
library(gridExtra)

#manual
#Load data
myfiles <- list.files(path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZU_files", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZU_files")
fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)
fs_comp_clean <- flow_auto_qc(fs_comp)
trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[,3:7]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)

#Visualise file
fs_comp_clean_trans[[1]]
autoplot(fs_comp_clean_trans[[1]])

#Basic gating
ggcyto(fs_comp_clean_trans[[1]], aes(x="FSC-A", y="SSC-A"))+geom_hex(bins=256)

#create the empty gating set
gs<-GatingSet(fs_comp_clean_trans)

#Cells - FSC SSC
rg1<-rectangleGate("FSC-A"=c(15000, Inf), filterId = "NoneDebris")
gs_pop_add(gs, rg1, parent="root")
recompute(gs)
gs_get_pop_paths(gs)
ggcyto(fs_comp_clean_trans[[1]], aes(x="FSC-A", y="SSC-A"))+geom_hex(bins=256)+geom_gate(gs_pop_get_gate(gs, "NoneDebris"))
gs_pop_get_stats(gs)

#Singlet gating
ggcyto(fs_comp_clean_trans[[1]], aes(x = "FSC-H", y = 'FSC-W'))+ geom_hex(bins = 256)
rg2 <- rectangleGate("FSC-H"=c(3.6, 4.2),"FSC-W"=c(50, 240))
gs_pop_add(gs, rg2, parent = "NoneDebris", name = "singlets")
gs_get_pop_paths(gs)
recompute(gs)
ggcyto(fs_comp_clean_trans, aes(x = "FSC-H", y = 'FSC-W'))+ geom_hex(bins = 256)+ geom_gate(gs_pop_get_gate(gs, "singlets"))

#exploring the gatingset
plot(gs)
gs_pop_get_stats(gs)
gs_pop_get_stats(gs, "NoneDebris", "percent")


#automatic
#Load data
myfiles <- list.files(path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZV_files", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles[1:2], path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZV_files", alter.names=TRUE)
fs_comp <-compensate(fs,spillover(fs[[1]])$SPILL)
#fix compensation matrix
matrix<-spillover(fs[[1]])$SPILL
matrix
colnames(matrix)<-c("X.FITC.A.", "X.Pacific.Blue.A.", "X.Alexa.680.A.", "X.APC.A.", "X.PE.Cy7.A.", "X.PE.Cy55.A.", "X.PE.Tx.RD.A.", "X.PE.Green.laser.A.")
fs_comp <-compensate(fs,matrix)
#continue
fs_comp_clean <- flow_auto_qc(fs_comp)
trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[,c(4,6:12)]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)
autoplot(fs_comp_clean_trans[[1]])

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
gs_pop_add(auto_gs, BGquad_gate, parent = "singlets", names = c("1", "2", "3", "4"))
recompute(auto_gs)
gs_get_pop_paths(auto_gs[[1]])
plot(auto_gs)
autoplot(auto_gs, x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.', gs_get_pop_paths(auto_gs)[4:7], bins = 256)

#fix plot
p<-ggcyto(auto_gs[1:2],aes(x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.'), subset="singlets", arrange = FALSE)
p<- p + geom_hex(bins=256)
p<- p + geom_gate(gs_get_pop_paths(auto_gs)[4:7]) 
p<- p + geom_stats(gs_get_pop_paths(auto_gs)[4:7])
p<- p + theme(strip.text = element_text(size = 7))
myPars <- ggcyto_par_set(limits = list(y = c(3,5), x = c(3,5)))
p<- p  + myPars
p

#Removing stuff
gs_pop_remove(auto_gs, "singlets")

#statistics
gs_pop_get_stats(auto_gs)
gs_pop_get_stats(auto_gs, "noneDebris_gate", "percent")
gs_pop_get_stats(auto_gs, "noneDebris_gate", type = pop.MFI)

pop.quantiles <- function(fr){
  chnls <- colnames(fr)
  res <- matrixStats::colQuantiles(exprs(fr), probs = 0.75)
  names(res) <- chnls
  res
}
gs_pop_get_stats(auto_gs, gs_get_pop_paths(auto_gs), type = pop.quantiles)

pop.mean <- function(fr){
  chnls <- colnames(fr)
  res <- colMeans(exprs(fr))
  names(res) <- chnls
  res
}
gs_pop_get_stats(auto_gs, gs_get_pop_paths(auto_gs), type = pop.mean)
