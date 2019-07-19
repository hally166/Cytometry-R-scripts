# The idea here is to be able to reuse simple code to analyse your flow cytometry data. Various bits have been liberated from elsewhere; from documentation, from pipelines, and from graft.
# Christopher Hall, Wellcome Sanger Institute

#install packages from BioConductor
install.packages("BiocManager")
BiocManager::install("flowCore")
library(flowCore)

#or from Rbase
install.packages("shiny")
library(shiny)

#flow cytometry packages
BiocManager::install("flowCore")
BiocManager::install("flowViz")
BiocManager::install("ggcyto")
BiocManager::install("flowWorkspace")
library(ggcyto)
library(flowCore)
library(flowViz)
library(flowWorkspace)

#packages for data cleaning
BiocManager::install("flowAI")
BiocManager::install("flowClean")
install.packages("devtools")
devtools::install_github("jmeskas/flowCut")
library(flowAI)
library(flowClean)
library(flowCut)

#load a fingle fcs
myfile <- "C:/FCSfiles/8peak400v.fcs"
fcsfile <- read.FCS(myfile)

#load many fcs files into a flow set
myfiles <- list.files(path="C:/FCSfiles/", pattern=".fcs$")
fs <- read.flowSet(myfiles, path="C:/FCSfiles/") # be careful about using path here, it seems unnecessary, and probably is on Linux, but not windows

#compensation
#single file
fcsfile_comp <-compensate(fcsfile, spillover(fcsfile)[[1]])
#flowset
comp <-fsApply(fs,function(x)spillover(x)[[1]], simplify=FALSE) # different cytometers require different [[1]]
fs_comp <-compensate(fs, comp)

#transform data
#single file
tf <- estimateLogicle(fcsfile_comp, colnames(fcsfile_comp[,7:24]))
fcsfile_trans <- transform(fcsfile_comp, trans)
#flow set
colnames(fs_comp) # to ID your detectors
tf <- estimateLogicle(fs_comp[[1]], channels = colnames(fs_comp[[1]][,7:24])) # change colnames for flowset
fs_trans <- transform(fs_comp, tf)
fs_trans

#explore data
#single file
exprs(fcsfile_trans)[1:10,]
summary(fcsfile_trans)
str(keyword(fcsfile_trans))
#flowSet
summary(fs_trans[[1]])
head(fsApply(fs_trans, nrow))
fsApply(fs_trans, each_col, median)

#check consistency and labels
keyword(fs,c("$P1V", "$P2V", "$P3V", "$P4V", "$P5V", "$P6V", "$P7V", "$P8V", "$P9V", "$P10V", "$P11V", "$P12V", "$P13V", "$P14V", "$P15V", "$P16V")) # check for labels
keyword(fs,c("$P1N", "$P2N", "$P3N", "$P4N", "$P5N", "$P6N", "$P7N", "$P8N", "$P9N", "$P10N", "$P11N", "$P12N", "$P13N", "$P14N", "$P15N", "$P16N")) # check paramaters
param<-c("$P1V", "$P2V", "$P3V", "$P4V", "$P5V", "$P6V", "$P7V", "$P8V", "$P9V", "$P10V", "$P11V", "$P12V", "$P13V", "$P14V", "$P15V", "$P16V") # 
cyt<-keyword(fs,"$P1V")
#loop through and check them all
for (check in param) {
  cyt<-keyword(fs,check)
  print(paste0(check, " ", all(sapply(cyt, FUN = identical, cyt[[1]]))))
}

#cleaning flowset
#flowCut
QC_flowCut<-fsApply(fs_trans, function(x)flowCut(x))
QC_flowCut # see data
QC_flowCut[[1]]$frame # access flowframe
plot(fs_trans[[1]], c("610/20 (561)-A", "450/50 (355)-A"))
plot(QC_flowCut[[1]], c("610/20 (561)-A", "450/50 (355)-A"))
# convert to flowset
files<-list() 
for(x in 1:length(resQC)){
  files<-append(files,list(resQC[[x]]$frame))
}
fs_flowCut<-flowSet(files)
#flowClean - never got it to work
fs_flowClean<-fsApply(fs_trans, function(x)clean(x,vectMarkers=c(7:24),filePrefixWithDir="sample_out", ext="fcs", diagnostic=TRUE))
#flowAI
fs_flowAI <- flow_auto_qc(fs_trans)
fs_flowAI <- flow_auto_qc(fs_trans, remove_from = c("FR", "FS")) # FR= flow rate, FS = signal, FR = dynamic range.  FR can be problamatic
head(fsApply(fs_trans, nrow))
head(fsApply(fs_flowAI, nrow))
plot(fs_trans[[1]], c("610/20 (561)-A", "450/50 (355)-A"))
plot(fs_flowAI[[1]], c("610/20 (561)-A", "450/50 (355)-A"))

#plotting
library(ggcyto)
autoplot(fs_trans[[1]], x="610/20 (561)-A", y="450/50 (355)-A", bins = 256)
autoplot(fs_trans, x="610/20 (561)-A", y="450/50 (355)-A", bins = 256)
p <- ggcyto(fs_trans, aes(x = "610/20 (561)-A", y =  "450/50 (355)-A"))
p <- p + geom_hex(bins = 128)
p

#gating
library(flowWorkspace)
gs <- GatingSet(fs_trans)
rg1 <- rectangleGate("FSC-H"=c(100000, Inf), filterId="NonDebris")
add(gs, rg1, parent = "root")
getNodes(gs)
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'SSC-H', "NonDebris", bins = 256)

rg2 <- rectangleGate("FSC-H"=c(100000, 150000),"FSC-W"=c(50000, 75000))
add(gs, rg2, parent = "NonDebris", name = "singlets")
getNodes(gs)
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'FSC-W', "singlets", bins = 256)

p <- ggcyto(fs_trans, aes(x = "FSC-H", y = 'FSC-W'))+ geom_hex(bins = 256)
g <- getGate(gs, "singlets")
p <- p + geom_gate(g)
p

#look at the gated data
plot(gs)
getStats(gs)
getStats(gs, "singlets", "percent")

autoplot(gs[[1]])
fs_singlets <- getData(gs, "/NonDebris/singlets")
fsApply(fs_singlets, each_col, median)

#automatic gating
library(flowCore)
library(flowWorkspace)
library(openCyto)

#gate the main population of events
gs <- GatingSet(fs_trans)
thisData <- getData(gs)
nonDebris_gate <- fsApply(thisData, function(fr) openCyto:::.flowClust.2d(fr, channels = c("FSC-H","SSC-H")))
add(gs, nonDebris_gate, parent = "root", name = "nonDebris")
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'SSC-H', "nonDebris", bins = 256)

#gate the singlets
thisData <- getData(gs, "nonDebris") #get parent data
singlet_gate <- fsApply(thisData, function(fr) openCyto:::.singletGate(fr, channels =c("FSC-H", "FSC-W")))
add(gs, singlet_gate, parent = "nonDebris", name = "singlets")
recompute(gs)
autoplot(gs,x = 'FSC-H', y = 'FSC-W', "singlets", bins = 256) 

#export data
library(gridExtra)
grid.arrange(as.ggplot(p1), as.ggplot(p2), nrow = 2)
g <- arrangeGrob(as.ggplot(p1), as.ggplot(p2), nrow = 2)
ggsave(file="plots.png", g)
write.csv(getPopStats(gs), "stats.csv")
write.csv(exprs(thisData[[1]]), "stats2.csv")

#clustering
#concatenating
concat_data<-fsApply(ds_fs_trans, function(x)rbind(exprs(x)))
#just data - I need a better way fo doing this
df1<-as.data.frame(exprs(ds_fs_trans[[1]]))
df2<-as.data.frame(exprs(ds_fs_trans[[2]]))
df1$source <- keyword(ds_fs_trans[[1]])$GUID.original
df2$source <- keyword(ds_fs_trans[[2]])$GUID.original
concat_data<-rbind(df1,df2)

#downsampling
#concatenated file
ds_concat<-concat_data[sample(1:nrow(concat_data), 3000, replace=FALSE),]
#single file
ds_trans <- fs_trans[[1]][sample(1:nrow(fs_trans[[1]]), 3000, replace=FALSE),]
#flowSet
ds_fs_trans <- fsApply(fs_trans, function(x)x[sample(1:nrow(x), 3000, replace=FALSE),])

#tSNE clustering
install.packages("Rtsne")
library(Rtsne)
install.packages("ggplot2")
set.seed(42) # Sets seed for reproducibility
tsne_out <- Rtsne(as.matrix(exprs(ds_fs_trans[[2]][,7:24])),perplexity=50) # Run TSNE, use check_duplicates = FALSE for similar data
plot(tsne_out$Y, col=exprs(ds_fs_trans[[2]][,7])) # Plot the result, colour by expression
#colour by source
plot(tsne_out$Y, col=1:length(concat_data$source))
#or
colors = rainbow(length(unique(concat_data$source)))
names(colors) = unique(concat_data$source)
plot(tsne_out$Y, main="tSNE", xlab="tSNE 1", ylab="tSNE 2", "cex.main"=2, "cex.lab"=1.5,col=colors)
legend(0,50,names(colors), text.col=rainbow(factor(colors)))

#flowSOM
BiocManager::install("FlowSOM")
library(FlowSOM)
out <- FlowSOM::ReadInput(fs_trans[[1]], transform = FALSE, scale = FALSE)
out <- FlowSOM::BuildSOM(out, colsToUse = 3:20)
out <- FlowSOM::BuildMST(out)
FlowSOM::PlotStars(out)

#export data
#gridextra
library(gridExtra)
grid.arrange(as.ggplot(p1), as.ggplot(p2), nrow = 2) # p1 and p2 are ggplot objects
g <- arrangeGrob(as.ggplot(p1), as.ggplot(p2), nrow = 2)
ggsave(file="plots.png", g)
write.csv(getPopStats(gs), "stats.csv")
write.csv(exprs(thisData[[1]]), "stats2.csv")
#simple png
to_plot <- paste("filename", ".png", sep = "")
png(file_plot, width = 600, height = 600)
plot(tsne_out)
dev.off()







