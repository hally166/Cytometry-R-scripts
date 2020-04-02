library(flowCore)
library(ggcyto)
library(flowClean)
library(flowCut)
library(flowAI)
library(gridExtra)

#artifical data
data(synPerturbed)

synPerturbed.flowclean <- clean(synPerturbed, vectMarkers=c(5:17), filePrefixWithDir="sampleName", ext="fcs")
synPerturbed.flowcut <- flowCut(synPerturbed, Channels=c(5:17) )
synPerturbed.flowai<-flowAI::flow_auto_qc(synPerturbed)

rg <- rectangleGate(filterId="gvb", list("GoodVsBad"=c(0, 9999)))
idx <- filter(synPerturbed.flowclean, rg)
synPerturbed.flowcleangated <- Subset(synPerturbed.flowclean, idx)

plot1<-ggcyto(synPerturbed, aes(x='Time', y='<V705-A>'))+ geom_hex(bins=512)+ coord_cartesian(xlim = c(0, 100000))+ coord_cartesian(ylim = c(-1000, 70000)) +labs(title =paste0("Raw:",nrow(synPerturbed)," cells"))+ expand_limits(x=0)
plot2<-ggcyto(synPerturbed.flowcleangated, aes(x='Time', y='<V705-A>'))+ geom_hex(bins=512)+ coord_cartesian(xlim = c(0, 100000))+ coord_cartesian(ylim = c(-1000, 70000)) +labs(title =paste0("flowClean:",nrow(synPerturbed.flowcleangated)," cells"))+ expand_limits(x=0)
plot3<-ggcyto(synPerturbed.flowcut$frame, aes(x='Time', y='<V705-A>'))+ geom_hex(bins=512)+ coord_cartesian(xlim = c(0, 100000))+ coord_cartesian(ylim = c(-1000, 70000)) +labs(title =paste0("flowCut:",nrow(synPerturbed.flowcut$frame)," cells"))+ expand_limits(x=0)
plot4<-ggcyto(synPerturbed.flowai, aes(x='Time', y='<V705-A>'))+ geom_hex(bins=512)+ coord_cartesian(xlim = c(0, 100000))+ coord_cartesian(ylim = c(-1000, 70000)) +labs(title =paste0("flowAI:",nrow(synPerturbed.flowai)," cells"))+ expand_limits(x=0)

#real data
realdata<-read.FCS("220662.fcs")

realdata.flowclean <- clean(realdata, vectMarkers=c(4:12), filePrefixWithDir="sampleName", ext="fcs")
realdata.flowcut <- flowCut(realdata, Channels=c(4:12) )
realdata.flowai<-flowAI::flow_auto_qc(realdata)

rg2 <- rectangleGate(filterId="gvb2", list("GoodVsBad"=c(0, 9999)))
idx2 <- filter(realdata.flowclean, rg2)
realdata.flowcleangated <- Subset(realdata.flowclean, idx2)

plot5<-ggcyto(realdata, aes(x='Time', y='<PE Cy55-A>'))+ geom_hex(bins=512)+ coord_cartesian(xlim = c(0, 4000))+ coord_cartesian(ylim = c(-1000, 3000)) +labs(title =paste0("Raw:",nrow(realdata)," cells"))+ expand_limits(x=0)
plot6<-ggcyto(realdata.flowcleangated, aes(x='Time', y='<PE Cy55-A>'))+ geom_hex(bins = 512)+ coord_cartesian(xlim = c(0, 4000))+ coord_cartesian(ylim = c(-1000, 3000))+labs(title =paste0("flowClean:",nrow(realdata.flowcleangated)," cells"))+ expand_limits(x=0)
plot7<-ggcyto(realdata.flowcut$frame, aes(x='Time', y='<PE Cy55-A>'))+ geom_hex(bins = 512)+ coord_cartesian(xlim = c(0, 4000))+ coord_cartesian(ylim = c(-1000, 3000))+labs(title =paste0("flowCut:",nrow(realdata.flowcut$frame)," cells"))+ expand_limits(x=0)
plot8<-ggcyto(realdata.flowai, aes(x='Time', y='<PE Cy55-A>'))+ geom_hex(bins = 512)+ coord_cartesian(xlim = c(0, 4000))+ coord_cartesian(ylim = c(-1000, 3000))+labs(title =paste0("flowAI:",nrow(realdata.flowai)," cells"))+ expand_limits(x=0)

#plot it all out
grid.arrange(as.ggplot(plot1), as.ggplot(plot2), as.ggplot(plot3), as.ggplot(plot4), ncol=2, top = "Artifical data")
grid.arrange(as.ggplot(plot5), as.ggplot(plot6), as.ggplot(plot7), as.ggplot(plot8), ncol=2, top = "Real data")
