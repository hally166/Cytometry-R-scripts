# Manual full spectrum flow cytometry unmxing
# Christopher Hall, Babraham Institute, Sept 2021
# This works better on Bigfoot data then on Aurora data, not sure why, but suspect it has to do with the unstained control baseline, I'm working on it

#load required packages
library(flowCore)

#load data
control_fs<-read.flowSet(list.files('C:/file_location',full.names=TRUE)) # flowset of positive control, you will want to pregate them first so they are 100% positive, I used FlowJo and exported the gate
file2unmix<-read.FCS('C:/file_location.fcs')

#Get the expression data for the controls
Control_Spectrums<-fsApply(control_fs,each_col,median)
Control_Spectrums<-as.data.frame(Control_Spectrums)
Control_Spectrums<-Control_Spectrums[,-grep("SC", names(Control_Spectrums))]
Control_Spectrums<-Control_Spectrums[,grep("-A", names(Control_Spectrums))]

#Get the expression data for the file to unmix
expresionData<-exprs(file2unmix)
expresionData<-as.data.frame(expresionData)
expresionData<-expresionData[,-grep("SC", names(expresionData))]
expresionData<-expresionData[,grep("-A", names(expresionData))]

#use lsfit() to do the linear least squares unmixing - taken from specUnmixCoFunction() in https://github.com/jtheorell/flowSpecs 
ls_corr <- lsfit(x = t(Control_Spectrums), y = t(expresionData), intercept = FALSE)
unmixResult <- t(ls_corr$coefficients)
unmixResult<-as.data.frame(unmixResult)

#plot the result - change the setting to your file and appearance preference. compare it with the instrument unmixed file for reference
plot(unmixResult$`PE_CD4.fcs`,unmixResult$`FITC_CD8.fcs`, log="xy", main="Unweighted", xlab="PE CD4", ylab="FITC CD8", pch=19, cex=0.5, cex.lab=1,cex.axis=1,cex.main=2)
