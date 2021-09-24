# Manual full spectrum flow cytometry unmxing
# Christopher Hall, Babraham Institute, Sept 2021

#load required packages
library(flowCore)

#load data - aurora data often needs truncate_max_range = TRUE
control_fs<-read.flowSet(list.files('',full.names=TRUE)) # flowset of positive control, you will want to pre-gate them first so they are 100% positive, I used FlowJo and exported the gate
file2unmix<-read.FCS('') 

#Get the expression data for the file to unmix
expresionData<-exprs(file2unmix)
expresionData<-as.data.frame(expresionData)
expresionData<-expresionData[,-grep("SC", names(expresionData))]
expresionData<-expresionData[,grep("-A", names(expresionData))]

#Get the expression data for the controls
Control_Spectrums<-fsApply(control_fs,each_col,median)
Control_Spectrums<-as.data.frame(Control_Spectrums)
Control_Spectrums<-Control_Spectrums[,-grep("SC", names(Control_Spectrums))]
Control_Spectrums<-Control_Spectrums[,grep("-A", names(Control_Spectrums))]
Control_Spectrums<-Control_Spectrums[names(expresionData)] #Reorder data to match the expression data

#Ignoring the unstained control

#use lsfit() to do the linear least squares unmixing - taken from specUnmixCoFunction() in https://github.com/jtheorell/flowSpecs 
ls_corr <- lsfit(x = t(Control_Spectrums), y = t(expresionData), intercept = FALSE)
unmixResult <- t(ls_corr$coefficients)
unmixResult<-as.data.frame(unmixResult)

#plot the result - change the setting to your file and appearance preference. compare it with the instrument unmixed file for reference
plot(unmixResult$`export_4 PE 20201029122042_PE.fcs`,unmixResult$`export_8 FITC 20201029121953_FITC.fcs`, log="xy", main="Unweighted", xlab="PE CD4", ylab="FITC CD8", pch=19, cex=0.1, cex.lab=1,cex.axis=1,cex.main=2)


#Using the the unstained control

negative_control<-read.FCS('') #aurora data often needs truncate_max_range = TRUE

#Get the expression data for the unstained
unstainedData<-exprs(negative_control)
unstainedData<-apply(unstainedData,2,median)
unstainedData<-as.data.frame(t(unstainedData))
unstainedData<-unstainedData[,-grep("SC", names(unstainedData))]
unstainedData<-unstainedData[,grep("-A", names(unstainedData))]
unstainedData<-unstainedData[names(expresionData)] #Reorder data to match the expression data

#Remove unstained medians from control medians
Control_Spectrums2<-mapply('-', Control_Spectrums, unstainedData, SIMPLIFY = TRUE)
rownames(Control_Spectrums2)<-rownames(Control_Spectrums)

#use lsfit() to do the linear least squares unmixing - taken from specUnmixCoFunction() in https://github.com/jtheorell/flowSpecs 
ls_corr2 <- lsfit(x = t(Control_Spectrums2), y = t(expresionData), intercept = FALSE)
unmixResult2 <- t(ls_corr2$coefficients)
unmixResult2<-as.data.frame(unmixResult2)

#plot the result - change the setting to your file and appearance preference. compare it with the instrument unmixed file for reference
plot(unmixResult2$`export_4 PE 20201029122042_PE.fcs`,unmixResult2$`export_8 FITC 20201029121953_FITC.fcs`, log="xy", main="weighted", xlab="PE CD4", ylab="FITC CD8", pch=19, cex=0.1, cex.lab=1,cex.axis=1,cex.main=2)
