#This is the R code from my January 2021 RMS Flow Facilities Meeting presentation
#Christopher Hall, Babraham Institute, UK

#Change the files locations to your own, remember to use / in file locations (NOT \).  It's an R thing. 

#Install packages
install.packages("BiocManager")
BiocManager::install("flowCore")
library(flowCore)

#load single file
myfile <- "C:/Users/hallc/Desktop/Database/140600.fcs"
fcsfile <- read.FCS(myfile)
fcsfile

#load multiple files into a flowset
myfiles <- list.files(path="C:/Users/hallc/Desktop/Snapshot", pattern=".fcs", full.names = TRUE)
flowSet <- read.flowSet(myfiles)
flowSet
flowSet[[2]]

#install flowAI
BiocManager::install("flowAI")
library(flowAI)
flow_auto_qc(flowSet)
flow_auto_qc(fs,folder_results="D:/cleaned_files")

#explore keywords
x<-keyword(flowSet[[2]])
keyword(flowSet[[1]])$"EXPORT USER NAME"
keyword(flowSet[[1]])$"TUBE NAME"

keyword(flowSet)$"EXPORT USER NAME"
fsApply(flowSet, keyword(flowSet)$"EXPORT USER NAME")
?fsApply
fsApply(flowSet, function(x) keyword(x)$"EXPORT USER NAME")
fsN<- fsApply(flowSet, function(x) keyword(x)$"EXPORT USER NAME")

#load csv file
read.csv("ResultsQC/QCmini.txt")
?read.csv
read.csv("ResultsQC/QCmini.txt", sep = "\t")
csv<-read.csv("ResultsQC/QCmini.txt", sep = "\t")
fsN
csv

#cahnge row names and merge csv
row.names(fsN)<-gsub('.fcs', '',row.names(x))
newCSV <- merge(x=csv, y=fsN, by.x="Name.file", by.y=0)
view(newCSV)
names(newCSV)[names(newCSV)=="V1"] <- "User"
View(newCSV)
write.csv(newCSV, "QCResults.csv")

#try on an actual database
myfiles <- list.files(path="C:/Users/hallc/Desktop/Database", pattern=".fcs", full.names = TRUE)
flowSet <- read.flowSet(myfiles)
for (file in myfiles[1:5]){
  flow_auto_qc(file, fcs_QC = FALSE, output = 0, remove_from = "FS_FM")
}

#cytoflex files
myfiles <- list.files(path="C:/Users/hallc/Desktop/Cytoflex", pattern=".fcs", full.names = TRUE, recursive = TRUE)

#attune files
myfiles <- list.files(path="C:/Users/hallc/Desktop/Attune", pattern=".fcs", full.names = TRUE, recursive = TRUE)
flowSet <- read.flowSet(myfiles, emptyValue = TRUE)

#flowSOM
install.packages("devtools")
devtools::install_github("saeyslab/FlowSOM", build_vignettes = TRUE)
library("FlowSOM")
?FlowSOM

# Read from file
fileName <- system.file("extdata", "68983.fcs", package="FlowSOM")
flowSOM.res <- FlowSOM(fileName, compensate=TRUE,transform=TRUE,
                       scale=TRUE,colsToUse=c(9,12,14:18),nClus=10)

# Plot results
PlotStars(flowSOM.res$FlowSOM,
          backgroundValues = flowSOM.res$metaclustering)
