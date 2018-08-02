#Christopher Hall, Wellcome Sanger Institute, ch15@sanger.ac.uk
#R3.5.1
#Bulk export of ARIA index information
#Change setwd to the directory of your files

source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")
library(flowCore)
setwd("")

files <- list.files(path=".", pattern=".fcs")
for (fileName in files) {
  flowfile<-read.FCS(fileName)
  editme<- data.frame(getIndexSort(flowfile))
  editme$WellID<-chartr("01234567", "ABCDEFGH", editme$XLoc)
  editme$WellID<-paste(editme$WellID,editme$YLoc +1, sep="")
  write.csv(editme, file=paste(fileName, "_index.csv", sep=""))
}
