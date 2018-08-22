#Christopher Hall, Wellcome Sanger Institute, ch15@sanger.ac.uk
#Bulk export of ARIA index information V3, now with added flowsetieness and 384 well plate support
#Change setwd to the directory of your folder

#load the required packages from BioConductor
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")
library(flowCore)
setwd("C:/Users/ch15/Desktop/files")

files <- list.files(path=".", pattern=".fcs$") #list files in directory the $ means it MUST end with .fcs
fs<-read.flowSet(files) #puts all the .fcs files into a flowset
fsApply(fs,function(frame){ #loop through the flowset
  comp <- keyword(frame)$SPILL #take the compensation matrix
  new_frame <- compensate(frame,comp) #apply the comepensation matrix
  editme<- data.frame(getIndexSort(new_frame)) #get the index data and convert to a data.frame to help with the WellID
  editme$XLoc<-chartr("0", "A", editme$XLoc)
  platemap<-setNames(c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","N","P"), c("A","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
  editme$WellID<-platemap[editme$XLoc]
  editme$XLoc<-chartr("A", "0", editme$XLoc)
  editme$WellID<-paste(editme$WellID,editme$YLoc +1, sep="") #add one to each YLoc to match the plate designations
  write.csv(editme, file=paste(keyword(new_frame,"GUID"), "_index.csv", sep="")) #write a csv using the GUID keyword as the filename
})
