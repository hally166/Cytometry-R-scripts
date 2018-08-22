#Christopher Hall, Wellcome Sanger Institute, ch15@sanger.ac.uk
#Bulk export of ARIA index information V2, now with added flowsetieness
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
  editme$WellID<-chartr("01234567", "ABCDEFGH", editme$XLoc) #Repalce the XLoc numbers with well letters
  editme$WellID<-paste(editme$WellID,editme$YLoc +1, sep="") #add one to each YLoc to match the plate designations
  write.csv(editme, file=paste(keyword(new_frame,"GUID"), "_index.csv", sep="")) #write a csv using the GUID keyword as the filename
})
