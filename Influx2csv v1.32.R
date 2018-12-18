#InfluxIndex2CSV
#v1.32 Jun 2017
#R 3.3.3
#Author : Christopher Hall, Wellcome Trust Sanger Institute, christopher.hall@sanger.ac.uk
#License : GPLv3 https://www.gnu.org/licenses/gpl-3.0.html

#Extracts the index data from an Influx .fcs file and saves it as a .csv (Sortware v1.2.0.142 only, but you can test it on other versions if you like).  


#Install required packages from bioconductor and CRAN.
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")
library(flowCore)
install.packages("fuzzyjoin")  
library(fuzzyjoin)

setwd("")  #Point the script to your directory, Remember to swap \ with / on windows machines

#This section loops through the index files and creates the .csv files
files <- list.files(path=".", pattern=".fcs")
for (fileName in files) {
  flowfile<-read.FCS(fileName,alter.names=TRUE) #you may or may not want to use data transformation here transformation=FALSE
  if (keyword(flowfile, "INDEXSORTPOSITIONS") %in% keyword(flowfile) ==TRUE){
    
    #Creates a datafrrame of the sorted events
    ofintrest<-exprs(flowfile)
    dataframe<-as.data.frame(ofintrest)
    indexed<-dataframe[dataframe$'Sort.Result.Bits' >0,]
    
    #Creates a dataframe of the well positions
    wellorder<-strsplit(keyword(flowfile,"INDEXSORTPOSITIONS")[[1]], ",")
    WellID <- wellorder[[1]][seq(1, length(wellorder[[1]]), 3)]
    x<-wellorder[[1]][seq(2, length(wellorder[[1]]), 3)]
    y<-wellorder[[1]][seq(3, length(wellorder[[1]]), 3)]
    Well_df<-data.frame(WellID,x,y)
    Well_df[, 2] <- as.numeric(as.character( Well_df[, 2] ))
    Well_df[, 3] <- as.numeric(as.character( Well_df[, 3] ))
    colnames(Well_df) <- c("Well","Tray.X", "Tray.Y")
    
    #Sortware does not store the index positions correctly so we use thr dplyr and fuzzylogic packages to correct for this
    indexdata<-difference_inner_join(Well_df,indexed,by=c("Tray.X", "Tray.Y"),max_dist = 3)
    colnames(indexdata)<-c('Well','Tray X (kw)','Tray Y (kw)',colnames(read.FCS(fileName))) #Reading the fcs file again to get the column names slows the script down by 1/3.  To avoid this pre populate this with the col names.  Doing it this way makes the script more resistant to changes in the .fcs files
    
    #indexdata<-subset(indexdata,select=c("Well", "FSC", "SSC")) #use this line to choose your paramaters if you do not want all the columns
    
    #creates the csv file
    filename=paste(substr(keyword(flowfile, "GUID"),1,(nchar(keyword(flowfile, "GUID"))-4))) #You could change 'GUID' to 'FIL' if you wish, i.e filename to actual saved name
    write.csv(indexdata, file=paste(filename, "_index.csv", sep=""))}
  }
  
