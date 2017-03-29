#Copyright (c) 2017 Genome Research Ltd.

#Influx2CSV
#v1.11 Mar 2017
#R 3.3.3
#Author : Christopher Hall, Wellcome Trust Sanger Institute, christopher.hall@sanger.ac.uk

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.


#Extracts the index data from an Influx .fcs file and saves it as a .csv (Sortware v1.2.0.142 only, but you can test it if you like).  


#Install required packages from bioconductor and CRAN.  'Robustbase' may be required, or may not.
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore", dependencies = TRUE)
#install.packages('robustbase')
library(flowCore)

#Use this section to look at an example file and to choose your parameters.  
setwd("")  #Remember to swap \ with / on windows machines
flowfile<-read.FCS('') #Choose your example file
flowfile
params<-c("FSC", "SSC", "*530/30 (488)") #Insert you parameters of interest here, add as many as you like

#This section loops through the index files and creates the .csv files
files <- list.files(path=".", pattern=".fcs")
for (fileName in files) {
  flowfile<-read.FCS(fileName) #transformation remains on for influx index files
  if (keyword(flowfile, "INDEXSORTPOSITIONS") %in% keyword(flowfile) ==TRUE){
    ofintrest<-exprs(flowfile)[c(1:nrow(flowfile)), c("Sort Result Bits", "Tray X", "Tray Y", params)]
    dataframe<-as.data.frame(ofintrest)
    indexed<-dataframe[dataframe$'Sort Result Bits' >0,]
    
    wellorder<-strsplit(keyword(flowfile,"INDEXSORTPOSITIONS")[[1]], ",")
    WellID <- wellorder[[1]][seq(1, length(wellorder[[1]]), 3)]
    x<-wellorder[[1]][seq(2, length(wellorder[[1]]), 3)]
    y<-wellorder[[1]][seq(3, length(wellorder[[1]]), 3)]
    Well_df<-data.frame(WellID,x,y)
    Well_df[, 2] <- as.numeric(as.character( Well_df[, 2] ))
    Well_df[, 3] <- as.numeric(as.character( Well_df[, 3] ))
    colnames(Well_df) <- c("Well","Tray X", "Tray Y")
    
    #Sortware does not store the index positions  correctly so we have to create a multiple flowframes to find them
    Well_df_x1<-Well_df
    Well_df_x1[,2]<-Well_df[,2] +1
    
    Well_df_x0<-Well_df
    Well_df_x0[,2]<-Well_df[,2] -1
    
    Well_df_y1<-Well_df
    Well_df_y1[,3]<-Well_df[,3] +1
    
    Well_df_y0<-Well_df
    Well_df_y0[,3]<-Well_df[,3] -1
    
    Well_df_xy2<-Well_df
    Well_df_xy2[,2]<-Well_df[,2] +1
    Well_df_xy2[,3]<-Well_df[,3] +1 
    
    Well_df_xy0<-Well_df
    Well_df_xy0[,2]<-Well_df[,2] -1
    Well_df_xy0[,3]<-Well_df[,3] -1
    
    Well_df_xy02<-Well_df
    Well_df_xy02[,2]<-Well_df[,2] -1
    Well_df_xy02[,3]<-Well_df[,3] +1
    
    Well_df_xy20<-Well_df
    Well_df_xy20[,2]<-Well_df[,2] +1
    Well_df_xy20[,3]<-Well_df[,3] -1
    ############
    Well_df_x2<-Well_df
    Well_df_x2[,2]<-Well_df[,2] +2
    
    Well_df_xm2<-Well_df
    Well_df_xm2[,2]<-Well_df[,2] -2
    
    Well_df_y2<-Well_df
    Well_df_y2[,3]<-Well_df[,3] +2
    
    Well_df_ym2<-Well_df
    Well_df_ym2[,3]<-Well_df[,3] -2
    
    Well_df_xy22<-Well_df
    Well_df_xy22[,2]<-Well_df[,2] +2
    Well_df_xy22[,3]<-Well_df[,3] +2 
    
    Well_df_xym2<-Well_df
    Well_df_xym2[,2]<-Well_df[,2] -2
    Well_df_xym2[,3]<-Well_df[,3] -2
    
    Well_df_xy04<-Well_df
    Well_df_xy04[,2]<-Well_df[,2] -2
    Well_df_xy04[,3]<-Well_df[,3] +2
    
    Well_df_xy40<-Well_df
    Well_df_xy40[,2]<-Well_df[,2] +2
    Well_df_xy40[,3]<-Well_df[,3] -2
    
    Well_df_xy24<-Well_df
    Well_df_xy24[,2]<-Well_df[,2] -1
    Well_df_xy24[,3]<-Well_df[,3] +2
    
    Well_df_xy42<-Well_df
    Well_df_xy42[,2]<-Well_df[,2] +2
    Well_df_xy42[,3]<-Well_df[,3] -1 
    ##
    Well_df_xy12<-Well_df
    Well_df_xy12[,2]<-Well_df[,2] +1
    Well_df_xy12[,3]<-Well_df[,3] +2
    
    Well_df_xy21<-Well_df
    Well_df_xy21[,2]<-Well_df[,2] +2
    Well_df_xy21[,3]<-Well_df[,3] +1 
    
    Well_df_xym12<-Well_df
    Well_df_xym12[,2]<-Well_df[,2] -1
    Well_df_xym12[,3]<-Well_df[,3] -2
    
    Well_df_xym21<-Well_df
    Well_df_xym21[,2]<-Well_df[,2] -2
    Well_df_xym21[,3]<-Well_df[,3] -1 
    
    Well_df_xyo21<-Well_df
    Well_df_xyo21[,2]<-Well_df[,2] -2
    Well_df_xyo21[,3]<-Well_df[,3] +1 
    
    Well_df_xyx21<-Well_df
    Well_df_xyx21[,2]<-Well_df[,2] +1
    Well_df_xyx21[,3]<-Well_df[,3] -2     
    
    a<-merge(indexed,Well_df_x1, c("Tray X", "Tray Y"))
    b<-merge(indexed,Well_df_x0, c("Tray X", "Tray Y"))
    c<-merge(indexed,Well_df_y1, c("Tray X", "Tray Y"))
    d<-merge(indexed,Well_df_y0, c("Tray X", "Tray Y"))
    e<-merge(indexed,Well_df, c("Tray X", "Tray Y"))
    
    f<-merge(indexed,Well_df_xy2, c("Tray X", "Tray Y"))
    g<-merge(indexed,Well_df_xy0, c("Tray X", "Tray Y"))
    h<-merge(indexed,Well_df_xy02, c("Tray X", "Tray Y"))
    i<-merge(indexed,Well_df_xy20, c("Tray X", "Tray Y"))
    ###
    j<-merge(indexed,Well_df_x2, c("Tray X", "Tray Y"))
    k<-merge(indexed,Well_df_xm2, c("Tray X", "Tray Y"))
    l<-merge(indexed,Well_df_y2, c("Tray X", "Tray Y"))
    m<-merge(indexed,Well_df_ym2, c("Tray X", "Tray Y"))
    n<-merge(indexed,Well_df_xy22, c("Tray X", "Tray Y"))
    
    o<-merge(indexed,Well_df_xym2, c("Tray X", "Tray Y"))
    p<-merge(indexed,Well_df_xy04, c("Tray X", "Tray Y"))
    q<-merge(indexed,Well_df_xy40, c("Tray X", "Tray Y"))
    r<-merge(indexed,Well_df_xy24, c("Tray X", "Tray Y"))
    
    s<-merge(indexed,Well_df_xy42, c("Tray X", "Tray Y"))
    t<-merge(indexed,Well_df_xy12, c("Tray X", "Tray Y"))
    u<-merge(indexed,Well_df_xy21, c("Tray X", "Tray Y"))
    v<-merge(indexed,Well_df_xym12, c("Tray X", "Tray Y"))
    w<-merge(indexed,Well_df_xym21, c("Tray X", "Tray Y"))
    aa<-merge(indexed,Well_df_xyo21, c("Tray X", "Tray Y"))
    bb<-merge(indexed,Well_df_xyx21, c("Tray X", "Tray Y"))   
    
    
    indexdata<-rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,aa,bb)
    indexdata<-subset(indexdata,select=c('Well', params))
    filename=paste(substr(keyword(flowfile, "GUID"),1,(nchar(keyword(flowfile, "GUID"))-4))) #You could change 'GUID' to 'FIL' if you wish, i.e filename to actual saved name
    write.csv(indexdata, file=paste(filename, "_index.csv", sep=""))}
}

