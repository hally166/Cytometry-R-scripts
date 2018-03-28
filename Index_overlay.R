#index overlay
#v1.0 Sept 2017
#R 3.4.1
#Author : Christopher Hall, Wellcome Trust Sanger Institute, christopher.hall@sanger.ac.uk

#This script creates an overlay dot plot from the index data extracted from the influx2csv script

#Use this line on the influx index 2 csv script to make a csv of the fcs file
write.csv(dataframe, file=paste(filename, "_index.csv", sep=""))

#install packages
install.packages("ggplot2")
library(ggplot2)

#read csv
fullfile<-read.csv("D:/sanger/fcsfile.csv")
indexfile<-read.csv("D:/sanger/indexfile.csv")

#what are the parameter names?
colnames(indexfile)
colnames(fullfile)

#Well numbers
indexfile[2]

#plot the data. Replace x & y with your parameters 
fullplot<-ggplot(fullfile, aes(x = X610.20..561., y = X530.30..488.)) + geom_point()
indexplot<-ggplot(indexfile, aes(x = X610.20..561., y = X530.30..488.)) + geom_point()

#format the data
fp<-fullplot+scale_x_log10(name="x",limits=c(1, 10000))+scale_y_log10(name="y",limits=c(1, 10000))
ip<-indexplot+scale_x_log10(name="x",limits=c(1, 10000))+scale_y_log10(name="y",limits=c(1, 10000))

#look at the data
fp
ip

#overlay the data
#If you want specific cells add [c(2,4,15),] to indexfile below (replace 2,4,15 with whatever you like)
#Replace x & y with your parameters
fp+geom_point(data=indexfile, aes(x = X610.20..561., y = X530.30..488.),colour="red1")
