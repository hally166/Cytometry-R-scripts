# This script is split into 3 sections, twice. The last bit is the same as the first, just with less comments. 
# Each section details how to use one of the flow cytometry FCS data cleaning packages.
# flowClean
# flowCut
# flowAI 
# I will show how to use them for single flowFrames and flowSets, as well as showing the most common variables.

#Author: Christopher Hall, EMBL:Rome

library(flowCore)

#FlowClean - https://pubmed.ncbi.nlm.nih.gov/27153628/
BiocManager::install("flowClean")
library(flowClean)

#using your own data
realdata<-read.FCS("220662.fcs") #this is your data file
realdata #display the metadata to decide which parameters to analyse, i.e. not FSC,SSC,Time
realdata.flowclean <- clean(realdata, vectMarkers=c(4:12), filePrefixWithDir="sampleName", ext="fcs")#vectMarkers are the pramameters you want to analyse. This also outputs diagnostics to your working directory folder
#realdata.flowclean is a new flowFrame with a GoodVsBad paparmeter added to use this do:
rg <- rectangleGate(filterId="gvb", list("GoodVsBad"=c(0, 9999))) #new gate of everyting under 10000 on the GoodVsBad parameter
idx <- filter(realdata.flowclean, rg) #create a filter of the data
realdata.flowcleangated <- Subset(realdata.flowclean, idx) #make a new flowFrame of the "good" data. You use this flowFram for your analysis

#Using flowClean in a flowSet
fs<- read.flowSet(files) #create a flowset from a list of files
#The following function loops through all the files in the flowSet and produces an outpt file using the name of the sample name as the template
fsclean<-fsApply(fs, function(fr) clean(fr, vectMarkers=c(4:12),filePrefixWithDir = paste0(keyword(fr)$GUID," --QCOutput"), ext = "fcs", diagnostic = TRUE))
#if you are doing this over many experiments you could create a function that will work out the vectMarkers automatically
autovect<- function(ff){
  c<- data.frame(ff@parameters@data)
  d<- grep("FSC|SSC|Time|-W|-H", c$name, invert = TRUE) #change the text within the "" to look like the parameters you do not want to analyse
  return(as.vector(d))
}
fsclean<-fsApply(fs, function(fr) clean(fr, vectMarkers = autovect(fr),filePrefixWithDir = paste0(keyword(fr)$GUID," --QCOutput"), ext = "fcs", diagnostic = TRUE))


#flowCut - https://github.com/jmeskas/flowCut
devtools::install_github("jmeskas/flowCut")
library(flowCut)

realdata<-read.FCS("220662.fcs") #this is your data file
realdata.flowcut <- flowCut(realdata, Channels=c(4:12)) #this runs flowCut and outputs diagnostics to /flowCut folder
#realdata.flowcut is now a flowCut object containing multiple pieces of data.  TO access each bit use:
realdata.flowcut #to see an overview
realdata.flowcut$frame # for the flowFrame to analyse
realdata.flowcut$data #for the metadata
realdata.flowcut$worstChan #for just the worst channel, useful to visalise any problems by using a plotting tool

#Using flowCut in a flowSet
fs<- read.flowSet(files) #create a flowset from a list of files
QC_flowCut<-flowCore::fsApply(fs, function(x)flowCut(x)) #a funcion to run all the files in a flowset 
#to access the data use:
QC_flowCut[[1]]$frame #for individual files
#or to make a new flowset of the data
files<-list() 
for(x in 1:length(QC_flowCut)){
  files<-append(files,list(QC_flowCut[[x]]$frame))
}
fs_flowCut<-flowSet(files) #fs_flowcut is your new flowset

#you could use the same function as in flowClean above called autovect to work out the parameters to analyse over multiple experiemnts
autovect<- function(ff){
  c<- data.frame(ff@parameters@data)
  d<- grep("FSC|SSC|Time|-W|-H", c$name, invert = TRUE) #change the text within the "" to look like the parameters you do not want to analyse
  return(as.vector(d))
}
QC_flowCut<-flowCore::fsApply(fs, function(x)flowCut(x, Channels = autovect(fr)))


#flowAI - https://pubmed.ncbi.nlm.nih.gov/26990501/
#By default running FlowAI will produce new fcs files in your working directory with the bad events highlighted. 
BiocManager::install("flowAI")
library(flowAI)

realdata<-read.FCS("220662.fcs") #this is your data file
resQC <- flowAI::flow_auto_qc(realdata) 
#resQC is a new flowFrame with your good data

#using flowAI in a flowSet
fs<- read.flowSet(files, emptyValue = FALSE) #create a flowset from a list of files
resQCfs <- flowAI::flow_auto_qc(fs) 
#resQCfs is a new flowSet with your good data. It will also procuce new FCS files in your working directory.

#read the flowAI docuemntation about how it works.  IT has one very important variable that can effect your results: "remove_from ="
resQC <- flow_auto_qc(realdata, remove_from = c("FR", "FS")) # FR= flow rate, FS = signal, FR = dynamic range.



# Without most of the comments



library(flowCore)

#FlowClean - https://pubmed.ncbi.nlm.nih.gov/27153628/
BiocManager::install("flowClean")
library(flowClean)

#using your own data
realdata<-read.FCS("220662.fcs")
realdata
realdata.flowclean <- clean(realdata, vectMarkers=c(4:12), filePrefixWithDir="sampleName", ext="fcs")

rg <- rectangleGate(filterId="gvb", list("GoodVsBad"=c(0, 9999)))
idx <- filter(realdata.flowclean, rg)
realdata.flowcleangated <- Subset(realdata.flowclean, idx)

#Using flowClean in a flowSet
fs<- read.flowSet(files)
fsclean<-fsApply(fs, function(fr) clean(fr, vectMarkers=c(4:12),filePrefixWithDir = paste0(keyword(fr)$GUID," --QCOutput"), ext = "fcs", diagnostic = TRUE))
autovect<- function(ff){
  c<- data.frame(ff@parameters@data)
  d<- grep("FSC|SSC|Time|-W|-H", c$name, invert = TRUE)
  return(as.vector(d))
}
fsclean<-fsApply(fs, function(fr) clean(fr, vectMarkers = autovect(fr),filePrefixWithDir = paste0(keyword(fr)$GUID," --QCOutput"), ext = "fcs", diagnostic = TRUE))


#flowCut - https://github.com/jmeskas/flowCut
devtools::install_github("jmeskas/flowCut")
library(flowCut)

realdata<-read.FCS("220662.fcs")
realdata.flowcut <- flowCut(realdata, Channels=c(4:12))
realdata.flowcut #to see an overview
realdata.flowcut$frame # for the flowFrame to analyse
realdata.flowcut$data #for the metadata
realdata.flowcut$worstChan #for just the worst channel

#Using flowCut in a flowSet
fs<- read.flowSet(files)
QC_flowCut<-flowCore::fsApply(fs, function(x)flowCut(x)) 
QC_flowCut[[1]]$frame
#to make a new flowset of the data
files<-list() 
for(x in 1:length(QC_flowCut)){
  files<-append(files,list(QC_flowCut[[x]]$frame))
}
fs_flowCut<-flowSet(files)

autovect<- function(ff){
  c<- data.frame(ff@parameters@data)
  d<- grep("FSC|SSC|Time|-W|-H", c$name, invert = TRUE)
  return(as.vector(d))
}
QC_flowCut<-flowCore::fsApply(fs, function(x)flowCut(x, Channels = autovect(fr)))


#flowAI - https://pubmed.ncbi.nlm.nih.gov/26990501/
BiocManager::install("flowAI")
library(flowAI)

realdata<-read.FCS("220662.fcs")
resQC <- flowAI::flow_auto_qc(realdata) 

#using flowAI in a flowSet
fs<- read.flowSet(files, emptyValue = FALSE)
resQCfs <- flowAI::flow_auto_qc(fs) 

#read the flowAI docuemntation about how it works.  IT has one very important variable that can effect your results: "remove_from ="
resQC <- flow_auto_qc(realdata, remove_from = c("FR", "FS")) # FR= flow rate, FS = signal, FR = dynamic range.
