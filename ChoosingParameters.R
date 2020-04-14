#You need to input which parameters you need to analyse many times when dealing with FCS data in R.
#Functions either require a list of names of parameters or an index of parameters (numbers).

#Following are two functions that will help and can be used standalone or within fsApply.

#Create a list of names

#Function - change the names within"" to exclude more (e.g. -W, FL1 etc..)
autovect_verbose<- function(ff){
  c<- data.frame(ff@parameters@data)
  d<- grep("FSC|SSC|Time", c$name, invert = TRUE, value = TRUE)
  return(unname(d))

#Example use - transformation
x<-autovect_verbose(fs[[1]])
tf<-estimateLogicle(fs[[1]], channels =  x)
fs_trans<-transform(fs,tf)
#Use in function
x<-fsApply(fs, function(fr) autovect_verbose(fr))


#Create an index

#Function - change the names within"" to exclude more (e.g. -W, FL1 etc..)
autovect<- function(ff){
  c<- data.frame(ff@parameters@data)
  d<- grep("FSC|SSC|Time", c$name, invert = TRUE)
  return(as.vector(d))
  
#Example use - in flowClean
fsclean<-fsApply(fs, function(fr) clean(fr, vectMarkers = autovect(fr),filePrefixWithDir = paste0(keyword(fr)$GUID," --QCOutput"), ext = "fcs", diagnostic = TRUE))
#Standard use
x<-autovect(fs_comp[[1]])
x<-fsApply(fs, function(fr) autovect(fr))
