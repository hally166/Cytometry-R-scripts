# load files (for Aurora)
file1<-read.FCS("C:/Users/hallc/Downloads/xxx/Reference Group/A2 Unstained (Cells).fcs",truncate_max_range = FALSE, transform=FALSE)
file2<-read.FCS("C:/Users/hallc/Downloads/xxx/Reference Group/B10 CD45 APC-Fire 810 (Cells).fcs",truncate_max_range = FALSE, transform=FALSE)

# works for Aurora
file4<-rbind(exprs(file1),exprs(file2))
exprs(file2)<-file4
#or if downsampling
exprs(file2)<-file4[sample(nrow(file4),5000),]

write.FCS(file2,"CD45_fire810_5000.fcs")
