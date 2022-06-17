# Installing the "basic" flow cytometry packages.  They all live on a repository called Bioconductor.
install.packages("BiocManager")

#You can either do library(BiocManager) and install("flowCore) or what I have done below.
BiocManager::install("flowCore") #interpret .fcs files
BiocManager::install("flowViz") #basic visulisation
BiocManager::install("ggcyto") #advanced visulisation using the ggPlot nomenclature
BiocManager::install("openCyto") #Used to link various analysis methodologies
BiocManager::install("flowWorkspace") #used to build anaysis templates
BiocManager::install("CytoML") #imports FlowJo and DiVA workspaces

#These packages largly require each other to work (except flowCore which is the "base package) 
#so will often load each other without my help.  For simplicty I have loaded them all.

#You will need to "clean" your data.  flowAI and flowCut are my recomendations.  
#flowClean is the original, but is supeceeded by flowCut
BiocManager::install("flowAI")
BiocManager::install("flowClean")
BiocManager::install("flowCut")

#The CytoExploreR package is not available on bioconductor and needs to be loaded straight from GitHub.  To do this you need the package devtools.
install.packages("devtools")

#An intresting project is CytoExploreR that trys to blend the power of R with the ease of use of a mouse.
devtools::install_github("DillonHammill/CytoExploreRData")
devtools::install_github("DillonHammill/CytoExploreR",build_vignettes = TRUE)

#Load a single fcs
myfile <- "C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZU_files/0001.FCS"
fcsfile <- flowCore::read.FCS(myfile)
library(flowCore)
fcsfile1 <- read.FCS(myfile)
fcsfile

#Load many fcs files into a flow set
myfiles <- list.files(path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZU_files", pattern=".FCS$")
fs <- flowCore::read.flowSet(myfiles, path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZU_files/")
fs
fs[[1]]
