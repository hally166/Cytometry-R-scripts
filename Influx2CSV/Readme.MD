#Influx2CSV
The BD Influx stores the index data (well positions) of a sort in the .fcs file.  This script batch exports this information to .csv file.

Copyright (c) 2017 Genome Research Ltd.

Author : Christopher Hall, Wellcome Trust Sanger Institute, christopher.hall@sanger.ac.uk

http://www.sanger.ac.uk/science/groups/cytometry-core-facility

This script only exports the index data from BD Influx files created with Sortware v1.2.0.142.  It has not been tested with any other versions of Sortware.  It does not know your gating scheme so no gates will be shown.  

This script is not pretty, or efficient, but it works.  If you want to analyse your data in R you should reseach 'flowCore' and use this script to understand how the index data is saved in the fcs file.

###Instructions
Learn the basics of how to use an R script.

Place your files into a folder.

Change the 'setwd()' line to that folder.

Read a single example file by adding it to 'flowfile<-read.FCS()'.

Choose yor parameters and add them to 'params<-c("FSC", "SSC")'.

Run the rest of the script.

Your files are saved in the same folder with '_index.csv' appended to the filename.
