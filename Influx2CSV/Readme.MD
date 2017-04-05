# Influx2CSV
The BD Influx flow cytometer sorter stores the index data (well positions) of a sort in the .fcs file.  This script batch exports this information to a .csv file.

Copyright (c) 2017 Genome Research Ltd.

Author : Christopher Hall, Wellcome Trust Sanger Institute, christopher.hall@sanger.ac.uk

http://www.sanger.ac.uk/science/groups/cytometry-core-facility

This script exports the index data from BD Influx files created with Sortware v1.2.0.142.  It has not been tested with any other versions of Sortware.  It does not know your gating scheme so no gates will be exported.  

### Instructions
Place your files into a folder

Change the 'setwd()' to that folder

Run the script

Your files are saved in the same folder with '_index.csv' appended to the filename

You can choose which parameters to export on the 'indexdata<-subset(indexdata.....' line
