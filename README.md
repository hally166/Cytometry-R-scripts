# R-Scripts
#### R scripts to help with your flow cytometry analysis
Copyright (c) 2017 Genome Research Ltd.

http://www.sanger.ac.uk/science/groups/cytometry-core-facility

## Influx2CSV
The BD Influx flow cytometer sorter stores the index data (well positions) of a sort in the .fcs file.  This script batch exports this information to a .csv file.

## Aria2csv
An R script to batch export data from Aria .fcs files.  NOTE: it won't work with 384 well plates because it is lazily coded.

## index_overlay
Overlays selected or all index data over a dot plot of all the data in the fcs file
