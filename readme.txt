This repository contains all the necessary code to reproduce the main figures and central data/statistics

Everything necessary to generate table and figures should already by found in the "data" folder upon download

Following Qiime2 processing, 16S sequence data have been arranged into a phyloseq R-Data object (from the "phyloseq" R package)
	- The phyloseq data represent the "initial" state of the data input into R scripts beginning with 1_qsip_rate_calculations.R
	- Other initial data are in csv or tab-separated text files
	- Subsequent scripts proceed with the calculations outlined in the methods section of the manuscript
	- To save on space, the phyloseq data object, and all intermediate data products have been saved as RData files and may be easily loaded by typing load('data_name.RData') into the R console

Several functions have been created for data processing and visualization and may be found in the "figure_creation_functions.R" script
