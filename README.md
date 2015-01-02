phrynomics-data
===============

![](http://barbbanbury.info/barbbanbury/Research_Projects_files/phrynoHead.jpg)

This is a repository for the workflows and data from the phrynomics project. The phrynomics package is a separate suite of tools for SNP data. See [phrynomics](https://github.com/bbanbury/phrynomics)

## Files

We present all of the files in this repository to maximize open-source reproducibility. 

The phrynomicsFunctions.R file includes many custom functions that do not belong in the general phrynomics package. Many of these are hard coded for our particular datasets. These functions range from setting up analyses to scraping the results files for info. Some of these functions are deprecated and were only used during data exploration. 

The file phrynoSetup.R is the R script that we used to manipulate the phrynosoma SNP datasets we got out of pyRad.  This file includes (commented out) shell and perl scripts for running all our analyses on the UW cluster.  

The files within Rscripts directory and files master.sh and sim.ctl were used for simulation. sim.ctl is the control file for the MCCOAL program in BPP, and contains all of the parameters for the simulated datasets.  The master.sh file is a hacky shell script that pieces together all of the components of the simulations, only extracting off the important pieces (such as total tree lengths). This script calls on MCCOAL and the Rscripts.  

The phrynoAnalyses.R file was used for all post-analyses scraping of results files.  This file contains scripts to remake many of the tables and figures found in our paper.  

Finally, the phrynoResults.Rdata file is an R-formatted data file that can be loaded into an R workspace and used to create the tables and figures.  It contains the full and snp only data files for all 14 minimum individual datasets, the RAxML trees, and branch length comparison tables. 