# NgonorrhoeaeLinkageGenomics
## This repository contains a collection of Perl scripts for population genetic analysis of genomic data and C++ programs used to simulate various evolutionary scenarios.

Instructions to replicate analyses:
Download or clone contents into a directory, move into this directory and type the following command:

>gunzip ./*.gz

to uncompress all gzip’d files in the current directory.

To uncompress the FASTA files used in the analysis, type

>tar zxf UK_FilteredFASTAS.tar.gz UK_FilteredFASTAS

To uncompress files only for the UK dataset. This may take a minute or so. The current script is set up
to analyze the UK dataset (using metadata specific to this population), but this may be changed by modifying
variables at the beginning of the NgonorrhoeaeLinkageGenomics.pl Perl script to analyse the US and NZ datasets
that are also mentioned in the manuscript.

To run the Perl script, type: 

>perl NgonorrhoeaeLinkageGenomics.pl tmp

and output files will be deposited into the new directory named “tmp”. To change parameters in the analysis, see
the code on lines 22-32. For instance, the number of bootstrap replicates may be modified by changing the $Bootreps
variable, which is currently set to 1 to speed computation for testing. Also, the variable $wd_fa contains the 
directory containing the FASTA files, which would need to be changed if you want to analyze the other datasets.

NOTE: some of these analyses may take a bit of memory, certianly more than ~7Gb and perhaps more than ~15Gb.

