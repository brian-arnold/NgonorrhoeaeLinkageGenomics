NgonorrhoeaeLinkageGenomics.pl uses the following input file: 
>NgonorrhoeaeFA1090_NmeningitidisAlpha14_AlleleTable.txt

that is essentially a reformatted version of the output from the program 
progressiveMauve, but is easier to parse for downstream analyses. This subdirectory
contains the necessary script, named progressiveMauveAlignment_2_AlleleTable.pl,
to recompute this table, using a reference genome named:

>NgonorrhoeaeFA1090_genomic.fna

along with the output of progressiveMauve, named:

>NgonorrhoeaeFA1090_NmeningitidisAlpha14_progMauve.alignment

The progressiveMauve alignment was produced by aligning the N. gonorrhoeae FA1090 
reference sequence with the N. meningitidis Alpha14 reference sequence.
