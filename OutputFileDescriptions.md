#The following are notes about each of the output files. Note that D stands for the letter used to calculate
pairwise associations between SNPs (D = Pij - Pi*Pj, where Pij is the frequency with which we observe the 
derived mutations on the same haplotype minus the expected frequency based on random chance). Also, rN-rS
stands for the difference between r calculated for nonsynonymous SNPs and for synonymous SNPs, where
r = D/sqrt(Pi(1-Pi)Pj(1-Pj))

###QCfile.txt
	This file contains a lot of output, mostly for sanity checks, but also general summary statistics
	that may be useful. For instance, it contains how many core genes were detected, how many
	synonymous and nonsynonymous SNPs were detected, but it also contains a list of genes in which 
	position information could not be accurately resolved.

###SFS.txt
	Contains the mutation site-frequency spectrum, broken down into synonymous, nonsynonymous, and
	premature stop codons (a type of nonsynonymous mutation). These data were used in Figure S2.

###PositionsOfNonsenseMutationsPerGene.txt
	For each gene that contains a premature stop codon polymorphisms, this file records at what position
	in the gene this polymorphism was found, along with the total gene length. This is done separately
	for biallelic sites and for multiallelic sites (which are exceedingly rare in N. gonorrhoeae).

###NumNonsenseMutsPerGene.txt
	Contains the number of nonsense mutations per gene.

###AlignmentStats_NumGappedSites.txt
	Contains information about the gene alignments for assessing quality, such as the total length of
	the gene alignment and how many of these positions were ungapped (i.e. A,G,C, or T) or gapped 
	(i.e. "-").

###ZerofoldDividedByFourfoldDiversity.txt
	Contains the overall ratio of zerofold-degenerate diversity divided by fourfold-degenerate
	diversity for each quantile (index starting at 0). The number of lines in this file depend
	on how many quantiles you specify to analyze in the main Perl script. The last column
	specifies how many genes were in that quantile and used to calculate the statistic. These data
	were used for Figure 3B in the manuscript.

###Polymorphism_and_SummaryStatistics.txt
	This file contains additional summary statistics for each gene, including the number of Bi-allelic
	and Multi-allelic sites, the number of synonymous and nonsynonymous polymorphisms and fixed 
	differences (with respect to the reference), as well as Tajima's D and Fay and Wu's H, which are
	both summaries of the mutation site-frequency spectrum. These data were used in Table S4 that 
	showed many of the genes that had the largest excess of NSyn couplings also had an excess of
	intermediate frequency alleles (highly positive Tajima's D).

###DensityNonSynPoly_meanByQuantile.txt
	Contains the mean density of nonsynonymous polymorphisms per quantile, along with the number of genes
	contributing to that quantile. For instance, if one analyzes genes by categorizing them by
	their nonsynonymous SNP density (as I have in my example script), this file will contain the mean 
	density of nonsynonymous polymorphisms per site for each quantile.

###D_R_Dprime_Rsq_vs_Dist_SYN_win*.txt
	There may be several of these files, where the asterisk is replaced by a number starting at 0.
	the number of genomic windows (win) depends on the window size specified in the Perl script.
	If a very long distance is specified (greater than whole genome), there will only be one window
	(win0). This file contains linkage information for SYNONYMOUS SNPs as a function of the distance 
	between SNPs. For each inter-SNP distance category (single bp resolution), this file shows the 
	number of synonymous SNPs that were found at exactly that distance apart, along with the expected
	values for three commonly used metrics to quantify linkage: D, r, D', and r-squared. These data
	were used to make Figure 1 in the manuscript.

###D_R_Dprime_Rsq_vs_Dist_NONSYN_win*.txt
	Same as file immediately above except for NONSYNONYMOUS SNPs.
	
###D_r_Dprime_vsQuantSYNWithinGene.txt
	This file contains the mean values of various linkage metrics for all genes within a SNP density
	quantile. These values were calculated only for synonymous SNPs. These data were used in Figure 4D
	in the manuscript.

###D_r_Dprime_vsQuantNONSYNWithinGene.txt
	This file is the same as the one above except linkage metrics are computed only for nonsynonymous
	SNPs. These data were used in Figure 4D in the manuscript.

###SynResamplingBinsNumObservations.txt
	For each inter-SNP distance bin, this file contains how many SYNONYMOUS SNP pairs were detected
	for each bin (NumObs_WG). It then breaks these observations down into the number of pairs that 
	were found within genes (NumObs_Intragenic) or between genes (NumObs_Intergenic).

###NonSynResamplingBinsNumObservations.txt
	Same as file immediately above except for NONSYNONYMOUS SNPs.

###NonsynRvSynR_WholeGenome_Quant*.txt
	There may be several of these files, where the asterisk is replaced by a number starting at 0.
	These numbers represent the quantile into which analyses were partitioned (e.g. by only
	analyzing genes within a certain nonsynonymous SNP density quantile, as I do in my example).
	For each of these files, each row represents an inter-SNP distance bin in which we've calculated
	rN-rS for that bin, downsampling whatever category (synonymous/nonsynonymous) is more frequent 
	and bootstrap resampling this process, randomly downsampling each time. This file shows the
	5%, 50%, and 95% quantile from the bootstrapped distribution, along with the number of 
	synonymous and nonsynonymous SNP pairs that contributed to that bin. These data were used in
	Figure 2 in the manuscript.

###NonsynRvSynR_Intragenic_Intergenic_BySNPdensityQuantile.txt
	This file also contains information about rN-rS, which is calculated as described above with
	boostrap resampling. However, instead of calculating rN-rS by inter-SNP distance bin, it is
	calculated for all genes in a particular quantile, where here genes are quantiled by nonsynonymous
	SNP density. Again, this value is calculated for many bootstrap replicates, downsampling 
	whatever category of SNPs is more frequent (synonymous or nonsynonymous). Reported in this
	file is the 5%, 50%, and 95% quantile of this bootstrapped distribution. These data were used 
	in Figure 3C in the manuscript.

###NonsynDprSynDpr_Intragenic_Intergenic_BySNPdensityQuantile.txt
	This file is the same as the one above except linkage is computed using D', instead of r. These
	data were used in Figure S12 in the manuscript.

###NonsynRvSynR_ByGene_Overall.txt
	Contains rN-rS for each individual gene, but may contain "NA" if there were too few polymorphisms.
	While the rN-rS value is reported for the entire gene, it is actually calculated by computing
	rN-rS for each inter-SNP distance bin within the gene, then taking the weighted average
	across all bins. This was done to build in a way always compare nonsynonymous and synonymous
	SNP pairs that have roughly similar linkage dynamics (for instance, genes that are very long could
	show that linkage decays within the gene, and it would be nice to compare rN and rS only between
	SNP pairs that have roughly similar linkage dynamics, such as those within 500bp are compared
	separately with one another from those pairs that are greater than 1kb apart). For my analyses,
	I just select an inter-SNP distance window that is larger than any gene (e.g. 20kb), such that 
	rN-rS is calculated across the entire gene, irrespective of distance. Values of rN-rS calculated 
	in this much simpler way were very similar to those in which genes were binned into smaller SNP 
	distance intervals, and a weighted average across the entire gene was used to calculate a final 
	rN-rS value. These data were used to make Figure 5 in the manuscript.
