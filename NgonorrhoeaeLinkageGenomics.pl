#!usr/bin/perl
use warnings ;
use strict ;
use lib "./" ; # This should contain the file path that contains the Perl module loaded next
use BacterialGenomics_GitHub5 ;

################################################################################################
#	Prerequisites
#	Bioinformatic pipeline
#		1.) assemble raw reads into de novo assemblies, and "polish" these assemblies with a program like Pilon
#		2.) annotate assemblies with Prokka, throw a completed reference genome in there to keep track of relative positions of genes in the more fragmented de novo assemblies 
#		3.) find orthologs within these annotated assemblies using Roary (which current pipeline is for, but I'm looking to explore alternatives)
#		3a.) I realigned the gene alignments output from Roary with a program like PAGAN, something that's codon-aware. Roary support PRANK but it's horrendously slow
#
#	In addition, I construct a multispecies allele table to determine ancestral state of nucleotides. I essentially take N. gonorrhoeae FA1090 reference genome and
#	align it to another genome (or set of genomes) with progressiveMauve. I use a Perl script to convert this progressiveMauve output to an allele table that is
#	easier to parse. This script is included, along with an example, in another directory "GenerateMultispeciesAlleletable"

# NOTE ABOUT CODE BELOW: datastructures with names ending in "_HR" are hash references, since many of these are passed by reference to downstream Perl functions
################################################################################################

my $subdir = $ARGV[0] ;			# subdirectory to put results in
my $Quantiles = 1 ;			# number of quantiles to bin genes into, e.g. by Nonsynonymous SNP density
my $min_AF = 0.0 ;			# minimum allele frequency to use for various analyses, such as calculating linkage or associations between derived mutations
my $max_AF = 1.0 ;			# maximum allele frequency to use for various analyses, such as calculating linkage or associations between derived mutations
my $window_size = 5000000 ; 		# Added to calculate pairwise linkage stats in genomic windows (bp), goes faster, for testing
my $max_dist = 5000000 ; 		# Maximum distance between SNPs for calculating pairwise linkage, lower values make script faster, for testing
my $min_dist = 0 ; 			# minimum distance between SNPs for calculating pairwise linkage
my $MinNonsynSNPperGene_forQuants = 0 ; # Minimum number of SNPs in a gene for it to be included in analysis; e.g., exclude genes with no polymorphism
my $Bootreps = 1 ;			# For analyses that have bootstrap replicates, this is set here
my $interSNPdist_binSize = 100 ;	# For linkage analyses that bin SNPs by distance (e.g. all SNPs within 10bp of each other as bin 0), this is set here
my $interSNPdist_numBins = 20000 ;	# the total number of distance bins to use; if you wish to analyze the whole genome (~2Mb for gonorrhoeae) this number times bin size above should roughly equal the genome size 

my $wd_fa = "./UK_FilteredFASTAS" ;					# name of directory that stores gene alignment fastas to be analyzed, e.g. from Roary
my @fasta_core = `ls ${wd_fa}/*` ;					# store names of all fasta files (1 per gene) in array
my $gene_presence_absence_file = "./UK_gene_presence_absence.csv" ;	# This file comes from Roary, it contains information about the number of copies per gene, i.e. whether a gene is core or accessory
my $Ref_gff_file = "./NgonorrhoeaeFA1090_modified.gff" ;		# Gff file of reference sequence, this is used to get position of genes with respect to a known, closed reference
my $Ref_multi_species_alignment = "./NgonorrhoeaeFA1090_NmeningitidisAlpha14_AlleleTable.txt" ; # multispecies allele table file used to infer ancestral state, mentioned above
my $Reference_ind = "NgonorrhoeaeFA1090_modified" ;			# Name of the reference individual used that is processed along with other samples (i.e with Prokka,Roary,etc..)

# Create subdirectory to store results
my $wd ;
my $base_dir = "./" ;
if($subdir){
	$wd = "./${subdir}" ;
	system("mkdir ${wd}") ;				
}else{
	$wd = $base_dir ;
}
open QC, ">>${wd}/QCfile.txt" ;		# file with some output for sanity-checking results

################################################
### SUBSAMPLE INDIVIDUALS ACCORDING TO METADATA 
### Here according to collection date
### This code was not made into subroutine since it needs to be modified for different datasets, and it's currently hard-coded for the UK data.
################################################
my %Individuals_total ;			# list of all individuals in the dataset
my %Individuals_subsample ;		# list of a subset of individuals, sampled accordion to metadata below
my $min_yr = 2013 ;			# earliest collection date
my $max_yr = 2013 ;			# latest collection date

# Get the total list of individuals
open IN, "<./UK_TotalListOfIsolates.txt" ; 
while(<IN>){
    chomp $_ ;
	my $ID = $_ ;
	$Individuals_total{$ID} ++ ;
}
close IN ;

# Subsample from total list of isolates using the metadata
open IN, "./UK_metadata_Table.txt" ; 
while(<IN>){
	chomp $_ ;
	my @line = split(/\t/, $_) ;
	my $date_tmp = $line[3] ;
	my @date_tmp2 = split(/-/, $date_tmp) ;
	my $sampleYear = $date_tmp2[0] ;
	my $ind = $line[5] ;
	if( (exists $Individuals_total{$ind}) && ($sampleYear >= $min_yr) && ($sampleYear <= $max_yr) ){
		# individual names have a "_modified" tag appended, since I had to modify the gff files output from PROKKA
		# so that they contained names according to the KEGG database, as opposed to the NCBI database (default)
		$ind = $ind."_modified" ;
		$Individuals_subsample{$ind}++ ;
	}
}
close IN ;

$Individuals_subsample{$Reference_ind} ++ ;						# Add in reference individual if it hasn't already been done
my $num_ind_NoRef = (scalar keys %Individuals_subsample)-1 ;	# the number of individuals NOT including the reference sequence
print QC "Total number of Individuals\/Isolates, including the reference sequence: ", scalar keys %Individuals_total, "\n" ;
print QC "Number of subsampled Individuals\/Isolates, including the reference sequence: ", scalar keys %Individuals_subsample, "\n" ;

################################################
### The following subroutines collect necessary information into data structures and load the sequence data.
################################################

################################################
### CREATE CODON TABLE HASH, MAPPING 3-NUCLEOTIDE CODONS TO AMINO ACIDS
################################################
my %Codon_Hash ; # $Codon_Hash{codon} = amino_acid ;
%Codon_Hash = Construct_Codon_hash () ;

################################################
### LOAD MULTIPLE-SPECIES ALIGNMENT, FOR ASSIGNING ANCESTRAL ALLELE
################################################
my $AlignFile_index_ForReference = 1 ; # 0 for N. gonorrhoeae; 1 for outgroup species e.g. N. meningitidis
my $MultiSpeciesAlign_HR ; # @{$MultiSpeciesAlign_HR{ pos }} = ( allele0, allele1) ; where allele0 is gono, allele 1 is outgroup 
$MultiSpeciesAlign_HR = MultipleSpeciesAlignment_2RefGenomes( \$Ref_multi_species_alignment, \$wd) ;

################################################
### OBTAIN GENOMIC POSITIONS (IN FA1090 REFERENCE) OF EACH GENE
################################################
my $Ref_gene_info_HR ;  # $Ref_gene_info_HR{$gene_name}{ "POSITION_REVERSE"/"POSITION_FORWARD" }{ "START"/"END" } = start/end
my $Ref_gene_positions_HR ;  # $Ref_gene_positions_HR{$start}{$end} = $gene_name ;
($Ref_gene_info_HR, $Ref_gene_positions_HR) = GetReferenceGenePositions( \$Ref_gff_file, \$wd ) ;

################################################
#### ONLY LOOK AT GENES WITH INFORMATION FROM REFERENCE (i.e. core genes that are also found in the FA1090 reference)
################################################
my @fasta_core_ref_overlap ;
@fasta_core_ref_overlap = FastaCoreOverlapWithReference ( \@fasta_core, $Ref_gene_info_HR, \$wd ) ;

################################################
# PROKKA ID -> ASSEMBLY ID
# Prokka creates names for fasta sequences that differ from the original assembly name, Roary keeps these names
# use the gene_presence_absence.csv output by Roary to map these Prokka name back to assembly names (or the names of the isolates).
################################################
my $ProkkaID_2_AssemblyName_HR ;# $ProkkaID_2_AssemblyName_HR{Prokka ID in fasta file} = assembly name
my $ProkkaID_2_GeneName_HR ;	# $ProkkaID_2_GeneName_HR{Prokka ID in fasta file} = gene name
($ProkkaID_2_AssemblyName_HR, $ProkkaID_2_GeneName_HR ) = ProkkaID_2_AssemblyID( \$gene_presence_absence_file, \%Individuals_subsample, \$wd ) ;

#############################################
# LOAD GENE SEQUENCES
# This function also checks for reading frame violations, likely due to alignment errors,
# which can be problematic when trying to infer functional effects of mutations (Nonsyn/Syn)
#############################################
my $core_loci_seqs_HR ; # $core_loci_seqs{ gene name }{ individual } = sequence
$core_loci_seqs_HR = LoadGeneSeqsForEachIndividual(\@fasta_core_ref_overlap, $ProkkaID_2_AssemblyName_HR, \%Individuals_subsample, \$Reference_ind, \$wd ) ;

#############################################
# MAP ROARY COG ALIGNMENT POSITION TO MULTISPECIES TABLE ALIGNMENT POSITION
# CREATE HASH WITH ANCESTRAL ALLELE
#############################################
my $COG_AlnPos_2_WGref_AlnPos_HR ;	# COG_AlnPos_2_WGref_AlnPos{gene}{Roary cog alignment pos in Reference_ind} = multispecies ref alignment pos
my $RefBase_HR ; # $RefBase_HR{ gene }{ site } = base ;
$COG_AlnPos_2_WGref_AlnPos_HR = RoaryCOGpos_2_MultispeciesAlleleTablePos( $core_loci_seqs_HR, $Ref_gene_info_HR, \$Reference_ind, $MultiSpeciesAlign_HR, \$wd ) ;
$RefBase_HR = MakeOutgroupBaseHash( $core_loci_seqs_HR, $COG_AlnPos_2_WGref_AlnPos_HR, $MultiSpeciesAlign_HR, $Ref_gene_info_HR, \$Reference_ind, \$AlignFile_index_ForReference) ;

#############################################
# GET CONSENSUS BASE PER POSITION FOR CALCULATING DEGENERACY IN REFERENCE SEQ
#############################################
my $Degeneracy_HR ; # $Degeneracy{ gene }{ site } = n-fold degeneracy, 0D,2D,3D,4D
$Degeneracy_HR = CalculatePositionDegeneracy($core_loci_seqs_HR, \$Reference_ind, \%Codon_Hash) ;


#############################################
# INTERLUDE
# The previous functions have gathered basic information from the reference individual, 
# the outgroup individual, and loaded the DNA sequences per gene, per individual.
# The following functions finally analyze these DNA sequences
#############################################

#############################################
# COLLECT INFORMATION ON POLYMOPHISM AND DIVERGENCE 
# A lot happens within this subroutine. Data structures are created that store the locations and frequencies of polymorphic sites within gene alignments,
# as well as sites with differences that are fixed with respect to the reference sequence. Additionally, the functional effect of these mutations
# (Synonymous, Nonsynonymous) is also stored. Lastly, a data structure stores information about 0-fold and 4-fold degenerate sites.
# NOTE: The following sites are ignored: Sites with gaps; sites in which functional effect (Syn/Nonsyn) cannot be inferred because of flanking gaps; sites in which no reference base is available
#############################################
my $Core_Biallelic_Segsites_HR ; 		# $Core_Biallelic_Segsites_HR{ gene }{ seg_site_index } = $site ; where seg_site_index ranges from 0 to total number of SNPs
my $Core_FixedDiffs_HR ;				# $Core_FixedDiffs_HR{ gene }{ seg_site_index } = $site ; where seg_site_index goes from 0 to to total number of fixed differences
my $Core_NumUngapped_Sites_HR ; 		# $Core_NumUngapped_Sites_HR{ gene } = number ; or the number of ungapped positions for the gene alignment
my $FunctionalEffect_BiAllelic_HR ; 	# $FunctionalEffect_BiAllelic_HR{ gene }{ site } =  "SYNONYMOUS" or "NONSYNONYMOUS:MISSENSE" or "NONSYNONYMOUS:NONSENSE"  ; 
my $FunctionalEffect_MultiAllelic_HR ; 	# $FunctionalEffect_MultiAllelic_HR{ gene }{ site }{ base1 }{ base2 } =  "SYNONYMOUS" or "NONSYNONYMOUS:MISSENSE" or "NONSYNONYMOUS:NONSENSE"; $base1 and $base2 include all permutations, e.g. b1=A/b2=T and b1=T/b2=A, plug in ReferenceBase for $base1
my $FunctionalEffect_Fixed_HR ; 		# $FunctionalEffect_Fixed_HR{ gene }{ site } =  "SYNONYMOUS" or "NONSYNONYMOUS:MISSENSE" or "NONSYNONYMOUS:NONSENSE"  ;
my $Biallelic_Sites_Freqs_HR  ; 		# $Biallelic_Sites_Freqs_HR{$gene}{ site }{ base } = frequency
my $Multiallelic_Sites_Freqs_HR ; 		# $Multiallelic_Sites_Freqs_HR{ gene }{ site }{ base } = frequency
my $Monomorphic_Sites_HR ;				# $Monomorphic_Sites_HR{ gene }{ site } = base (A,G,C,T) ;
my $ZeroFourfold_NumSegSites_HR ; 		# $ZeroFourfold_NumSegSites_HR{ gene }{"0D" or "4D"} = Number of sites ;
my $ZeroFourfold_NumUnGappedSites_HR ; 	# $ZeroFourfold_NumUnGappedSites_HR{ gene }{"0D" or "4D"} = Number of sites ;
($Core_Biallelic_Segsites_HR, $Core_FixedDiffs_HR, $Core_NumUngapped_Sites_HR, $FunctionalEffect_BiAllelic_HR, $FunctionalEffect_MultiAllelic_HR, $FunctionalEffect_Fixed_HR, $Biallelic_Sites_Freqs_HR, $Multiallelic_Sites_Freqs_HR, $Monomorphic_Sites_HR, $ZeroFourfold_NumSegSites_HR, $ZeroFourfold_NumUnGappedSites_HR) = CollectSegSite_FixedSite_FunctionalEffects($core_loci_seqs_HR, \$Reference_ind, \%Codon_Hash, $RefBase_HR, $Degeneracy_HR, \$wd) ;

#############################################
## CHARACTERIZE NONSENSE MUTATIONS
#############################################
# This collects information about Nonsense mutations within genes.
my $NonsenseMutsPerGene_HR ; # NonsenseMutsPerGene_HR{gene} = number
$NonsenseMutsPerGene_HR = CharacterizeNonsenseMutations($core_loci_seqs_HR, $FunctionalEffect_BiAllelic_HR, $FunctionalEffect_MultiAllelic_HR, \$Reference_ind, \$wd) ;

#############################################
# CALCULATE MUTATION SITE-FREQUENCY SPECTRUM
#############################################
my $SFS_polarized_w_RefSeq_HR ; # $SFS_polarized_w_RefSeq_HR{"SYN"/"NONSYN"}{ freq in sample } = number of occurrences in genome ;
$SFS_polarized_w_RefSeq_HR = CalculateSFS($Biallelic_Sites_Freqs_HR, $FunctionalEffect_BiAllelic_HR, $RefBase_HR, \$num_ind_NoRef, \$wd) ;

#############################################
# GET POSITIONS OF SNPS WITH RESPECT TO REFERENCE
# This allows us to think about linkage disequilibrium beyond a single gene, by using the positions of the genes in the reference sequence
# This information is put in a data structure with an index/counter from 0 to the number of SNPs, ordered with respect to the reference
# which facilitates (1) iterating through all SNPs and (2) grouping them into genomic windows (e.g. 50kb windows)
#############################################
my $Biallelic_Syn_SNP_Ref_Pos_HR ; 				# Biallelic_Syn_SNP_Ref_Pos_HR{counter} = reference position ; counter ranges from 0 to number of Syn SNPs
my $Biallelic_Syn_SNP_Alignment_Coord_HR ;		# @{$Biallelic_Syn_SNP_Alignment_Coord_HR{counter}} = (gene, seg_site_index) ; counter ranges from 0 to number of Syn SNPs, seg_site_index ranges from 0 to total number of SNPs in gene
my $Biallelic_NonSyn_SNP_Ref_Pos_HR ; 			# Biallelic_NonSyn_SNP_Ref_Pos_HR{counter} = reference position ; counter ranges from 0 to number of NonSyn SNPs
my $Biallelic_NonSyn_SNP_Alignment_Coord_HR ; 	# @{$Biallelic_NonSyn_SNP_Alignment_Coord_HR{counter}} = (gene, seg_site_index) ; counter ranges from 0 to number of NonSyn SNPs, seg_site_index ranges from 0 to total number of SNPs in gene
## NOTE: counter is different for Syn and Nonsyn data structures!
($Biallelic_Syn_SNP_Ref_Pos_HR, $Biallelic_Syn_SNP_Alignment_Coord_HR, $Biallelic_NonSyn_SNP_Ref_Pos_HR, $Biallelic_NonSyn_SNP_Alignment_Coord_HR) = GetPositionsSynNonsynBiallelicSites($Ref_gene_positions_HR, $Ref_gene_info_HR, $RefBase_HR, $Core_Biallelic_Segsites_HR, $Biallelic_Sites_Freqs_HR, \$num_ind_NoRef, $FunctionalEffect_BiAllelic_HR, \$min_AF, \$max_AF, \$wd) ;

#############################################
## CALCULATE DIVERSITY AND SFS SUMMARY STATISTICS
#############################################
## This calculates the proportion of Bi-allelic and Multi-allelic sites for each category of codon degeneracy (i.e. 0-fold degenerate sites, 2-fold, 3-fold, 4-fold)
my $Num_Variable_Sites_By_Degeneracy_HR ;	#			$Num_Variable_Sites_By_Degeneracy_HR{"total"/"biallelic"/"multiallelic"}{ 0D,2D,3D,4D } = num occurrences;
($Num_Variable_Sites_By_Degeneracy_HR) = CalculateNumberVariableSites($Degeneracy_HR, $FunctionalEffect_BiAllelic_HR, $Multiallelic_Sites_Freqs_HR, \$wd) ;

# Calculate Diversity statistics PER GENE (i.e. statistics ARE NOT divided by gene lengths), which are later used for SFS summary statistics
my $Wattersons_Theta_Bi_HR ;
my $Pi_Theta_Bi_HR ;
my $L_theta_Bi_HR ;
my $SegSites_Bi_HR ;
($Wattersons_Theta_Bi_HR, $Pi_Theta_Bi_HR, $L_theta_Bi_HR, $SegSites_Bi_HR) = CalculateThetasBiAllelicByGene($Biallelic_Sites_Freqs_HR, $FunctionalEffect_BiAllelic_HR, $RefBase_HR, \$num_ind_NoRef) ;

# SFS SUMMARY STATISTICS
my $TajD_PerGene_HR ;
my $FayWuH_PerGene_HR ;
($TajD_PerGene_HR, $FayWuH_PerGene_HR) = SFS_SummStats($Wattersons_Theta_Bi_HR, $Pi_Theta_Bi_HR, $L_theta_Bi_HR, $SegSites_Bi_HR, \$num_ind_NoRef) ;

############################################
# CONSTRUCT WINDOWS
# If you wish to bin reference positions into windows of e.g. 50kb, done separately for Synonymous and Nonsynonymous mutations
#############################################
my $Windows_Syn ; # @{$Windows_Syn{$window}} = (all the counts)
my $Windows_NonSyn ;
($Windows_Syn, $Windows_NonSyn) = ConstructWindows($Biallelic_Syn_SNP_Ref_Pos_HR, $Biallelic_NonSyn_SNP_Ref_Pos_HR, \$window_size) ;

#############################################
###  CATEGORIZE GENES BY SNP DENSITY, PRINT POLYMORPHISM INFORMATION AND SFS SUMMARIES
#############################################
my $DensityNonSynPoly_byGene_Quantiled_HR ; # $$DensityNonSynPoly_byGene_Quantiled_HR{$gene} = quantile, [0, $Quantiles-1]
($DensityNonSynPoly_byGene_Quantiled_HR) = CategorizeBySNPdensity_PrintSFSsummaries($Core_NumUngapped_Sites_HR, $Biallelic_Sites_Freqs_HR, $Multiallelic_Sites_Freqs_HR, $FunctionalEffect_BiAllelic_HR, $FunctionalEffect_MultiAllelic_HR, $FunctionalEffect_Fixed_HR, $Ref_gene_info_HR, $RefBase_HR, \$Quantiles, $TajD_PerGene_HR, $FayWuH_PerGene_HR, $NonsenseMutsPerGene_HR, \$MinNonsynSNPperGene_forQuants, \$wd) ;
# Print the number of genes per quantile; general function in case other metrics are used to quantile the data
PrintNumGenesPerQuantile($DensityNonSynPoly_byGene_Quantiled_HR, \$wd) ;

#############################################
## 0FOLD VS 4FOLD DIVERSITY BY QUANTILE
#############################################
ZerofoldFourfoldDivByQuantile($ZeroFourfold_NumSegSites_HR, $ZeroFourfold_NumUnGappedSites_HR, $DensityNonSynPoly_byGene_Quantiled_HR, \$wd) ;

#############################################
### CALCULATE PAIRWISE LINKAGE STATS BY DISTANCE WITHIN WINDOWS
#############################################
my $PairWise_LD_vs_dist_SYN ; # @{$PairWise_LD_vs_dist_SYN{$win}{$dist}}, ($LD) ;
my $PairWise_PC_vs_dist_SYN ; # @{$PairWise_PC_vs_dist_SYN{$win}{$dist}}, $compatibility ;
my $PairWise_LD_vs_dist_NONSYN ; # @{$PairWise_LD_vs_dist_SYN{$win}{$dist}}, ($LD) ;
my $PairWise_PC_vs_dist_NONSYN ; # @{$PairWise_PC_vs_dist_SYN{$win}{$dist}}, $compatibility ;

# This calculates linkage disequilibrium as a function of inter-SNP distance
CalculatePairwiseLinkageStats_ByDistance($Windows_Syn, $Windows_NonSyn, $Biallelic_Syn_SNP_Alignment_Coord_HR, $Biallelic_Syn_SNP_Ref_Pos_HR, $Biallelic_NonSyn_SNP_Alignment_Coord_HR, $Biallelic_NonSyn_SNP_Ref_Pos_HR, $Core_Biallelic_Segsites_HR, $RefBase_HR, $core_loci_seqs_HR, \$Reference_ind, \$max_dist, \$wd) ;
CalculatePairwiseLinkageStats_ByGeneCategory($Biallelic_Syn_SNP_Alignment_Coord_HR, $Biallelic_Syn_SNP_Ref_Pos_HR, $Biallelic_NonSyn_SNP_Alignment_Coord_HR, $Biallelic_NonSyn_SNP_Ref_Pos_HR, $Core_Biallelic_Segsites_HR, $RefBase_HR, $core_loci_seqs_HR, \$Reference_ind, \$min_dist, \$wd, $DensityNonSynPoly_byGene_Quantiled_HR) ;

#############################################
#### CALCULATE R BY INTER-SNP DISTANCE BIN 
##############################################
## Min and Max allele frequencies imposed when constructing "Biallelic_Syn_SNP_Alignment_Coord_HR" hash
BinLDbyDistForBoostrapResampling_WholeGenome(\$interSNPdist_binSize, \$interSNPdist_numBins, $Biallelic_Syn_SNP_Alignment_Coord_HR, $Biallelic_Syn_SNP_Ref_Pos_HR, $Biallelic_NonSyn_SNP_Alignment_Coord_HR, $Biallelic_NonSyn_SNP_Ref_Pos_HR, $Core_Biallelic_Segsites_HR, $RefBase_HR, $core_loci_seqs_HR, \$Reference_ind, $DensityNonSynPoly_byGene_Quantiled_HR, \$Quantiles, \$Bootreps, \$wd) ; 
BinLDbyDistForBoostrapResampling_IntraInterGenic(\$interSNPdist_binSize, \$interSNPdist_numBins, $Biallelic_Syn_SNP_Alignment_Coord_HR, $Biallelic_Syn_SNP_Ref_Pos_HR, $Biallelic_NonSyn_SNP_Alignment_Coord_HR, $Biallelic_NonSyn_SNP_Ref_Pos_HR, $Core_Biallelic_Segsites_HR, $RefBase_HR, $core_loci_seqs_HR, \$Reference_ind, $DensityNonSynPoly_byGene_Quantiled_HR, \$Quantiles, \$Bootreps, \$wd) ;
BinLDbyDistForBoostrapResampling_ByGene(\$interSNPdist_binSize, \$interSNPdist_numBins, $Core_Biallelic_Segsites_HR, $RefBase_HR, $core_loci_seqs_HR, \$Reference_ind, $FunctionalEffect_BiAllelic_HR, \$min_AF, \$max_AF, \$Bootreps, $DensityNonSynPoly_byGene_Quantiled_HR, \$wd) ;

close QC ;





exit ;

