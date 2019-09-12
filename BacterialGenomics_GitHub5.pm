package BacterialGenomics_GitHub5 ;
use strict ;
use warnings ;

# Perl Exporter allows export of functions/variables of modules to user’s namespace using
# the standard import method.
# @EXPORT contains list of functions to be exported into the caller namespace

use Exporter qw(import) ;
our @EXPORT  ; 
push @EXPORT, qw(MultipleSpeciesAlignment_2RefGenomes) ;
push @EXPORT, qw(Construct_Codon_hash) ;
push @EXPORT, qw(GetReferenceGenePositions) ;
push @EXPORT, qw(FastaCoreOverlapWithReference) ;
push @EXPORT, qw(ProkkaID_2_AssemblyID) ;
push @EXPORT, qw(LoadGeneSeqsForEachIndividual) ;
push @EXPORT, qw(RoaryCOGpos_2_MultispeciesAlleleTablePos) ;
push @EXPORT, qw(MakeOutgroupBaseHash) ;
push @EXPORT, qw(CalculatePositionDegeneracy) ;
push @EXPORT, qw(CollectSegSite_FixedSite_FunctionalEffects) ;
push @EXPORT, qw(CalculateSFS) ;
push @EXPORT, qw(GetPositionsSynNonsynBiallelicSites) ;
push @EXPORT, qw(CalculateNumberVariableSites) ;
push @EXPORT, qw(ConstructWindows) ;
push @EXPORT, qw(CalculatePairwiseLinkageStats_ByDistance);
push @EXPORT, qw(CalculatePairwiseLinkageStats_ByGeneCategory);
push @EXPORT, qw(CategorizeBySNPdensity_PrintSFSsummaries) ;
push @EXPORT, qw(SFS_SummStats) ;
push @EXPORT, qw(CalculateThetasBiAllelicByGene) ;
push @EXPORT, qw(BinLDbyDistForBoostrapResampling_WholeGenome) ;
push @EXPORT, qw(BinLDbyDistForBoostrapResampling_IntraInterGenic) ;
push @EXPORT, qw(BinLDbyDistForBoostrapResampling_ByGene) ;
push @EXPORT, qw(CharacterizeNonsenseMutations) ;
push @EXPORT, qw(ZerofoldFourfoldDivByQuantile); 
push @EXPORT, qw(PrintNumGenesPerQuantile) ;

sub MultipleSpeciesAlignment_2RefGenomes{
	# Creates a data structure containing the 2 different outgroup alleles for each 
	# position in reference, using the Allele Table created from a progressiveMauve 
	# alignment; this function can be easily modified to tolerate any file structure 
	
	my $Ref_multi_species_alignment = $_[0] ;
	my $wd = $_[1] ;
	
	my %MultiSpeciesAlign ; # @{$MultiSpeciesAlign{ reference position }} = ( allele1, allele2 ) ;
	open IN, "<$${Ref_multi_species_alignment}" ;
	while(<IN>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		my $pos = $line[0] ;
		my $allele1 = $line[2] ;
		my $allele2 = $line[3] ;
		if($_ !~ m/^RefPos/){
			if($pos !~ m/NA/){ # eq "NA" if there's DNA present in outgroup not present in Reference
				@{$MultiSpeciesAlign{$pos}} = ( $allele1, $allele2 ) ;
			}
		}else{
			# Double check you're reading in the right file; this part should be modified 
			# if using different file structures!
			if(scalar @line != 4){
				print "YOU MAY NOT HAVE THE RIGHT MULTISPECIESALIGNMENTFILE! EXITING...\n" ;
				exit ;
			}
		}
	}
	close IN ;
	open OUT, ">>${$wd}/QCfile.txt" ;
	print OUT "Positions in MultipleSpeciesAlignment: ", scalar keys %MultiSpeciesAlign, "\n" ;
	close OUT ;
	return(\%MultiSpeciesAlign) ;

}

sub GetReferenceGenePositions {
	# Creates 2 data structures that store positions of genes found in a gff file of a 
	# reference sequence; the 2 data structures allow accessing position information in
	# different ways
	
	my $Ref_gff_file = $_[0] ;
	my $wd = $_[1] ;
	
	# %Ref_gene_info is 3D hash: gene name, orientation of gene, start or end of gene 
	my %Ref_gene_info ; # $Ref_gene_info{gene_name}{"POSITION_FORWARD"}{"START"} = start
	# %Ref_gene_positions is 2D hash: start position of gene, end position of gene
	my %Ref_gene_positions ; 

	open OUT, ">>${$wd}/QCfile.txt" ;
	#my %New_Prokka_Gene_Boundaries ;
	open IN, "<$${Ref_gff_file}" ;
	while(<IN>){
		chomp $_ ;
		if($_ =~ m/CDS/){	# only look at lines in gff file describing CDS 
			my @line = split(/\t/, $_) ; 
			my $start = $line[3] ;
			my $end = $line[4] ;
			my $direction = $line[6] ;
			my @info = split(";", $line[8]) ;
				
			# $info[1] field contains "gene=NGO0002", where NGO0002 is name of gene	
			my @gene = split("=", $info[1]) ;
			if($gene[0] =~ m/gene/){
				my $gene_name = $gene[1] ;
				if($direction eq "+"){
					$Ref_gene_info{$gene_name}{"POSITION_FORWARD"}{"START"} = $start ;
					$Ref_gene_info{$gene_name}{"POSITION_FORWARD"}{"END"} = $end ;
				}elsif($direction eq "-"){ # coding sequence is on opposite strand
					$Ref_gene_info{$gene_name}{"POSITION_REVERSE"}{"START"} = $start ;
					$Ref_gene_info{$gene_name}{"POSITION_REVERSE"}{"END"} = $end ;
				}
				$Ref_gene_positions{$start}{$end} = $gene_name ;
			}
		}
	}
	print OUT "Number of genes in PROKKA annotation of reference: ", scalar keys %Ref_gene_info, "\n" ;
	close IN ;
	close OUT ;
	return( \%Ref_gene_info, \%Ref_gene_positions ) ;
}

sub FastaCoreOverlapWithReference {
	# Filters list of core genes to only those with information from reference
	
	my $fasta_core = $_[0] ;
	my $Ref_gene_info = $_[1] ;
	my $wd = $_[2] ;
	
	my @fasta_core_ref_overlap ;
	
	foreach my $file (@$fasta_core){
		chomp $file ;
		my @split = split(/\//, $file) ;
		my $gene = $split[$#split] ;
		$gene =~ s/\.fa\.aln// ;
		# Only accept core genes that also have information from reference
		if( exists $$Ref_gene_info{$gene} ){
				push @fasta_core_ref_overlap, $file ;
		}
	}
	open OUT, ">>${$wd}/QCfile.txt" ;	
	print OUT "Total Core Genes: ", scalar @$fasta_core, "\n" ;
	undef @$fasta_core ;
	print OUT "Total Core Genes overlapped with reference: ", scalar @fasta_core_ref_overlap, "\n" ;
	close OUT ;
	return(@fasta_core_ref_overlap) ;
}

sub ProkkaID_2_AssemblyID{
	# Creates 2 data structures that map Prokka identifiers back to the gene name and the 
	# assembly from which they came

	my $gene_presence_absence_file = $_[0] ;
	my $Individuals_to_include = $_[1] ;
	my $wd = $_[2] ;
	
	my %ProkkaID_2_AssemblyName ;
	my %ProkkaID_2_GeneName ;
	open IN, "<$${gene_presence_absence_file}" ;
	my %Assembly_names ; # Maps column number to assembly name, $Assembly_names{column} = name
	my $Sanity_checker ;
	while(<IN>){
		chomp $_ ;
		$_ =~ s/\r//g ; 				# remove hidden formatting at end of line
		my @line = split(/\",\"/, $_) ; # "," field seperator
		if($_ =~ m/^\"Gene\"/){			# Header line, get assembly names
			foreach my $index (14..$#line){
				$line[$index] =~ s/\"//g ;
				my $assembly_name = $line[$index] ;
				$Assembly_names{$index} = $assembly_name ;
				if(exists $$Individuals_to_include{$Assembly_names{$index}}){
					$Sanity_checker++ ;
				}
			}
		}else{
			my $gene = $line[0] ;
			$gene =~ s/\"//g ;
			foreach my $index (14..$#line){
				$line[$index] =~ s/\"//g ;
				my $ProkkaID = $line[$index] ;
				$ProkkaID_2_AssemblyName{ $ProkkaID } = $Assembly_names{$index} ;
				$ProkkaID_2_GeneName{ $ProkkaID } = $gene ;
			}
		}
	}
	close IN ;
	open OUT, ">>${$wd}/QCfile.txt" ;	
	print OUT "Note: I have assumed individual names start on column 15 of gene_presence_absence file\n" ;
	print OUT "Number individuals identified in gene_presence_absence.csv: ", $Sanity_checker, "\n" ;
	close OUT ;
	return(\%ProkkaID_2_AssemblyName, \%ProkkaID_2_GeneName) ;
}

sub LoadGeneSeqsForEachIndividual{
	# Load gene sequences from fasta files

	my $fasta_core_ref_overlap= $_[0] ;
	my $ProkkaID_2_AssemblyName = $_[1] ;
	my $Individuals_to_include = $_[2] ;
	my $Reference_ind = $_[3] ;
	my $wd = $_[4] ;

	my %core_loci_seqs  ;
	open OUT, ">>${$wd}/QCfile.txt" ;	

	foreach my $file (@$fasta_core_ref_overlap){
		chomp $file ;
		my @split = split(/\//, $file) ;
		my $gene = $split[$#split] ;
		$gene =~ s/\.fa\.aln// ;
		open IN, "<${file}" ;	# These should be standard fasta files
		my $id ;
		my $ind ;
		while(<IN>){
			chomp $_ ;
			if($_ =~ m/^>/){
				$_ =~ s/>//g ;
				$id = $_ ; # id is unique Prokka ID for this gene,
				$ind = $$ProkkaID_2_AssemblyName{$id} ;
				if( exists $$Individuals_to_include{ $ind } ){ # subsample individuals
					my $x = <IN> ; chomp $x ;
					$core_loci_seqs{$gene}{ $ind } = $x ;
					next ;
				}
			}else{
				if( exists $$Individuals_to_include{ $ind } ){
					$core_loci_seqs{$gene}{ $ind } = $core_loci_seqs{$gene}{ $ind }.${_} ;
				}
			}
		}
		close IN ;
	}
	print OUT "Number of core genes with sequencess: ", scalar keys %core_loci_seqs, "\n" ;
	my $sum1 = 0;
	foreach my $gene (keys %core_loci_seqs){
   		$sum1 +=  length($core_loci_seqs{$gene}{$$Reference_ind}) ;
	}
	print OUT "Total DNA analyzed (bp): ", $sum1, "\n" ;
	
	## See which gene have reading frame violations
	## This could occur from frameshift mutations or, more likely, alignment errors
	my %RF_viol ; # tallies the number of reading frame violations per gene
	foreach my $gene (sort {$a cmp $b} keys %core_loci_seqs){
		foreach my $ind (sort {$a cmp $b} keys %{$core_loci_seqs{$gene}}){
			if($ind ne $$Reference_ind){
				my $tmp = $core_loci_seqs{$gene}{$ind} ;
				$tmp =~ s/-//g ;
				for (my $site = 0; $site <=length($core_loci_seqs{$gene}{$ind})-3; $site +=3 ){
					my $triplet = substr($core_loci_seqs{$gene}{$ind}, $site, 3) ;
					if( ($triplet =~ tr/[AGCTN]//) == 3 || ($triplet =~ tr/-//) == 3){ 
						## do nothing; reading frame intact, no frameshift mutations/gaps
					}else{
						$RF_viol{$gene}++ ;
					}
				}
			}
		}
	}
	print OUT "GENES WITH READING FRAME VIOLATIONS (codons containing nucleotides AND gaps) : ", scalar keys %RF_viol, "\n" ;
	foreach my $gene (sort {$a cmp $b} keys %core_loci_seqs){
		if(scalar keys %{$core_loci_seqs{$gene}} != scalar keys %$Individuals_to_include){
			print "WARNING: gene ${gene} isn't a core gene! May happen b/c gene duplicate: ind assigned seq name compesed of 2 file names\n" ;
		}
	}
	close OUT ;
	return(\%core_loci_seqs) ;
}

sub RoaryCOGpos_2_MultispeciesAlleleTablePos{

	my $core_loci_seqs = $_[0] ;
	my $Ref_gene_info = $_[1] ;
	my $Reference_ind = $_[2] ;
	my $MultiSpeciesAlign = $_[3] ;
	my $wd = $_[4] ;
	open OUT, ">>${$wd}/QCfile.txt" ;	
	
	my %COG_AlnPos_2_WGref_AlnPos ; # this hash maps positions in COG alignments to positions in the reference genome

	foreach my $gene (keys %$core_loci_seqs){
		my $MultiSpeciesAlign_pos_tracker = 0 ; # Keeps track of gaps in reference seq from alignment procedure
		foreach my $site (0 .. length($$core_loci_seqs{$gene}{$$Reference_ind})-1){
			if(substr($$core_loci_seqs{$gene}{$$Reference_ind}, $site, 1) eq "-"){
				$COG_AlnPos_2_WGref_AlnPos{$gene}{$site} = "NA" ; ## NA if there's a gap in Reference_ind in COG alignment
				$MultiSpeciesAlign_pos_tracker++ ;
			}else{
				# base corresponds to position in reference genome
				# use MultiSpeciesAlign_pos_tracker to know which reference gene site you're at
				if( exists $$Ref_gene_info{$gene}{"POSITION_FORWARD"} ){
					$COG_AlnPos_2_WGref_AlnPos{$gene}{$site} = $$Ref_gene_info{$gene}{"POSITION_FORWARD"}{"START"} + ($site-$MultiSpeciesAlign_pos_tracker) ;
				}elsif( exists $$Ref_gene_info{$gene}{"POSITION_REVERSE"} ){
					$COG_AlnPos_2_WGref_AlnPos{$gene}{$site} = $$Ref_gene_info{$gene}{"POSITION_REVERSE"}{"END"} - ($site-$MultiSpeciesAlign_pos_tracker) ;
				}
			}
		}
	}
	
	
	## Check if reference base in COG alignment is the same as base in MultiSpecies allele table, else delete
	my %Bad_genes ;
	my %Good_genes ;
	foreach my $gene (keys %COG_AlnPos_2_WGref_AlnPos ){
		my $flag = 0 ;
		foreach my $site ( sort {$a <=> $b} keys %{$COG_AlnPos_2_WGref_AlnPos{$gene}} ){
			if( !exists $$MultiSpeciesAlign{ $COG_AlnPos_2_WGref_AlnPos{$gene}{$site} } && $COG_AlnPos_2_WGref_AlnPos{$gene}{$site} eq "NA"){ 
				# exists fails probably because $COG_AlnPos_2_WGref_AlnPos{$gene}{$site} == "NA"
					next ;
			}else{
				my $base = ${$$MultiSpeciesAlign{ $COG_AlnPos_2_WGref_AlnPos{$gene}{$site} }}[0] ;
				if( exists $$Ref_gene_info{$gene}{"POSITION_REVERSE"} ){ # take reverse complement
					$base =~ tr/AGCT/TCGA/ ;
				}
				if( substr($$core_loci_seqs{$gene}{$$Reference_ind}, $site, 1) ne $base ){
					$flag++ ;
				}
			}
		}
		if($flag){
			$Bad_genes{$gene} = $flag ;
		}else{
			$Good_genes{$gene}++ ;
		}
	}
	print OUT "\nCOGS where Reference gene didn't match nucleotide in MultiSpeciesAlignment (the orientation of these genes may have been misspecified)\n" ;
	print OUT "Gene\tNumMismatchSites\tRefPos's\tOrientation\n" ;
	foreach my $gene ( keys %Bad_genes ){
		print OUT $gene, "\t", $Bad_genes{$gene}, "\t" ;
		if(exists $$Ref_gene_info{$gene}{"POSITION_FORWARD"}){
			print OUT $$Ref_gene_info{$gene}{"POSITION_FORWARD"}{"START"}, ",", $$Ref_gene_info{$gene}{"POSITION_FORWARD"}{"END"}, "\tFORWARD\n" ; 
		}elsif(exists $$Ref_gene_info{$gene}{"POSITION_REVERSE"}){
			print OUT $$Ref_gene_info{$gene}{"POSITION_REVERSE"}{"START"}, ",", $$Ref_gene_info{$gene}{"POSITION_REVERSE"}{"END"}, "\tREVERSE\n" ; 
		}
	}
	print OUT "\n" ;
	print OUT "good genes with reference info across COG and MultiSpeciesAlgiments:\t", scalar keys %Good_genes, "\n" ;
	print OUT "bad genes with BAD reference info across COG and MultiSpeciesAlgiments:\t", scalar keys %Bad_genes, "\n" ;
	print OUT "deleting these bad ones...\n" ;
	foreach my $gene (keys %Bad_genes){
		delete $$core_loci_seqs{$gene} ;
	}
	close OUT ;
	return(\%COG_AlnPos_2_WGref_AlnPos) ;
}

sub MakeOutgroupBaseHash {

	my $core_loci_seqs = $_[0] ;
	my $COG_AlnPos_2_WGref_AlnPos = $_[1] ;
	my $MultiSpeciesAlign = $_[2] ;
	my $Ref_gene_info = $_[3] ;
	my $Reference_ind = $_[4] ;
	my $AlignFile_index_ForReference = $_[5] ;	

	my %RefBase ;
	# Refbase for genes
	foreach my $gene (keys %$core_loci_seqs){
		foreach my $site (0 .. length($$core_loci_seqs{$gene}{$$Reference_ind})-1){
			if( $$COG_AlnPos_2_WGref_AlnPos{$gene}{$site} ne "NA" ){## NA if there's a gap in Reference_ind in COG alignment; insertion has no reference 
				my $base = ${$$MultiSpeciesAlign{ $$COG_AlnPos_2_WGref_AlnPos{$gene}{$site} }}[ $$AlignFile_index_ForReference ] ;
				$base = uc($base) ;
				if($base =~ m/[AGCT]/){ # this excludes outgroup sites that are N's or -'s
					if(exists $$Ref_gene_info{$gene}{"POSITION_REVERSE"}){
						$base =~ tr/AGCT/TCGA/ ; # Ref seq may have reverse complemented seq. of gene
					}
					$RefBase{$gene}{$site} = $base ;
				}
			}
		}
	}	
	return( \%RefBase ) ;
}

sub CalculatePositionDegeneracy{
	# Using the reference sequence, we can calculate the n-fold degeneracy of each position
	# n-fold degeneracy -> n nucleotides specify same amino acid (a.a.)
	my $core_loci_seqs = $_[0] ;
	my $Reference_ind = $_[1] ;
	my $Codon_Hash = $_[2] ;

	my %Degeneracy ;
	foreach my $gene (keys %$core_loci_seqs){
		foreach my $site (0 .. length($$core_loci_seqs{$gene}{$$Reference_ind})-1){
			if(substr($$core_loci_seqs{$gene}{$$Reference_ind}, $site, 1) eq "-"){
				# skip if at gap, current algorithm aligned with codon-aware PRANK or PAGAN
				# makes gaps usually occur in triplets
				next;
			}else{
				my $codon_position = ($site+1)%3 ; # add 1 to $site because it starts at 0
				my $flank_site_1 ;
				my $flank_site_2 ;
				if($codon_position == 1){ 		# site is in 1st position of codon
					$flank_site_1 = $site+1 ;
					$flank_site_2 = $site+2 ;
				}elsif($codon_position == 2){ 	# site is in 2nd position of codon
					$flank_site_1 = $site-1 ;
					$flank_site_2 = $site+1 ;
				}elsif($codon_position == 0){	 # site is in 3rd position of codon
					$flank_site_1 = $site-2 ;
					$flank_site_2 = $site-1 ;
				}
				# store reference bases that flank the current position
				my $fb1 = substr($$core_loci_seqs{$gene}{$$Reference_ind}, $flank_site_1, 1) ;
				my $fb2 = substr($$core_loci_seqs{$gene}{$$Reference_ind}, $flank_site_2, 1) ;
			
				# with flanking reference bases, cycle through 4 possible nucleotides
				# record how many a.a.'s specified using %Codon_Hash
				my %codon_degeneracy ; 
				foreach my $base ( "A", "G", "C", "T" ){
					if($codon_position == 1){ #site in 1st position
						$codon_degeneracy{ $$Codon_Hash{$base.$fb1.$fb2} } ++ ; 
						# where $$Codon_Hash{$base.$fb1.$fb2} is an a.a.
					}elsif($codon_position == 2){ #site in 2nd position
						$codon_degeneracy{ $$Codon_Hash{$fb1.$base.$fb2} } ++ ;
					}elsif($codon_position == 0){ #site in 3rd position
						$codon_degeneracy{ $$Codon_Hash{$fb1.$fb2.$base} } ++ ;
					}
				}
				
				if((scalar keys %codon_degeneracy) == 1){
					$Degeneracy{$gene}{$site} = "4D" ; # or 4-fold Degenerate
				}elsif((scalar keys %codon_degeneracy) == 2){
					 # two a.a.'s encoded, could be 2-fold or 3-fold degenerate
					 # temporarily record lower frequency one to check
					(my $LowerFreqAA) = sort {$codon_degeneracy{$a} <=> $codon_degeneracy{$b}} keys %codon_degeneracy ;
					if( $codon_degeneracy{$LowerFreqAA} == 2){ # a.a. specified twice, therefore other also twice for four total
						$Degeneracy{$gene}{$site} = "2D" ;
					}elsif($codon_degeneracy{$LowerFreqAA} == 1){ # a.a. specified once, therefore other three times for four total
						$Degeneracy{$gene}{$site} = "3D" ;
					}
				}elsif((scalar keys %codon_degeneracy) == 3){ # always 2-fold degenerate
					$Degeneracy{$gene}{$site} = "2D" ;
				}elsif((scalar keys %codon_degeneracy) == 4){ # always 0-fold degenerate
					$Degeneracy{$gene}{$site} = "0D" ;
				}
			}
		}
	}
	return ( \%Degeneracy ) ;
}

sub CollectSegSite_FixedSite_FunctionalEffects{

	my $core_loci_seqs = $_[0] ;
	my $Reference_ind = $_[1] ;
	my $Codon_Hash = $_[2] ;
	my $RefBase = $_[3] ;
	my $Degeneracy = $_[4] ;
	my $wd = $_[5] ;
	open OUT, ">>${$wd}/QCfile.txt" ;	

	my %core_ref_biallelic_segsites ; # $core_ref_biallelic_segsites_HR{ gene }{ seg_site_index } = $site ; where seg_site_index ranges from 0 to total number of SNPs
	my %core_ref_fixed_diffs ; # $core_ref_fixed_diffs_HR{ gene }{ seg_site_index } = $site ; where seg_site_index goes from 0 to to total number of fixed differences
	my %core_ref_fixed_same ;
	my %core_ref_NumUngapped_sites ; # $core_ref_NumUngapped_sites_HR{ gene } = number ; or the number of ungapped positions for the gene alignment
	my %core_ref_NumGapped_sites ;
	my %Functional_Effect_BiAllelic ; # $Functional_Effect_BiAllelic_HR{ gene }{ site } =  "SYNONYMOUS" or "NONSYNONYMOUS:MISSENSE" or "NONSYNONYMOUS:NONSENSE"  ; 
	my %Functional_Effect_MultiAllelic ; # $Functional_Effect_MultiAllelic_HR{ gene }{ site }{ base1 }{ base2 } =  "SYNONYMOUS" or "NONSYNONYMOUS:MISSENSE" or "NONSYNONYMOUS:NONSENSE"; $base1 and $base2 include all permutations, e.g. b1=A/b2=T and b1=T/b2=A, plug in ReferenceBase for $base1
	my %Functional_Effect_Fixed ; # $Functional_Effect_Fixed_HR{ gene }{ site } =  "SYNONYMOUS" or "NONSYNONYMOUS:MISSENSE" or "NONSYNONYMOUS:NONSENSE"  ;
	my %Biallelic_sites_freqs  ; # $Biallelic_sites_freqs_HR{$gene}{ site }{ base } = frequency
	my %Multiallelic_sites_freqs ; # $Multiallelic_sites_freqs_HR{ gene }{ site }{ base } = frequency
	my %Monomorphic_sites ; # $Monomorphic_sites{ gene }{ site } = base (A,G,C,T) ;
	my %ZeroFourfoldSegSites ; # $ZeroFourfoldSegSites_HR{ gene }{"0D" or "4D"} = Number of sites ;
	my %ZeroFourfoldSites_UnGapped ; # $ZeroFourfoldSites_UnGapped_HR{ gene }{"0D" or "4D"} = Number of sites ;

	my $SitesNotConsideredCuzIndelsOrNs = 0 ;
	my $SitesNotConsideredCuzNoRefBase = 0 ;
	# Go through each site of each gene, count number of bases/alleles present
	foreach my $gene (keys %$core_loci_seqs){
		my $seg_site_index = 0 ;
		my $fix_site_index = 0 ;
		$core_ref_NumUngapped_sites{$gene} = 0 ;
		$core_ref_NumGapped_sites{$gene} = 0 ;
		$ZeroFourfoldSites_UnGapped{$gene}{"0D"} = 0 ;
		$ZeroFourfoldSites_UnGapped{$gene}{"4D"} = 0 ;
		$ZeroFourfoldSegSites{$gene}{"0D"} = 0 ;
		$ZeroFourfoldSegSites{$gene}{"4D"} = 0 ;
		foreach my $site (0 .. length($$core_loci_seqs{$gene}{$$Reference_ind})-1){
			# Excludes sites with no reference info, can occur because
				# 1.) reference sequence had gap in COG alignment with other sequences
				# 2.) MultiSpecies Allele table had gap, perhaps because sequence wasn't present in outgroup
			if(exists $$RefBase{$gene}{$site}){
				my %poly_bases ; 	# Keeps track of bases across individuals
				my $indel = 0 ;		# Checks for presence of indel/gap
				my $N = 0 ;			# Checks for presence of ambiguous nucleotide "N"
				foreach my $ind (keys %{$$core_loci_seqs{$gene}}){
					if($ind ne $$Reference_ind){ # Ignore reference individual, we only care about non-reference sequences
						my $base = uc(substr($$core_loci_seqs{$gene}{$ind}, $site, 1)) ;
						if( $base =~ m/[AGCT]/){
							$poly_bases{$base} ++ ;
						}elsif( $base =~ m/-/ ){
							$indel ++ ;
						}elsif( $base =~ m/N/ ){
							$N ++ ;
						}
					}
				}
				if( $indel==0 && $N==0){	# Only consider sites without gaps or "N"
					$core_ref_NumUngapped_sites{$gene}++ ;
					if(exists $$Degeneracy{$gene}{$site}){
						if ($$Degeneracy{$gene}{$site} eq "0D"){
							$ZeroFourfoldSites_UnGapped{$gene}{"0D"}++ ;	
						}elsif($$Degeneracy{$gene}{$site} eq "4D"){
							$ZeroFourfoldSites_UnGapped{$gene}{"4D"}++ ;	
						}
					}
					if((scalar keys %poly_bases > 1)){
						if(exists $$Degeneracy{$gene}{$site}){
							if ($$Degeneracy{$gene}{$site} eq "0D"){
								$ZeroFourfoldSegSites{$gene}{"0D"}++ ;	
							}elsif($$Degeneracy{$gene}{$site} eq "4D"){
								$ZeroFourfoldSegSites{$gene}{"4D"}++ ;	
							}
						}	
					}
					
					my $codon_position ; 	# position of $site within the codon
					my $fb1 ;				# flanking base 1
					my $fb2 ;				# flanking base 2
					# $codon_position, $fb1 and $fb2 are references! dereference with extra $
					($codon_position, $fb1, $fb2) = GetCodonPosFlankingBasesInReference($core_loci_seqs, \$gene, \$site, $Reference_ind) ;
					
					# We need to deal with Invariant, BiAllelic, and MultiAllelic sites differently 
					##############
					# IS THIS SITE MULTIALLELIC?
					##############
					if( (scalar keys %poly_bases >= 3) ){
						my @MultiAllelicBases ;
						foreach my $base ( keys %poly_bases ){
							$Multiallelic_sites_freqs{$gene}{$site}{$base} = $poly_bases{$base} ;
							push @MultiAllelicBases, $base ;
						}
						if($$fb1 ne "-" && $$fb2 ne "-"){
							# cycle through all ordered permutations of bases
							foreach my $i1 ( 0..$#MultiAllelicBases ){
								foreach my $i2 ( 0..$#MultiAllelicBases ){
									if($i1 != $i2){
										#construct biallelic codon to infer functional effect
										my $base1 = $MultiAllelicBases[$i1] ;
										my $base2 = $MultiAllelicBases[$i2] ;
										my $c1 ; my $c2 ;								
										if($$codon_position == 1){ #site in 1st position
											$c1 = $base1.$$fb1.$$fb2 ;
											$c2 = $base2.$$fb1.$$fb2 ;
										}elsif($$codon_position == 2){ #site in 2nd position
											$c1 = $$fb1.$base1.$$fb2 ;
											$c2 = $$fb1.$base2.$$fb2 ;						
										}elsif($$codon_position == 0){ #site in 3rd position
											$c1 = $$fb1.$$fb2.$base1 ;
											$c2 = $$fb1.$$fb2.$base2 ;						
										}
										
										#if biallilic codon produces same a.a. -> synonymous
										if(!exists $$Codon_Hash{$c1} || !exists $$Codon_Hash{$c2}){
											print OUT "WARNING: ERROR IN INTERPRETING CODONS: ", $c1, " " ,$c2, "\t", "in gene: ", $gene, " ", $site+1, "\n" ;
										}else{
											if( $$Codon_Hash{$c1} eq $$Codon_Hash{$c2} ){
												$Functional_Effect_MultiAllelic{$gene}{$site}{$base1}{$base2} = "SYNONYMOUS" ;
											}else{
												$Functional_Effect_MultiAllelic{$gene}{$site}{$base1}{$base2} = "NONSYNONYMOUS:MISSENSE" ;
												## is this polymorphic for a LossOfFunction mutation?
												if( $$Codon_Hash{$c1} eq "Stop" || $$Codon_Hash{$c2} eq "Stop"){
													$Functional_Effect_MultiAllelic{$gene}{$site}{$base1}{$base2} = "NONSYNONYMOUS:NONSENSE" ;
												}
											}
										}								
									}
								}
							}
						}														
					}
					##############
					# IS THIS SITE BIALLELIC?
					##############
					if( (scalar keys %poly_bases == 2) ){
						if($$fb1 ne "-" && $$fb2 ne "-"){
							$core_ref_biallelic_segsites{$gene}{$seg_site_index} = $site  ;
							$seg_site_index++ ;
							#construct biallelic codons, only 2 according to if(condition) above
							my @biallelic_codon ;
							foreach my $base ( keys %poly_bases ){
								$Biallelic_sites_freqs{$gene}{$site}{$base} = $poly_bases{$base} ;
								if($$codon_position == 1){ #site in 1st position
									push @biallelic_codon, $base.$$fb1.$$fb2 ;
								}elsif($$codon_position == 2){ #site in 2nd position
									push @biallelic_codon, $$fb1.$base.$$fb2 ;
								}elsif($$codon_position == 0){ #site in 3rd position
									push @biallelic_codon, $$fb1.$$fb2.$base ;
								}
							}
							#if biallilic codon produces same a.a. -> synonymous
							if(!exists $$Codon_Hash{$biallelic_codon[0]} || !exists $$Codon_Hash{$biallelic_codon[1]}){
								print OUT "WARNING: ERROR IN INTERPRETING CODONS: ", "@biallelic_codon", , "\t", "in gene: ", $gene, " ", $site+1, "\n" ;
							}else{
								if( $$Codon_Hash{$biallelic_codon[0]} eq $$Codon_Hash{$biallelic_codon[1]} ){
									$Functional_Effect_BiAllelic{$gene}{$site} = "SYNONYMOUS" ;
								}else{
									$Functional_Effect_BiAllelic{$gene}{$site} = "NONSYNONYMOUS:MISSENSE" ;
									## is this polymorphic for a LossOfFunction mutation?, codons don't equal
									if( $$Codon_Hash{$biallelic_codon[0]} eq "Stop" || $$Codon_Hash{$biallelic_codon[1]} eq "Stop"){
										$Functional_Effect_BiAllelic{$gene}{$site} = "NONSYNONYMOUS:NONSENSE" ;
									}
								}
							}
						}
					}
					##############
					# IS THIS SITE FIXED?
					##############
					if( (scalar keys %poly_bases == 1) ){
						(my $base) = keys %poly_bases ;
						$Monomorphic_sites{$gene}{$site} = $base ;
						if( ($base ne $$RefBase{$gene}{$site}) ){
							if($$fb1 ne "-" && $$fb2 ne "-"){
								$core_ref_fixed_diffs{$gene}{$fix_site_index} = $site  ;
								$fix_site_index++ ;
								#construct biallelic codons, only 2 according to if(condition) above
								my @biallelic_codon ;				
								if($$codon_position == 1){ #site in 1st position
									@biallelic_codon = ($base.$$fb1.$$fb2, $$RefBase{$gene}{$site}.$$fb1.$$fb2) ;
								}elsif($$codon_position == 2){ #site in 2nd position
									@biallelic_codon = ($$fb1.$base.$$fb2, $$fb1.$$RefBase{$gene}{$site}.$$fb2) ;
								}elsif($$codon_position == 0){ #site in 3rd position
									@biallelic_codon = ($$fb1.$$fb2.$base, $$fb1.$$fb2.$$RefBase{$gene}{$site}) ;
								}
								#if biallilic codon produces same a.a. -> synonymous
								if(!exists $$Codon_Hash{$biallelic_codon[0]} || !exists $$Codon_Hash{$biallelic_codon[1]}){
									print OUT "WARNING: ERROR IN INTERPRETING CODONS: ", "@biallelic_codon", , "\t", "in gene: ", $gene, " ", $site+1, "\n" ;
								}else{
									if( $$Codon_Hash{$biallelic_codon[0]} eq $$Codon_Hash{$biallelic_codon[1]} ){
										$Functional_Effect_Fixed{$gene}{$site} = "SYNONYMOUS" ;
									}else{
										$Functional_Effect_Fixed{$gene}{$site} = "NONSYNONYMOUS:MISSENSE" ;
										## is this polymorphic for a LossOfFunction mutation?, codons don't equal
										if( $$Codon_Hash{$biallelic_codon[0]} eq "Stop" || $$Codon_Hash{$biallelic_codon[1]} eq "Stop"){
											$Functional_Effect_Fixed{$gene}{$site} = "NONSYNONYMOUS:NONSENSE" ;
										}									
									}
								}
							}
						}
					}
				}else{
					$core_ref_NumGapped_sites{$gene}++ ;
					$SitesNotConsideredCuzIndelsOrNs ++ ;
				}
			}else{
				$SitesNotConsideredCuzNoRefBase++ ;
			}
		}	
	}
	print OUT "Sites ignored because of Indels or Ns: ", $SitesNotConsideredCuzIndelsOrNs, "\n" ;
	print OUT "Sites ignored because no Reference Base: ", $SitesNotConsideredCuzNoRefBase, "\n" ;
	my $fixsum_syn = 0 ;
	my $fixsum_nonsyn = 0 ;
	foreach my $gene (keys %Functional_Effect_Fixed){
		foreach my $site (keys %{$Functional_Effect_Fixed{$gene}}){
			if( $Functional_Effect_Fixed{$gene}{$site} eq "SYNONYMOUS" ){
				$fixsum_syn++ ;
			}else{
				$fixsum_nonsyn++ ;
			}
		
		}
	}
	print OUT "Number Synonymous Fixed mutations: ", $fixsum_syn, "\n" ;
	print OUT "Number Nonsynonymous Fixed mutations: ", $fixsum_nonsyn, "\n" ;
	close OUT ;
	
	open OUT, ">${$wd}/AlignmentStats_NumGappedSites.txt" ;	
	print OUT "Gene", "\t", "TotalNumSites", "\t", "NumUngappedSites", "\t", "NumGappedSites", "\n" ;
	foreach my $gene (keys %$core_loci_seqs){
		print OUT $gene, "\t", length($$core_loci_seqs{$gene}{$$Reference_ind}), "\t", $core_ref_NumUngapped_sites{$gene}, "\t", $core_ref_NumGapped_sites{$gene}, "\n" ;
	}
	close OUT ;
	return( \%core_ref_biallelic_segsites, \%core_ref_fixed_diffs, \%core_ref_NumUngapped_sites, \%Functional_Effect_BiAllelic, \%Functional_Effect_MultiAllelic, \%Functional_Effect_Fixed, \%Biallelic_sites_freqs, \%Multiallelic_sites_freqs, \%Monomorphic_sites, \%ZeroFourfoldSegSites, \%ZeroFourfoldSites_UnGapped ) ;

}

sub GetCodonPosFlankingBasesInReference{
	my $core_loci_seqs = $_[0] ;
	my $gene = $_[1] ;
	my $site = $_[2] ;
	my $Reference_ind = $_[3] ;

	# find codon position and positions of other sites in same codon
	my $codon_position = ($$site+1)%3 ; #add 1 to $site b/c starts at 0
	my $flank_site_1 ;
	my $flank_site_2 ;
	if($codon_position == 1){ #site in 1st position
		$flank_site_1 = $$site+1 ;
		$flank_site_2 = $$site+2 ;
	}elsif($codon_position == 2){ #site in 2nd position
		$flank_site_1 = $$site-1 ;
		$flank_site_2 = $$site+1 ;
	}elsif($codon_position == 0){ #site in 3rd position
		$flank_site_1 = $$site-2 ;
		$flank_site_2 = $$site-1 ;
	}
	#USE THE REFERENCE SEQUENCE FOR FLANKING BASES
	my $fb1 = substr($$core_loci_seqs{$$gene}{$$Reference_ind}, $flank_site_1, 1) ;
	my $fb2 = substr($$core_loci_seqs{$$gene}{$$Reference_ind}, $flank_site_2, 1) ;

	return(\$codon_position, \$fb1, \$fb2) ;
}

sub CharacterizeNonsenseMutations{
	# This function quantifies how many genes have nonsense mutations and if so, where they occur within the gene

	my $core_loci_seqs = $_[0] ;
	my $Functional_Effect_BiAllelic = $_[1] ; # $Functional_Effect_BiAllelic_HR{$gene}{$site} =  "SYNONYMOUS"/"NONSYNONYMOUS:MISSENSE"/"NONSYNONYMOUS:NONSENSE"  ; 
	my $Functional_Effect_MultiAllelic = $_[2] ; # $Functional_Effect_MultiAllelic_HR{$gene}{$site}{$base1}{$base2} =  "SYNONYMOUS"/"NONSYNONYMOUS:MISSENSE"/"NONSYNONYMOUS:NONSENSE"  ; $base1 and $base2 include all permutations, e.g. $b1=A/$b2=T and $b1=T/$b2=A, plug in RefBase for $base1
	my $Reference_ind = $_[3] ;
	my $wd = $_[4] ;

	my %NonsenseMutsPerGene ;
	open OUT, ">${$wd}/PositionsOfNonsenseMutationsPerGene.txt" ;
	print OUT "BIALLELIC SITES\n" ;
	print OUT "Gene", "\t", "Site", "\t", "GeneLength", "\n" ;

	### THIS WAS FOR PRINTING POSITION OF NONSENSE MUTATIONS
	foreach my $gene ( sort{$a cmp $b} keys %{$Functional_Effect_BiAllelic} ){
		foreach my $site(sort {$a <=> $b} keys %{$$Functional_Effect_BiAllelic{$gene}} ){
			if( $$Functional_Effect_BiAllelic{$gene}{$site} eq "NONSYNONYMOUS:NONSENSE"){
				$NonsenseMutsPerGene{$gene}++ ;
				print OUT $gene, "\t", $site, "\t", length($$core_loci_seqs{$gene}{$$Reference_ind}), "\n" ;
			}
		}
	}
	print OUT "\n\n" ;
	print OUT "MULTIALLELIC SITES\n" ;
	print OUT "Gene", "\t", "Site", "\t", "GeneLength", "\n" ;
	foreach my $gene ( sort{$a cmp $b} keys %{$Functional_Effect_MultiAllelic} ){
		foreach my $site( sort{$a <=> $b} keys %{$$Functional_Effect_MultiAllelic{$gene}} ){
			my $tmp = 0 ;
			foreach my $b1 ( keys %{$$Functional_Effect_MultiAllelic{$gene}{$site}} ){
				foreach my $b2 ( keys %{$$Functional_Effect_MultiAllelic{$gene}{$site}{$b1}} ){
					if( $$Functional_Effect_MultiAllelic{$gene}{$site}{$b1}{$b2} eq "NONSYNONYMOUS:NONSENSE" ){
						$tmp++ ;
					}
				}
			}
			if( $tmp ){
				$NonsenseMutsPerGene{$gene}++ ;
				print OUT $gene, "\t", $site, "\t", length($$core_loci_seqs{$gene}{$$Reference_ind}), "\n" ;
			}
		}
	}
	close OUT; 
	open OUT, ">${$wd}/NumNonsenseMutsPerGene.txt" ;
	print OUT "Gene\tNumber\n" ;
	foreach my $gene (sort {$a cmp $b} keys %NonsenseMutsPerGene){
		print OUT $gene, "\t", $NonsenseMutsPerGene{$gene}, "\n" ;	
	}
	close OUT ;
	return(\%NonsenseMutsPerGene) ;
}

sub CalculateSFS{
	# This function calculates and prints the mutation site-frequency spectrum
	
	my $Biallelic_sites_freqs = $_[0] ;
	my $Functional_Effect = $_[1] ;
	my $RefBase_HASHREF = $_[2] ;
	my $num_ind_NoRef = $_[3] ;
	my $wd = $_[4] ;

	my %SFS_polarized_w_Ref ;
	foreach my $gene( keys %$Biallelic_sites_freqs ){
		foreach my $site ( keys %{$$Biallelic_sites_freqs{$gene}} ){
			my $ref_base = $$RefBase_HASHREF{$gene}{$site} ;
			## are one of the bases in the reference?
			if( exists $$Biallelic_sites_freqs{$gene}{$site}{$ref_base} ){
				foreach my $base ( keys %{$$Biallelic_sites_freqs{$gene}{$site}} ){
					if($base ne $ref_base){
						my $freq = $$Biallelic_sites_freqs{$gene}{$site}{$base} ;
						if( $$Functional_Effect{$gene}{$site} eq "SYNONYMOUS" ){
							$SFS_polarized_w_Ref{"SYN"}{ $freq }++ ;
						}else{
							$SFS_polarized_w_Ref{"NONSYN"}{ $freq }++ ;
							if($$Functional_Effect{$gene}{$site} eq "NONSYNONYMOUS:NONSENSE"){
								$SFS_polarized_w_Ref{"NONSENSE"}{ $freq }++ ;
							}
						}
					}
				}
			}
		}
	}

	open OUT, ">${$wd}/SFS.txt" ;
	print OUT "Freq", "\t", "Num_Syn", "\t", "Num_NonSyn", "\t", "Num_NonSyn_StopCodon", "\n" ;
	foreach my $freq (1 .. $$num_ind_NoRef-1){
		print OUT $freq, "\t" ; 
		
		if( exists $SFS_polarized_w_Ref{"SYN"}{$freq} ){
			print OUT $SFS_polarized_w_Ref{"SYN"}{$freq}, "\t" ;
		}else{
			print OUT "0", "\t"
		}
		if( exists $SFS_polarized_w_Ref{"NONSYN"}{$freq} ){
			print OUT $SFS_polarized_w_Ref{"NONSYN"}{$freq}, "\t" ;
		}else{
			print OUT "0", "\t"
		}
		if(exists $SFS_polarized_w_Ref{"NONSENSE"}{$freq}){
			print OUT $SFS_polarized_w_Ref{"NONSENSE"}{$freq}, "\n" ;
		}else{
			print OUT "0", "\n" ;
		}
	}
	close OUT ;


	return( \%SFS_polarized_w_Ref ) ;
}

sub GetPositionsSynNonsynBiallelicSites{

	my $Ref_gene_positions_HR = $_[0] ;
	my $Ref_gene_info_HR = $_[1] ;
	my $RefBase_HR = $_[2] ;
	my $core_ref_biallelic_segsites = $_[3] ;
	my $Biallelic_sites_freqs = $_[4] ;
	my $num_ind_NoRef = $_[5] ;
	my $Functional_Effect = $_[6] ;
	my $min_AF = $_[7] ;
	my $max_AF = $_[8] ;
	my $wd = $_[9] ;
	open OUT, ">>${$wd}/QCfile.txt" ;	
	
	my %Biallelic_Syn_SNP_Ref_Pos ;
	my %Biallelic_Syn_SNP_Alignment_Coord ;
	my %Biallelic_NonSyn_SNP_Ref_Pos ;
	my %Biallelic_NonSyn_SNP_Alignment_Coord ;
	
	my $count_syn = 0 ;
	my $count_nonsyn = 0 ;
	foreach my $start (sort{$a <=> $b} keys %$Ref_gene_positions_HR){
		foreach my $end (sort{$a <=> $b} keys %{$$Ref_gene_positions_HR{$start}}){
			my $gene = $$Ref_gene_positions_HR{$start}{$end} ;
			if( exists $$core_ref_biallelic_segsites{$gene} ){ # defined only for sites in which (exists $$RefBase_HR{$gene}{$site})
				if(exists $$Ref_gene_info_HR{$gene}{"POSITION_FORWARD"}){
					#go through in ascending order
					foreach my $seg_site_index (sort{$a <=> $b} keys %{$$core_ref_biallelic_segsites{$gene}}){
						my $site = $$core_ref_biallelic_segsites{$gene}{$seg_site_index} ;
						## are one of the bases in the reference?
						my $refbase = $$RefBase_HR{$gene}{$site} ;
						if( exists $$Biallelic_sites_freqs{$gene}{$site}{$refbase} && $refbase ne "-" ){
							my $alt_base ;
							foreach my $b (keys %{$$Biallelic_sites_freqs{$gene}{$site}} ){
								if($b ne $refbase){
									$alt_base = $b ;
								}
							}
						
							my $af = $$Biallelic_sites_freqs{$gene}{$site}{$alt_base}/$$num_ind_NoRef ;
						
							if( $af>=$$min_AF && $af<=$$max_AF ){
								if( $$Functional_Effect{$gene}{$site} eq "SYNONYMOUS" ){
									$Biallelic_Syn_SNP_Ref_Pos{$count_syn} = $start+$site ;
									@{$Biallelic_Syn_SNP_Alignment_Coord{$count_syn}} = ($gene, $seg_site_index) ;
									$count_syn++ ;
								}
								if( $$Functional_Effect{$gene}{$site} =~ m/NONSYNONYMOUS/ ){
									$Biallelic_NonSyn_SNP_Ref_Pos{$count_nonsyn} = $start+$site ;
									@{$Biallelic_NonSyn_SNP_Alignment_Coord{$count_nonsyn}} = ($gene, $seg_site_index) ;
									$count_nonsyn++ ;
								}
							}
						}
					}
				}
				if(exists $$Ref_gene_info_HR{$gene}{"POSITION_REVERSE"}){
					#go through in descending order to cycle through SNPs along the genome from LEFT (3') to RIGHT (5')
					foreach my $seg_site_index (sort{$b <=> $a} keys %{$$core_ref_biallelic_segsites{$gene}}){
						my $site = $$core_ref_biallelic_segsites{$gene}{$seg_site_index} ;
						## are one of the bases in the reference?
						my $refbase = $$RefBase_HR{$gene}{$site} ;
						if( exists $$Biallelic_sites_freqs{$gene}{$site}{$refbase} && $refbase ne "-" ){
							my $alt_base ;
							foreach my $b (keys %{$$Biallelic_sites_freqs{$gene}{$site}} ){
								if($b ne $refbase){
									$alt_base = $b ;
								}
							}
						
							my $af = $$Biallelic_sites_freqs{$gene}{$site}{$alt_base}/$$num_ind_NoRef ;
						
							if( $af>=$$min_AF && $af<=$$max_AF ){
								if( $$Functional_Effect{$gene}{$site} eq "SYNONYMOUS" ){
									$Biallelic_Syn_SNP_Ref_Pos{$count_syn} = $end-$site ;
									@{$Biallelic_Syn_SNP_Alignment_Coord{$count_syn}} = ($gene, $seg_site_index) ;
									$count_syn++ ;
								}
								if( $$Functional_Effect{$gene}{$site} =~ m/NONSYNONYMOUS/ ){
									$Biallelic_NonSyn_SNP_Ref_Pos{$count_nonsyn} = $end-$site ;
									@{$Biallelic_NonSyn_SNP_Alignment_Coord{$count_nonsyn}} = ($gene, $seg_site_index) ;
									$count_nonsyn++ ;
								}
							}
						}
					}
				}
			}
		}
	}

	## NONSYNONYMOUS

	print OUT "min_AF: ", $$min_AF, "\n" ;
	print OUT "max_AF: ", $$max_AF, "\n" ;
	print OUT "NUM SYN BIALLELIC SITES, FILTERED: ", scalar keys %Biallelic_Syn_SNP_Ref_Pos, "\n" ;
	print OUT "NUM NONSYN BIALLELIC SITES, FILTERED: ", scalar keys %Biallelic_NonSyn_SNP_Ref_Pos, "\n" ;
	close OUT ;
	return(\%Biallelic_Syn_SNP_Ref_Pos, \%Biallelic_Syn_SNP_Alignment_Coord, \%Biallelic_NonSyn_SNP_Ref_Pos, \%Biallelic_NonSyn_SNP_Alignment_Coord ) ;
}

sub CalculateNumberVariableSites{
	# This function calculates the proportion of sites that are biallelic or multiallelic.
	# Sites are categorized by their n-fold degeneracy (0D, 2D, 3D, 4D)

	my $Degeneracy_HR = $_[0] ;
	my $Functional_Effect = $_[1] ;
	my $Multiallelic_sites_freqs = $_[2] ;
	my $wd = $_[3] ;
	open OUT, ">>${$wd}/QCfile.txt" ;	

	my %Num_Variable_Sites_By_Degeneracy ;
	foreach my $gene (sort{$a cmp $b} keys %$Degeneracy_HR){
		foreach my $site (sort{$a <=> $b} keys %{$$Degeneracy_HR{$gene}}){
			$Num_Variable_Sites_By_Degeneracy{"total"}{ $$Degeneracy_HR{$gene}{$site} } ++ ;
			if( exists $$Functional_Effect{$gene}{$site} ){
				$Num_Variable_Sites_By_Degeneracy{"biallelic"}{ $$Degeneracy_HR{$gene}{$site} } ++ ;
			}elsif(exists $$Multiallelic_sites_freqs{$gene}{$site}){
				$Num_Variable_Sites_By_Degeneracy{"multiallelic"}{ $$Degeneracy_HR{$gene}{$site} } ++ ;
			}
		}
	}
	print OUT "\nDegeneracy\tTotalNumSites\tProportion_BiAllelic\tProportion_MultiAllelic\n" ;
	foreach my $key (sort{$a cmp $b} keys %{$Num_Variable_Sites_By_Degeneracy{"total"}} ){
		print OUT $key, "\t", $Num_Variable_Sites_By_Degeneracy{"total"}{$key}, "\t", $Num_Variable_Sites_By_Degeneracy{"biallelic"}{$key}/$Num_Variable_Sites_By_Degeneracy{"total"}{$key}, "\t", $Num_Variable_Sites_By_Degeneracy{"multiallelic"}{$key}/$Num_Variable_Sites_By_Degeneracy{"total"}{$key},"\n" ;
	}
	close OUT ;
	return (\%Num_Variable_Sites_By_Degeneracy) ;
}

sub CalculateThetasBiAllelicByGene{
	# This function calculates diversity statistics per gene, which are subsequently used
	# for site-frequency spectrum summaries such as Tajima's D and Fay and Wu's H
	
	my $Biallelic_sites_freqs = $_[0] ; # $Biallelic_sites_freqs{$gene}{$site}{$base} = frequency
	my $Functional_Effect_BiAllelic = $_[1] ; # $Functional_Effect_BiAllelic_HR{$gene}{$site} =  "SYNONYMOUS"/"NONSYNONYMOUS:MISSENSE"/"NONSYNONYMOUS:NONSENSE"  ; 
	my $RefBase = $_[2] ;  # $RefBase_HR{$gene}{$site} = $base ;
	my $n = $_[3] ;
	
	my %Wattersons_Theta_Bi ;
	my %Pi_Theta_Bi ;
	my %L_theta_Bi ;
	my %SegSites_Bi ;
	
	my $a1 = 0 ;
	foreach ( 1 .. ($$n-1) ){ # sample size may vary depending on $max_miss_data
		$a1 += (1/$_) ;
	}
	foreach my $gene( keys %$Biallelic_sites_freqs ){
		foreach my $site ( keys %{$$Biallelic_sites_freqs{$gene}} ){
			my $ref_base = $$RefBase{$gene}{$site} ;
			## are one of the bases in the reference?
			if( exists $$Biallelic_sites_freqs{$gene}{$site}{$ref_base} ){
				$SegSites_Bi{$gene}++ ;
				$Wattersons_Theta_Bi{$gene} += 1/$a1 ;
				foreach my $base ( keys %{$$Biallelic_sites_freqs{$gene}{$site}} ){
					if($base ne $ref_base){
						my $af = $$Biallelic_sites_freqs{$gene}{$site}{$base} ;
						$Pi_Theta_Bi{$gene} += ($af*($$n-$af))/bc($$n,2) ;
						$L_theta_Bi{$gene} += ($af)/($$n-1) ;						
					}
				}
			}
		}
	}
	return(\%Wattersons_Theta_Bi, \%Pi_Theta_Bi, \%L_theta_Bi, \%SegSites_Bi) ;
}

sub SFS_SummStats{
	# This function calculates Tajima's D and Fay and Wu's H per gene!

	my $Wattersons_Theta_Bi = $_[0] ;
	my $Pi_Theta_Bi = $_[1] ;
	my $L_theta_Bi = $_[2] ;
	my $SegSites_Bi = $_[3] ;
	my $n = $_[4] ;
	
	my %TajD_PerGene ;
	my %FayWuH_PerGene ;
	######## CONSTANTS FOR TAJ_D + FAYWUH
	my $a1 = 0 ; 
	my $a2 = 0 ;
	foreach ( 1 .. ($$n - 1) ) { 
		$a1 = $a1 + 1/$_ ; 
		$a2 = $a2 + 1/($_)**2 ;
	}
	my $e1 = (1/$a1)*(($$n+1)/(3*($$n-1)) - (1/$a1)) ;
	my $e2 = (1/($a1**2 + $a2))*((2*(($$n**2) + $$n + 3))/(9*$$n*($$n-1)) - ($$n+2)/($$n*$a1) + ($a2)/($a1**2)) ;
	####### CONSTANTS FOR FAYWUH
	my $b2 = 0 ;
	foreach ( 1 .. ($$n) ) { 
		$b2 = $b2 + 1/($_)**2 ;
	}
	#####################################
	foreach my $gene ( keys %$SegSites_Bi ){
		my $S = $$SegSites_Bi{$gene} ;
		## TAJ_D
		if( (sqrt($S*$e1 + $S*($S-1)*$e2)) > 0 ) { 
			$TajD_PerGene{$gene} = ($$Pi_Theta_Bi{$gene} - $$Wattersons_Theta_Bi{$gene}) / (sqrt($S*$e1 + $S*($S-1)*$e2)) ; 
		}else{
			$TajD_PerGene{$gene} = "NA" ;
		}
		## FayWuH
		my $theta_square = $S*($S-1)/(($a1**2) + $a2) ;
		my $var_1 = ($$n-2)*$$Wattersons_Theta_Bi{$gene}/(6*($$n-1)) ;
		my $var_2 = (18*($$n**2))*((3*$$n) + 2)*$b2 - ((88*($$n**3))+(9*($$n**2))-(13*$$n)+6) ;
		$var_2 = $var_2*($theta_square)/(9*$$n*(($$n-1)**2));
		my $var = $var_1 + $var_2 ;
		if (sqrt($var) > 0){
			$FayWuH_PerGene{$gene} = ($$Pi_Theta_Bi{$gene} - $$L_theta_Bi{$gene})/(sqrt($var)) ;
		}else{
			$FayWuH_PerGene{$gene} = "NA" ;
		}		
		
	}	
	return(\%TajD_PerGene, \%FayWuH_PerGene) ;
}

sub ConstructWindows{
	# This function bins Syn and Nonsyn SNPs into genomic windows, using their relative position in the reference sequence

	my $Biallelic_Syn_SNP_Ref_Pos = $_[0] ;
	my $Biallelic_NonSyn_SNP_Ref_Pos = $_[1] ;
	my $window_size = $_[2] ;

	my %Windows_Syn ;
	my %Windows_NonSyn ;

	foreach my $cnt ( sort{$a <=> $b} keys %$Biallelic_Syn_SNP_Ref_Pos ){
		my $site = $$Biallelic_Syn_SNP_Ref_Pos{$cnt} ;
		my $win = int($site/$$window_size) ;
		push @{$Windows_Syn{$win}}, $cnt  ;
	}
	foreach my $cnt ( sort{$a <=> $b} keys %$Biallelic_NonSyn_SNP_Ref_Pos ){
		my $site = $$Biallelic_NonSyn_SNP_Ref_Pos{$cnt} ;
		my $win = int($site/$$window_size) ;
		push @{$Windows_NonSyn{$win}}, $cnt  ;
	}
	
	return(\%Windows_Syn, \%Windows_NonSyn) ;
}

sub CalculatePairwiseLinkageStats_ByDistance{
	# This function calculates linkage statistics as a function of distance between SNPs

	my $Windows_Syn = $_[0] ;
	my $Windows_NonSyn = $_[1] ;
	my $Biallelic_Syn_SNP_Alignment_Coord = $_[2] ;
	my $Biallelic_Syn_SNP_Ref_Pos = $_[3] ;
	my $Biallelic_NonSyn_SNP_Alignment_Coord = $_[4] ;
	my $Biallelic_NonSyn_SNP_Ref_Pos = $_[5] ;	
	my $core_ref_biallelic_segsites = $_[6] ;
	my $RefBase = $_[7] ;
	my $core_loci_seqs = $_[8] ;
	my $Reference_ind = $_[9] ;
	my $max_dist = $_[10] ;
	my $wd = $_[11] ;

	open OUT, ">>${$wd}/QCfile.txt" ;	

	##### MAKE THIS USE DERIVED ALLELE, NOT ANCESTRAL ONE
	#right now
	#my $site1_MAF_base = $RefBase{$gene1}{$site1} ;
	#makes it such that you're measuring freq of ancestral base?

	my %PairWise_D_vs_dist_SYN ; # @{$PairWise_D_vs_dist_SYN{$win}{$dist}}, ($LD) ;
	my %PairWise_R_vs_dist_SYN ; # @{$PairWise_R_vs_dist_SYN{$win}{$dist}}, ($LD) ;
	my %PairWise_Dprime_vs_dist_SYN ; # @{$PairWise_Dprime_vs_dist_SYN{$win}{$dist}}, $compatibility ;
	my %PairWise_Rsq_vs_dist_SYN ;
	my %PairWise_D_vs_dist_NONSYN ; # @{$PairWise_D_vs_dist_SYN{$win}{$dist}}, ($LD) ;
	my %PairWise_R_vs_dist_NONSYN ; # @{$PairWise_R_vs_dist_SYN{$win}{$dist}}, ($LD) ;
	my %PairWise_Dprime_vs_dist_NONSYN ; # @{$PairWise_Dprime_vs_dist_SYN{$win}{$dist}}, $compatibility ;
	my %PairWise_Rsq_vs_dist_NONSYN ;

	my @IsNonsyn = (1,0) ; # 1 for NonSyn, 0 for Syn; used for making code below more efficient
	foreach my $nonsyn (@IsNonsyn){	
		my ($lower, $upper) ; 	# temporary variables to store indices of 2 SNPs
		my %Windows ;			# we make a copy of the genomic windows data structure
		if($nonsyn){
			%Windows = %$Windows_NonSyn ;
		}else{
			%Windows = %$Windows_Syn ;
		}
		foreach my $win (sort{$a <=> $b} keys %Windows){
			my @counts = @{$Windows{$win}} ; # these are $$Biallelic_Syn_SNP_Ref_Pos count indices
			# Go through each pair of indices in @counts, calculate LD each time
			for ($lower = 0; $lower<((scalar @counts)-1); $lower++){
				my $gene1 ; 
				my $seg_site_index1 ;
				if($nonsyn){
					$gene1 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$counts[$lower]}}[0] ;
					$seg_site_index1 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$counts[$lower]}}[1] ;
				}else{
					$gene1 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$counts[$lower]}}[0] ;
					$seg_site_index1 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$counts[$lower]}}[1] ;
				}				
				my $site1 = $$core_ref_biallelic_segsites{$gene1}{$seg_site_index1} ;
				my $site1_ancestral_base = $$RefBase{$gene1}{$site1} ;
		
				for ($upper = ($lower+1); $upper<(scalar @counts); $upper++){ # note $upper > $lower
					my $dist ;
					if($nonsyn){
						$dist = ($$Biallelic_NonSyn_SNP_Ref_Pos{$counts[$upper]}-$$Biallelic_NonSyn_SNP_Ref_Pos{$counts[$lower]}) ;
					}else{
						$dist = ($$Biallelic_Syn_SNP_Ref_Pos{$counts[$upper]}-$$Biallelic_Syn_SNP_Ref_Pos{$counts[$lower]}) ;
					}
					if( $dist <= $$max_dist ){ # only calculate LD between SNPs within certain distance
						my $gene2 ;
						my $seg_site_index2 ;
						if($nonsyn){
							$gene2 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$counts[$upper]}}[0] ;
							$seg_site_index2 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$counts[$upper]}}[1] ;
						}else{
							$gene2 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$counts[$upper]}}[0] ;
							$seg_site_index2 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$counts[$upper]}}[1] ;
						}
						my $site2 = $$core_ref_biallelic_segsites{$gene2}{$seg_site_index2} ;
						my $site2_ancestral_base = $$RefBase{$gene2}{$site2} ;
				
						my $site1_site2_dblDAF_haplo = 0 ;
						my $site1_DAF = 0 ;
						my $site2_DAF = 0 ;
						my $samp_size_for_pair = 0 ; # this should be same for all SNP pairs, but added to extend analyses to accessory genome later
						my $compatibility ;
						foreach my $ind (keys %{$$core_loci_seqs{$gene1}}){
							if($ind ne $$Reference_ind){
								if( exists $$core_loci_seqs{$gene1}{$ind} && exists $$core_loci_seqs{$gene2}{$ind} ){
									$samp_size_for_pair++ ;
									my $b1 = uc(substr($$core_loci_seqs{$gene1}{$ind}, $site1, 1)) ;
									my $b2 = uc(substr($$core_loci_seqs{$gene2}{$ind}, $site2, 1)) ;
									if( $b1 ne $site1_ancestral_base){
										$site1_DAF++ ;
										if( $b2 ne $site2_ancestral_base ){
											$site1_site2_dblDAF_haplo++ ;
										}
									}
										if( $b2 ne $site2_ancestral_base ){
										$site2_DAF++ ;
									}
								}
							}
						}
						$site1_DAF = $site1_DAF/$samp_size_for_pair ;
						$site2_DAF = $site2_DAF/$samp_size_for_pair ;
						$site1_site2_dblDAF_haplo = $site1_site2_dblDAF_haplo/$samp_size_for_pair ;
						if($site1_DAF > 0 && $site2_DAF >0){				
							my $LD = $site1_site2_dblDAF_haplo - ($site1_DAF*$site2_DAF) ;
							# Normalize D according to Lewontin (1964)
							# if Dab positive, Dmax = min( Pa*(1-Pb), Pb*(1-Pa) )
							# if Dab negative, Dmax = min( Pa*(Pb), (1-Pa)*(1-Pb) )
							my $Dmax = 0 ;
							if ($LD>=0){
								if( $site1_DAF*(1-$site2_DAF) < $site2_DAF*(1-$site1_DAF) ){
									$Dmax = $site1_DAF*(1-$site2_DAF) ;
								}else{
									$Dmax = $site2_DAF*(1-$site1_DAF) ;
								}
							}elsif($LD<0){
								if( ($site1_DAF*$site2_DAF) < (1-$site1_DAF)*(1-$site2_DAF) ){
									$Dmax = $site1_DAF*$site2_DAF ;
								}else{
									$Dmax = (1-$site1_DAF)*(1-$site2_DAF) ;
								}
							}
							if($nonsyn){
								push @{$PairWise_D_vs_dist_NONSYN{$win}{$dist}}, ($LD) ;
								push @{$PairWise_R_vs_dist_NONSYN{$win}{$dist}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
								push @{$PairWise_Dprime_vs_dist_NONSYN{$win}{$dist}}, ($LD)/$Dmax ;
								push @{$PairWise_Rsq_vs_dist_NONSYN{$win}{$dist}}, ($LD**2)/($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
							}else{
								push @{$PairWise_D_vs_dist_SYN{$win}{$dist}}, ($LD) ;
								push @{$PairWise_R_vs_dist_SYN{$win}{$dist}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
								push @{$PairWise_Dprime_vs_dist_SYN{$win}{$dist}}, ($LD)/$Dmax ;
								push @{$PairWise_Rsq_vs_dist_SYN{$win}{$dist}}, ($LD**2)/($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
							}
						}else{
							if($site1_DAF==0){
								print OUT "WARNING: NOT A SEGREGATING SITE: ", $gene1, "\t", $site1, "\n";
							}elsif($site2_DAF==0){
								print OUT "WARNING: NOT A SEGREGATING SITE: ", $gene2, "\t", $site2, "\n" ;
							}
							next ;
						}	
					}
				}
			}

		}
	}
	# Print results for Synonymous SNPs
	my $total_D_SYN = 0 ;
	my $total_D_SYN_counts = 0 ;
	foreach my $win (sort{$a <=> $b} keys %PairWise_D_vs_dist_SYN){
		open OUT1, ">${$wd}/D_R_Dprime_Rsq_vs_Dist_SYN_win${win}.txt" ;
		print OUT1 "Distance", "\t", "Number_SNP_pairs", "\t", "D", "\t", "R", "\t", "Dprime", "\t", "Rsq", "\n" ;
		# Foreach inter-SNP distance category, calculate the mean value for each linkage statistic
		foreach my $dist (sort{$a <=> $b} keys %{$PairWise_D_vs_dist_SYN{$win}}){
			my $sum_D = 0 ;
			my $sum_R = 0 ;
			my $sum_Dpr = 0 ;
			my $sum_Rsq = 0 ;			
			foreach my $x ( @{$PairWise_D_vs_dist_SYN{$win}{$dist}} ){
				$sum_D += $x ;
				$total_D_SYN += $x ;
				$total_D_SYN_counts++ ;
			}
			foreach my $x ( @{$PairWise_R_vs_dist_SYN{$win}{$dist}} ){
				$sum_R += $x ;
			}			
			foreach my $x ( @{$PairWise_Dprime_vs_dist_SYN{$win}{$dist}} ){
				$sum_Dpr += $x ;
			}
			foreach my $x ( @{$PairWise_Rsq_vs_dist_SYN{$win}{$dist}} ){
				$sum_Rsq += $x ;
			}			
			print OUT1 $dist, "\t", scalar @{$PairWise_D_vs_dist_SYN{$win}{$dist}}, "\t" ;
			print OUT1 $sum_D/(scalar @{$PairWise_D_vs_dist_SYN{$win}{$dist}}), "\t" ;
			print OUT1 $sum_R/(scalar @{$PairWise_R_vs_dist_SYN{$win}{$dist}}), "\t" ;
			print OUT1 $sum_Dpr/(scalar @{$PairWise_Dprime_vs_dist_SYN{$win}{$dist}}), "\t" ;
			print OUT1 $sum_Rsq/(scalar @{$PairWise_Rsq_vs_dist_SYN{$win}{$dist}}), "\n" ;
		}
		close OUT1 ;
	}
	print OUT "AVG LD SYN: ", $total_D_SYN/$total_D_SYN_counts, "\n" ;
	
	# Print results for Nonsynonymous SNPs
	my $total_D_NONSYN = 0 ;
	my $total_D_NONSYN_counts = 0 ;
	foreach my $win (sort{$a <=> $b} keys %PairWise_D_vs_dist_NONSYN){
		open OUT1, ">${$wd}/D_R_Dprime_Rsq_vs_Dist_NONSYN_win${win}.txt" ;
		print OUT1 "Distance", "\t", "Number_SNP_pairs", "\t", "D", "\t", "R", "\t", "Dprime", "\t", "Rsq", "\n" ;
		foreach my $dist (sort{$a <=> $b} keys %{$PairWise_D_vs_dist_NONSYN{$win}}){
			my $sum_D = 0 ;
			my $sum_R = 0 ;
			my $sum_Dpr = 0 ;
			my $sum_Rsq = 0 ;			
			foreach my $x ( @{$PairWise_D_vs_dist_NONSYN{$win}{$dist}} ){
				$sum_D += $x ;
				$total_D_NONSYN += $x ;
				$total_D_NONSYN_counts++ ;
			}
			foreach my $x ( @{$PairWise_R_vs_dist_NONSYN{$win}{$dist}} ){
				$sum_R += $x ;
			}			
			foreach my $x ( @{$PairWise_Dprime_vs_dist_NONSYN{$win}{$dist}} ){
				$sum_Dpr += $x ;
			}
			foreach my $x ( @{$PairWise_Rsq_vs_dist_NONSYN{$win}{$dist}} ){
				$sum_Rsq += $x ;
			}			
			print OUT1 $dist, "\t", scalar @{$PairWise_D_vs_dist_NONSYN{$win}{$dist}}, "\t" ;
			print OUT1 $sum_D/(scalar @{$PairWise_D_vs_dist_NONSYN{$win}{$dist}}), "\t" ;
			print OUT1 $sum_R/(scalar @{$PairWise_R_vs_dist_NONSYN{$win}{$dist}}), "\t" ;
			print OUT1 $sum_Dpr/(scalar @{$PairWise_Dprime_vs_dist_NONSYN{$win}{$dist}}), "\t" ;
			print OUT1 $sum_Rsq/(scalar @{$PairWise_Rsq_vs_dist_NONSYN{$win}{$dist}}), "\n" ;
		}
		close OUT1 ;
	}
	print OUT "AVG LD NONSYN: ", $total_D_NONSYN/$total_D_NONSYN_counts, "\n" ;
	close OUT ;
}

sub CalculatePairwiseLinkageStats_ByGeneCategory{

	my $Biallelic_Syn_SNP_Alignment_Coord = $_[0] ;
	my $Biallelic_Syn_SNP_Ref_Pos = $_[1] ;
	my $Biallelic_NonSyn_SNP_Alignment_Coord = $_[2] ;
	my $Biallelic_NonSyn_SNP_Ref_Pos = $_[3] ;	
	my $core_ref_biallelic_segsites = $_[4] ;
	my $RefBase = $_[5] ;
	my $core_loci_seqs = $_[6] ;
	my $Reference_ind = $_[7] ;
	my $min_dist = $_[8] ;
	my $wd = $_[9] ;
	my $Genes_Quantiled_HR = $_[10] ;

	open OUT, ">>${$wd}/QCfile.txt" ;	

	my %PairWise_D_SYN_WithinGene ; # @{$PairWise_D_SYN{$dist}}, ($LD) ;
	my %PairWise_Dprime_SYN_WithinGene ; # @{$PairWise_Dprime_SYN{$dist}}, $compatibility ;
	my %PairWise_r_SYN_WithinGene ; # @{$PairWise_Dprime_SYN{$dist}}, $compatibility ;
	my %PairWise_D_NONSYN_WithinGene ; # @{$PairWise_D_vs_dist_SYN{$dist}}, ($LD) ;
	my %PairWise_Dprime_NONSYN_WithinGene ; # @{$PairWise_Dprime_SYN{$dist}}, $compatibility ;
	my %PairWise_r_NONSYN_WithinGene ; # @{$PairWise_Dprime_SYN{$dist}}, $compatibility ;

	my @count_Syn ;
	my @count_NonSyn ;
	foreach my $cnt ( sort{$a <=> $b} keys %$Biallelic_Syn_SNP_Ref_Pos ){
		push @count_Syn, $cnt ;
	}
	foreach my $cnt ( sort{$a <=> $b} keys %$Biallelic_NonSyn_SNP_Ref_Pos ){
		push @count_NonSyn, $cnt ;
	}

	### SYNONYMOUS
	my ($lower, $upper) ;
	for ($lower = 0; $lower<((scalar @count_Syn)-1); $lower++){
		my $gene1 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$count_Syn[$lower]}}[0] ;
		my $seg_site_index1 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$count_Syn[$lower]}}[1] ;
		my $site1 = $$core_ref_biallelic_segsites{$gene1}{$seg_site_index1} ;
		my $site1_ancestral_base = $$RefBase{$gene1}{$site1} ;
	
		for ($upper = ($lower+1); $upper<(scalar @count_Syn); $upper++){
			my $gene2 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$count_Syn[$upper]}}[0] ;
			my $dist = ($$Biallelic_Syn_SNP_Ref_Pos{$count_Syn[$upper]}-$$Biallelic_Syn_SNP_Ref_Pos{$count_Syn[$lower]}) ;
			## ARE GENES IN THE SAME QUANTILE??
			if(exists $$Genes_Quantiled_HR{$gene1} && exists $$Genes_Quantiled_HR{$gene2} ){
				if( $$Genes_Quantiled_HR{$gene1} == $$Genes_Quantiled_HR{$gene2} ){
					my $quant = $$Genes_Quantiled_HR{$gene1} ;
					my $seg_site_index2 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$count_Syn[$upper]}}[1] ;
					my $site2 = $$core_ref_biallelic_segsites{$gene2}{$seg_site_index2} ;
					my $site2_ancestral_base = $$RefBase{$gene2}{$site2} ;
			
					my $site1_site2_dblDAF_haplo = 0 ;
					my $site1_DAF = 0 ;
					my $site2_DAF = 0 ;
					# go through individuals of lower samp size gene, if they differ
					my $samp_size_for_pair = 0 ;
					foreach my $ind (keys %{$$core_loci_seqs{$gene1}}){
						if($ind ne $$Reference_ind){
							if( exists $$core_loci_seqs{$gene1}{$ind} && exists $$core_loci_seqs{$gene2}{$ind} ){
								$samp_size_for_pair++ ;
								my $b1 = uc(substr($$core_loci_seqs{$gene1}{$ind}, $site1, 1)) ;
								my $b2 = uc(substr($$core_loci_seqs{$gene2}{$ind}, $site2, 1)) ;
								if( $b1 ne $site1_ancestral_base){
									$site1_DAF++ ;
									if( $b2 ne $site2_ancestral_base ){
										$site1_site2_dblDAF_haplo++ ;
									}
								}
								if( $b2 ne $site2_ancestral_base ){
									$site2_DAF++ ;
								}
							}
						}
					}
				
					$site1_DAF = $site1_DAF/$samp_size_for_pair ;
					$site2_DAF = $site2_DAF/$samp_size_for_pair ;
					$site1_site2_dblDAF_haplo = $site1_site2_dblDAF_haplo/$samp_size_for_pair ;
			
					my $LD = $site1_site2_dblDAF_haplo - ($site1_DAF*$site2_DAF) ;
					# if Dab positive, Dmax = min( Pa*(1-Pb), Pb*(1-Pa) )
					# if Dab negative, Dmax = min( Pa*(Pb), (1-Pa)*(1-Pb) )
					my $Dmax = 0 ;
					if ($LD>=0){
						if( $site1_DAF*(1-$site2_DAF) < $site2_DAF*(1-$site1_DAF) ){
							$Dmax = $site1_DAF*(1-$site2_DAF) ;
						}else{
							$Dmax = $site2_DAF*(1-$site1_DAF) ;
						}
					}elsif($LD<0){
						if( ($site1_DAF*$site2_DAF) < (1-$site1_DAF)*(1-$site2_DAF) ){
							$Dmax = $site1_DAF*$site2_DAF ;
						}else{
							$Dmax = (1-$site1_DAF)*(1-$site2_DAF) ;
						}
					}
					if($gene1 eq $gene2){
						push @{$PairWise_D_SYN_WithinGene{$quant}}, ($LD) ;
						push @{$PairWise_Dprime_SYN_WithinGene{$quant}}, ($LD)/$Dmax ;
						push @{$PairWise_r_SYN_WithinGene{$quant}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF)) ;
					}
				}
			}
		}
	}
	
	### NONSYNONYMOUS
	for ($lower = 0; $lower<((scalar @count_NonSyn)-1); $lower++){
		my $gene1 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$count_NonSyn[$lower]}}[0] ;
		my $seg_site_index1 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$count_NonSyn[$lower]}}[1] ;
		my $site1 = $$core_ref_biallelic_segsites{$gene1}{$seg_site_index1} ;
		my $site1_ancestral_base = $$RefBase{$gene1}{$site1} ;
	
		for ($upper = ($lower+1); $upper<(scalar @count_NonSyn); $upper++){
			my $gene2 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$count_NonSyn[$upper]}}[0] ;
			my $dist = ($$Biallelic_NonSyn_SNP_Ref_Pos{$count_NonSyn[$upper]}-$$Biallelic_NonSyn_SNP_Ref_Pos{$count_NonSyn[$lower]}) ;
			if(exists $$Genes_Quantiled_HR{$gene1} && exists $$Genes_Quantiled_HR{$gene2} ){
				if( $$Genes_Quantiled_HR{$gene1} == $$Genes_Quantiled_HR{$gene2} ){
					my $quant = $$Genes_Quantiled_HR{$gene1} ;
					my $seg_site_index2 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$count_NonSyn[$upper]}}[1] ;
					my $site2 = $$core_ref_biallelic_segsites{$gene2}{$seg_site_index2} ;
					my $site2_ancestral_base = $$RefBase{$gene2}{$site2} ;
			
					my $site1_site2_dblDAF_haplo = 0 ;
					my $site1_DAF = 0 ;
					my $site2_DAF = 0 ;
					# go through individuals of lower samp size gene, if they differ
					my $samp_size_for_pair = 0 ;
					foreach my $ind (keys %{$$core_loci_seqs{$gene1}}){
						if($ind ne $$Reference_ind){
							if( exists $$core_loci_seqs{$gene1}{$ind} && exists $$core_loci_seqs{$gene2}{$ind} ){
								$samp_size_for_pair++ ;
								my $b1 = uc(substr($$core_loci_seqs{$gene1}{$ind}, $site1, 1)) ;
								my $b2 = uc(substr($$core_loci_seqs{$gene2}{$ind}, $site2, 1)) ;
								if( $b1 ne $site1_ancestral_base){
									$site1_DAF++ ;
									if( $b2 ne $site2_ancestral_base ){
										$site1_site2_dblDAF_haplo++ ;
									}
								}
								if( $b2 ne $site2_ancestral_base ){
									$site2_DAF++ ;
								}
							}
						}
					}
				
					$site1_DAF = $site1_DAF/$samp_size_for_pair ;
					$site2_DAF = $site2_DAF/$samp_size_for_pair ;
					$site1_site2_dblDAF_haplo = $site1_site2_dblDAF_haplo/$samp_size_for_pair ;
			
					my $LD = $site1_site2_dblDAF_haplo - ($site1_DAF*$site2_DAF) ;
					# if Dab positive, Dmax = min( Pa*(1-Pb), Pb*(1-Pa) )
					# if Dab negative, Dmax = min( Pa*(Pb), (1-Pa)*(1-Pb) )
					my $Dmax = 0 ;
					if ($LD>=0){
						if( $site1_DAF*(1-$site2_DAF) < $site2_DAF*(1-$site1_DAF) ){
							$Dmax = $site1_DAF*(1-$site2_DAF) ;
						}else{
							$Dmax = $site2_DAF*(1-$site1_DAF) ;
						}
					}elsif($LD<0){
						if( ($site1_DAF*$site2_DAF) < (1-$site1_DAF)*(1-$site2_DAF) ){
							$Dmax = $site1_DAF*$site2_DAF ;
						}else{
							$Dmax = (1-$site1_DAF)*(1-$site2_DAF) ;
						}
					}
					if($gene1 eq $gene2){
						push @{$PairWise_D_NONSYN_WithinGene{$quant}}, ($LD) ;
						push @{$PairWise_Dprime_NONSYN_WithinGene{$quant}}, ($LD)/$Dmax ;
						push @{$PairWise_r_NONSYN_WithinGene{$quant}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF)) ;
					}
				}
			}
		}
	}
	
	#################
	## PRINT WITHIN GENE RESULTS
	#################
	my $total_D_SYN = 0 ;
	my $total_D_SYN_tally = 0 ;
	open OUT1, ">${$wd}/D_R_Dprime_vsQuantSYNWithinGene.txt" ;
	print OUT1 "Quantile\tNumObservations\tD\tr\tDprime\n" ;
	foreach my $quant (sort{$a <=> $b} keys %PairWise_D_SYN_WithinGene){
		my $sum_D = 0 ;
		my $sum_Dpr = 0 ;
		my $sum_R = 0 ;			
		foreach my $x ( @{$PairWise_D_SYN_WithinGene{$quant}} ){
			$sum_D += $x ;
			$total_D_SYN += $x ;
			$total_D_SYN_tally++ ;
		}
		foreach my $x ( @{$PairWise_Dprime_SYN_WithinGene{$quant}} ){
			$sum_Dpr += $x ;
		}
		foreach my $x ( @{$PairWise_r_SYN_WithinGene{$quant}} ){
			$sum_R += $x ;
		}
		print OUT1 $quant, "\t", scalar @{$PairWise_D_SYN_WithinGene{$quant}}, "\t" ;
		print OUT1 $sum_D/(scalar @{$PairWise_D_SYN_WithinGene{$quant}}), "\t", $sum_R/(scalar @{$PairWise_r_SYN_WithinGene{$quant}}), "\t", $sum_Dpr/(scalar @{$PairWise_Dprime_SYN_WithinGene{$quant}}), "\n" ;
	}
	close OUT1 ;
	print OUT "AVG LD SYN WITHIN GENE: ", $total_D_SYN/$total_D_SYN_tally, "\n" ;
	#######
	my $total_D_NONSYN = 0 ;
	my $total_D_NONSYN_tally = 0 ;
	open OUT1, ">${$wd}/D_R_Dprime_vsQuantNONSYNWithinGene.txt" ;
	print OUT1 "Quantile\tNumObservations\tD\tr\tDprime\n" ;
	foreach my $quant (sort{$a <=> $b} keys %PairWise_D_NONSYN_WithinGene){
		my $sum_D = 0 ;
		my $sum_Dpr = 0 ;
		my $sum_R = 0 ;			
		foreach my $x ( @{$PairWise_D_NONSYN_WithinGene{$quant}} ){
			$sum_D += $x ;
			$total_D_NONSYN += $x ;
			$total_D_NONSYN_tally++ ;
		}
		foreach my $x ( @{$PairWise_Dprime_NONSYN_WithinGene{$quant}} ){
			$sum_Dpr += $x ;
		}
		foreach my $x ( @{$PairWise_r_NONSYN_WithinGene{$quant}} ){
			$sum_R += $x ;
		}
		print OUT1 $quant, "\t", scalar @{$PairWise_D_NONSYN_WithinGene{$quant}}, "\t" ;
		print OUT1 $sum_D/(scalar @{$PairWise_D_NONSYN_WithinGene{$quant}}), "\t", $sum_R/(scalar @{$PairWise_r_NONSYN_WithinGene{$quant}}), "\t", $sum_Dpr/(scalar @{$PairWise_Dprime_NONSYN_WithinGene{$quant}}), "\n" ;
	}
	close OUT1 ;

	print OUT "AVG LD NONSYN WITHIN GENE: ", $total_D_NONSYN/$total_D_NONSYN_tally, "\n" ;
	
	close OUT ;

}

sub CategorizeBySNPdensity_PrintSFSsummaries{
	# This function prints a lot of information per gene, including the number of ungapped sites in the alignment,
	# the number of Syn and Nonsyn polymorphisms, and site-frequency spectrum summary statistics.
	# However, while we're at it, since we have to cylce through this info we can also categorize genes
	# by SNP density, which is used later on.

	my $core_ref_NumUngapped_sites = $_[0] ;
	my $Biallelic_sites_freqs = $_[1] ;
	my $Multiallelic_sites_freqs = $_[2] ;
	my $Functional_Effect_BiAllelic = $_[3] ;
	my $Functional_Effect_MultiAllelic = $_[4] ;
	my $Functional_Effect_Fixed = $_[5] ;
	my $Ref_gene_info = $_[6] ;
	my $RefBase = $_[7] ;
	my $dNdS_quantiles = $_[8] ;
	my $TajD_PerGene = $_[9] ;
	my $FayWuH_PerGene = $_[10] ;
	my $NonsenseMutsPerGene_HR = $_[11] ;
	my $MinNonsynSNPperGene_forQuants = $_[12] ;
	my $wd = $_[13] ;

	my %DensityNonSynPoly_byGene ;
	my %DensityNonSynPoly_byGene_Quantiled ; 
	
	open OUT, ">${$wd}/Polymorphism_and_SummaryStatistics.txt" ;
	print OUT "Gene", "\t", "PositionInReference", "\t", "Num_ungappedSites", "\t", "Num_BiAl", "\t", "Num_MultAl", "\t", "Num_SynPoly", "\t", "Num_NonSynPoly", "\t", "Num_SynFixed", "\t", "Num_NonSynFixed", "\t", "TajimasD", "\t", "FayWuH", "\t", "Num_NonsenseMuts\n" ; 
	foreach my $gene (keys %$core_ref_NumUngapped_sites){
		print OUT $gene, "\t" ;
		# print position
		if( exists $$Ref_gene_info{$gene}{"POSITION_FORWARD"} ){
			print OUT $$Ref_gene_info{$gene}{"POSITION_FORWARD"}{"START"}, "\t" ;
		}else{
			print OUT $$Ref_gene_info{$gene}{"POSITION_REVERSE"}{"START"}, "\t" ;
		}
	
		print OUT $$core_ref_NumUngapped_sites{$gene}, "\t" ;
		if(exists $$Biallelic_sites_freqs{$gene}){
			print OUT scalar keys %{$$Biallelic_sites_freqs{$gene}}, "\t" ;
		}else{
			print OUT "0\t" ;
		}
		if(exists $$Multiallelic_sites_freqs{$gene}){
			print OUT scalar keys %{$$Multiallelic_sites_freqs{$gene}}, "\t" ;
		}else{
			print OUT "0\t" ;
		}

		my $num_syn = 0 ;
		my $num_nonsyn = 0 ;	
		if(exists $$Functional_Effect_BiAllelic{$gene}){
			foreach my $site  (keys %{$$Functional_Effect_BiAllelic{$gene}} ){
				if( $$Functional_Effect_BiAllelic{$gene}{$site} eq "SYNONYMOUS" ){
					$num_syn++ ;
				}elsif( $$Functional_Effect_BiAllelic{$gene}{$site} =~ m/NONSYNONYMOUS/ ){
					$num_nonsyn++ ;
				}
			}
		}
		if(exists $$Functional_Effect_MultiAllelic{$gene}){
			foreach my $site ( keys %{$$Functional_Effect_MultiAllelic{$gene}} ){
				my $refbase = $$RefBase{$gene}{$site} ;
				foreach my $b2 (keys %{$$Functional_Effect_MultiAllelic{$gene}{$site}{$refbase}} ){
					if( $$Functional_Effect_MultiAllelic{$gene}{$site}{$refbase}{$b2} eq "SYNONYMOUS" ){
						$num_syn++ ;
					}elsif( $$Functional_Effect_MultiAllelic{$gene}{$site}{$refbase}{$b2} =~ m/NONSYNONYMOUS/ ){
						$num_nonsyn++ ;
					}
				}
			}
		}
		if($num_syn){
			print OUT $num_syn, "\t" ;
		}else{
			print OUT "0", "\t" ;
		}
		## BIN GENES INTO QUANTILES BY NONSYN SNP DENSITY
		if($$core_ref_NumUngapped_sites{$gene} && $num_nonsyn >= $$MinNonsynSNPperGene_forQuants ){
			$DensityNonSynPoly_byGene{$gene} = $num_nonsyn/$$core_ref_NumUngapped_sites{$gene} ;
		}
		if($num_nonsyn){
			print OUT $num_nonsyn, "\t" ;
			#requires at least 1 nonsyn polymorphism
			#$NumNonSynPoly_byGene{$gene} = $num_nonsyn ;
			#$DensityNonSynPoly_byGene{$gene} = $num_nonsyn/$$core_ref_NumUngapped_sites{$gene} ;
			
		}else{
			print OUT "0", "\t" ;
		}			
		if(exists $$Functional_Effect_Fixed{$gene}){
			my $num_syn = 0 ;
			my $num_nonsyn = 0 ;
			foreach my $site  (keys %{$$Functional_Effect_Fixed{$gene}} ){
				if( $$Functional_Effect_Fixed{$gene}{$site} eq "SYNONYMOUS" ){
					$num_syn++ ;
				}elsif( $$Functional_Effect_Fixed{$gene}{$site} =~ m/NONSYNONYMOUS/ ){
					$num_nonsyn++ ;
				}
			}
			print OUT $num_syn, "\t", $num_nonsyn, "\t" ;
		}else{
			print OUT "0", "\t", "0", "\t" ;
		}		
		if(exists $$TajD_PerGene{$gene}){
			print OUT $$TajD_PerGene{$gene}, "\t" ;
		}else{	
			print OUT "NA", "\t" ;
		}
		
		if(exists $$FayWuH_PerGene{$gene}){
			print OUT $$FayWuH_PerGene{$gene}, "\t" ;
		}else{	
			print OUT "NA", "\t" ;
		}
		if(exists $$NonsenseMutsPerGene_HR{$gene}){
			print OUT $$NonsenseMutsPerGene_HR{$gene}, "\n" ;
		}else{	
			print OUT "0", "\n" ;
		}
		
	}
	close OUT ;
	
	# Make DensityNonsynPoly into quantiles
	my $numGenesPerQuantile = (scalar keys %DensityNonSynPoly_byGene)/$$dNdS_quantiles ;
	my $counter = 0 ;
	my %Quant_means ;
	foreach my $gene ( sort{ $DensityNonSynPoly_byGene{$a} <=> $DensityNonSynPoly_byGene{$b} } keys %DensityNonSynPoly_byGene){
		my $quant = int($counter/$numGenesPerQuantile) ;
		$DensityNonSynPoly_byGene_Quantiled{$gene} = $quant ;
		$Quant_means{$quant} += $DensityNonSynPoly_byGene{$gene}/$numGenesPerQuantile ;
		$counter++ ;
	}
	open OUT, ">${$wd}/DensityNonSynPoly_meanByQuantile.txt" ;
	print OUT "Quantile", "\t", "Mean", "\t", "NumGenes", "\n" ;
	foreach my $quant( sort{$a<=>$b} keys %Quant_means ){
		print OUT $quant, "\t", $Quant_means{$quant}, "\t", $numGenesPerQuantile, "\n" ;
	}
	close OUT ;	
	
	return(\%DensityNonSynPoly_byGene_Quantiled) ;
}

sub PrintNumGenesPerQuantile{
	# This function prints the number of genes per quantile
		
	my $Genes_Quantiled = $_[0] ;
	my $wd = $_[1] ;

 	open OUT, ">>${$wd}/QCfile.txt" ;

    my %NQ ; # temporary data structure to hold the quantile values

    #Calculate number of genes per quantile
    foreach my $value (values %$Genes_Quantiled ){
        $NQ{$value}++ ;
    }
    my $NumQuantiles = scalar keys %NQ ;
    print OUT "Number of categories/quantiles: ", $NumQuantiles, "\n" ;
    print OUT "Quantile\tNumberOfGenes\n" ;
	foreach my $quant ( keys %NQ ){
        print OUT $quant, "\t", $NQ{$quant}, "\n" ;
	}


}

sub ZerofoldFourfoldDivByQuantile{
	# This function calculates the relative amounts of 0D site diversity to 4D site diversity

	my $ZeroFourfoldSegSites_HR = $_[0] ;
	my $ZeroFourfoldSites_UnGapped_HR = $_[1] ;
	my $Gene_Quantiled_HR = $_[2] ;
	my $wd = $_[3] ;

	my %ZeroFourDiv ;
	my %tally ;

	foreach my $gene (sort {$a cmp $b} keys %$Gene_Quantiled_HR ){
		my $quant = $$Gene_Quantiled_HR{$gene} ;
		if($$ZeroFourfoldSegSites_HR{$gene}{"4D"}){
			$tally{$quant}++ ;
			$ZeroFourDiv{$quant}+= $$ZeroFourfoldSegSites_HR{$gene}{"0D"}/($$ZeroFourfoldSegSites_HR{$gene}{"4D"}) ;
		}
	}
	open OUT, ">${$wd}/ZerofoldDividedByFourfoldDiversity.txt" ;
	print OUT "Quantile\tZeroDivFour_diversity\tNumGenes\n" ;
	foreach my $quant (sort {$a <=> $b} keys %ZeroFourDiv){
		print OUT $quant, "\t", $ZeroFourDiv{$quant}/$tally{$quant}, "\t", $tally{$quant}, "\n" ;
	}
	close OUT ;
}

sub BinLDbyDistForBoostrapResampling_WholeGenome{
	# Computes rN-rS for SNPs within the same distance bin
	# This function also outputs the number of Syn and NSyn observations per bin

	my $interSNPdist_binSize = $_[0] ;
	my $interSNPdist_numBins = $_[1] ;
	my $Biallelic_Syn_SNP_Alignment_Coord = $_[2] ;
	my $Biallelic_Syn_SNP_Ref_Pos = $_[3] ;
	my $Biallelic_NonSyn_SNP_Alignment_Coord = $_[4] ;
	my $Biallelic_NonSyn_SNP_Ref_Pos = $_[5] ;	
	my $core_ref_biallelic_segsites = $_[6] ;
	my $RefBase = $_[7] ;
	my $core_loci_seqs = $_[8] ;
	my $Reference_ind = $_[9] ;
	my $DensityNonSynPoly_byGene_Quantiled = $_[10] ; # $$NumNonSynPoly_byGene_Quantiled_HR{$gene} = quantile, [0, $dNdS_quantiles-1]
	my $NumQuantiles = $_[11] ;
	my $Bootreps = $_[12] ;
	my $wd = $_[13] ;

	open OUT, ">>${$wd}/QCfile.txt" ;	
	print OUT "BinSize: ", $$interSNPdist_binSize, "\n" ;
	print OUT "NumBins: ", $$interSNPdist_numBins, "\n" ;

	# Synonynomous LD binned by interSNP distance, for bootstrap resampling
	my %SynR_WholeGenome_ResamplingBins ; # @{$SynR_WholeGenome_ResamplingBins{$quant}{$bin}} = (values)
	my %SynR_Intragenic_ResamplingBins ;  # Same structure as %SynR_WholeGenome_ResamplingBins, but contains rS values only for SNP pairs within the same gene
	my %SynR_Intergenic_ResamplingBins ;  # Same structure as %SynR_WholeGenome_ResamplingBins, but contains rS values only for SNP pairs from different genes
	# Nonsynonymous LD binned by interSNP distance
	my %NonsynR_WholeGenome_BinByDist ;
	my %NonsynR_Intragenic_BinByDist ;
	my %NonsynR_Intergenic_BinByDist ;

	
	# INITIALIZE
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			@{$SynR_WholeGenome_ResamplingBins{$quant}{$bin}} = () ;
			@{$SynR_Intragenic_ResamplingBins{$quant}{$bin}} = () ;
			@{$SynR_Intergenic_ResamplingBins{$quant}{$bin}} = () ;
			@{$NonsynR_WholeGenome_BinByDist{$quant}{$bin}} = () ;
			@{$NonsynR_Intragenic_BinByDist{$quant}{$bin}} = () ;
			@{$NonsynR_Intergenic_BinByDist{$quant}{$bin}} = () ;
		}
	}

	# Gather rS values per bin
	my $firstIndex ;
	my $lastIndex ;
	($firstIndex) = sort {$a <=> $b} keys %$Biallelic_Syn_SNP_Ref_Pos ; # Indexed from 0 to (Num Syn SNPs)-1
	($lastIndex) = sort {$b <=> $a} keys %$Biallelic_Syn_SNP_Ref_Pos ; # Indexed from 0 to (Num Syn SNPs)-1

	my ($lower, $upper) ;
	for ($lower = $firstIndex; $lower<=($lastIndex-1); $lower++){
		my $gene1 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$lower}}[0] ;
		my $seg_site_index1 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$lower}}[1] ;
		my $site1 = $$core_ref_biallelic_segsites{$gene1}{$seg_site_index1} ;
		my $site1_ancestral_base = $$RefBase{$gene1}{$site1} ;
	
		for ($upper = ($lower+1); $upper<$lastIndex; $upper++){
			my $gene2 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$upper}}[0] ;
			if( $$DensityNonSynPoly_byGene_Quantiled{$gene1} == $$DensityNonSynPoly_byGene_Quantiled{$gene2} ){
				my $dist = ($$Biallelic_Syn_SNP_Ref_Pos{$upper}-$$Biallelic_Syn_SNP_Ref_Pos{$lower}) ;
				my $seg_site_index2 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$upper}}[1] ;
				my $site2 = $$core_ref_biallelic_segsites{$gene2}{$seg_site_index2} ;
				my $site2_ancestral_base = $$RefBase{$gene2}{$site2} ;
				my $quant = $$DensityNonSynPoly_byGene_Quantiled{$gene1} ; # should be same as gene2
				my $bin = int($dist/$$interSNPdist_binSize) ;
				
				if( $bin <= $$interSNPdist_numBins-1 && $dist>0 ){
					my $site1_site2_dblDAF_haplo = 0 ;
					my $site1_DAF = 0 ;
					my $site2_DAF = 0 ;
					my $samp_size_for_pair = 0 ;
					foreach my $ind (keys %{$$core_loci_seqs{$gene1}}){
						if($ind ne $$Reference_ind){
							if( exists $$core_loci_seqs{$gene1}{$ind} && exists $$core_loci_seqs{$gene2}{$ind} ){
								$samp_size_for_pair++ ;
								my $b1 = uc(substr($$core_loci_seqs{$gene1}{$ind}, $site1, 1)) ;
								my $b2 = uc(substr($$core_loci_seqs{$gene2}{$ind}, $site2, 1)) ;
								if( $b1 ne $site1_ancestral_base){
									$site1_DAF++ ;
									if( $b2 ne $site2_ancestral_base ){
										$site1_site2_dblDAF_haplo++ ;
									}
								}
								if( $b2 ne $site2_ancestral_base ){
									$site2_DAF++ ;
								}
							}
						}
					}
					$site1_DAF = $site1_DAF/$samp_size_for_pair ;
					$site2_DAF = $site2_DAF/$samp_size_for_pair ;
					$site1_site2_dblDAF_haplo = $site1_site2_dblDAF_haplo/$samp_size_for_pair ;
					if($site1_DAF > 0 && $site2_DAF >0){				
						my $LD = $site1_site2_dblDAF_haplo - ($site1_DAF*$site2_DAF) ;
						push @{$SynR_WholeGenome_ResamplingBins{$quant}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
						if( $gene1 eq $gene2 ){
							push @{$SynR_Intragenic_ResamplingBins{$quant}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
						}else{
							push @{$SynR_Intergenic_ResamplingBins{$quant}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
						}
					}else{
						if($site1_DAF==0){
							print OUT "NOT A SEG SITE: ", $gene1, "\t", $site1, "\n";
						}elsif($site2_DAF==0){
							print OUT "NOT A SEG SITE: ", $gene2, "\t", $site2, "\n" ;
						}
						next ;
					}	
				}
			}
		}
	}

	open OUT1, ">${$wd}/SynResamplingBinsNumObservations.txt" ;
	print OUT1 "Quantile\tBinNum\tNumObs_WG\tNumObs_Intragenic\tNumObs_Intergenic\n" ;
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			print OUT1 $quant, "\t", $bin, "\t", scalar(@{$SynR_WholeGenome_ResamplingBins{$quant}{$bin}}), "\t", scalar(@{$SynR_Intragenic_ResamplingBins{$quant}{$bin}}), "\t", scalar(@{$SynR_Intergenic_ResamplingBins{$quant}{$bin}}),  "\n" ;	
		}
	}
 	close OUT1 ;
 
	# Gather rN values per bin
	($firstIndex) = sort {$a <=> $b} keys %$Biallelic_NonSyn_SNP_Ref_Pos ; # Indexed from 0 to (Num Syn SNPs)-1
	($lastIndex) = sort {$b <=> $a} keys %$Biallelic_NonSyn_SNP_Ref_Pos ; # Indexed from 0 to (Num Syn SNPs)-1
	for ($lower = $firstIndex; $lower<=($lastIndex-1); $lower++){
		my $gene1 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$lower}}[0] ;
		my $seg_site_index1 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$lower}}[1] ;
		my $site1 = $$core_ref_biallelic_segsites{$gene1}{$seg_site_index1} ;
		my $site1_ancestral_base = $$RefBase{$gene1}{$site1} ;
	
		for ($upper = ($lower+1); $upper<$lastIndex; $upper++){
			my $gene2 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$upper}}[0] ;
			if( $$DensityNonSynPoly_byGene_Quantiled{$gene1} == $$DensityNonSynPoly_byGene_Quantiled{$gene2} ){
				my $dist = ($$Biallelic_NonSyn_SNP_Ref_Pos{$upper}-$$Biallelic_NonSyn_SNP_Ref_Pos{$lower}) ;
				my $seg_site_index2 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$upper}}[1] ;
				my $site2 = $$core_ref_biallelic_segsites{$gene2}{$seg_site_index2} ;
				my $site2_ancestral_base = $$RefBase{$gene2}{$site2} ;
				my $quant = $$DensityNonSynPoly_byGene_Quantiled{$gene1} ; # should be same as gene2
				my $bin = int($dist/$$interSNPdist_binSize) ;
				
				if( $bin <= $$interSNPdist_numBins-1 && $dist>0 ){
					my $site1_site2_dblDAF_haplo = 0 ;
					my $site1_DAF = 0 ;
					my $site2_DAF = 0 ;
					my $samp_size_for_pair = 0 ;
					foreach my $ind (keys %{$$core_loci_seqs{$gene1}}){
						if($ind ne $$Reference_ind){
							if( exists $$core_loci_seqs{$gene1}{$ind} && exists $$core_loci_seqs{$gene2}{$ind} ){
								$samp_size_for_pair++ ;
								my $b1 = uc(substr($$core_loci_seqs{$gene1}{$ind}, $site1, 1)) ;
								my $b2 = uc(substr($$core_loci_seqs{$gene2}{$ind}, $site2, 1)) ;
								if( $b1 ne $site1_ancestral_base){
									$site1_DAF++ ;
									if( $b2 ne $site2_ancestral_base ){
										$site1_site2_dblDAF_haplo++ ;
									}
								}
								if( $b2 ne $site2_ancestral_base ){
									$site2_DAF++ ;
								}
							}
						}
					}
					$site1_DAF = $site1_DAF/$samp_size_for_pair ;
					$site2_DAF = $site2_DAF/$samp_size_for_pair ;
					$site1_site2_dblDAF_haplo = $site1_site2_dblDAF_haplo/$samp_size_for_pair ;
					if($site1_DAF > 0 && $site2_DAF >0){				
						my $LD = $site1_site2_dblDAF_haplo - ($site1_DAF*$site2_DAF) ;
						push @{$NonsynR_WholeGenome_BinByDist{$quant}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
						if( $gene1 eq $gene2 ){
							push @{$NonsynR_Intragenic_BinByDist{$quant}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
						}else{
							push @{$NonsynR_Intergenic_BinByDist{$quant}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
						}
					}else{
						if($site1_DAF==0){
							print OUT "NOT A SEG SITE: ", $gene1, "\t", $site1, "\n";
						}elsif($site2_DAF==0){
							print OUT "NOT A SEG SITE: ", $gene2, "\t", $site2, "\n" ;
						}
						next ;
					}	
				}
			}
		}
	}
	open OUT1, ">${$wd}/NonSynResamplingBinsNumObservations.txt" ;
	print OUT1 "Quantile\tBinNum\tNumObs_WG\tNumObs_Intragenic\tNumObs_Intergenic\n" ;
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			print OUT1 $quant, "\t", $bin, "\t", scalar(@{$NonsynR_WholeGenome_BinByDist{$quant}{$bin}}), "\t", scalar(@{$NonsynR_Intragenic_BinByDist{$quant}{$bin}}), "\t", scalar(@{$NonsynR_Intergenic_BinByDist{$quant}{$bin}}),  "\n" ;      
		}
	}
 	close OUT1 ;

	# RANDOMIZE BINS SO RESAMPLING RANDOM 
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			if( @{$SynR_WholeGenome_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$SynR_WholeGenome_ResamplingBins{$quant}{$bin}}) };
			if( @{$SynR_Intragenic_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$SynR_Intragenic_ResamplingBins{$quant}{$bin}}) };
			if( @{$SynR_Intergenic_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$SynR_Intergenic_ResamplingBins{$quant}{$bin}}) };
			if( @{$NonsynR_WholeGenome_BinByDist{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$NonsynR_WholeGenome_BinByDist{$quant}{$bin}}) };
			if( @{$NonsynR_Intragenic_BinByDist{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$NonsynR_Intragenic_BinByDist{$quant}{$bin}}) };
			if( @{$NonsynR_Intergenic_BinByDist{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$NonsynR_Intergenic_BinByDist{$quant}{$bin}}) };
		}
	}

	# Compare rN and rS by bin for entire genome
	my %NonsynRvSynR_WholeGenome_ByBin ; # @{$NonsynRvSynR_WholeGenome_ByBin{quant}{bin}} = (NonsynR - SynR, ... )
	# Compare rN and rS within and between genes
	my %NonsynRvSynR_Intragenic ; # @{$NonsynRvSynR_Intragenic{quant}} = (NonsynR - SynR, ... )
	my %NonsynRvSynR_Intergenic ; # @{$NonsynRvSynR_Intergenic{quant}} = (NonsynR - SynR, ... )
	my %NonsynRvSynR_Intragenic_ByBin ; # @{$NonsynRvSynR_Intragenic_ByBin{quant}{bin}} = (NonsynR - SynR, ... )
	my %NonsynRvSynR_Intergenic_ByBin ; # @{$NonsynRvSynR_Intergenic_ByBin{quant}{bin}} = (NonsynR - SynR, ... )
	# Initialize structures
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			@{$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}} = () ;
			@{$NonsynRvSynR_Intragenic_ByBin{$quant}{$bin}} = () ;
			@{$NonsynRvSynR_Intergenic_ByBin{$quant}{$bin}} = () ;
		}
	}
	# Foreach quantile, go through all bins, and for each bootstrap replicate, sample without replacement from whatever bin (Syn or NSyn) has more observations
	# Calculate mean rN-rS for each bin
	# First do analysis for whole genome
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			# Go through bootreps for each bin, so you don't uneccesarily recalculate the mean of whatever category (SYN/NONSYN) that has fewer observations
			# which you would need to determine every bootstrap replicate
			my $SynSNPpairs = scalar @{$SynR_WholeGenome_ResamplingBins{$quant}{$bin}} ;
			my $NonsynSNPpairs = scalar @{$NonsynR_WholeGenome_BinByDist{$quant}{$bin}} ;
			# Only compare Syn and NSyn when there are observations for both categories
			if( $SynSNPpairs && $NonsynSNPpairs ){ 
				if( $SynSNPpairs >= $NonsynSNPpairs ){ # Synonymous SNPs have more observations; sample them w/out replacement
					my $NonsynMean = 0 ;
					foreach ( @{$NonsynR_WholeGenome_BinByDist{$quant}{$bin}} ){
						$NonsynMean += $_ ;
					}
					$NonsynMean = $NonsynMean/$NonsynSNPpairs ;
					foreach my $rep ( 1 .. $$Bootreps ){
						my @indices_left_syn = ( 0 .. $SynSNPpairs-1 ) ;
						my @Syn_SampWOreplacementIndices ;
						foreach ( 1 .. $NonsynSNPpairs ){ # Subsample same number of Syn Observations as Nonsyn observations you just took mean of
							my $tmp = int(rand( scalar @indices_left_syn )) ;
							push @Syn_SampWOreplacementIndices, $indices_left_syn[$tmp] ;
							splice @indices_left_syn, $tmp, 1 ;
						}
						my $SynMean = 0 ;
						foreach my $index ( @Syn_SampWOreplacementIndices ){
							$SynMean += ${$SynR_WholeGenome_ResamplingBins{$quant}{$bin}}[$index] ;
						}
						$SynMean = $SynMean/$NonsynSNPpairs ; # $NonsynSNPpairs in denom since you sampled from Syn SNP pairs the same number of Nonsyn SNP pairs
						push @{$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}}, ($NonsynMean-$SynMean) ;
					}
				}else{ # Nonsynonymous SNPs have more observations; sample w/out replacement from these
					my $SynMean = 0 ;
					foreach ( @{$SynR_WholeGenome_ResamplingBins{$quant}{$bin}} ){
						$SynMean += $_ ;
					}
					$SynMean = $SynMean/$SynSNPpairs ;
					foreach my $rep ( 1 .. $$Bootreps ){
						my @indices_left_nonsyn = ( 0 .. $NonsynSNPpairs-1 ) ;
						my @Nonsyn_SampWOreplacementIndices ;
						foreach ( 1 .. $SynSNPpairs ){  # Subsample same number of NonSyn Observations as Syn observations you just took mean of
							my $tmp = int(rand( scalar @indices_left_nonsyn )) ;
							push @Nonsyn_SampWOreplacementIndices, $indices_left_nonsyn[$tmp] ;
							splice @indices_left_nonsyn, $tmp, 1 ;
						}
						my $NonsynMean = 0 ;
						foreach my $index ( @Nonsyn_SampWOreplacementIndices ){
							$NonsynMean += ${$NonsynR_WholeGenome_BinByDist{$quant}{$bin}}[$index] ;
						}
						$NonsynMean = $NonsynMean/$SynSNPpairs ; # $SynSNPpairs since you sampled from Nonsyn SNP pairs the same number of Syn SNP pairs
						push @{$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}}, ($NonsynMean-$SynMean) ;
					}
				}
			}else{
				$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin} = "NA" ;
			}
		}
	}

	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		open OUT1, ">${$wd}/NonsynRvSynR_WholeGenome_Quant${quant}.txt" ;
		print OUT1 "BinNum\t0.05pQuant\t0.5pQuant\t0.95pQuant\tNumSynSNPpairs\tNumNonsynSNPpairs\n" ;
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			print OUT1 $bin, "\t" ;
			if($NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin} ne "NA"){	
				@{$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}} = sort {$a <=> $b} @{$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}} ;
				my $lower90CI = int(scalar(@{$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}})*0.05) - 1 ;
				my $median = int(scalar(@{$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}})*0.5) - 1 ;
				my $upper90CI = int(scalar(@{$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}})*0.95) - 1 ;
				print OUT1 ${$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}}[$lower90CI], "\t" ;
				print OUT1 ${$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}}[$median], "\t" ;
				print OUT1 ${$NonsynRvSynR_WholeGenome_ByBin{$quant}{$bin}}[$upper90CI], "\t" ;
				print OUT1 scalar(@{$SynR_WholeGenome_ResamplingBins{$quant}{$bin}}), "\t" ;
				print OUT1 scalar(@{$NonsynR_WholeGenome_BinByDist{$quant}{$bin}}), "\n" ;
			}else{
				print OUT1 "NA\tNA\tNA\tNA\tNA\n" ;
			}
		}
		close OUT1 ;
	}
}

sub BinLDbyDistForBoostrapResampling_IntraInterGenic{

	my $interSNPdist_binSize = $_[0] ;
	my $interSNPdist_numBins = $_[1] ;
	my $Biallelic_Syn_SNP_Alignment_Coord = $_[2] ;
	my $Biallelic_Syn_SNP_Ref_Pos = $_[3] ;
	my $Biallelic_NonSyn_SNP_Alignment_Coord = $_[4] ;
	my $Biallelic_NonSyn_SNP_Ref_Pos = $_[5] ;	
	my $core_ref_biallelic_segsites = $_[6] ;
	my $RefBase = $_[7] ;
	my $core_loci_seqs = $_[8] ;
	my $Reference_ind = $_[9] ;
	my $DensityNonSynPoly_byGene_Quantiled = $_[10] ; # $$NumNonSynPoly_byGene_Quantiled_HR{$gene} = quantile, [0, $dNdS_quantiles-1]
	my $NumQuantiles = $_[11] ;
	my $Bootreps = $_[12] ;
	my $wd = $_[13] ;

	open OUT, ">>${$wd}/QCfile.txt" ;	
	print OUT "BinSize: ", $$interSNPdist_binSize, "\n" ;
	print OUT "NumBins: ", $$interSNPdist_numBins, "\n" ;

	# Synonynomous LD binned by interSNP distance, for bootstrap resampling
	my %SynR_Intragenic_ResamplingBins ;
	my %SynR_Intergenic_ResamplingBins ;
	# Nonsynonymous LD binned by interSNP distance
	my %NonsynR_Intragenic_ResamplingBins ;
	my %NonsynR_Intergenic_ResamplingBins ;

	# Synonynomous LD binned by interSNP distance, for bootstrap resampling
	my %SynDpr_Intragenic_ResamplingBins ;
	my %SynDpr_Intergenic_ResamplingBins ;
	# Nonsynonymous LD binned by interSNP distance
	my %NonsynDpr_Intragenic_ResamplingBins ;
	my %NonsynDpr_Intergenic_ResamplingBins ;

	# INITIALIZE
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			@{$SynR_Intragenic_ResamplingBins{$quant}{$bin}} = () ;
			@{$SynR_Intergenic_ResamplingBins{$quant}{$bin}} = () ;
			@{$NonsynR_Intragenic_ResamplingBins{$quant}{$bin}} = () ;
			@{$NonsynR_Intergenic_ResamplingBins{$quant}{$bin}} = () ;

			@{$SynDpr_Intragenic_ResamplingBins{$quant}{$bin}} = () ;
			@{$SynDpr_Intergenic_ResamplingBins{$quant}{$bin}} = () ;
			@{$NonsynDpr_Intragenic_ResamplingBins{$quant}{$bin}} = () ;
			@{$NonsynDpr_Intergenic_ResamplingBins{$quant}{$bin}} = () ;
		}
	}

	### SYNONYMOUS
	my $firstIndex ;
	my $lastIndex ;
	($firstIndex) = sort {$a <=> $b} keys %$Biallelic_Syn_SNP_Ref_Pos ; # Indexed from 0 to (Num Syn SNPs)-1
	($lastIndex) = sort {$b <=> $a} keys %$Biallelic_Syn_SNP_Ref_Pos ; # Indexed from 0 to (Num Syn SNPs)-1

	my ($lower, $upper) ;
	for ($lower = $firstIndex; $lower<=($lastIndex-1); $lower++){
		my $gene1 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$lower}}[0] ;
		my $seg_site_index1 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$lower}}[1] ;
		my $site1 = $$core_ref_biallelic_segsites{$gene1}{$seg_site_index1} ;
		my $site1_ancestral_base = $$RefBase{$gene1}{$site1} ;
	
		for ($upper = ($lower+1); $upper<$lastIndex; $upper++){
			my $gene2 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$upper}}[0] ;
			if( $$DensityNonSynPoly_byGene_Quantiled{$gene1} == $$DensityNonSynPoly_byGene_Quantiled{$gene2} ){
				my $dist = ($$Biallelic_Syn_SNP_Ref_Pos{$upper}-$$Biallelic_Syn_SNP_Ref_Pos{$lower}) ;
				my $seg_site_index2 = ${$$Biallelic_Syn_SNP_Alignment_Coord{$upper}}[1] ;
				my $site2 = $$core_ref_biallelic_segsites{$gene2}{$seg_site_index2} ;
				my $site2_ancestral_base = $$RefBase{$gene2}{$site2} ;
				my $quant = $$DensityNonSynPoly_byGene_Quantiled{$gene1} ; # should be same as gene2
				my $bin = int($dist/$$interSNPdist_binSize) ;
				
				if( $bin <= $$interSNPdist_numBins-1 && $dist>0 ){
					my $site1_site2_dblDAF_haplo = 0 ;
					my $site1_DAF = 0 ;
					my $site2_DAF = 0 ;
					my $samp_size_for_pair = 0 ;
					foreach my $ind (keys %{$$core_loci_seqs{$gene1}}){
						if($ind ne $$Reference_ind){
							if( exists $$core_loci_seqs{$gene1}{$ind} && exists $$core_loci_seqs{$gene2}{$ind} ){
								$samp_size_for_pair++ ;
								my $b1 = uc(substr($$core_loci_seqs{$gene1}{$ind}, $site1, 1)) ;
								my $b2 = uc(substr($$core_loci_seqs{$gene2}{$ind}, $site2, 1)) ;
								if( $b1 ne $site1_ancestral_base){
									$site1_DAF++ ;
									if( $b2 ne $site2_ancestral_base ){
										$site1_site2_dblDAF_haplo++ ;
									}
								}
								if( $b2 ne $site2_ancestral_base ){
									$site2_DAF++ ;
								}
							}
						}
					}
					$site1_DAF = $site1_DAF/$samp_size_for_pair ;
					$site2_DAF = $site2_DAF/$samp_size_for_pair ;
					$site1_site2_dblDAF_haplo = $site1_site2_dblDAF_haplo/$samp_size_for_pair ;
					if($site1_DAF > 0 && $site2_DAF >0){				
						my $LD = $site1_site2_dblDAF_haplo - ($site1_DAF*$site2_DAF) ;
						# if Dab positive, Dmax = min( Pa*(1-Pb), Pb*(1-Pa) )
						# if Dab negative, Dmax = min( Pa*(Pb), (1-Pa)*(1-Pb) )
						my $Dmax = 0 ;
						if ($LD>=0){
							if( $site1_DAF*(1-$site2_DAF) < $site2_DAF*(1-$site1_DAF) ){
								$Dmax = $site1_DAF*(1-$site2_DAF) ;
							}else{
								$Dmax = $site2_DAF*(1-$site1_DAF) ;
							}
						}elsif($LD<0){
							if( ($site1_DAF*$site2_DAF) < (1-$site1_DAF)*(1-$site2_DAF) ){
								$Dmax = $site1_DAF*$site2_DAF ;
							}else{
								$Dmax = (1-$site1_DAF)*(1-$site2_DAF) ;
							}
						}
						if( $gene1 eq $gene2 ){
							push @{$SynR_Intragenic_ResamplingBins{$quant}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
							push @{$SynDpr_Intragenic_ResamplingBins{$quant}{$bin}}, ($LD)/$Dmax  ;
						}else{
							push @{$SynR_Intergenic_ResamplingBins{$quant}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
							push @{$SynDpr_Intergenic_ResamplingBins{$quant}{$bin}}, ($LD)/$Dmax  ;
						}
					}else{
						if($site1_DAF==0){
							print OUT "NOT A SEG SITE: ", $gene1, "\t", $site1, "\n";
						}elsif($site2_DAF==0){
							print OUT "NOT A SEG SITE: ", $gene2, "\t", $site2, "\n" ;
						}
						next ;
					}	
				}
			}
		}
	}

	($firstIndex) = sort {$a <=> $b} keys %$Biallelic_NonSyn_SNP_Ref_Pos ; # Indexed from 0 to (Num Syn SNPs)-1
	($lastIndex) = sort {$b <=> $a} keys %$Biallelic_NonSyn_SNP_Ref_Pos ; # Indexed from 0 to (Num Syn SNPs)-1
	for ($lower = $firstIndex; $lower<=($lastIndex-1); $lower++){
		my $gene1 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$lower}}[0] ;
		my $seg_site_index1 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$lower}}[1] ;
		my $site1 = $$core_ref_biallelic_segsites{$gene1}{$seg_site_index1} ;
		my $site1_ancestral_base = $$RefBase{$gene1}{$site1} ;
	
		for ($upper = ($lower+1); $upper<$lastIndex; $upper++){
			my $gene2 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$upper}}[0] ;
			if( $$DensityNonSynPoly_byGene_Quantiled{$gene1} == $$DensityNonSynPoly_byGene_Quantiled{$gene2} ){
				my $dist = ($$Biallelic_NonSyn_SNP_Ref_Pos{$upper}-$$Biallelic_NonSyn_SNP_Ref_Pos{$lower}) ;
				my $seg_site_index2 = ${$$Biallelic_NonSyn_SNP_Alignment_Coord{$upper}}[1] ;
				my $site2 = $$core_ref_biallelic_segsites{$gene2}{$seg_site_index2} ;
				my $site2_ancestral_base = $$RefBase{$gene2}{$site2} ;
				my $quant = $$DensityNonSynPoly_byGene_Quantiled{$gene1} ; # should be same as gene2
				my $bin = int($dist/$$interSNPdist_binSize) ;
				
				if( $bin <= $$interSNPdist_numBins-1 && $dist>0 ){
					my $site1_site2_dblDAF_haplo = 0 ;
					my $site1_DAF = 0 ;
					my $site2_DAF = 0 ;
					my $samp_size_for_pair = 0 ;
					foreach my $ind (keys %{$$core_loci_seqs{$gene1}}){
						if($ind ne $$Reference_ind){
							if( exists $$core_loci_seqs{$gene1}{$ind} && exists $$core_loci_seqs{$gene2}{$ind} ){
								$samp_size_for_pair++ ;
								my $b1 = uc(substr($$core_loci_seqs{$gene1}{$ind}, $site1, 1)) ;
								my $b2 = uc(substr($$core_loci_seqs{$gene2}{$ind}, $site2, 1)) ;
								if( $b1 ne $site1_ancestral_base){
									$site1_DAF++ ;
									if( $b2 ne $site2_ancestral_base ){
										$site1_site2_dblDAF_haplo++ ;
									}
								}
								if( $b2 ne $site2_ancestral_base ){
									$site2_DAF++ ;
								}
							}
						}
					}
					$site1_DAF = $site1_DAF/$samp_size_for_pair ;
					$site2_DAF = $site2_DAF/$samp_size_for_pair ;
					$site1_site2_dblDAF_haplo = $site1_site2_dblDAF_haplo/$samp_size_for_pair ;
					if($site1_DAF > 0 && $site2_DAF >0){				
						my $LD = $site1_site2_dblDAF_haplo - ($site1_DAF*$site2_DAF) ;
						# if Dab positive, Dmax = min( Pa*(1-Pb), Pb*(1-Pa) )
						# if Dab negative, Dmax = min( Pa*(Pb), (1-Pa)*(1-Pb) )
						my $Dmax = 0 ;
						if ($LD>=0){
							if( $site1_DAF*(1-$site2_DAF) < $site2_DAF*(1-$site1_DAF) ){
								$Dmax = $site1_DAF*(1-$site2_DAF) ;
							}else{
								$Dmax = $site2_DAF*(1-$site1_DAF) ;
							}
						}elsif($LD<0){
							if( ($site1_DAF*$site2_DAF) < (1-$site1_DAF)*(1-$site2_DAF) ){
								$Dmax = $site1_DAF*$site2_DAF ;
							}else{
								$Dmax = (1-$site1_DAF)*(1-$site2_DAF) ;
							}
						}
						if( $gene1 eq $gene2 ){
							push @{$NonsynR_Intragenic_ResamplingBins{$quant}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
							push @{$NonsynDpr_Intragenic_ResamplingBins{$quant}{$bin}}, ($LD)/$Dmax  ;
						}else{
							push @{$NonsynR_Intergenic_ResamplingBins{$quant}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
							push @{$NonsynDpr_Intergenic_ResamplingBins{$quant}{$bin}}, ($LD)/$Dmax  ;
						}
					}else{
						if($site1_DAF==0){
							print OUT "NOT A SEG SITE: ", $gene1, "\t", $site1, "\n";
						}elsif($site2_DAF==0){
							print OUT "NOT A SEG SITE: ", $gene2, "\t", $site2, "\n" ;
						}
						next ;
					}	
				}
			}
		}
	}

	# RANDOMIZE BINS SO RESAMPLING RANDOM 
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			if( @{$SynR_Intragenic_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$SynR_Intragenic_ResamplingBins{$quant}{$bin}}) };
			if( @{$SynR_Intergenic_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$SynR_Intergenic_ResamplingBins{$quant}{$bin}}) };
			if( @{$NonsynR_Intragenic_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$NonsynR_Intragenic_ResamplingBins{$quant}{$bin}}) };
			if( @{$NonsynR_Intergenic_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$NonsynR_Intergenic_ResamplingBins{$quant}{$bin}}) };

			if( @{$SynDpr_Intragenic_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$SynDpr_Intragenic_ResamplingBins{$quant}{$bin}}) };
			if( @{$SynDpr_Intergenic_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$SynDpr_Intergenic_ResamplingBins{$quant}{$bin}}) };
			if( @{$NonsynDpr_Intragenic_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$NonsynDpr_Intragenic_ResamplingBins{$quant}{$bin}}) };
			if( @{$NonsynDpr_Intergenic_ResamplingBins{$quant}{$bin}} ){ fisher_yates_shuffle(\@{$NonsynDpr_Intergenic_ResamplingBins{$quant}{$bin}}) };
		}
	}
	# COMPARE NONSYN_R AND SYN_R WITHIN AND BETWEEN GENES
	my %NonsynRvSynR_Intragenic ; # @{$NonsynRvSynR_Intragenic{quant}} = (NonsynR - SynR, ... )
	my %NonsynRvSynR_Intergenic ; # @{$NonsynRvSynR_Intergenic{quant}} = (NonsynR - SynR, ... )

	my %NonsynDprSynDpr_Intragenic ; # @{$NonsynRvSynR_Intragenic{quant}} = (NonsynR - SynR, ... )
	my %NonsynDprSynDpr_Intergenic ; # @{$NonsynRvSynR_Intergenic{quant}} = (NonsynR - SynR, ... )
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		# CALCULATE MEAN NonsynR-SynR FOR BOOTREP, COMBINING RESULTS ACROSS BINS 
		foreach my $rep ( 1 .. $$Bootreps ){
			my %NonsynR_minus_SynR_Intragenic  ; # Temporary structure for calculating weighted avg (NonsynR-SynR) across all bins, for each rep
			my %NonsynR_minus_SynR_Intergenic  ;
			my %NonsynDpr_minus_SynDpr_Intragenic  ; # Temporary structure for calculating weighted avg (NonsynR-SynR) across all bins, for each rep
			my %NonsynDpr_minus_SynDpr_Intergenic  ;
			my $Intragenic_TotalObs = 0 ;
			my $Intergenic_TotalObs = 0 ;
			foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
				# FIRST CALCULATE INTRAGENIC METRICS
				@{$NonsynR_minus_SynR_Intragenic{$bin}} = () ; # tuple containing mean(NonsynR-SynR) and num observations for that bin 
				@{$NonsynDpr_minus_SynDpr_Intragenic{$bin}} = () ; # tuple containing mean(NonsynR-SynR) and num observations for that bin 
				my $SynSNPpairs = scalar @{$SynR_Intragenic_ResamplingBins{$quant}{$bin}} ;
				my $NonsynSNPpairs = scalar @{$NonsynR_Intragenic_ResamplingBins{$quant}{$bin}} ;
				if( $SynSNPpairs && $NonsynSNPpairs ){ # Can't compare Syn with Nonsyn if one category has 0 observations
					if( $SynSNPpairs >= $NonsynSNPpairs ){ # Synonymous SNPs have more observations; sample them w/out replacement
						my @indices_left_syn = ( 0 .. $SynSNPpairs-1 ) ;
						my @Syn_SampWOreplacementIndices ;
						foreach ( 1 .. $NonsynSNPpairs ){ # Subsample same number of Syn Observations as Nonsyn observations you just took mean of
							my $tmp = int(rand( scalar @indices_left_syn )) ;
							push @Syn_SampWOreplacementIndices, $indices_left_syn[$tmp] ;
							splice @indices_left_syn, $tmp, 1 ;
						}
						my $SynMeanR = 0 ;
						my $SynMeanDpr = 0 ;
						foreach my $index ( @Syn_SampWOreplacementIndices ){
							$SynMeanR += ${$SynR_Intragenic_ResamplingBins{$quant}{$bin}}[$index] ;
							$SynMeanDpr += ${$SynDpr_Intragenic_ResamplingBins{$quant}{$bin}}[$index] ;
						}
						$SynMeanR = $SynMeanR/$NonsynSNPpairs ;
						$SynMeanDpr = $SynMeanDpr/$NonsynSNPpairs ;

						my $NonsynMeanR = 0 ;
						my $NonsynMeanDpr = 0 ;
						foreach ( @{$NonsynR_Intragenic_ResamplingBins{$quant}{$bin}} ){
							$NonsynMeanR += $_ ;
						}
						foreach ( @{$NonsynDpr_Intragenic_ResamplingBins{$quant}{$bin}} ){
							$NonsynMeanDpr += $_ ;
						}
						$NonsynMeanR = $NonsynMeanR/$NonsynSNPpairs ;
						$NonsynMeanDpr = $NonsynMeanDpr/$NonsynSNPpairs ;

						@{$NonsynR_minus_SynR_Intragenic{$bin}} = (($NonsynMeanR - $SynMeanR), $NonsynSNPpairs) ;
						@{$NonsynDpr_minus_SynDpr_Intragenic{$bin}} = (($NonsynMeanDpr - $SynMeanDpr), $NonsynSNPpairs) ;
						$Intragenic_TotalObs += $NonsynSNPpairs ;	
					}else{  # NonSynonymous SNPs have more observations; sample them w/out replacement
						my @indices_left_nonsyn = ( 0 .. $NonsynSNPpairs-1 ) ;
						my @Nonsyn_SampWOreplacementIndices ;
						foreach ( 1 .. $SynSNPpairs ){ # Subsample same number of NonSyn Observations as Syn observations you just took mean of
							my $tmp = int(rand( scalar @indices_left_nonsyn )) ;
							push @Nonsyn_SampWOreplacementIndices, $indices_left_nonsyn[$tmp] ;
							splice @indices_left_nonsyn, $tmp, 1 ;
						}
						my $NonsynMeanR = 0 ;
						my $NonsynMeanDpr = 0 ;
						foreach my $index ( @Nonsyn_SampWOreplacementIndices ){
							$NonsynMeanR += ${$NonsynR_Intragenic_ResamplingBins{$quant}{$bin}}[$index] ;
						}
						foreach my $index ( @Nonsyn_SampWOreplacementIndices ){
							$NonsynMeanDpr += ${$NonsynDpr_Intragenic_ResamplingBins{$quant}{$bin}}[$index] ;
						}
						$NonsynMeanR = $NonsynMeanR/$SynSNPpairs ;
						$NonsynMeanDpr = $NonsynMeanDpr/$SynSNPpairs ;

						my $SynMeanR = 0 ;
						my $SynMeanDpr = 0 ;
						foreach ( @{$SynR_Intragenic_ResamplingBins{$quant}{$bin}} ){
							$SynMeanR += $_ ;
						}
						foreach ( @{$SynDpr_Intragenic_ResamplingBins{$quant}{$bin}} ){
							$SynMeanDpr += $_ ;
						}
						$SynMeanR = $SynMeanR/$SynSNPpairs ;
						$SynMeanDpr = $SynMeanDpr/$SynSNPpairs ;
						@{$NonsynR_minus_SynR_Intragenic{$bin}} = (($NonsynMeanR - $SynMeanR), $SynSNPpairs) ;
						@{$NonsynDpr_minus_SynDpr_Intragenic{$bin}} = (($NonsynMeanDpr - $SynMeanDpr), $SynSNPpairs) ;
						$Intragenic_TotalObs += $SynSNPpairs ;	
					}
				}else{
					@{$NonsynR_minus_SynR_Intragenic{$bin}} = (0,0) ;
					@{$NonsynDpr_minus_SynDpr_Intragenic{$bin}} = (0,0) ;
				}
				# SECOND CALCULATE INTERGENIC METRICS
				@{$NonsynR_minus_SynR_Intergenic{$bin}} = () ; # tuple containing mean(NonsynR-SynR) and num observations for that bin 
				@{$NonsynDpr_minus_SynDpr_Intergenic{$bin}} = () ; # tuple containing mean(NonsynR-SynR) and num observations for that bin 
				$SynSNPpairs = scalar @{$SynR_Intergenic_ResamplingBins{$quant}{$bin}} ;
				$NonsynSNPpairs = scalar @{$NonsynR_Intergenic_ResamplingBins{$quant}{$bin}} ;
				if( $SynSNPpairs && $NonsynSNPpairs ){ # Can't compare Syn with Nonsyn if one category has 0 observations
					if( $SynSNPpairs >= $NonsynSNPpairs ){ # Synonymous SNPs have more observations; sample them w/out replacement
						my @indices_left_syn = ( 0 .. $SynSNPpairs-1 ) ;
						my @Syn_SampWOreplacementIndices ;
						foreach ( 1 .. $NonsynSNPpairs ){ # Subsample same number of Syn Observations as Nonsyn observations you just took mean of
							my $tmp = int(rand( scalar @indices_left_syn )) ;
							push @Syn_SampWOreplacementIndices, $indices_left_syn[$tmp] ;
							splice @indices_left_syn, $tmp, 1 ;
						}
						my $SynMeanR = 0 ;
						my $SynMeanDpr = 0 ;
						foreach my $index ( @Syn_SampWOreplacementIndices ){
							$SynMeanR += ${$SynR_Intergenic_ResamplingBins{$quant}{$bin}}[$index] ;
							$SynMeanDpr += ${$SynDpr_Intergenic_ResamplingBins{$quant}{$bin}}[$index] ;
						}
						$SynMeanR = $SynMeanR/$NonsynSNPpairs ;
						$SynMeanDpr = $SynMeanDpr/$NonsynSNPpairs ;
						
						my $NonsynMeanR = 0 ;
						my $NonsynMeanDpr = 0 ;
						foreach ( @{$NonsynR_Intergenic_ResamplingBins{$quant}{$bin}} ){
							$NonsynMeanR += $_ ;
						}
						foreach ( @{$NonsynDpr_Intergenic_ResamplingBins{$quant}{$bin}} ){
							$NonsynMeanDpr += $_ ;
						}
						$NonsynMeanR = $NonsynMeanR/$NonsynSNPpairs ;
						$NonsynMeanDpr = $NonsynMeanDpr/$NonsynSNPpairs ;
						@{$NonsynR_minus_SynR_Intergenic{$bin}} = (($NonsynMeanR - $SynMeanR), $NonsynSNPpairs) ;
						@{$NonsynDpr_minus_SynDpr_Intergenic{$bin}} = (($NonsynMeanDpr - $SynMeanDpr), $NonsynSNPpairs) ;
						$Intergenic_TotalObs += $NonsynSNPpairs ;	
					}else{  # NonSynonymous SNPs have more observations; sample them w/out replacement
						my @indices_left_nonsyn = ( 0 .. $NonsynSNPpairs-1 ) ;
						my @Nonsyn_SampWOreplacementIndices ;
						foreach ( 1 .. $SynSNPpairs ){ # Subsample same number of NonSyn Observations as Syn observations you just took mean of
							my $tmp = int(rand( scalar @indices_left_nonsyn )) ;
							push @Nonsyn_SampWOreplacementIndices, $indices_left_nonsyn[$tmp] ;
							splice @indices_left_nonsyn, $tmp, 1 ;
						}
						my $NonsynMeanR = 0 ;
						my $NonsynMeanDpr = 0 ;
						foreach my $index ( @Nonsyn_SampWOreplacementIndices ){
							$NonsynMeanR += ${$NonsynR_Intergenic_ResamplingBins{$quant}{$bin}}[$index] ;
						}
						foreach my $index ( @Nonsyn_SampWOreplacementIndices ){
							$NonsynMeanDpr += ${$NonsynDpr_Intergenic_ResamplingBins{$quant}{$bin}}[$index] ;
						}
						$NonsynMeanR = $NonsynMeanR/$SynSNPpairs ;
						$NonsynMeanDpr = $NonsynMeanDpr/$SynSNPpairs ;

						my $SynMeanR = 0 ;
						my $SynMeanDpr = 0 ;
						foreach ( @{$SynR_Intergenic_ResamplingBins{$quant}{$bin}} ){
							$SynMeanR += $_ ;
						}
						foreach ( @{$SynDpr_Intergenic_ResamplingBins{$quant}{$bin}} ){
							$SynMeanDpr += $_ ;
						}
						$SynMeanR = $SynMeanR/$SynSNPpairs ;
						$SynMeanDpr = $SynMeanDpr/$SynSNPpairs ;
						@{$NonsynR_minus_SynR_Intergenic{$bin}} = (($NonsynMeanR - $SynMeanR), $SynSNPpairs) ;
						@{$NonsynDpr_minus_SynDpr_Intergenic{$bin}} = (($NonsynMeanDpr - $SynMeanDpr), $SynSNPpairs) ;
						$Intergenic_TotalObs += $SynSNPpairs ;	
					}
				}else{
					@{$NonsynR_minus_SynR_Intergenic{$bin}} = (0,0) ;
					@{$NonsynDpr_minus_SynDpr_Intergenic{$bin}} = (0,0) ;
				}
				
			}
			# THEN TAKE WEIGHTED AVG OF MEANS
			my $IntragenicMeanR = 0 ;
			my $IntragenicMeanDpr = 0 ;
			foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
				$IntragenicMeanR += ${$NonsynR_minus_SynR_Intragenic{$bin}}[0]*${$NonsynR_minus_SynR_Intragenic{$bin}}[1] ;
				$IntragenicMeanDpr += ${$NonsynDpr_minus_SynDpr_Intragenic{$bin}}[0]*${$NonsynDpr_minus_SynDpr_Intragenic{$bin}}[1] ;
			}
			if($Intragenic_TotalObs){
				$IntragenicMeanR = $IntragenicMeanR/$Intragenic_TotalObs ;
				$IntragenicMeanDpr = $IntragenicMeanDpr/$Intragenic_TotalObs ;
			}else{
				$IntragenicMeanR = "NA" ;
				$IntragenicMeanDpr = "NA" ;
			}
			my $IntergenicMeanR = 0 ;
			my $IntergenicMeanDpr = 0 ;
			foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
				$IntergenicMeanR += ${$NonsynR_minus_SynR_Intergenic{$bin}}[0]*${$NonsynR_minus_SynR_Intergenic{$bin}}[1] ;
				$IntergenicMeanDpr += ${$NonsynDpr_minus_SynDpr_Intergenic{$bin}}[0]*${$NonsynDpr_minus_SynDpr_Intergenic{$bin}}[1] ;
			}
			if($Intergenic_TotalObs){
				$IntergenicMeanR = $IntergenicMeanR/$Intergenic_TotalObs ;
				$IntergenicMeanDpr = $IntergenicMeanDpr/$Intergenic_TotalObs ;
			}else{
				$IntergenicMeanR = "NA" ;
				$IntergenicMeanDpr = "NA" ;
			}

			push @{$NonsynRvSynR_Intragenic{$quant}}, $IntragenicMeanR ;
			push @{$NonsynRvSynR_Intergenic{$quant}}, $IntergenicMeanR ;
			push @{$NonsynDprSynDpr_Intragenic{$quant}}, $IntragenicMeanDpr ;
			push @{$NonsynDprSynDpr_Intergenic{$quant}}, $IntergenicMeanDpr ;
		}
	}
	open OUT1, ">${$wd}/NonsynRvSynR_Intragenic_Intergenic_BySNPdensityQuantile.txt" ;
	print OUT1 "Intra\/Inter\tQuantile\t0.05pQuant\t0.5pQuant\t0.95Quant\n" ;
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		@{$NonsynRvSynR_Intragenic{$quant}} = sort {$a <=> $b} @{$NonsynRvSynR_Intragenic{$quant}} ;
		my $lower90CI = int(scalar(@{$NonsynRvSynR_Intragenic{$quant}})*0.05) - 1 ;
		my $median = int(scalar(@{$NonsynRvSynR_Intragenic{$quant}})*0.5) - 1 ;
		my $upper90CI = int(scalar(@{$NonsynRvSynR_Intragenic{$quant}})*0.95) - 1 ;
		print OUT1 "Intragenic", "\t" ;
		print OUT1 $quant, "\t" ;
		print OUT1 ${$NonsynRvSynR_Intragenic{$quant}}[$lower90CI], "\t" ;
		print OUT1 ${$NonsynRvSynR_Intragenic{$quant}}[$median], "\t" ;
		print OUT1 ${$NonsynRvSynR_Intragenic{$quant}}[$upper90CI], "\n" ;

		@{$NonsynRvSynR_Intergenic{$quant}} = sort {$a <=> $b} @{$NonsynRvSynR_Intergenic{$quant}} ;
		$lower90CI = int(scalar(@{$NonsynRvSynR_Intergenic{$quant}})*0.05) - 1 ;
		$median = int(scalar(@{$NonsynRvSynR_Intergenic{$quant}})*0.5) - 1 ;
		$upper90CI = int(scalar(@{$NonsynRvSynR_Intergenic{$quant}})*0.95) - 1 ;
		print OUT1 "Intergenic", "\t" ;
		print OUT1 $quant, "\t" ;
		print OUT1 ${$NonsynRvSynR_Intergenic{$quant}}[$lower90CI], "\t" ;
		print OUT1 ${$NonsynRvSynR_Intergenic{$quant}}[$median], "\t" ;
		print OUT1 ${$NonsynRvSynR_Intergenic{$quant}}[$upper90CI], "\n" ;

	}
	close OUT1 ;

	open OUT1, ">${$wd}/NonsynDprSynDpr_Intragenic_Intergenic_BySNPdensityQuantile.txt" ;
	print OUT1 "Intra\/Inter\tQuantile\t0.05pQuant\t0.5pQuant\t0.95Quant\n" ;
	foreach my $quant ( 0 .. $$NumQuantiles-1 ){
		@{$NonsynDprSynDpr_Intragenic{$quant}} = sort {$a <=> $b} @{$NonsynDprSynDpr_Intragenic{$quant}} ;
		my $lower90CI = int(scalar(@{$NonsynDprSynDpr_Intragenic{$quant}})*0.05) - 1 ;
		my $median = int(scalar(@{$NonsynDprSynDpr_Intragenic{$quant}})*0.5) - 1 ;
		my $upper90CI = int(scalar(@{$NonsynDprSynDpr_Intragenic{$quant}})*0.95) - 1 ;
		print OUT1 "Intragenic", "\t" ;
		print OUT1 $quant, "\t" ;
		print OUT1 ${$NonsynDprSynDpr_Intragenic{$quant}}[$lower90CI], "\t" ;
		print OUT1 ${$NonsynDprSynDpr_Intragenic{$quant}}[$median], "\t" ;
		print OUT1 ${$NonsynDprSynDpr_Intragenic{$quant}}[$upper90CI], "\n" ;

		@{$NonsynDprSynDpr_Intergenic{$quant}} = sort {$a <=> $b} @{$NonsynDprSynDpr_Intergenic{$quant}} ;
		$lower90CI = int(scalar(@{$NonsynDprSynDpr_Intergenic{$quant}})*0.05) - 1 ;
		$median = int(scalar(@{$NonsynDprSynDpr_Intergenic{$quant}})*0.5) - 1 ;
		$upper90CI = int(scalar(@{$NonsynDprSynDpr_Intergenic{$quant}})*0.95) - 1 ;
		print OUT1 "Intergenic", "\t" ;
		print OUT1 $quant, "\t" ;
		print OUT1 ${$NonsynDprSynDpr_Intergenic{$quant}}[$lower90CI], "\t" ;
		print OUT1 ${$NonsynDprSynDpr_Intergenic{$quant}}[$median], "\t" ;
		print OUT1 ${$NonsynDprSynDpr_Intergenic{$quant}}[$upper90CI], "\n" ;

	}
	close OUT1 ;
	close OUT ;
}

sub BinLDbyDistForBoostrapResampling_ByGene{
	# This function is very similar to BinLDbyDistForBoostrapResampling_WholeGenome, with the exception that rN-rS is calculated by gene, not by bin

	my $interSNPdist_binSize = $_[0] ;
	my $interSNPdist_numBins = $_[1] ;
	my $core_ref_biallelic_segsites = $_[2] ;
	my $RefBase = $_[3] ;
	my $core_loci_seqs = $_[4] ;
	my $Reference_ind = $_[5] ;
	my $Functional_Effect_BiAllelic_HR = $_[6] ;
	my $min_AF = $_[7] ;
	my $max_AF = $_[8] ;
 	my $Bootreps = $_[9] ;
	my $Genes_Quantiled = $_[10] ;
    	my $wd = $_[11] ;
     
	open OUT, ">>${$wd}/QCfile.txt" ;     

	# Synonynomous LD binned by interSNP distance, for bootstrap resampling
	my %SynR_ByGene_ResamplingBins ; # @{$SynR_ByGene_ResamplingBins{$gene}{$bin}}
	# Nonsynonymous LD binned by interSNP distance
	my %NonsynR_ByGene_BinByDist ;
	my %TotalPairwiseComparisons ; # $TotalPairwiseComparisons{$gene} = $Intragenic_TotalObs/$$Bootreps ;
	my %PairwiseComparisonsSynNonsyn ;

	# INITIALIZE
	foreach my $gene ( keys %$core_ref_biallelic_segsites ){
		$TotalPairwiseComparisons{$gene} = 0 ;
		$PairwiseComparisonsSynNonsyn{$gene}{"SYN"} = 0 ;
		$PairwiseComparisonsSynNonsyn{$gene}{"NONSYN"} = 0 ;
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			@{$SynR_ByGene_ResamplingBins{$gene}{$bin}} = () ; 
			@{$NonsynR_ByGene_BinByDist{$gene}{$bin}} = () ;  
		}    
	}

	# MAKE RESAMPLING BINS
    	foreach my $gene (sort {$a cmp $b} keys %$core_ref_biallelic_segsites){
		if(scalar keys %{$$core_ref_biallelic_segsites{$gene}} >=2 ){
			(my $last_Segsite_index) = sort{$b <=> $a} keys %{$$core_ref_biallelic_segsites{$gene}} ;
			for (my $i=0; $i < $last_Segsite_index; $i++){
				my $site1 = $$core_ref_biallelic_segsites{$gene}{$i} ;
				my $site1_ancestral_base = $$RefBase{$gene}{$site1} ;
				for (my $j = $i+1; $j<=$last_Segsite_index; $j++){
					my $site2 = $$core_ref_biallelic_segsites{$gene}{$j} ;
					my $site2_ancestral_base = $$RefBase{$gene}{$site2} ;
					my $dist = ($site2-$site1) ;
					if($dist <=0){
						print "Warning: Detected negative distance\n" ;
					}
					my $bin = int($dist/$$interSNPdist_binSize) ;

					if( $bin <= $$interSNPdist_numBins-1 && $dist>0 ){
						my $site1_site2_dblDAF_haplo = 0 ;
						my $site1_DAF = 0 ;
						my $site2_DAF = 0 ;
						my $samp_size_for_pair = 0 ;
						foreach my $ind (keys %{$$core_loci_seqs{$gene}}){
							if($ind ne $$Reference_ind){
								$samp_size_for_pair++ ;
								my $b1 = uc(substr($$core_loci_seqs{$gene}{$ind}, $site1, 1)) ;
								my $b2 = uc(substr($$core_loci_seqs{$gene}{$ind}, $site2, 1)) ;
								if( $b1 ne $site1_ancestral_base){
									$site1_DAF++ ;
									if( $b2 ne $site2_ancestral_base ){
											$site1_site2_dblDAF_haplo++ ;
									}
								}
								if( $b2 ne $site2_ancestral_base ){
									$site2_DAF++ ;
								}
							}
						}
						$site1_DAF = $site1_DAF/$samp_size_for_pair ;
						$site2_DAF = $site2_DAF/$samp_size_for_pair ;
						$site1_site2_dblDAF_haplo = $site1_site2_dblDAF_haplo/$samp_size_for_pair ;

						if( ($site1_DAF > $$min_AF && $site1_DAF < $$max_AF) && ($site2_DAF > $$min_AF && $site2_DAF < $$max_AF)){
							my $LD = $site1_site2_dblDAF_haplo - ($site1_DAF*$site2_DAF) ;
							if($$Functional_Effect_BiAllelic_HR{$gene}{$site1} eq "SYNONYMOUS" && $$Functional_Effect_BiAllelic_HR{$gene}{$site2} eq "SYNONYMOUS"){
								push @{$SynR_ByGene_ResamplingBins{$gene}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF)) ;
							}elsif($$Functional_Effect_BiAllelic_HR{$gene}{$site1} =~ m/NONSYNONYMOUS/ && $$Functional_Effect_BiAllelic_HR{$gene}{$site2} =~ m/NONSYNONYMOUS/){
								push @{$NonsynR_ByGene_BinByDist{$gene}{$bin}}, ($LD)/sqrt($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF)) ;
							}
						}
					}
				}
			}
		}
        }

	# RANDOMIZE BINS SO RESAMPLING RANDOM 
	foreach my $gene (sort {$a cmp $b} keys %$core_ref_biallelic_segsites){
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			if( @{$SynR_ByGene_ResamplingBins{$gene}{$bin}} ){ fisher_yates_shuffle(\@{$SynR_ByGene_ResamplingBins{$gene}{$bin}}) };
			if( @{$NonsynR_ByGene_BinByDist{$gene}{$bin}} ){ fisher_yates_shuffle(\@{$NonsynR_ByGene_BinByDist{$gene}{$bin}}) };
		}
	}
                              
      # COMPARE NONSYN_R AND SYN_R WITHIN AND BETWEEN GENES
        my %NonsynRvSynR_ByGene_Overall ; # @{$NonsynRvSynR_ByGene_Overall{gene}} = (NonsynR - SynR, ... )
	foreach my $gene (sort {$a cmp $b} keys %$core_ref_biallelic_segsites){
			@{$NonsynRvSynR_ByGene_Overall{$gene}} = () ;
        }
	foreach my $gene (sort {$a cmp $b} keys %$core_ref_biallelic_segsites){
		# CALCULATE MEAN NonsynR-SynR FOR BOOTREP, COMBINING RESULTS ACROSS BINS 
		foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
			$PairwiseComparisonsSynNonsyn{$gene}{"SYN"} += scalar @{$SynR_ByGene_ResamplingBins{$gene}{$bin}} ;
			$PairwiseComparisonsSynNonsyn{$gene}{"NONSYN"} += scalar @{$NonsynR_ByGene_BinByDist{$gene}{$bin}} ;
		}	
		foreach my $rep ( 1 .. $$Bootreps ){
			my %NonsynR_minus_SynR_ByBin  ; # Temporary structure for calculating weighted avg (NonsynR-SynR) across all bins, for each rep
			my $Intragenic_TotalObs = 0 ;
			foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
				# FIRST CALCULATE INTRAGENIC METRICS
				@{$NonsynR_minus_SynR_ByBin{$bin}} = () ; # tuple containing mean(NonsynR-SynR) and num observations for that bin 
				my $SynSNPpairs = scalar @{$SynR_ByGene_ResamplingBins{$gene}{$bin}} ;
				my $NonsynSNPpairs = scalar @{$NonsynR_ByGene_BinByDist{$gene}{$bin}} ;
				if( $SynSNPpairs && $NonsynSNPpairs ){ # Can't compare Syn with Nonsyn if one category has 0 observations
					if( $SynSNPpairs >= $NonsynSNPpairs ){ # Synonymous SNPs have more observations; sample them w/out replacement
						my @indices_left_syn = ( 0 .. $SynSNPpairs-1 ) ;
						my @Syn_SampWOreplacementIndices ;
						foreach ( 1 .. $NonsynSNPpairs ){ # Subsample same number of Syn Observations as Nonsyn observations you just took mean of
							my $tmp = int(rand( scalar @indices_left_syn )) ;
							push @Syn_SampWOreplacementIndices, $indices_left_syn[$tmp] ;
							splice @indices_left_syn, $tmp, 1 ;
						}
						my $SynMean = 0 ;
						foreach my $index ( @Syn_SampWOreplacementIndices ){
							$SynMean += ${$SynR_ByGene_ResamplingBins{$gene}{$bin}}[$index] ;
						}
						$SynMean = $SynMean/$NonsynSNPpairs ;
						my $NonsynMean = 0 ;
						foreach ( @{$NonsynR_ByGene_BinByDist{$gene}{$bin}} ){
							$NonsynMean += $_ ;
						}
						$NonsynMean = $NonsynMean/$NonsynSNPpairs ;
						@{$NonsynR_minus_SynR_ByBin{$bin}} = (($NonsynMean - $SynMean), $NonsynSNPpairs) ;
						$Intragenic_TotalObs += $NonsynSNPpairs ;
					}else{  # NonSynonymous SNPs have more observations; sample them w/out replacement
						my @indices_left_nonsyn = ( 0 .. $NonsynSNPpairs-1 ) ;
						my @Nonsyn_SampWOreplacementIndices ;
						foreach ( 1 .. $SynSNPpairs ){ # Subsample same number of NonSyn Observations as Syn observations you just took mean of
							my $tmp = int(rand( scalar @indices_left_nonsyn )) ;
							push @Nonsyn_SampWOreplacementIndices, $indices_left_nonsyn[$tmp] ;
							splice @indices_left_nonsyn, $tmp, 1 ;
						}
						my $NonsynMean = 0 ;
						foreach my $index ( @Nonsyn_SampWOreplacementIndices ){
							$NonsynMean += ${$NonsynR_ByGene_BinByDist{$gene}{$bin}}[$index] ;
						}
						$NonsynMean = $NonsynMean/$SynSNPpairs ;
						my $SynMean = 0 ;
						foreach ( @{$SynR_ByGene_ResamplingBins{$gene}{$bin}} ){
							$SynMean += $_ ;
						}
						$SynMean = $SynMean/$SynSNPpairs ;
						@{$NonsynR_minus_SynR_ByBin{$bin}} = (($NonsynMean - $SynMean), $SynSNPpairs) ;
						$Intragenic_TotalObs += $SynSNPpairs ;
					}
				}else{
					@{$NonsynR_minus_SynR_ByBin{$bin}} = (0,0) ;
				}
				$TotalPairwiseComparisons{$gene} += $Intragenic_TotalObs ;
			}
			# THEN TAKE WEIGHTED AVG OF MEANS
			my $IntragenicMean = 0 ;
			foreach my $bin ( 0 .. $$interSNPdist_numBins-1 ){
				$IntragenicMean += ${$NonsynR_minus_SynR_ByBin{$bin}}[0]*${$NonsynR_minus_SynR_ByBin{$bin}}[1] ;
			}
			if($Intragenic_TotalObs){
				$IntragenicMean = $IntragenicMean/$Intragenic_TotalObs ;
				push @{$NonsynRvSynR_ByGene_Overall{$gene}}, $IntragenicMean ;
			}
		}
	}

        open OUT1, ">${$wd}/NonsynRvSynR_ByGene_Overall.txt" ;
	print OUT1 "Gene\tQuantile\tSynPairs\tNonsynPairs\t0.05pQuant\t0.5Quant\t0.95Quant\n" ;
        foreach my $gene (sort {$a cmp $b} keys %TotalPairwiseComparisons){
		@{$NonsynRvSynR_ByGene_Overall{$gene}} = sort {$a <=> $b} @{$NonsynRvSynR_ByGene_Overall{$gene}} ;
		my $lower90CI = int(scalar(@{$NonsynRvSynR_ByGene_Overall{$gene}})*0.05) - 1 ;
		my $median = int(scalar(@{$NonsynRvSynR_ByGene_Overall{$gene}})*0.5) - 1 ;
		my $upper90CI = int(scalar(@{$NonsynRvSynR_ByGene_Overall{$gene}})*0.95) - 1 ;
        	print OUT1 $gene, "\t", $$Genes_Quantiled{$gene}, "\t"  ;
		print OUT1 $PairwiseComparisonsSynNonsyn{$gene}{"SYN"}, "\t", $PairwiseComparisonsSynNonsyn{$gene}{"NONSYN"}, "\t" ;
		if(${$NonsynRvSynR_ByGene_Overall{$gene}}[$lower90CI]){
			print OUT1 ${$NonsynRvSynR_ByGene_Overall{$gene}}[$lower90CI], "\t" ;
		}else{
			print OUT1"NA", "\t" ;
		}
		if(${$NonsynRvSynR_ByGene_Overall{$gene}}[$median]){
			print OUT1 ${$NonsynRvSynR_ByGene_Overall{$gene}}[$median], "\t" ;
		}else{
			print OUT1 "NA", "\t" ;
		}
		if(${$NonsynRvSynR_ByGene_Overall{$gene}}[$upper90CI]){
			print OUT1 ${$NonsynRvSynR_ByGene_Overall{$gene}}[$upper90CI], "\n" ;
		}else{
			print OUT1 "NA", "\n" ;
		}
        }    
	close OUT1 ;
        close OUT ;    
}          

sub Construct_Codon_hash{
    # Loads the genetic code into a Hash for checking the functional effect of mutations

    my %hash ;
    $hash{"ATT"} = "I" ;
    $hash{"ATC"} = "I" ;
    $hash{"ATA"} = "I" ;
    
    $hash{"CTT"} = "L" ;
    $hash{"CTC"} = "L" ;
    $hash{"CTA"} = "L" ;
    $hash{"CTG"} = "L" ;
    $hash{"TTA"} = "L" ;
    $hash{"TTG"} = "L" ;
    
    $hash{"GTT"} = "V" ;
    $hash{"GTC"} = "V" ;
    $hash{"GTA"} = "V" ;
    $hash{"GTG"} = "V" ;
    
    $hash{"TTT"} = "F" ;
    $hash{"TTC"} = "F" ;
    
    $hash{"ATG"} = "M" ;
    
    $hash{"TGT"} = "C" ;
    $hash{"TGC"} = "C" ;
    
    $hash{"GCT"} = "A" ;
    $hash{"GCC"} = "A" ;
    $hash{"GCA"} = "A" ;
    $hash{"GCG"} = "A" ;
    
    $hash{"GGT"} = "G" ;
    $hash{"GGC"} = "G" ;
    $hash{"GGA"} = "G" ;
    $hash{"GGG"} = "G" ;
    
    $hash{"CCT"} = "P" ;
    $hash{"CCC"} = "P" ;
    $hash{"CCA"} = "P" ;
    $hash{"CCG"} = "P" ;
    
    $hash{"ACT"} = "T" ;
    $hash{"ACC"} = "T" ;
    $hash{"ACA"} = "T" ;
    $hash{"ACG"} = "T" ;
    
    $hash{"TCT"} = "S" ;
    $hash{"TCC"} = "S" ;
    $hash{"TCA"} = "S" ;
    $hash{"TCG"} = "S" ;
    $hash{"AGT"} = "S" ;
    $hash{"AGC"} = "S" ;
    
    $hash{"TAT"} = "Y" ;
    $hash{"TAC"} = "Y" ;
    
    $hash{"TGG"} = "W" ;
    
    $hash{"CAA"} = "Q" ;
    $hash{"CAG"} = "Q" ;
    
    $hash{"AAT"} = "N" ;
    $hash{"AAC"} = "N" ;
    
    $hash{"CAT"} = "H" ;
    $hash{"CAC"} = "H" ;
    
    $hash{"GAA"} = "E" ;
    $hash{"GAG"} = "E" ;
    
    $hash{"GAT"} = "D" ;
    $hash{"GAC"} = "D" ;
    
    $hash{"AAA"} = "K" ;
    $hash{"AAG"} = "K" ;
    
    $hash{"CGT"} = "R" ;
    $hash{"CGC"} = "R" ;
    $hash{"CGA"} = "R" ;
    $hash{"CGG"} = "R" ;
    $hash{"AGA"} = "R" ;
    $hash{"AGG"} = "R" ;
    
    $hash{"TAA"} = "Stop" ;
    $hash{"TAG"} = "Stop" ;
    $hash{"TGA"} = "Stop" ;
    
    return (%hash) ;
}

sub bc {
    # Calculates a binomial coefficient, "n choose k"

    my ($n,$k) = @_;
    my $r=1;
    $r*=$n/($n-$k),$n--while$n>$k;
    $r;
}

sub fisher_yates_shuffle {
    # Shuffles an array, randomizing the order of the elements

    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }    
}


1;

