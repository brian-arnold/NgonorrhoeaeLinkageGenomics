#!usr/bin/perl
use warnings ;
use strict ;

########################
# This script takes a progressive mauve alignment for 2 species and converts the information into an allele table
# that is considerably more easy to parse in downstream scripts. One of the 2 genomes used in the ProgMauve alignment
# is treated as a "reference", such that the output is indexed with respect to the positions in that sequence. 
# For instance, here I use the N. gonorrhoeae FA1090 sequence as a reference, and the output is arranged such that 
# positions correspond to those in the FA1090 reference genome. If the outgroup species contains DNA not present in 
# the reference, inducing gaps,these positions are output as "NA".
#
# The output consist of the position in the reference genome (column 1), the base in the reference genome used as input
# (column 2), the base of the reference genome in the progMauve alignment (column 3), and the base of the outgorup
# genome used in the progMauve alignment (column 4)
########################

########################
### LOAD FULL REFERENCE GENOME
### This is the same N. gonorrhoeae FA1090 reference used in the progressive Mauve alignment
### It is included in the output file as a sanity check, 2nd column of output, should equal base
### present in 3rd column
########################
open IN, "<./NgonorrhoeaeFA1090_genomic.fna" ;
my $Ref_genome_full = "" ;
while(<IN>){
	chomp $_ ;
	if($_ !~ m/^>/){
		$Ref_genome_full = $Ref_genome_full.${_} ;
	}
}
close IN ;

my $genome_OutGroup = "NmeningitidisAlpha14" ;
my $genome_Reference = "NgonorrhoeaeFA1090" ;
my $contig = 1 ;
my $current_genome ;
my %MultiSpeciesAlign_Seq ;
my %MultiSpeciesAlign_Pos ;
my %MultiSpeciesAlign_Pos2 ;

########################
### LOAD MULTI SPECIES ALIGNMENT
########################
open IN, "<./NgonorrhoeaeFA1090_NmeningitidisAlpha14_progMauve.alignment" ;
while(<IN>){
	chomp $_ ;
	if($_ =~ m/^#/){
		next ;
	}elsif($_ =~ m/=/){
		$contig++ ; 
		next ;
	}elsif($_ =~ m/^>/){
		if($_ =~ m/${genome_OutGroup}/){
			$current_genome = $genome_OutGroup ;
		}elsif($_ =~ m/${genome_Reference}/){
			$current_genome = $genome_Reference ;
		}
		my @line1 = split(/\s/, $_) ;
		my @line2 = split(/:/, $line1[1]) ;# file structure == contig:coord1-coord2
		my @line3 = split(/-/, $line2[1]) ;
		@{$MultiSpeciesAlign_Pos{$contig}{$current_genome}} = @line3 ; #lower and upper coordinates
		@{$MultiSpeciesAlign_Pos2{$current_genome}{$line3[0]}} = ($contig, $line3[1]) ;
		
		next ;
	}else{
		if(exists $MultiSpeciesAlign_Seq{$contig}{$current_genome} ){
			$MultiSpeciesAlign_Seq{$contig}{$current_genome} = $MultiSpeciesAlign_Seq{$contig}{$current_genome}.${_} ;
		}else{
			$MultiSpeciesAlign_Seq{$contig}{$current_genome} = $_ ;
		}
	}
}
close IN ;

open OUT, ">AlleleTable.txt" ;
print OUT "RefPosition", "\t", "RefGenomeAllele", "\t", "RefGenomeProgMauveAllele", "\t", "OutgroupGenomeProgMauveAllele", "\n" ;
foreach my $start (sort {$a <=> $b} keys %{$MultiSpeciesAlign_Pos2{$genome_Reference}} ){
	my $con = ${$MultiSpeciesAlign_Pos2{$genome_Reference}{$start}}[0] ;
	my $end = ${$MultiSpeciesAlign_Pos2{$genome_Reference}{$start}}[1] ;
	my $len = length($MultiSpeciesAlign_Seq{$con}{$genome_Reference}) ;
	my $num_dashes = 0 ;
	my $curr_position_in_reference ;
	foreach my $x ( 0..$len-1 ){
		$curr_position_in_reference = ($start+$x)-$num_dashes ; # keep track of position in NgonorrhoeaeFA1090 reference genome, so skip dashes
		my $refNotDash = 0 ;
		if( substr($MultiSpeciesAlign_Seq{$con}{$genome_Reference}, $x, 1) !~ m/-/ ){
			$refNotDash = 1 ;
		}

		if( $refNotDash ){
			print OUT $curr_position_in_reference, "\t" ;
		}else{
			$num_dashes++ ;
			print OUT "NA\t" ;
		}
		$curr_position_in_reference = ($start+$x)-$num_dashes ;
		if( $refNotDash ){
			print OUT substr($Ref_genome_full, $curr_position_in_reference-1, 1 ), "\t" ;
		}else{
			print OUT "NA\t" ;
		}
		print OUT substr($MultiSpeciesAlign_Seq{$con}{$genome_Reference}, $x, 1), "\t" ;
		if( exists $MultiSpeciesAlign_Seq{$con}{$genome_OutGroup} ){
			print OUT substr($MultiSpeciesAlign_Seq{$con}{$genome_OutGroup}, $x, 1), "\n" ;
		}else{
			print OUT "-\n" ;
		}

	}
}
exit ;

