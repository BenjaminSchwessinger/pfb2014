#! /usr/bin/perl
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;

#This will be a script that parses cuffdiff data into a hash starting from a tab delimited file 
#the user can specify different parameter  as shown in the options. The programs returns files containing 
# the gene ids of all genes matching these parameters individually and in combination. A second file containing the coordinates of these gene



my $qchr = 0;
my $qstartposition = 0;
my $qendposition = 0;
my $p_value= 0;
my $fold_change= 0;


GetOptions(	"chr=i" => \$qchr,
		"start|s=i" => \$qstartposition,
		"end|e=i" => \$qendposition,
		"p_value|p=f" => \$p_value,
		"fold_change|f=f" => \$fold_change,
	);

my $usage = "$0 [--chr] [--start] [--end] [--p_value] [--fold_change]";

unless($qchr and $qstartposition and $qendposition and $p_value and $fold_change){
	warn "$usage\n";
	exit(1);
	}

my $filename = shift;

open my $IN, '<', $filename;

my $cuffdiffhash ={}; #reference for the hash the data will be stort in
my $firstline = <$IN>; #handles the first line and makes the keys for the cuffdiff hash
chomp $firstline;
my @keys =split /\t/,$firstline;
#while loop slurps up the data and puts it into the hash with gene as key with an hash reference as value
while (my $line = <$IN>){
	chomp $line;
	my @elements = split /\t/, $line;
	my %hash;
	@hash{@keys}=@elements;
	${$cuffdiffhash}{$elements[2]}= \%hash;
};

my @genesinregion = genesinregion($qchr, $qstartposition, $qendposition, $cuffdiffhash); #getgenes in a region specificed on the command line with -chr -s -e
my @siggenes = siggenes($p_value, $fold_change, $cuffdiffhash); #get genes with a certain p_value and foldchange specificed on the command line with -p -f

my @combined_array = array_overlap (\@genesinregion, \@siggenes);#gets the genes that fullfill both criteria 

open my $OUT5, '>', "$filename.$qchr:$qstartposition-$qendposition-and-p$p_value.f$fold_change.txt"; #prints the list of genes in a region to file
print_IDs_to_file($OUT5, @combined_array);
close $OUT5;

open my $OUT, '>', "$filename.$qchr:$qstartposition-$qendposition.txt"; #prints the list of genes in a region to file
print_IDs_to_file($OUT, @genesinregion);
close $OUT;

open my $OUT2, '>', "$filename.p$p_value.f$fold_change.txt"; #prints differentially regulated genes to a file
print_IDs_to_file($OUT2, @siggenes);
close $OUT2;

my @combined_array_coord = coordinates($cuffdiffhash, @combined_array);

my @genesinregion_coord = coordinates($cuffdiffhash, @genesinregion);#get the coordinates for the genes within the selected region

my @siggenes_coord = coordinates($cuffdiffhash, @siggenes); #gets the coordinates for the gene signficantly altered 

open my $OUT4, '>', "$filename.$qchr:$qstartposition-$qendposition.coordinates.txt"; #prints the list of genes in a region to file
print_IDs_to_file($OUT4, @genesinregion_coord);
close $OUT4;


open my $OUT3, '>', "$filename.p$p_value.f$fold_change.coordinates.txt"; #prints differentially regulated genes coordinates to a file
print_IDs_to_file($OUT3, @siggenes_coord);
close $OUT3;


say "out of loop";
say 'next control step';
#####subs##########

#sub genes in region returns the genes within a specific region on a chromosom, the input is the chromosome, start, endposition and the hash to search. Returns an array of ids of genes in the region

sub genesinregion {
	my @input = @_;
	my @genesinregion;	
	foreach my $genes (keys $input[3]){
		${$input[3]}{$genes}{locus} =~ /chr([0-9]*):([0-9]*)-([0-9]*)/; #extracts the chromosome, start and end position
		my $chr = $1;
		my $start = $2;
		my $end =$3;
		
if ($chr==$input[0] and ($start >= $input[1] or $end <= $input[2])){ #compars the position of the gene with the query position and adds genes within the region to an array
		push @genesinregion, $genes;	
				}
	 }
	return @genesinregion;
}	


###sub siggenes returns the significant genes within a hash after conisdering fold change and p_value

sub siggenes {
	my @input = @_;
	my @siggens;
	foreach my $genes (keys $input[2]){
		if ($input[0] >= ${$input[2]}{$genes}{p_value} and $input[1] >= ${$input[2]}{$genes}{'log2(fold_change)'}){
			push @siggens, $genes;
		}
	}
	return @siggens;
}

#sub that prints ids to a file taking a file handle and an array
sub print_IDs_to_file {
	my $OUT = shift;
	my @genes = @_; 
	foreach my $gene (@genes){
	print $OUT "$gene \n";
	}
}


#sub my coordinates returns the coorinates in a parired array in which the first value of the pair contains the cchromosome number and the second the range. It takes a  list of genes and a hash containing the coordinates as locus key value pair

sub coordinates {
	my $hash = shift;
	my @genes = @_;
	my @coordinates;
	foreach my $genes(keys $hash){
		${$hash}{$genes}{locus} =~ /chr([0-9]*):([0-9]*)-([0-9]*)/; #extracts the chromosome, start and end position
		
my $coordinates = "$1\t$2-$3";
		push @coordinates, $coordinates;
		}	
	return @coordinates;	
			}
###sub array_overlap checks for the overlap of two arrays that are passed into the sub as a reference

sub array_overlap{
	my $array1 = shift;
	my $array2 = shift;
	my @combined_array;
	foreach my $id1(@{$array1}){
		foreach my $id2(@{$array2}){ 
		if ($id1 eq $id2){
		push @combined_array, $id1;
		}
		}
	}
	return @combined_array;
}

