package utility;
#filename: utitlity.pm
use strict;
use warnings;
use autodie;
use feature 'say';
use base 'Exporter';
use Data::Dumper;
our @EXPORT = qw(cuffdiff_parser genesinregion siggenes print_IDs_to_file cuffdiff_parser_HTML array_overlap print_array print_array_HTML) ;

####sub cuffdiff_parser parses cuffdiff files into a hash and returns reference to this hash 
sub cuffdiff_parser{
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
	return $cuffdiffhash;
}

####sub cuffdiff_parser parses cuffdiff files that are passed as an array into a hash and returns reference to this hash 
sub cuffdiff_parser_HTML{
	my @filecontent = @_;
	my $cuffdiffhash ={}; #reference for the hash the data will be stort in
	my $firstline = shift @filecontent; #handles the first line and makes the keys for the cuffdiff hash
	chomp $firstline;
	my @keys =split /\t/,$firstline;
#while loop slurps up the data and puts it into the hash with gene as key with an hash reference as value
	foreach my $line (@filecontent){
	chomp $line;
	my @elements = split /\t/, $line;
	my %hash;
	@hash{@keys}=@elements;
	${$cuffdiffhash}{$elements[2]}= \%hash;
	};
	return $cuffdiffhash;
}

#sub genes in region returns the genes within a specific region on a chromosom, the input is the chromosome, start, endposition and the hash to search. Returns an array of ids of genes in the region

sub genesinregion {
	my @input = @_;
	my @genesinregion;	
	foreach my $genes (keys $input[3]){
		${$input[3]}{$genes}{locus} =~ /chr([0-9]*):([0-9]*)-([0-9]*)/; #extracts the chromosome, start and end position
		my $chr = $1;
		my $start = $2;
		my $end =$3;
		
if ($chr==$input[0] and ($start >= $input[1] and $end <= $input[2])){ #compars the position of the gene with the query position and adds genes within the region to an array
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
		if ($input[0] >= ${$input[2]}{$genes}{p_value} and $input[1] <= ${$input[2]}{$genes}{'log2(fold_change)'}){
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

#sub that prints ids to STDOUT taking a file handle and an array
sub print_array {
	my @genes = @_; 
	foreach my $gene (@genes){
	print "$gene\n";
	}
}
#sub prints HTML output with each id on a new line
sub print_array_HTML {
	my @genes = @_; 
	foreach my $gene (@genes){
	print "$gene \n";
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
1;
