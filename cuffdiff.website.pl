#! /usr/bin/perl
#Filename: cuffdiffwebsit.pl
use strict;
use warnings;
use autodie;
use feature 'say';
use CGI ':standard';
use CGI::Carp 'fatalsToBrowser';
use utility;

use Data::Dumper;

print header;
print start_html('\\The variants\\'),
h1('Expression Analysis'),
p,
"This webform will enable you to search for genes in a certain region and if they meet certain expression criteria.",
p,
hr,
start_multipart_form,
"Upload your cuffdiff file here as tab delimited file:",br,
filefield(-name => 'uploaded_file'),hr,
"Enter chromosome number, start and end position here:",br,
textarea(-name =>'chromosome', -rows=>1, -cols=>1),
textarea(-name =>'startposition', -rows=>1, -cols=>10),
textarea(-name => 'endposition', -rows=>1, -cols=>10),
hr,
"Enter p_value here:",br,
textarea(-name => 'p_value',-rows =>1, -cols=>10),br,
"Enter log2 fold-change here:",br,
textarea(-name => 'fold_change', -rows =>1, -cols=>10),br,
submit('Submit'),


endform,
hr;

if ( param ){
	my $hash;
	if (my $filename = param('uploaded_file')){
		 $hash = cuffdiff_parser($filename);
	}else{
	print "Please enter filename to upload cuffdiff file",hr;
	}
	my $chr = param('chromosome');
	my $start = param('startposition');
	my $end = param('endposition');
	my @genesinregion;
	@genesinregion = genesinregion($chr, $start, $end, $hash);
	if (scalar @genesinregion != 0){
	print h2("Genes with the region chr$chr:$start-$end!");
	br,
	p,
 	print_array_HTML(@genesinregion);
	p,	
 	}else{
	say "No genes in the incated region";
	br,
	}
	my $p_value = param('p_value');
	my $fold_change = param('fold_change');
	my @siggenes= siggenes($p_value, $fold_change, $hash);
	if (@siggenes){
	print h2("Genes that are differentially express: p_value $p_value, log2_fold_change $fold_change");
	br,
	p,
	print_array_HTML(@siggenes);
	}else{
	say "No genes fullfill this criteria";
	br,
	}
}
print end_html;






 
