#! /usr/bin/perl
#Filename: cuffdiffwebsit.pl
use strict;
use warnings;
use autodie;
use feature 'say';
use CGI ':standard';
use CGI::Carp 'fatalsToBrowser';
use utility;

#use Data::Dumper;

print header;#starts the HTML script
print start_html('\\The variants\\'),#generates the header of the site
h1('Expression Analysis'),
p,
"This webform will enable you to search for genes in a certain region and if they meet certain expression criteria.",
p,
hr,
start_multipart_form, #start the form requesting all the infromation to be searched for
"Upload your cuffdiff file here as tab delimited file:",br,
filefield(-name => 'uploaded_file'),hr, #call for the file that content will be read into an array an passed into a subroutinethe rest of the fields queries all the pramaters
"Enter chromosome number, start and end position here:",br,
textfield(-name =>'chromosome', -rows=>1, -cols=>2),
textfield(-name =>'startposition', -rows=>1, -cols=>10),
textfield(-name => 'endposition', -rows=>1, -cols=>10),
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
	my $filename = upload('uploaded_file');
	my @filecontent = <$filename>;# reads the file content into an array
	if ( $filename){
		 $hash = cuffdiff_parser_HTML(@filecontent);# passes the file content as an array into the subroutine
	}else{
	print "Please enter filename to upload cuffdiff file",hr;
	}
	my $chr = param('chromosome');
	my $start = param('startposition');
	my $end = param('endposition');
	my @genesinregion;
	@genesinregion = genesinregion($chr, $start, $end, $hash);#generates the an array of genes 
	if (scalar @genesinregion != 0){
	print h2("Genes with the region chr$chr:$start-$end!");
	br,
	p,
 	print_array_HTML(@genesinregion);
	p;	
 	}else{
	say "No genes in the incated region";
	br,
	}
	my $p_value = param('p_value');
	my $fold_change = param('fold_change');
	my @siggenes= siggenes($p_value, $fold_change, $hash);#generates an array of genes
	if (@siggenes){
	print h2("Genes that are differentially express: p_value $p_value, log2_fold_change $fold_change");
	br,
	p,
	print_array_HTML(@siggenes);
	p,
	print hr;
	}else{
	say "No genes fullfill this criteria";
	br,
	}
	my @combined_genes = array_overlap(\@siggenes, \@genesinregion);
	if(@combined_genes){
	print h2("These are the differentailly regulated genes in the selected region");
	br,
	p,
	print_array_HTML(@combined_genes);
	p,
	}else{
	hr,
	say "No genes in the region chr$chr:$start-$end are differentially expressed";
	hr;
	}		
}
print end_html;






 
