#!/usr/bin/perl
# Written by Dr. Gulum Kosova
use strict;
use warnings;
use ProcessVCF;
use Data::Dumper;
use Getopt::Long;

open(OUT, '>', 'cleaned_variants.txt');
my $infile = shift;
my ($qual_cutoff, $dp_cutoff, $hwe_cutoff, $help); 
GetOptions("q|qual=f", \$qual_cutoff,
           "d|depth=f", \$dp_cutoff,
           "h|hwe=f", \$hwe_cutoff,
           "help" , \$help,);
my $usage = "script options: q=QUAL score (default=50), d=Read Depth (default=500), h=HWE-Pval (default=0.001)";
die $usage if $help;           


##Start Parser##
my ($vcfdata, $ids) = VCF_parse($infile);
print OUT '#__There were ', scalar @{$vcfdata}, " variants in the input VCF file\n";
print OUT '#__Following QC criteria were applied:', "\n", '#______Quality score cutoff: ', $qual_cutoff, "\n",  
      '#______Total read depth cutoff: ', $dp_cutoff, "\n", '#______HWE P-value cutoff: ', $hwe_cutoff, "\n" ; 


##Start genotype counter & HWE calculator##
my $vcf_hwe = geno_assess($vcfdata, $ids);
#print   Dumper  @{$vcf_hwe}[1];

##Start filtration##
my $vcf_cleaned = $vcf_hwe;
$vcf_cleaned  = filter_qual($vcf_cleaned,$qual_cutoff);
$vcf_cleaned = filter_dp($vcf_cleaned,$dp_cutoff);
$vcf_cleaned = filter_hwe($vcf_cleaned, $hwe_cutoff);

print OUT '#__', scalar @{$vcf_hwe}-scalar @{$vcf_cleaned}, " variants have been filtered based on the given QC criteria","\n", 
'#__',scalar @{$vcf_cleaned}, " variants listed below passed the QC filtering","\n",'#', "\n";

print OUT "#Chr:position \t SNP_ID \t Ref/Alt \t VAF \ ref/ref \t ref/alt \t alt/alt \t NoCallRate \t QUAL \t Total Depth \t HWE-Pval \n";


foreach my $var (@{$vcf_cleaned}){
    print OUT "$var->{chr}:$var->{pos} \t $var->{ID} \t $var->{REF}/$var->{ALT} \t $var->{INFO}->{AF} \t $var->{GENO_COUNT}->{'ref/ref'} \t $var->{GENO_COUNT}->{'ref/alt'} \t $var->{GENO_COUNT}->{'alt/alt'} \t $var->{NoCallRate} \t  $var->{QUAL} \t $var->{INFO}->{DP} \t $var->{HWE_Pval} \t $var->{INFO}->{AF} \n";
}
print "Done! \n\n";

