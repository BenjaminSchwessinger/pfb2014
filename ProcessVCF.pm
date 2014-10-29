package ProcessVCF;
#file: ProcessVCF.pm, written by Dr. Gulum Kosova
use strict;
use warnings;
use base 'Exporter';
use Statistics::Distributions;
use Getopt::Long;

our @EXPORT = qw (VCF_parse);
our @EXPORT_OK = qw (geno_assess filter_qual filter_dp filter_hwe);

sub VCF_parse {
    my $in = shift;

    open( IN,     '<', $in );
    open( HEADER, '>', 'header.txt' );

    my @ID;
    my @vcf;
    my @genotypes;
    while ( my $line = <IN> ) {
        chomp($line);
        if ( $line =~ /^##/ ) {
            print HEADER $line;
        }
        elsif ( $line =~ /^#/ ) {
            my $header2 = $line;
            print HEADER $line;
            @ID = split( /\t/, $line );
            @ID = splice( @ID, 9 );
        }
        else {
            push( @vcf, $line );
        }
    }

    my @content;
    foreach my $line (@vcf) {
        my @array = split( /\t/, $line );
        my %hash;
        $hash{'chr'}    = shift @array;
        $hash{'pos'}    = shift @array;
        $hash{'ID'}     = shift @array;
        $hash{'REF'}    = shift @array;
        $hash{'ALT'}    = shift @array;
        $hash{'QUAL'}   = shift @array;
        $hash{'FILTER'} = shift @array;

        ##INFO  COLUMN##
        my @infoarray = split( /;/, shift @array );
        foreach my $var (@infoarray) {
            next unless ( $var =~ /\=/ );
            $var =~ /(.*)=(.*)/;
            my $newkey   = $1;
            my $newvalue = $2;
            $hash{'INFO'}{$newkey} = $newvalue;
        }

        ## GENOTYPE COLUMN
        @genotypes = @array;
        @genotypes = splice( @genotypes, 1 );
        my $i = 0;
        foreach my $geno (@genotypes) {
            $geno =~
              /^([\d.]\/[\d.]) :* ([\d+,]*) :* (\d*) :* (\d*) :* ([\d+,]*)/x;
            $hash{'GENOTYPES'}{'GT'}{ $ID[$i] } = $1;
            $hash{'GENOTYPES'}{'AD'}{ $ID[$i] } = $2;
            $hash{'GENOTYPES'}{'DP'}{ $ID[$i] } = $3;
            $hash{'GENOTYPES'}{'GQ'}{ $ID[$i] } = $4;
            $hash{'GENOTYPES'}{'PL'}{ $ID[$i] } = $5;
            $i++;
        }
        push( @content, \%hash );
    }
    return ( \@content, \@ID );
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##this subroutine is to count genotypes##

sub geno_assess {
    my $vcfdata = shift;
    my $ids   = shift;
    open (OUT, '>', 'genotype_counts_HWE.txt');
    print OUT 
"Var \t r/r \t r/a \t a/a \t nc \t total \t nc_Rate \t Ref_Ct \t Alt_ct \t Alt_Freq \t HW_chi2 \t HW_Pval \n";

    foreach my $var ( @{$vcfdata} ) {
        my $rrct = 0;
        my $ract = 0;
        my $aact = 0;
        for ( my $i = 0 ; $i < scalar @{$ids} ; $i++ ) {
            if ( $var->{GENOTYPES}->{GT}->{ @{$ids}[$i] } eq '0/0' ) {
                $rrct++;
            }
            elsif ( $var->{GENOTYPES}->{GT}->{ @{$ids}[$i] } eq '0/1' ) {
                $ract++;
            }
            elsif ( $var->{GENOTYPES}->{GT}->{ @{$ids}[$i] } eq '1/1' ) {
                $aact++;
            }

           #elsif ($var->{GENOTYPES}->{GT}->{@{$ids}[$i]} eq './.') {$ncct++;}
        }
        my $total      = $rrct + $ract + $aact;
        my $ncct       = scalar @{$ids} - ($total);
        my $nc_rate    = $ncct / scalar @{$ids};
        my $ref_count  = ( ( 2 * $rrct ) + $ract );
        my $alt_count  = ( ( 2 * $aact ) + $ract );
        my $alt_freq   = ( ( 2 * $aact ) + $ract ) / ( 2 * $total );
        my $exp_homz   = ( ( 1 - $alt_freq )**2 ) * $total;
        my $exp_het    = ( ( 2 * $alt_freq * ( 1 - $alt_freq ) ) * $total );
        my $exp_althom = ( $alt_freq**2 ) * $total;
        my $hw_chi2;

        if ( $exp_homz == 0 || $exp_het == 0 || $exp_althom == 0 ) {
            $hw_chi2 = 0;
        }
        else {
            $hw_chi2 =
              ( ( ( $rrct - $exp_homz )**2 ) / $exp_homz ) +
              ( ( ( $ract - $exp_het )**2 ) / $exp_het ) +
              ( ( ( $aact - $exp_althom )**2 ) / $exp_althom );
        }
        my $chi2Pval = Statistics::Distributions::chisqrprob( 1, $hw_chi2 );
        
        $var->{'GENO_COUNT'}->{'ref/ref'} = $rrct;
        $var->{'GENO_COUNT'}->{'ref/alt'} = $ract;
        $var->{'GENO_COUNT'}->{'alt/alt'} = $aact;
        $var->{'GENO_COUNT'}->{'n/n'} = $ncct;
        $var->{'HWE_Pval'} = $chi2Pval; 
        $var->{NoCallRate} = $nc_rate;
        
        print OUT "$var->{chr}:$var->{pos}", "\t", $rrct, "\t", $ract, "\t", $aact,
          "\t", $ncct, "\t", $total, "\t", $nc_rate, "\t",
          , $ref_count, "\t", $alt_count, "\t", $alt_freq, "\t", $hw_chi2, "\t",
          $chi2Pval, "\n";
    }
return ($vcfdata);
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## These next subroutines are for filtering variants for desired QC metrics ##
## Filter: QUAL < 50##

sub filter_qual {
    my $vcf = shift;
    my $cutoff = shift || 50;
    GetOptions("c|cutoff=f", \$cutoff);
    my @vcf_clean;
    foreach my $var ( @{$vcf} ) {
        if ( $var->{'QUAL'} >= $cutoff ) {
            push (@vcf_clean, $var);
        }
     }
     return \@vcf_clean;       
}        

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Filter: DP <500 ##
sub filter_dp {
    my $vcf = shift;
    my $cutoff = shift || 500;
    GetOptions("c|cutoff=f", \$cutoff);
    my @vcf_clean;
    foreach my $var ( @{$vcf} ) {
        if ( $var->{'INFO'}->{'DP'} >= $cutoff ) {
            push (@vcf_clean, $var);
        }
    }
    return \@vcf_clean;
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Filter: hwe > 0.001 ##
sub filter_hwe {
    my $vcf = shift;
    my $cutoff = shift || 0.001;
    GetOptions("c|cutoff=f", \$cutoff);
    my @vcf_clean;
    foreach my $var ( @{$vcf} ) {
        if ( $var->{'HWE_Pval'} >= $cutoff ) {
            push (@vcf_clean, $var);
        }
    }
return \@vcf_clean;
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

1;
