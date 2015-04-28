#!usr/bin/perl

##################################################
# qual_filter_vcf.pl v1.0
# filters out stuff from VCF files for further analysis
# 
# 2012 - cristel thomas
##################################################

use strict; use warnings;

die "usage <*.filt.vcf><min qual><min DP><max DP>\n" unless @ARGV==4;

my ($file,$qual,$dpmin,$dpmax)=@ARGV;
my $dp=0;
my $i =0;
my $j=0;

my ($vcf)=$file=~/(.{1,})\Q.\Efilt\Q.\Evcf/;

open(IN, "$vcf.filt.vcf") or die "can't open input\n";
open (OUT, ">$vcf.filtered.vcf") or die "can't open output\n";
open (OUT1, ">$vcf.multiallelic_sites") or die "can't open output 2\n";

while (<IN>){
	chomp;
	$j++;
	my ($bl, $qu, $reste) = (split /\t/) [4,5,9];
	if ($bl =~/,/){$i++;print OUT1 $_}
	if ($bl =~/[CTAGN]{1,}/){
		$dp=(split /:/, $reste)[2]; #die "ohhh$dp";
	}
	else {
		$dp=(split /:/, $reste)[1]; #die "$dp";
	}

	if ($qu > $qual) {
                if ($dp>$dpmin){
                        if ($dp<$dpmax) {print OUT $_,"\n"}
                }
        }
}
close OUT;
close IN;
close OUT1;

print "$j records, $i multi\n";

if ($i==0){system("rm $vcf.multiallelic_sites")}

__END__
version note:
version 2: merged with count_multi.pl

