#!/usr/bin/perl

################################################################
# gets total bp represented by different features and detected w/ confidence, SNPs and MAFs 
# 
# detected_genome.pl  ---  v1.1
# 2012 - cris 
################################################################
use strict;use warnings;

use Getopt::Std;
getopts('hvu:t:g:T:p:m:c:r:f:');
our($opt_h,$opt_v,$opt_g,$opt_r,$opt_f,$opt_c,$opt_u,$opt_t,$opt_T, $opt_p, $opt_m);

my $vers="version 1.2";
my $usage ="
usage: [options] -g <trimmed genome table> -u <UTR table> -t <TSS table> -T <fa.out from repeat masker> -m INT -c <contig.list> -f <filtered.comp>
options:
   -h help
   -v version
   
	INPUT FILES
   -c contigs list (req)
   -f filtered.comp (req)
   -g trimmed genome table (req)
   -t TSS table (req)
   -T RepeatMasker output .fa.out (req)
   -u UTR table (req)

	PARAMETERS
   -m minimum threshold for detected heterozygosity. 0 <= m <= 1 (by default, .2)
   -r reference a sample? yes/no (by default, yes)

	OUTPUT
   -p prefix (req)
";

my $help="
Gets distance detected with confidence in each type of feature, also get SNPs and MAFs.

Note: the minimum threshold for detected heterozygosity corresponds to the minimum proportion of total reads detected at that position supporting the least frequent allele.
Version note:
1.1 - takes a list of contigs (text file, ie remWS230_contig.lengths.txt.0) to look at. similar to option -c in briggsae version.
1.2 - added reference as a potential sample.
";

if ($opt_h) {
	print $usage;
	print $help;
	exit;
}

if ($opt_v) {
	print $vers;
	exit;
}
my $th=0.2;
if ($opt_m){
	if ($opt_m>1){die "$usage\n$help\nm should be between 0 and 1\n"}
	$th=$opt_m;
}
unless ($opt_p) {die "$usage\nplease indicate prefix"}
unless ($opt_u) {die "$usage\nplease indicate UTR.tbl"}
unless ($opt_t) {die "$usage\nplease indicate TSS.tbl"}
unless ($opt_c) {die "$usage\nplease indicate list of contigs"}
unless ($opt_g) {die "$usage\nplease indicate trimmed_genome.tbl"}
unless ($opt_T) {die "$usage\nplease indicate RepeatMasker output"}
unless ($opt_f) {die "$usage\nplease indicate filtered.comp file"}

my ($addnm) = $opt_c =~/\Q.\E(\d+)$/;

my $reff=1;
if ($opt_r){
	unless ($opt_r eq 'N' or $opt_r eq 'no' or $opt_r eq 'Y' or $opt_r eq 'yes'){die "$usage\nplease indicate Y/yes (default) or no/N for -r"}
	if ($opt_r eq 'no' or $opt_r eq 'N'){$reff=0}
}
die "$usage" unless @ARGV==0;

my %contigs;
open(IN,$opt_c) or die "can't open $opt_c";
while(<IN>){
	chomp;
	my ($cc) = (split /\t/)[0];
	$contigs{$cc}++;
}
close IN;

my %pos;
# 0)get TE positions
my %te;
my $TE=0;
open (IN, $opt_T) or die "can't open TE file"; ## TE file is output .fa.out of repeatmasker
while (my $lin=<IN>){
	#next if ~ /matching|score|Low_complexity|Simple_repeat|Satellite/);
	if ($lin=~ /^\s{1,}\d/){
		my ($Ctg, $St, $End)=(split /\s{1,}/, $lin)[5,6,7];
		next unless (exists $contigs{$Ctg});
		if ( ($End-$St) > 50) {
			$te{$Ctg}{$St}=$End;
			$TE++;
		}
	}
}
close IN;

foreach my $ctgte (keys %te){
	foreach my $stte (keys %{$te{$ctgte}}){
		for (my $i=$stte; $i<=$te{$ctgte}{$stte};$i++){
			$pos{$ctgte}{$i}="TE";
		}
	}
}
undef(%te);
print "TE done\n";

# 1) get tRNA, exons, introns etc except if they're in TEs
my %trna;
my %ex;
my %pse;
my %int;

my ($TRNA,$EX,$PSE,$INT)=(0,0,0,0);

open (IN, $opt_g) or die "can't open trimmed genome table\n";
while (<IN>){
	next if /ContigID/;
	chomp;
	my ($ctgID,$st,$ed,$type,undef)=split; 
	next unless (exists $contigs{$ctgID});
	if ($type =~ /tRNA/) {$trna{$ctgID}{$st}=$ed;$TRNA++}
	if ($type =~ /exon/){
		if ($type =~ /pseudo/){$pse{$ctgID}{$st}=$ed;$PSE++}
		else {$ex{$ctgID}{$st}=$ed;$EX++}
	}
	if ($type =~ /intron/){
		if ($type =~ /pseudo/){$pse{$ctgID}{$st}=$ed}
		else {$int{$ctgID}{$st}=$ed;$INT++}
	}
}
close IN;

foreach my $ctgtrna (keys %trna){
	foreach my $sttrna (keys %{$trna{$ctgtrna}}){
		for (my $j=$sttrna; $j<=$trna{$ctgtrna}{$sttrna};$j++){
			unless (exists $pos{$ctgtrna}{$j}){
				$pos{$ctgtrna}{$j}="tRNA";
			}
		}
	}
}
undef(%trna);

foreach my $ctgpse (keys %pse){
	foreach my $stpse (keys %{$pse{$ctgpse}}){
		for (my $k=$stpse; $k<=$pse{$ctgpse}{$stpse};$k++){
			unless (exists $pos{$ctgpse}{$k}){
				$pos{$ctgpse}{$k}="pseudogene";
			}
		}
	}
}
undef(%pse);

foreach my $ctgex (keys %ex){
	foreach my $stex (keys %{$ex{$ctgex}}){
		for (my $l=$stex; $l<=$ex{$ctgex}{$stex};$l++){
			unless (exists $pos{$ctgex}{$l}){
				$pos{$ctgex}{$l}="exon";
			}
		}
	}
}
undef(%ex);

foreach my $ctgint (keys %int){
	foreach my $stint (keys %{$int{$ctgint}}){
		for (my $m=$stint; $m<=$int{$ctgint}{$stint};$m++){
			unless (exists $pos{$ctgint}{$m}){
				$pos{$ctgint}{$m}="intron";
			}
		}
	}
}
undef(%int);

print "Got tRNA, exons, introns and pseudogenes\n";

my %utr;
# 2) get UTRs except if they're in TEs
my ($CI,$TS)=(0,0);

open (IN, $opt_u) or die "peux pas ouvrir UTR table\n";
while (<IN>){
	next if /ContigID/;
	chomp;
	my ($ID3,$st3,$ed3,$type)=(split /\t/)[0,1,2,3];
	next unless (exists $contigs{$ID3}); 
	$utr{$ID3}{$st3}={ed=>$ed3,ty=>$type};
	if ($type =~ /3/){$TS++}
	if ($type =~ /5/){$CI++}
}
close IN;

foreach my $ctgu (keys %utr){
	foreach my $stu (keys %{$utr{$ctgu}}){
		for (my $n=$stu; $n<=$utr{$ctgu}{$stu}{ed};$n++){
			unless (exists $pos{$ctgu}{$n}){
				$pos{$ctgu}{$n}=$utr{$ctgu}{$stu}{ty};
			}
		}
	}
}
undef(%utr);
my $TSS;
my %tss;
# 3) get TSSs except if they're in TEs, and attribute "TSS/UTR" if in 3'UTR
open (IN, $opt_t) or die "peux pas ouvrir tss table\n";
while (<IN>){
	next if /ContigID/;
	chomp;
	my ($ID,$st2,$ed2)=(split /\t/)[0,1,2];
	next unless (exists $contigs{$ID});
	$tss{$ID}{$st2}=$ed2;
	$TSS++;
}
close IN;

foreach my $ctgt (keys %tss){
	foreach my $stt (keys %{$tss{$ctgt}}){
		for (my $o=$stt; $o<=$tss{$ctgt}{$stt};$o++){
			if (exists $pos{$ctgt}{$o}){
				if ($pos{$ctgt}{$o} =~ /3/){
					$pos{$ctgt}{$o}="TSS-UTR";
				}
			}
			else { 
				$pos{$ctgt}{$o}="TSS";
			}
		}
	}
}

undef(%tss);

print "got info for all positions\n";

# 4) check actual coverage
# count nb of SNPs at the same time

my %det;
my %snps;

my $p=0;
my $q=0;
my $rf=0;
my $r=0;

open(IN,$opt_f) or die "can't open comp file";
while(<IN>){
	chomp;
	my ($ctgvcf,$ps,$rec)=split;
	next unless (exists $contigs{$ctgvcf});
	$det{$ctgvcf}{$ps}{ct}++;	#b/c may have several lines for same position
	my $tpr="";
	my %tpa;
	my %tp;
	my %mf;
	my %sple;
	my @var = split(/-/, $rec);
	foreach my $stu (@var){
		my ($nm,$ref,$at) = (split /;/, $stu)[0,1,2]; 
		$sple{$ctgvcf}{$ps}{$nm}++;
		if ($at =~/^[ACTG]{1}$/ and $ref =~/^.{1}$/) { 
			my ($AD)= (split /;/, $stu)[4];	
			push(@{$tp{$ctgvcf}{$ps}{$at}{AD}},$AD);
			$tp{$ctgvcf}{$ps}{$at}{ct}++;
			$tpr=$ref;
			$tpa{$at}++;
			if ($rf>0){$p++}
		}
		if ($at =~ /\Q.\E/){
			$rf=$rf+2;
			if (exists $tp{$ctgvcf}{$ps}) {$p++}
		}
	}
	# get nb of samples
	foreach my $sam (keys %{$sple{$ctgvcf}{$ps}}){
		$det{$ctgvcf}{$ps}{spl}++;
		if ($sple{$ctgvcf}{$ps}{$sam}>1){$r++}
	}
	$q=$det{$ctgvcf}{$ps}{spl}*2;

	foreach my $alt (keys %{$tp{$ctgvcf}{$ps}} ) {
		foreach my $ads (@{$tp{$ctgvcf}{$ps}{$alt}{AD}}){
			my ($adr,$ada)= split ("," , $ads);
	
			if ($adr >= (($adr+$ada)*$th)) { 
				$det{$ctgvcf}{$ps}{het}++;
				$rf++;
				$mf{$alt}=$tp{$ctgvcf}{$ps}{$alt}{ct};
			}
			else {			
				$tp{$ctgvcf}{$ps}{$alt}{ct}++;
			}
		}
	}
	unless (exists $det{$ctgvcf}{$ps}{het}) {
               $det{$ctgvcf}{$ps}{het}=0;
 
	}
	# correct for reference counting as a homozygote sample if needed
	if ($reff==1){
		$det{$ctgvcf}{$ps}{spl}++;
		$rf=$rf+2;
	}

        # if completely homozygote and {ct}==nb de samples => not a SNPs
	# else it's a snp...
	foreach my $alte (keys %{$tp{$ctgvcf}{$ps}}) {
		if ($det{$ctgvcf}{$ps}{het}>0) {$p++}
		else {
			unless ($tp{$ctgvcf}{$ps}{$alte}{ct}==($det{$ctgvcf}{$ps}{spl}*2)) {
				$p++;
				$mf{$alte}=$tp{$ctgvcf}{$ps}{$alte}{ct};
			}
		}
	}
	
	if ($p>=1){
		$snps{$ctgvcf}{$ps}{rc}=$rec;

		# get hom counts
		if ($r>=1){
			$snps{$ctgvcf}{$ps}{wr}="check_het"
		} 
		else {
			$snps{$ctgvcf}{$ps}{wr}="na"
		}

		$det{$ctgvcf}{$ps}{hom}=$det{$ctgvcf}{$ps}{spl}-$det{$ctgvcf}{$ps}{het};

		# get MAF
		foreach my $maf (keys %mf){
        		if ($mf{$maf}<$q){$q=$mf{$maf}}
       		}
		if ($q>$rf){$q=$rf}
		$snps{$ctgvcf}{$ps}{maf}=$q/($det{$ctgvcf}{$ps}{spl}*2);
		
		# get type of SNPs
		my $truc=0;
		foreach my $bll (keys %tpa) {$truc++}
		if ($truc > 1){
			$snps{$ctgvcf}{$ps}{type}="multi";
		}
		else { 
			foreach my $base (keys %tpa){
				if ($tpr=~/A|G/){
					if ($base=~/A|G/){
						$snps{$ctgvcf}{$ps}{type}="transition";
					}
					else {
						$snps{$ctgvcf}{$ps}{type}="transversion";
					}	
				}
				if ($tpr=~/C|T/){
					if ($base=~/C|T/){
						$snps{$ctgvcf}{$ps}{type}="transition";
					}
					else {
						$snps{$ctgvcf}{$ps}{type}="transversion";
					}
				}
			}
		}
	} 
	$p=0;
	$rf=0;
	$r=0;
	undef(%tp);
	undef(%tpa);
	undef(@var);
	undef(%mf);
	undef(%sple);
}
close IN;
print "got SNPs info\n";

# 5) counts (rest is intergenic ) and prints outputs

my $snp=0;
open (OUT, ">$opt_p.$addnm.snps.txt") or die "can't open SNP output";
my $varr="NOT considered"; 
if ($reff==1) {$varr="considered"}
print OUT "### Warning: Reference is ",$varr," a sample\n";
print OUT "contig\tposition\tnb_samples\tnb_het\tnb_hom\tMAF\ttype\tlocus_type\trecords\n";
foreach my $ctgsnp (keys %snps){
	foreach my $pn (keys %{$snps{$ctgsnp}}){
		$snp++;
		if (exists $pos{$ctgsnp}{$pn}){
			print OUT "$ctgsnp\t$pn\t$det{$ctgsnp}{$pn}{spl}\t$det{$ctgsnp}{$pn}{het}\t$det{$ctgsnp}{$pn}{hom}\t$snps{$ctgsnp}{$pn}{maf}\t$snps{$ctgsnp}{$pn}{type}\t$pos{$ctgsnp}{$pn}\t$snps{$ctgsnp}{$pn}{rc}\n";
		}
		else {
		print OUT "$ctgsnp\t$pn\t$det{$ctgsnp}{$pn}{spl}\t$det{$ctgsnp}{$pn}{het}\t$det{$ctgsnp}{$pn}{hom}\t$snps{$ctgsnp}{$pn}{maf}\t$snps{$ctgsnp}{$pn}{type}\tintergenic\t$snps{$ctgsnp}{$pn}{rc}\n";
		}
	}
}
close OUT;

open (OUT, ">$opt_p.$addnm.det_pos.txt") or die "can't open det_pos output";
print OUT "contig\tposition\tlocus_type\n";
foreach my $cng (keys %det){
	foreach my $pn (keys %{$det{$cng}}){
		if (exists $pos{$cng}{$pn}){
			print OUT "$cng\t$pn\t$pos{$cng}{$pn}\n";
		}
		else {
			print OUT "$cng\t$pn\tintergenic\n";
		}
	}
}
close OUT;

my ($tr,$exo,$intr,$pseudo,$tp,$fp,$tsss,$dbl,$inter,$det,$tee)=(0,0,0,0,0,0,0,0,0,0,0);
foreach my $contig (keys %det){
	foreach my $loc (keys %{$det{$contig}}){
		$det++;
		if (exists $pos{$contig}{$loc}) { 
			if ($pos{$contig}{$loc}=~/tRNA/){$tr++}
			if ($pos{$contig}{$loc}=~/TE/){$tee++}
			if ($pos{$contig}{$loc}=~/exon/){$exo++}
			if ($pos{$contig}{$loc}=~/intron/){$intr++}
			if ($pos{$contig}{$loc}=~/pseud/){$pseudo++}
			if ($pos{$contig}{$loc}=~/3/){$tp++}
			if ($pos{$contig}{$loc}=~/5/){$fp++;$tsss++}
			if ($pos{$contig}{$loc} eq "TSS"){$tsss++}
			if ($pos{$contig}{$loc}=~/TSS-UTR/){$dbl++;$tsss++;$tp++}
		}
		else {$inter++}		
	}	
}

my ($try,$exoy,$intry,$pseudoy,$tey,$tpy,$fpy,$tsssy,$dbly,$intery)=(0,0,0,0,0,0,0,0,0,0);
foreach my $contigs (keys %pos){
	foreach my $locs (keys %{$pos{$contigs}}){
		if (exists $pos{$contigs}{$locs}) { 
			if ($pos{$contigs}{$locs}=~/tRNA/){$try++}
			if ($pos{$contigs}{$locs}=~/TE/){$tey++}
			if ($pos{$contigs}{$locs}=~/exon/){$exoy++}
			if ($pos{$contigs}{$locs}=~/intron/){$intry++}
			if ($pos{$contigs}{$locs}=~/pseud/){$pseudoy++}
			if ($pos{$contigs}{$locs}=~/3/){$tpy++}
			if ($pos{$contigs}{$locs}=~/5/){$fpy++;$tsssy++}
			if ($pos{$contigs}{$locs} eq "TSS"){$tsssy++}
			if ($pos{$contigs}{$locs}=~/TSS-UTR/){$dbly++;$tsssy++;$tpy++}
		}
		else {$intery++}		
	}	
}


open (OUT, ">$opt_p.$addnm.coverage_report.txt") or die "can't open report output";
print OUT "
total nb of bases detected\t$det
distance in TEs\t$tey\tdetected\t$tee\tnb of features\t$TE
distance in tRNAs\t$try\tdetected\t$tr\tnb of features\t$TRNA
distance in exons\t$exoy\tdetected\t$exo\tnb of features\t$EX
distance in introns\t$intry\tdetected\t$intr\tnb of features\t$INT
distance in pseudogenes\t$pseudoy\tdetected\t$pseudo\tnb of features\t$PSE
distance in 3'UTR\t$tpy\tdetected\t$tp\tnb of features\t$TS
distance in 5'UTR\t$fpy\tdetected\t$fp\tnb of features\t$CI
distance in TSS\t$tsssy\tdetected\t$tsss\tnb of features\t$TSS
distance in 3'UTR-TSS overlap\t$dbly\tdetected\t$dbl
distance in intergenic\t$intery(should be 0)\tdetected\t$inter
"; 
close OUT;


