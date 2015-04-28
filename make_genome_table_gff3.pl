#!usr/bin/perl

################
# creates single table <contigID><start><end><type of locus> 
# types being tRNA / 3'UTR / 5'UTR / exon / intron / 1kb upstream TSS
# 2012 - cris
################
use strict;use warnings;

use Getopt::Std;
getopts('hvtp:g:i:m:');
our($opt_h,$opt_v,$opt_p,$opt_g, $opt_i, $opt_m);

my $vers="version 2.1 -- for gff3";
my $usage ="
usage: [options] -p <prefix> -g <gff3> -m <miRNA.gff3> -i <piRNA.gff3>
options:
   -h help
   -v version
	
	INPUT
   -p prefix (req)
   -g genome gff3 (req)
   -m miRNA gff3 (if available)
   -i piRNA gff3 (if available)
";

my $help="
outputs tables easier to work with than gff2...

Version note :
v2.0 added options and merged with other scripts (trim_gen_tbl.pl, getTSS.pl and utrs_coord.pl) to get all .tbl files.
v2.1 merged in previous steps of the pipeline.
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
unless ($opt_p) {die "$usage\nplease indicate prefix for output files"}
unless ($opt_g) {die "$usage\nplease indicate gff3 file"}

die $usage unless @ARGV==0;

my %gns;
my %starts;
my $i=0;
my %nb;
my @cds;
my %ctgs;
my %ctg_gns;
my %ckgn;
open (IN, $opt_g);
open (OUTG, ">exons_$opt_p.gff3") or die "cant open gns output";
open (OUTU, ">UTR_$opt_p.tbl") or die "cant open UTR tbl";
open (OUTCC, ">$opt_p.ctgs_length.txt") or die "cant open contig length";
while (<IN>){
	next unless (/^scaff|^linkage_/);
	my ($id, $met, $st, $ed)=(split "\t")[0,2,3,4];
	if ($met eq 'CDS'){print OUTG}
	elsif ($met =~ /_prime_UTR/) {
		my ($grpe)=(split /\t/)[8];
		print OUTU "$id\t$st\t$ed\t";
		if ($met=~/five/) {
			my ($nme)= $grpe =~ /ID\=(.+)\:five_pr/;
			print OUTU "5'UTR\t$nme\n";
		}
		elsif ($met =~/three/){
			my ($nme3)= $grpe =~ /ID\=(.+)\:three_pr/; 
			print OUTU "3'UTR\t$nme3\n";
		}
	}
	elsif ($met eq 'contig'){
		print OUTCC "$id\t$st\t$ed\n";
		$ctgs{$id}=$ed;		#stores contig sizes
	}
	elsif ($met eq 'gene'){
		my ($strand, $group)=(split /\t/)[6,8];
		my ($iso) = $group =~ /ID\=(.+)\;Name\=/; #die "$group\n\n$iso\n"; 
		$gns{$iso} = {ctgid=>$id, st=>$st, end=>$ed, strd=>$strand};
		$starts{$iso}=$st; 		#start position of genes
		$i++;
		push (@{$ctg_gns{$id}}, $iso);
		my $taille = $ed-$st;
		push(@{$ckgn{$taille}},$iso);
	}	
}
close OUTG;
close OUTU;
close OUTCC;
close IN;

## exon / intron
my %tmp;
my %chr;
my %strd;
open (IN, "exons_$opt_p.gff3") or die "can't open gns\n";
while (<IN>) {
	chomp;
	my($ID,$st,$end,$strd,$grp)=(split /\t/)[0,3,4,6,8]; 
	my ($nm)= $grp =~ /ID\=(.+)\:cds/; 
	#exons
	$strd{$nm}=$strd;
	$chr{$nm}=$ID;
	$tmp{$nm}{$st}=$end;
}
close IN;

my %tmpgs;
my %tmpge;
foreach my $geneID (keys %tmp){
	my $l=0;
	foreach my $start (sort {$a <=> $b} keys %{$tmp{$geneID}}) {
		$tmpgs{$geneID}[$l]=$start;
		$tmpge{$geneID}[$l]=$tmp{$geneID}{$start};
		$l++;
	}
}

open (OUT, ">trimmed_$opt_p.genome.tbl") or die "can't open output\n";
#print OUT "ContigID\tstart\tend\ttype\n";
foreach my $gid (keys %tmpgs){
	my $trm = scalar(@{$tmpgs{$gid}}) - 2;
	for my $m (0..$trm){
		print OUT "$chr{$gid}\t$tmpgs{$gid}[$m]\t$tmpge{$gid}[$m]\texon\t$gid\n";
		print OUT "$chr{$gid}\t", $tmpge{$gid}[$m] + 1 ,"\t", $tmpgs{$gid}[$m+1] -1 ,"\tintron\t$gid\n";
	}
	print OUT "$chr{$gid}\t$tmpgs{$gid}[$trm+1]\t$tmpge{$gid}[$trm+1]\texon\t$gid\n";
}
print "genes and utrs done.\n";

if ($opt_i){
	print " looking at pi_rnas\n";
	open (IN, $opt_i) or die "can't open $opt_i";
	while (<IN>){
		next unless (/^link|^scaff/);
		my($IDp,$stp,$endp,$grpi)=(split /\t/)[0,3,4,8]; 
		my ($nmp)= $grpi =~ /Name\=(.+)\;Target/; 
		print OUT "$IDp\t$stp\t$endp\tpiRNA\t$nmp\n"
	}
	close IN;

}
if ($opt_m){
	print " looking at mi_rnas\n";
	open (IN, $opt_m) or die "can't open $opt_m";
	while (<IN>){
		next unless (/^link|^scaff/);
		my($IDm,$stm,$endm,$grm)=(split /\t/)[0,3,4,8]; 
		my ($nmm)= $grm =~ /Name\=(.+)\;Target/; 
		print OUT "$IDm\t$stm\t$endm\tmiRNA\t$nmm\n"
	}
	close IN;
}

close OUT;
print "looking into tssssssss\n";
## TSSs
my %tss;
my ($j, $k)=(0,0);

my %tre;
my %dbb;
## find genes within genes.
foreach my $sz (sort {$b<=>$a} keys %ckgn){
	foreach my $nom (@{$ckgn{$sz}}){
		for my $r ($gns{$nom}{st}..$gns{$nom}{end}){
			my $pjk = join(";", $gns{$nom}{ctgid}, $r);
			if (exists $tre{$pjk}){$dbb{$nom}++}
			else {$tre{$pjk}++}
		}
	}
}

my %temp;
my %sorted_gn;
## gives genes index # per contig 
foreach my $bip (keys %ctgs) {
	foreach my $gnID (@{$ctg_gns{$bip}}) {
		unless (exists $dbb{$gnID}){
			$temp{$gnID} = $starts{$gnID};
			push(@cds,$gnID);
			$nb{$bip}++;
		}
	}
	foreach my $pos (sort {$temp{$a} <=> $temp{$b}} keys %temp){
		$k++;
		$sorted_gn{$bip}{$k}=$pos;	##keeps record of genes on contig in order
		$gns{$pos}{idx}=$k;
	}
	%temp=(); $k=0;
}

## gets 1kb or less from TSS coordinate
my $dis=0;
my $posplus=0;
my $posneg=0;

foreach my $gene (@cds){
	## if only one gene on contig
	if ($nb{$gns{$gene}{ctgid}}==1){
		if ($gns{$gene}{strd} =~/\Q+\E/){
			## if st > 1000 => 1000; note = contig edge;
			if ($gns{$gene}{st} > 1000){
				$tss{$gene}={stTSS=> $gns{$gene}{st}-1000, endTSS=> $gns{$gene}{st}, dis=> 1000, note=>"unique"};
			}	
			## else => st; note = contig edge; $j++; die "";
			else {
				$tss{$gene}={stTSS=> 1, endTSS=> $gns{$gene}{st}, dis=> $gns{$gene}{st}, note=>"unique"};
			}
		}
		## else 
		else {
			## if end contig - st > 1000 => 1000; 
			if (($ctgs{$gns{$gene}{ctgid}}-$gns{$gene}{end})>1000){
				$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $gns{$gene}{end}+1000, dis=> 1000, note=>"unique"};
			}
			## else => end contig - st;
			else {
				$dis = $ctgs{$gns{$gene}{ctgid}}-$gns{$gene}{end};
				$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $ctgs{$gns{$gene}{ctgid}}, dis=> $dis, note=>"unique"};
			}
		}
	}

	## if first gene of a contig
	elsif ($gns{$gene}{idx}==1){
		## if strand +
		if ($gns{$gene}{strd} =~/\Q+\E/){
			## if st > 1000 => 1000; note = contig edge;
			if ($gns{$gene}{st} > 1000){
				$tss{$gene}={stTSS=> $gns{$gene}{st}-1000, endTSS=> $gns{$gene}{st}, dis=> 1000, note=>"NA"};
			}	
			## else => st; note = contig edge; $j++; die "";
			else {
				$tss{$gene}={stTSS=> 1, endTSS=> $gns{$gene}{st}, dis=> $gns{$gene}{st}, note=>"ctg_edge"};
			}
		}
		## else 
		else {
			## if strand gene #2 +
			if ($gns{$sorted_gn{$gns{$gene}{ctgid}}{2}}{strd} =~ /\Q+\E/){
				## if distance entre les starts des genes > 2000 => 1000;  
				if (($gns{$sorted_gn{$gns{$gene}{ctgid}}{2}}{st}-$gns{$gene}{end}) > 2000){
					$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $gns{$gene}{end}+1000, dis=> 1000, note=>"NA"};
				}
				## else => dist / 2; note = split;
				else {
					$dis = ($gns{$sorted_gn{$gns{$gene}{ctgid}}{2}}{st}-$gns{$gene}{end}) / 2;
					$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $gns{$gene}{end}+$dis, dis=> $dis, note=>"split"};
				}
			}
			## else 
			else {
				## if distance entre start et fin du gn #2 > 1000 => 1000; 
				if (($gns{$sorted_gn{$gns{$gene}{ctgid}}{2}}{st}-$gns{$gene}{end}) > 1000){
					$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $gns{$gene}{end}+1000, dis=> 1000, note=>"NA"};
				}
				## else => dist = distance entre le start et la fin du gn; 
				else {
					$dis = $gns{$sorted_gn{$gns{$gene}{ctgid}}{2}}{st}-$gns{$gene}{end};
					$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $gns{$gene}{end}+$dis, dis=> $dis, note=>"short"};
				}
			}
		}
	}
	## if last gene of a contig
	elsif ( $gns{$gene}{idx} == $nb{$gns{$gene}{ctgid}}){
		## if strand -
		if ($gns{$gene}{strd} =~/-/){ #print OUT "$gene\n";
			## if end contig - st > 1000 => 1000; 
			if (($ctgs{$gns{$gene}{ctgid}}-$gns{$gene}{end})>1000){
				$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $gns{$gene}{end}+1000, dis=> 1000, note=>"NA"};
			}
			## else => end contig - st;
			else {
				$dis = $ctgs{$gns{$gene}{ctgid}}-$gns{$gene}{end};
				$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $ctgs{$gns{$gene}{ctgid}}, dis=> $dis, note=>"ctg_edge"};
			}
		}
		## else -- 
		else {
			$posneg = $gns{$gene}{idx} - 1; 
			## if strand gene n-1 -    #$sorted_gn{contigID}{position}=CRENNNN  #$gns{CRENNNN}{idx}=position
			if ($gns{$sorted_gn{$gns{$gene}{ctgid}}{$posneg}} {strd} =~ /-/){  
				## if distance entre les starts des genes > 2000 => 1000;
				if (($gns{$gene}{st} - $gns{$sorted_gn{$gns{$gene}{ctgid}}{$posneg}}{end}) > 2000){
					$tss{$gene}={stTSS=> $gns{$gene}{st}-1000, endTSS=> $gns{$gene}{st}, dis=> 1000, note=>"NA"};
				}
				## else => dist / 2; note = split;
				else {
					## else => dist / 2; note = split;
					$dis = ($gns{$gene}{st} - $gns{$sorted_gn{$gns{$gene}{ctgid}}{$posneg}}{end}) / 2;
					$tss{$gene}={stTSS=> $gns{$gene}{st}-$dis, endTSS=> $gns{$gene}{st}, dis=> $dis, note=>"split"};
				}
			}
			## else 
			else {
				## if distance entre start et fin du gn n-1 > 1000 => 1000;
				if (($gns{$gene}{st} - $gns{$sorted_gn{$gns{$gene}{ctgid}}{$posneg}}{end}) > 1000){  
					$tss{$gene}={stTSS=> $gns{$gene}{st}-1000, endTSS=> $gns{$gene}{st}, dis=> 1000, note=>"NA"};
				}
				## else => dist = distance entre le start et la fin du gn; note = short; $j++; die "";
				else {
					$dis = $gns{$gene}{st} - $gns{$sorted_gn{$gns{$gene}{ctgid}}{$posneg}}{end};
					$tss{$gene}={stTSS=> $gns{$gene}{st}-$dis, endTSS=> $gns{$gene}{st}, dis=> $dis, note=>"short"};
				}
			}
		}
	}
	## else any other position
	else {
		$posplus = $gns{$gene}{idx} + 1;
		$posneg = $gns{$gene}{idx} - 1;

		## if strand +
		if ($gns{$gene}{strd} =~/\Q+\E/){
			## if strand gene n-1 est -
			if ($gns{$sorted_gn{$gns{$gene}{ctgid}}{$posneg}}{strd} =~ /-/){
				## if distance entre les starts des genes > 2000 => 1000; 
				if (($gns{$gene}{st} - $gns{$sorted_gn{$gns{$gene}{ctgid}}{$posneg}}{end}) > 2000){
					$tss{$gene}={stTSS=> $gns{$gene}{st}-1000, endTSS=> $gns{$gene}{st}, dis=> 1000, note=>"NA"};
				}
				## else needs to check gn n-1 
				else {
					## if => dist / 2; note = split;
					$dis = ($gns{$gene}{st} - $gns{$sorted_gn{$gns{$gene}{ctgid}}{$posneg}}{end}) / 2;
					$tss{$gene}={stTSS=> $gns{$gene}{st}-$dis, endTSS=> $gns{$gene}{st}, dis=> $dis, note=>"split"};
				}
			}		
			## else 
			else {
				## if distance entre start et fin du gn n-1 > 1000 => 1000; 
				if (($gns{$gene}{st} - $gns{$sorted_gn{$gns{$gene}{ctgid}}{$posneg}}{end}) > 1000){
					$tss{$gene}={stTSS=> $gns{$gene}{st}-1000, endTSS=> $gns{$gene}{st}, dis=> 1000, note=>"NA"};
				}
				## else => dist = distance entre le start et la fin du gn; 
				else {
					$dis = $gns{$gene}{st} - $gns{$sorted_gn{$gns{$gene}{ctgid}}{$posneg}}{end};
					$tss{$gene}={stTSS=> $gns{$gene}{st}-$dis, endTSS=> $gns{$gene}{st}, dis=> $dis, note=>"short"};
				}
			}
		}
		## else -- 
		elsif ($gns{$gene}{strd} =~/-/) { 
			## if strand gene n+1 + 
			if ($gns{$sorted_gn{$gns{$gene}{ctgid}}{$posplus}}{strd} =~ /\Q+\E/){ 
				## if distance entre les starts des genes > 2000 => 1000; 
				if (($gns{$sorted_gn{$gns{$gene}{ctgid}}{$posplus}}{st} - $gns{$gene}{end}) > 2000){
					$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $gns{$gene}{end}+1000, dis=> 1000, note=>"NA"};
				}
				## else check gn + 1
				else {
					$dis = ($gns{$sorted_gn{$gns{$gene}{ctgid}}{$posplus}}{st} - $gns{$gene}{end}) / 2;
					$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $gns{$gene}{end}+$dis, dis=> $dis, note=>"split"};
				}
			}
			## else 
			else {
				## if distance entre start et fin du gn n+1 > 1000 => 1000; note = na; $j++; die "";
				if (($gns{$sorted_gn{$gns{$gene}{ctgid}}{$posplus}}{st} - $gns{$gene}{end}) > 1000){ 
					$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $gns{$gene}{end}+1000, dis=> 1000, note=>"NA"};
				}
				## else => dist = distance entre le start et la fin du gn; note = short; $j++; die "";
				else {
					$dis = $gns{$sorted_gn{$gns{$gene}{ctgid}}{$posplus}}{st} - $gns{$gene}{end};
					$tss{$gene}={stTSS=> $gns{$gene}{end}, endTSS=> $gns{$gene}{end}+$dis, dis=> $dis, note=>"short"};
				}
			}
		}
	}
}
open (OUT, ">TSS_$opt_p.tbl");

foreach my $TSS (keys %tss){
	print OUT "$gns{$TSS}{ctgid}\t$tss{$TSS}{stTSS}\t$tss{$TSS}{endTSS}\tTSS\t$TSS\t$tss{$TSS}{dis}\t$tss{$TSS}{note}\n"
}
close OUT;

__END__
