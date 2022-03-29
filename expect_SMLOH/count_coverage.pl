#!/usr/bin/perl
use warnings;
use strict;

my $pid=$ARGV[0];
my $ref="~/tcga_data/GRCh38.d1.vd1.fa";
if(defined $ARGV[1]){$ref=$ARGV[1];}

my @triplet = qw(aca acc acg act cca ccc ccg cct gca gcc gcg gct tca tcc tcg tct ata atc atg att cta ctc ctg ctt gta gtc gtg gtt tta ttc ttg ttt);
my %coverednum=();
foreach my $trip (@triplet){
	$trip=uc $trip;
	$coverednum{$trip}=0;
}

my %trip_rev=();
my @atgc=qw(A T G C);
foreach my $a (@atgc){
	foreach my $b (@atgc){
		foreach my $c (@atgc){
			my $trip="$a$b$c";
			my $out=$trip;
			if(!defined $coverednum{$trip}){
				$out=reverse($trip);
				$out=~tr/ATGC/TACG/;
				if(!defined $coverednum{$out}){die "why???\n"}
			}
			$trip_rev{$trip}=$out;
		}
	}
}



my $covfile="coveraged.txt";
if(!-s $covfile){die "ERROR::coveraged.txt have no data at $pid\n";}
open(COV,"$covfile")or die "ERROR::cannot open $covfile\n";
my($chr,$start,$end)=("",0,0);
while(<COV>){
	chomp;
	my ($chrom,$pos)=split(/\t/,);
	if($chrom !~ /^chr[\dXY]+$/){next;}
	if($end +1 != $pos){
		if($chr ne ""){&count_coverednum($chr,$start,$end);}
		$chr=$chrom;$start=$pos;$end=$pos;
	}else{
		$end=$pos;
	}
}
if($chr ne ""){&count_coverednum($chr,$start,$end);}
close COV;

open(OUT,">count_triplet_num.txt");
foreach my $trip (@triplet){
	print OUT "$pid\t$trip\t$coverednum{$trip}\n";
}

sub count_coverednum( $ $ $ ){
	my($chr,$star,$end)=@_;
	$start--;$end++;
	my $fasta=`/usr/local/package/samtools/0.1.20/bin/samtools faidx $ref $chr:$start-$end`;
	my $seq = &get_seq($fasta);
	for(my $i=1;$i<$end-$start;$i++){
		my $trip=substr($seq,$i-1,3);
		if(!defined$trip_rev{$trip}){
			print "WARNING::what sequece?? $chr:$start-$end pos:$i => $trip\n";
			next;
		}
		$coverednum{$trip_rev{$trip}}++;
	}
}

sub get_seq( $ ){
	my $seq="";
	my @lines=split(/\n/,$_[0]);
	foreach my $line (@lines){
		if($line =~/^>/){next;}
		$seq.=$line;
	}
	$seq=uc $seq;
	return($seq);
}


