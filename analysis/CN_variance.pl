#!/usr/bin/perl
use warnings;
use strict;

my $now_dir = "/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion";
my $pwd = `pwd`;chomp $pwd;
if($now_dir ne $pwd){die "ERROR::do on wrong dir\n";}

open(ASCN,"gunzip -c all_patient_ascat.tsv.gz|");
my $header=<ASCN>;chomp$header;
my %col=&header2hash($header);
my %cn_variance=();
my %cn_length=();
while(<ASCN>){
		chomp;
		my @line = split(/\t/,);
		if($line[$col{"chr"}] eq "X"){next;}
		my$cnv = abs($line[$col{nMajor}]-1)+abs($line[$col{nMinor}]-1);
		my$length =$line[$col{endpos}]-$line[$col{startpos}]+1;
		$cn_variance{$line[$col{sample}]}+=$cnv*$length;
		if(($line[$col{nMajor}]==1)&&($line[$col{nMinor}]==0)){$cn_length{$line[$col{sample}]}{"1_0"}+=$length;}
		elsif(($line[$col{nMajor}]==2)&&($line[$col{nMinor}]==0)){$cn_length{$line[$col{sample}]}{"2_0"}+=$length;}
		elsif(($line[$col{nMajor}]==1)&&($line[$col{nMinor}]==1)){$cn_length{$line[$col{sample}]}{"1_1"}+=$length;}
		else{$cn_length{$line[$col{sample}]}{other}+=$length;}
}
close ASCN;

open(OUT,">by_patient_CN_variance.tsv");
print OUT "patient_id\tcn_variance\t1_0\t2_0\t1_1\tother\n";
foreach my $sid(sort keys %cn_variance){
		my @cn=("1_0","2_0","1_1","other");
		foreach my $cn (@cn){if(!defined($cn_length{$sid}{$cn})){$cn_length{$sid}{$cn}=0;}}
		print OUT "$sid\t$cn_variance{$sid}\t$cn_length{$sid}{'1_0'}\t$cn_length{$sid}{'2_0'}\t$cn_length{$sid}{'1_1'}\t$cn_length{$sid}{other}\n";
}
close OUT;







#####################################################################################################3
sub header2hash ( $ ){
		my $header = $_[0];
		my @colm = split(/\t/,$header);
		my %out = ();
		for(my $i=0; $i < scalar(@colm); $i++){
				$out{$colm[$i]}=$i;
		}
		return(%out);
}
