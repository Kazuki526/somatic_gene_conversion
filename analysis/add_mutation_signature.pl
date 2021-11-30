#!/usr/bin/perl
use warnings;
use strict;

my $now_dir = "/Users/tkaz/work/SomaticGeneConversion/extract_raw_maf/somatic_conversion";
my $pwd = `pwd`;chomp $pwd;
if($now_dir ne $pwd){die "ERROR::do on wrong dir\n";}

my $all_tbl = "all_pass_with_dist_position_CPE.maf.gz";
-e $all_tbl or die "ERROR::$all_tbl is not exist\n";
my $fasta="/Users/tkaz/work/GRCh38.d1.vd1.fa";
-e $fasta or die "ERROR::$fasta is not exist\n";

my $patient_file = "sample_list.tsv";
-e $patient_file or die "ERROR::$patient_file is not exist\n";
my %patients = ();
open(CF,"$patient_file");
my $header = <CF>;chomp $header;
my %col = &header2hash($header);
while(<CF>){
		chomp;
		my @line = split(/\t/,);
		$patients{$line[$col{tumor_sample_id}]} = $line[$col{purity}];
}
close CF;

open(TBL,"gunzip -c $all_tbl|");
$header = <TBL>;chomp $header;
%col = &header2hash($header);
open(OUT, "|gzip -c >/Users/tkaz/Dropbox/work/somatic_gene_conversion/revise/all_pass_SBS.tsv.gz");
print OUT  "sample_id\tchr\tstart\tref\talt\tvariant_type\tt_depth\tt_ref\tt_alt\tgenotype\tpurity\tdcv_posi\tsbs_ref\tsbs_alt\tsbs_pattern\n";
while(<TBL>){
		chomp;
		my @line=split(/\t/,);
		if($line[$col{variant_type}] ne "SNP"){next;}
		my $sid = $line[$col{sample_id}];
# purity
		if(!defined $patients{$sid}){next;}
		if($patients{$sid} <0.8){next;}
# 0< allele num <=2
		my $genotype="";
		if(($line[$col{ascat_major}]==1) && ($line[$col{ascat_minor}]==1)){$genotype="AB";
		}elsif(($line[$col{ascat_major}]==1) && ($line[$col{ascat_minor}]==0)){$genotype="A";
		}elsif(($line[$col{ascat_major}]==2) && ($line[$col{ascat_minor}]==0)){$genotype="AA";
		}else{next;}
# depth correlation value
		my$dcv_posi = $line[$col{mutect_dcv_posi}]/$line[$col{mutect_mut_num}];
		if(($dcv_posi<0.1)||($dcv_posi>0.9)){next;}
		my $out  = "$sid\t$line[$col{chr}]\t$line[$col{start}]\t$line[$col{ref}]\t$line[$col{alt}]\t$line[$col{variant_type}]\t";
		$out    .= "$line[$col{t_depth}]\t$line[$col{t_ref}]\t$line[$col{t_alt}]\t$genotype\t$patients{$sid}\t$dcv_posi";
		my ($reverse,$ref,$alt)=(0,$line[$col{ref}],$line[$col{alt}]);
		if($line[$col{ref}] !~/[CT]/){
				$reverse=1;
				$ref =~ tr/gGaA/cCtT/;
				$alt =~ tr/gGaAcCtT/cCtTgGaA/;
		}
		my ($st,$en)=($line[$col{start}]-1,$line[$col{start}]+1);
		my $region ="$line[$col{chr}]:$st-$en";
		my $trinuc=`samtools faidx $fasta $region|awk 'NR>1'`;
		chomp $trinuc;
		if($reverse){
				$trinuc =~ tr/gGaAcCtT/cCtTgGaA/;
				$trinuc = reverse($trinuc);
		}
		print OUT "$out\t$ref\t$alt\t$trinuc\n";
}
close TBL;
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
