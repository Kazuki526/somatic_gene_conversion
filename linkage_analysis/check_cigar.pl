#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my $ref = "$ENV{HOME}/GRCh38.d1.vd1.fa";
-e $ref or die "ERROR::not exist $ref\n";

my $token_path=`ls $ENV{HOME}/|grep 'gdc-user-token'`;
if(!$token_path){die "ERROR::token file not exitst!!";}
chomp $token_path;
$token_path="$ENV{HOME}/$token_path";
my $token=`cat $token_path`;

my %gene_conversion=();
#open(CAND, "candidate_conversion_variants.tsv") or die  "ERROR::cannot open candidate_conversion_variants.tsv\n";
open(CAND, "candidate_conversion_variants_test.tsv") or die  "ERROR::cannot open candidate_conversion_variants_test.tsv\n";
my $candlist_header = <CAND>; chomp $candlist_header;
my %col_can=&header2hash($candlist_header);
while(<CAND>){
		chomp;
		my @line = split(/\t/,);
		my $sid=$line[$col_can{sample_id}];
		my ($start,$end) = ($line[$col_can{start}]-600,$line[$col_can{start}]+600);
		$gene_conversion{$sid}{site}.="$line[$col_can{chr}]:$line[$col_can{start}]:$line[$col_can{ref}]:$line[$col_can{alt}],";
		$gene_conversion{$sid}{purity}=$line[$col_can{purity}];
		$gene_conversion{$sid}{line}.="$_\n";
}
close CAND;

my %norm2tumor=();
open(TNID,"bam_check/sample_list/tumor_normal_sample_id.tsv") or die "ERROR::cannot tumor_normal_sample_id.tsv";
my $header=<TNID>;chomp $header;
my %col = &header2hash($header);
while(<TNID>){
		chomp;
		my @line = split(/\t/,);
		my $tid = $line[$col{tumor_sample_id}];
		my $nid = $line[$col{normal_sample_id}];
		if(defined $gene_conversion{$tid}){
				$gene_conversion{$tid}{nid}=$nid;
				$norm2tumor{$nid}=$tid
		}
}
close TNID;

open(RES,"nkf -Lu bam_check/sample_list/response.tsv|") or die "ERROR::cannot open response.tsv\n";
$header=<RES>;chomp$header;
%col = &header2hash($header);
while(<RES>){
		chomp;
		my @line = split(/\t/,);
		my $sid = $line[$col{'cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id'}];
		if((!defined $gene_conversion{$sid})&&(!defined $norm2tumor{$sid})){next;}
		my $stype = $line[$col{'cases.0.samples.0.sample_type'}];
		if($stype eq "Primary Tumor"){
				$gene_conversion{$sid}{uuid} = $line[$col{id}];
				$gene_conversion{$sid}{filename} = $line[$col{file_name}];
		}else{
				my $tid = $norm2tumor{$sid};
				$gene_conversion{$tid}{nuuid} = $line[$col{id}];
				$gene_conversion{$tid}{nfilename} = $line[$col{file_name}];
		}
}
close RES;


open(LOG,">log.txt");
my %cigar_check=();
foreach my $sid (sort keys %gene_conversion){
		print LOG "$sid\n";$|=1;
		my $dir="bam_check/$sid";
		my ($tfile,$nfile)=("$dir/$gene_conversion{$sid}{filename}","$dir/$gene_conversion{$sid}{nfilename}");
		my %tcigar = &check_cigar($tfile);
		foreach my $cig (keys%tcigar){
				$cigar_check{$cig}+=$tcigar{$cig};
		}
}
foreach my $cig (keys%cigar_check){
		print "$cig = $cigar_check{$cig}\n";
}
		

###################################################################################################
sub check_cigar( $ ){
		my $bam = $_[0];
		open(BAM,"samtools view $bam|");
		my %cigar_num=();
		while(<BAM>){
				chomp;
				my @bam = split(/\t/,);
				my $cigar = $bam[5];
				if($cigar eq "*"){next;}
				while($cigar){
						$cigar =~ s/^(\d+)(\w)//;
						my($a,$b)=($1,$2);
						$cigar_num{$b}+=$a;
				}
				if($bam[11] =~ /^MC:Z:(.+)$/){
						$cigar=$1;
						if($cigar eq "*"){next;}
						while($cigar){
								$cigar =~ s/^(\d+)(\w)//;
								my($a,$b)=($1,$2);
								$cigar_num{$b}+=$a;
						}
				}
		}
		return(%cigar_num);
}


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
#######################################################################################################
sub read_vcv ( $ ){
		my $vcf_name =$_[0];
		my %vcf=();
		open(SNP,"$vcf_name.snp.vcf");
		while(<SNP>){
				if($_=~ /^#/){next;}
				chomp;
				my @line = split(/\t/,);
				if($line[8] ne "GT:GQ:DP:RD:AD:FREQ:DP4"){die "ERROR::$vcf_name.snp.vcf have not FORMAT=GT:GQ:DP:RD:AD:FREQ:DP4\n$_\n";}
				if($line[9] !~ /^0\/1:/){next;}
				my $site="$line[0]:$line[1]:$line[3]:$line[4]";
				my @norm = split(/:/,$line[9]);
				my @tumor= split(/:/,$line[10]);
				$vcf{$site}="$norm[2]:$norm[3]$norm[4]\t$tumor[2]:$tumor[3]:$tumor[4]";
		}
		close SNP;
		open(IND,"$vcf_name.indel.vcf");
		while(<IND>){
				if($_=~ /^#/){next;}
				chomp;
				my @line = split(/\t/,);
				if($line[8] ne "GT:GQ:DP:RD:AD:FREQ:DP4"){die "ERROR::$vcf_name.snp.vcf have not FORMAT=GT:GQ:DP:RD:AD:FREQ:DP4\n$_\n";}
				if($line[9] !~ /^0\/1:/){next;}
				my $site="$line[0]:$line[1]:$line[3]:$line[4]";
				my @norm = split(/:/,$line[9]);
				my @tumor= split(/:/,$line[10]);
				$vcf{$site}="$norm[2]:$norm[3]$norm[4]\t$tumor[2]:$tumor[3]:$tumor[4]";
		}
		close IND;
		return(%vcf);
}


