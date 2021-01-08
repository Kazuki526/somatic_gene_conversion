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
open(CAND, "candidate_conversion_variants.tsv") or die  "ERROR::cannot open candidate_conversion_variants.tsv\n";
#open(CAND, "candidate_conversion_variants_test.tsv") or die  "ERROR::cannot open candidate_conversion_variants_test.tsv\n";
my $candlist_header = <CAND>; chomp $candlist_header;
my %col=&header2hash($candlist_header);
while(<CAND>){
		chomp;
		my @line = split(/\t/,);
		my $sid=$line[$col{sample_id}];
		my ($start,$end) = ($line[$col{start}]-500,$line[$col{start}]+500);
		$gene_conversion{$sid}{posi}.="$line[$col{chr}]:$start-$end,";
		$gene_conversion{$sid}{purity}=$line[$col{purity}];
		$gene_conversion{$sid}{line}.="$_\n";
}
close CAND;

my %norm2tumor=();
open(TNID,"bam_check/sample_list/tumor_normal_sample_id.tsv") or die "ERROR::cannot tumor_normal_sample_id.tsv";
my $header=<TNID>;chomp $header;
%col = &header2hash($header);
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

open(LOG,">download_error.log");
my $max_processes = 10;
my $pm = new Parallel::ForkManager($max_processes);
foreach my $sid (sort keys %gene_conversion){
		$pm->start and next;
		print "start $sid\n";
		my $dir="bam_check/$sid";
		mkdir $dir;
		open(OUT,">$dir/candidate.tsv");
		print OUT "$candlist_header\n$gene_conversion{$sid}{line}";
		close OUT;

		my $json_file="$dir/candidate_region.json";
		open(JSON,">$json_file");
		my @regions = split(/,/,$gene_conversion{$sid}{posi});
		my $region = join("\",\n\t\t\"",@regions);
		print JSON "{\n\t\"regions\":[\n\t\t\"$region\"\n\t]\n}";
		close JSON;

		#bam slicing
		my ($tuuid,$tfile)=($gene_conversion{$sid}{uuid},"$dir/$gene_conversion{$sid}{filename}");
		my ($tumor_download,$norm_download)=(0,0);
		while(($tumor_download < 10) && (&bam_check($tfile) ne "ok")){
				system("curl --header \"X-Auth-Token: $token\" --request POST https://api.gdc.cancer.gov/slicing/view/$tuuid --header \"Content-Type: application/json\" -d\@$json_file --output $tfile > /dev/null 2>&1");
				$tumor_download++;
		}
		my ($nuuid,$nfile)=($gene_conversion{$sid}{nuuid},"$dir/$gene_conversion{$sid}{nfilename}");
		while(($norm_download < 10) && (&bam_check($nfile) ne "ok")){
				system("curl --header \"X-Auth-Token: $token\" --request POST https://api.gdc.cancer.gov/slicing/view/$nuuid --header \"Content-Type: application/json\" -d\@$json_file --output $nfile > /dev/null 2>&1");
				$norm_download++;
		}
		if(($tumor_download==0)&&($norm_download==0)){$pm->finish;}
		if($tumor_download == 10){print LOG "$sid: $tfile is has download error\n";$pm->finish;}
		if($norm_download == 10){print LOG "$sid: $nfile is has download error\n";$pm->finish;}
		`samtools index $tfile`;
		`samtools index $nfile`;

		#varscan
		my $mpile = "samtools mpileup -q 10 -f $ref  $nfile $tfile";
		`zsh -c \"varscan somatic <\($mpile\) $dir/$sid --tumor-purity $gene_conversion{$sid}{purity} --p-value 0.1 --output-vcf 1 --mpileup 1\"`;

		print "end $sid\n";
		$pm->finish;
}
$pm->wait_all_children;
		


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
sub bam_check ( $ ){
		my $bam = $_[0];
		if(!-e $bam){
				return("not exist");
		}else{
				my $samtools_head = `samtools view $bam 2>&1|head -n 1`;
				if(($samtools_head !~ /EOF\smarker\sis\sabsent/) && ($samtools_head !~ /Parse\serror/) && ($samtools_head !~ /fail\sto\sread\sthe\sheader/)){
						return "ok";
				}else{
						return "download error";
				}
		}
}
