#!/usr/bin/perl
use warnings;
use strict;

my $now_dir = "/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion";
my $pwd = `pwd`;chomp $pwd;
if($now_dir ne $pwd){die "ERROR::do on wrong dir\n";}

my $raw_dir = "/Volumes/areca42TB2/gdc/somatic_maf/raw_maf";
my $manifest = "$raw_dir/gdc_manifest.2019-06-20.txt";
-e $manifest or die "ERROR::$manifest is not exist\n";

my $patient_file = "patient_list.tsv";
-e $patient_file or die "ERROR::$patient_file is not exist\n";
my %patients = ();
open(CF,"$patient_file");
my $header = <CF>;chomp $header;
my %col = &header2hash($header);
while(<CF>){
		chomp;
		my @line = split(/\t/,);
		$patients{$line[$col{tumor_sample_id}]}{pid} = $line[$col{patient_id}];
		$patients{$line[$col{tumor_sample_id}]}{canty} = $line[$col{cancer_type}];
		$patients{$line[$col{tumor_sample_id}]}{ascat} = $line[$col{ASCAT}];
		$patients{$line[$col{tumor_sample_id}]}{he} = $line[$col{HE_staining}];
		$patients{$line[$col{tumor_sample_id}]}{cpe} = $line[$col{CPE}];
		$patients{$line[$col{tumor_sample_id}]}{max} =	$line[$col{MAX}];
}
close CF;

open(IN,"$manifest");
$header = <IN>;chomp $header;
%col = &header2hash($header);
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		my ($ct,$soft);
		if($line[$col{filename}] =~ /^TCGA\.(\w+)\.(\w+)\./){($ct,$soft)=($1,$2);}else{die "ERROR::what file $line[$col{filename}]\n";}
		if($soft ne "mutect"){next;}
		print "start $ct =>";$|=1;
		&mutation_pickup("$raw_dir/$line[$col{id}]",$line[$col{filename}]);
		print "done $ct\n";
}
close IN;





#####################################################################################################3
sub mutation_pickup( $ ){
		my ($dir,$file) = @_;
		my ($ct,$soft);
		if($file =~ /^TCGA\.(\w+)\.(\w+)\./){($ct,$soft)=($1,$2);}else{die "ERROR::what file $file\n";}
		open(MAF,"gunzip -c $dir/$file|");
		my $header = <MAF>;
		while($header =~/^#/){$header=<MAF>;}
		chomp $header;
		my %col = &header2hash($header);
		my %mutation=();
		while(<MAF>){
				chomp;
				my @line = split(/\t/,);
				if($line[$col{Chromosome}] !~ /^chr\d+$/){next;}
				my $pid;
				if($line[$col{Tumor_Sample_Barcode}] =~ /^(TCGA-..-....)/){$pid = $1;}else{die "ERROR::sample name error\n$_\n";}
				my $sample = $line[$col{Tumor_Sample_Barcode}];
				if(!defined $patients{$sample}){next;}
				if(($line[$col{t_depth}] <10)||($line[$col{n_depth}] <10)){next;}
				my ($chr,$posi) = ($line[$col{Chromosome}], $line[$col{Start_Position}]);
				my $var ="$line[$col{Reference_Allele}]\t$line[$col{Tumor_Seq_Allele2}]\t$line[$col{Consequence}]\t$line[$col{FILTER}]\t";
				$var .= "$line[$col{t_depth}]\t$line[$col{n_depth}]";
				$mutation{$pid}{$chr}{$posi}{var}=$var;
				$mutation{$pid}{$chr}{$posi}{sample}=$line[$col{Tumor_Sample_Barcode}];
		}
		close MAF;
		open(OUT,"|gzip -c >patients_mutect/patient_$ct.maf.gz");
		print OUT "patient_id\tsample_id\tcancer_type\tchr\tstart\tref\talt\tconsequence\tfilter\t";
		print OUT "t_depth\tn_depth\tallele_num\tASCAT\tHE_staining\tCPE\tMAX\n";
		
		print "go to ascat =>";$|=1;

		open(ASCAT,"gunzip -c all_patient_ascat.tsv.gz|") or die "ERROR::cannot open all_patient ascat file\n";
		$header = <ASCAT>;chomp $header;
		if($header ne "sample\tchr\tstartpos\tendpos\tnMajor\tnMinor\tploidy\tpurity"){die "ERROR::ASCAR file colum changed ??\n";}
		while(<ASCAT>){
				chomp;
				my @line = split(/\t/,);
				my $chr = "chr$line[1]";
				if(!defined $mutation{$line[0]}{$chr}){next;}
				foreach my $posi (keys%{$mutation{$line[0]}{$chr}}){
						if(($posi >= $line[2])&&($posi <= $line[3])){
								my $sample = $mutation{$line[0]}{$chr}{$posi}{sample};
								my $allele_num = $line[4]+$line[5];
								print OUT "$line[0]\t$sample\t$ct\t$chr\t$posi\t$mutation{$line[0]}{$chr}{$posi}{var}\t$allele_num\t";
								print OUT "$patients{$sample}{ascat}\t$patients{$sample}{he}\t";
								print OUT "$patients{$sample}{cpe}\t$patients{$sample}{max}\n";
						}
				}
		}
		close ASCAT;
		close OUT;
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

