#!/usr/bin/perl
use warnings;
use strict;

my $now_dir = "/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf";
my $pwd = `pwd`;chomp $pwd;
if($now_dir ne $pwd){die "ERROR::do on wrong dir\n";}


my $patient_file = "somatic_conversion/patient_list.tsv";
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
		$patients{$line[$col{tumor_sample_id}]}{max} = $line[$col{MAX}];
}
close CF;

print "go to all PASS\n";$|=1;

my $file = "all_pass/extracted_all_pass.maf.gz";
open(MAF,"gunzip -c $file|");
$header = <MAF>;chomp $header;
%col = &header2hash($header);
my %mutation=();
my %count=();
while(<MAF>){
		chomp;
		my @line = split(/\t/,);
		my $sample=$line[$col{tumor_sample}];
		if(!defined $patients{$sample}){next;}
		my $var ="$line[$col{ref}]\t$line[$col{alt}]\t$line[$col{gene}]\t$line[$col{Variant_Classification}]\t";
		$var .= "$line[$col{Variant_Type}]\t$line[$col{Consequence}]\t$line[$col{IMPACT}]\t";
		$var .= join("\t",split(":",$line[$col{mutect_tdepth}]))."\t";
		$var .= join("\t",split(":",$line[$col{mutect_ndepth}]));
		$mutation{$patients{$sample}{pid}}{$line[$col{chr}]}{$line[$col{start}]}{var}=$var;
		$mutation{$patients{$sample}{pid}}{$line[$col{chr}]}{$line[$col{start}]}{sample}=$sample;
		$mutation{$patients{$sample}{pid}}{$line[$col{chr}]}{$line[$col{start}]}{vt}=$line[$col{Variant_Type}];
		$count{all}{$line[$col{Variant_Type}]}++;
}
close MAF;
open(OUT,"|gzip -c >somatic_conversion/all_pass_with_ascat.maf.gz");
print OUT "patient_id\tsample_id\tcancer_type\tchr\tstart\tref\talt\tgene\tvariant_classification\tvariant_type\tconsequence\timpact\t";
print OUT "t_depth\tt_ref\tt_alt\tn_depth\tn_ref\tn_alt\tascat_major\tascat_minor\tallele_num\tascat_region\tASCAT\tHE_staining\tCPE\tMAX\n";
print "go to ascat\n";$|=1;

open(ASCAT,"gunzip -c somatic_conversion/all_patient_ascat.tsv.gz|") or die "ERROR::cannot open all_patient ascat file\n";
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
						my $out = "$line[4]\t$line[5]\t$allele_num\t$line[2]-$line[3]";
						print OUT "$line[0]\t$sample\t$patients{$sample}{canty}\t$chr\t$posi\t$mutation{$line[0]}{$chr}{$posi}{var}\t$out\t";
						print OUT "$patients{$sample}{ascat}\t$patients{$sample}{he}\t";
						print OUT "$patients{$sample}{cpe}\t$patients{$sample}{max}\n";
						if($line[5]==0){$count{homo}{$mutation{$line[0]}{$chr}{$posi}{vt}}++;
						}else{$count{hetero}{$mutation{$line[0]}{$chr}{$posi}{vt}}++;}
				}
		}
}
close ASCAT;
close OUT;

open(COUNT, ">somatic_conversion/count_mutation.txt")or die "ERROR::cannot open count out file\n";
print COUNT "all mutation\nSNP: $count{all}{SNP}\nINS: $count{all}{INS}\nDEL: $count{all}{DEL}\n\n";
print COUNT "homo mutation\nSNP: $count{homo}{SNP}\nINS: $count{homo}{INS}\nDEL: $count{homo}{DEL}\n\n";
print COUNT "hetero mutation\nSNP: $count{hetero}{SNP}\nINS: $count{hetero}{INS}\nDEL: $count{hetero}{DEL}\n";
close COUNT;

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
