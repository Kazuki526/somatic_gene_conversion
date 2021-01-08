#!/usr/bin/perl
use warnings;
use strict;

my $now_dir = "/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/all_pass";
my $pwd = `pwd`;chomp $pwd;
if($now_dir ne $pwd){die "ERROR::do on wrong dir\n";}

open(IN,"gunzip -c extracted_all_pass.maf.gz|") or die "ERROR::cannot open extracted_all_pass.maf.gz\n";
my $header = <IN>;chomp $header;
my %col = &header2hash($header);
my %patient =();
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		$patient{$line[$col{tumor_sample}]}{pid} = $line[$col{patient_id}];
		$patient{$line[$col{tumor_sample}]}{canty} = $line[$col{cancer_type}];
		$patient{$line[$col{tumor_sample}]}{num}++;
}
close IN;

open(OUT,">maf_patient_list.tsv");
print OUT "patient_id\ttumor_sample_id\tcancer_type\tmutation_num\n";
foreach my $sample (sort keys %patient){
		print OUT "$patient{$sample}{pid}\t$sample\t$patient{$sample}{canty}\t$patient{$sample}{num}\n";
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
