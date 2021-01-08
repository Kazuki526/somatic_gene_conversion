#!/usr/bin/perl
use strict;
use warnings;

open(MAF,"gunzip -c extracted_all_pass.maf.gz|")or die "ERROR::cannot open maf file\n";
my $header = <MAF>;chomp $header;
my %col=&header2hash($header);
my %done_sample=();
open(OUT,">tumor_normal_sample_id.tsv");
print OUT "patient_id\tcancer_type\ttumor_sample_id\tnormal_sample_id\n";
while(<MAF>){
		chomp;
		my @line=split(/\t/,);
		if(defined $done_sample{$line[$col{tumor_sample}]}){next;}
		print OUT "$line[$col{patient_id}]\t$line[$col{cancer_type}]\t$line[$col{tumor_sample}]\t$line[$col{norm_sample}]\n";
		$done_sample{$line[$col{tumor_sample}]}="done";
}
close MAF;
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
