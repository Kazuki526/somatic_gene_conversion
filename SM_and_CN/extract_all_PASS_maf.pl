#!/usr/bin/perl
use warnings;
use strict;

my $now_dir = "/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/all_pass";
my $pwd = `pwd`;chomp $pwd;
if($now_dir ne $pwd){die "ERROR::do on wrong dir\n";}

my $raw_dir = "/Volumes/areca42TB2/gdc/somatic_maf/raw_maf";
my $manifest = "$raw_dir/gdc_manifest.2019-06-20.txt";
-e $manifest or die "ERROR::$manifest is not exist\n";

open(IN,"$manifest");
my $header = <IN>;chomp $header;
my %col = &header2hash($header);
open(OUT,"|gzip -c >extracted_all_pass.maf.gz");
print OUT "patient_id\tcancer_type\tgene\tVariant_Classification\tVariant_Type\tchr\tstart\tref\t";
print OUT "alt\ttumor_sample\tnorm_sample\tmutect_tdepth\tmutect_ndepth\tConsequence\tIMPACT\n";

my %mutations=();
my %indel=();
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		my ($ct,$soft);
		if($line[$col{filename}] =~ /^TCGA\.(\w+)\.(\w+)\./){($ct,$soft)=($1,$2);}else{die "ERROR::what file $line[$col{filename}]\n";}
		&imput_mutation("$raw_dir/$line[$col{id}]",$line[$col{filename}]);
		if(scalar(keys %{$mutations{$ct}}) ==4){
				foreach my $var (sort keys %{$mutations{$ct}{muse}}){
						if((defined $mutations{$ct}{mutect}{$var}{depth})&&
						   (defined $mutations{$ct}{somaticsniper}{$var}{depth})&&
						   (defined $mutations{$ct}{varscan}{$var}{depth})){
								print OUT "$mutations{$ct}{muse}{$var}{vinf}\t$var\t";
								print OUT "$mutations{$ct}{mutect}{$var}{depth}\t";
								print OUT "$mutations{$ct}{muse}{$var}{inf}\n";
						}
				}
				undef %{$mutations{$ct}};
				foreach my $var (sort keys %{$indel{$ct}{varscan}}){
						if(defined $indel{$ct}{mutect}{$var}{depth}){
								print OUT "$indel{$ct}{varscan}{$var}{vinf}\t$var\t";
								print OUT "$indel{$ct}{mutect}{$var}{depth}\t";
								print OUT "$indel{$ct}{varscan}{$var}{inf}\n";
						}
				}
				undef %{$indel{$ct}};
				print "done $ct\n";
		}
}
close IN;





#####################################################################################################3
sub imput_mutation( $ ){
		my ($dir,$file) = @_;
		my ($cancertype,$soft);
		if($file =~ /^TCGA\.(\w+)\.(\w+)\./){($cancertype,$soft)=($1,$2);}else{die "ERROR::what file $file\n";}
		open(MAF,"gunzip -c $dir/$file|");
		my $header = <MAF>;
		while($header =~/^#/){$header=<MAF>;}
		chomp $header;
		my %col = &header2hash($header);
		while(<MAF>){
				chomp;
				my @line = split(/\t/,);
				if($line[$col{Chromosome}] !~ /^chr\d+$/){next;}
				if($line[$col{FILTER}] ne "PASS"){next;}
				my $sample;
				if($line[$col{Tumor_Sample_Barcode}] =~ /^(TCGA-..-....)/){$sample = $1;}else{die "ERROR::sample name error\n$_\n";}
				my $var ="$line[$col{Chromosome}]\t$line[$col{Start_Position}]\t$line[$col{Reference_Allele}]\t$line[$col{Tumor_Seq_Allele2}]\t$line[$col{Tumor_Sample_Barcode}]\t$line[$col{Matched_Norm_Sample_Barcode}]";
				if($line[$col{Variant_Type}] eq "SNP"){
						if($soft eq "muse"){
								$mutations{$cancertype}{$soft}{$var}{vinf} = "$sample\t$cancertype\t$line[$col{Hugo_Symbol}]\t$line[$col{Variant_Classification}]\t$line[$col{Variant_Type}]";
								$mutations{$cancertype}{$soft}{$var}{inf} = "$line[$col{Consequence}]\t$line[$col{IMPACT}]";
						$mutations{$cancertype}{$soft}{$var}{depth} = "have";
						}elsif($soft eq "mutect"){
								$mutations{$cancertype}{$soft}{$var}{depth} = "$line[$col{t_depth}]:$line[$col{t_ref_count}]:$line[$col{t_alt_count}]\t$line[$col{n_depth}]:$line[$col{n_ref_count}]:$line[$col{n_alt_count}]";
						}else{
								$mutations{$cancertype}{$soft}{$var}{depth} = "have";
						}
				}else{
						if($soft eq "varscan"){
								$indel{$cancertype}{$soft}{$var}{vinf} = "$sample\t$cancertype\t$line[$col{Hugo_Symbol}]\t$line[$col{Variant_Classification}]\t$line[$col{Variant_Type}]";
								$indel{$cancertype}{$soft}{$var}{inf} = "$line[$col{Consequence}]\t$line[$col{IMPACT}]";
								$indel{$cancertype}{$soft}{$var}{depth} = "have";
						}else{
								$indel{$cancertype}{$soft}{$var}{depth} = "$line[$col{t_depth}]:$line[$col{t_ref_count}]:$line[$col{t_alt_count}]\t$line[$col{n_depth}]:$line[$col{n_ref_count}]:$line[$col{n_alt_count}]";
						}
				}

		}
		close MAF;
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

