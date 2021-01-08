#!/usr/bin/perl
use warnings;
use strict;

#depth correlation value (dcv) = beta/alpha

my $purity_class= $ARGV[0];
if(($purity_class ne "ASCAT")&&($purity_class ne "HE_staining")&&($purity_class ne "CPE")&&($purity_class ne "MAX")){
		die "ERROR::ARGV is wrong. input purity class\n";
}

my $now_dir = "/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion";
my $pwd = `pwd`;chomp $pwd;
if($now_dir ne $pwd){die "ERROR::do on wrong dir\n";}

#read all "PASS" mutations
open(MAF,"gunzip -c all_pass_with_ascat.maf.gz|") or die "ERROR::cannot open all_pass_with_ascat.maf.gz\n";
my $header = <MAF>;chomp$header;
my %col = &header2hash($header);
my %maf=();
while(<MAF>){
		chomp;
		my @line = split(/\t/,);
		if($line[$col{n_depth}]==0){next;}
		my $pur = $line[$col{$purity_class}];
		if($pur eq "NA"){$pur = $line[$col{HE_staining}];}
		if($pur eq "NA"){$pur = $line[$col{ASCAT}];}
		my ($td,$nd,$an)  =($line[$col{t_depth}],$line[$col{n_depth}],$line[$col{allele_num}]);
		if($an==0){$an=1;}
		if(($pur eq "NA")||($pur ==0)){next;}
		my $depth_correction_value;
		if((1-$pur)+$pur*$an/2 ==0){
				$depth_correction_value = 100;
		}else{
				$depth_correction_value = $td / $nd / ((1-$pur)+$pur*$an/2);
		}
		my $site = "$line[$col{chr}]\t$line[$col{start}]\t$line[$col{ref}]\t$line[$col{alt}]";
		$maf{$line[$col{cancer_type}]}{$line[$col{sample_id}]}{$site}{line}=join("\t",@line[$col{gene}..$#line])."\t$depth_correction_value";
		$maf{$line[$col{cancer_type}]}{$line[$col{sample_id}]}{$site}{dcv}=$depth_correction_value;
}
close MAF;

open(OUT,"|gzip -c >all_pass_with_dist_position_$purity_class.maf.gz") or die "ERROR::cannot open out file\n";
print OUT "$header\tdepth_correction_value\tmutect_mut_num\tmutect_dcv_posi\n";

foreach my $canty (sort keys %maf){
		print "start $canty => ";$|=1;
		open(MAF,"gunzip -c patients_mutect/patient_$canty.maf.gz|") or die "ERROR::cannnot open $canty patient_mutect file\n";
		my $head = <MAF>;chomp$head;
		my %c = &header2hash($head);
		my %correct_value =();
		while(<MAF>){
				chomp;
				my @line = split(/\t/,);
				if($line[$c{allele_num}]==0){next;}
				my $pur = $line[$c{$purity_class}];
				if($pur eq "NA"){$pur = $line[$c{HE_staining}];}
				if($pur eq "NA"){$pur = $line[$c{ASCAT}];}
				my ($td,$nd,$an)  = ($line[$c{t_depth}],$line[$c{n_depth}],$line[$c{allele_num}]);
				if(($pur eq "NA")||($pur ==0)){next;}
				my $depth_correction_value = $td / $nd / ((1-$pur)+$pur*$an/2);
				$correct_value{$line[$c{sample_id}]}.="$depth_correction_value,";
		}
		close MAF;
		foreach my $sample (keys %{$maf{$canty}}){
				if(!defined$correct_value{$sample}){next;}
				my $pid; 
				if($sample =~ /^(TCGA-..-....)/){$pid=$1;}else{die "ERROR::this sample is not include patient id $1\n";}
				my @dcv = sort {$a <=> $b} split(/,/,$correct_value{$sample});
				my $munum = scalar(@dcv);
				foreach my $site (keys%{$maf{$canty}{$sample}}){
						if($dcv[$#dcv] <= $maf{$canty}{$sample}{$site}{dcv}){
								print OUT "$pid\t$sample\t$canty\t$site\t$maf{$canty}{$sample}{$site}{line}\t$munum\t$munum\n";
						}else{
								for(my $i=0;$i < $munum;$i++){
										if($dcv[$i] >= $maf{$canty}{$sample}{$site}{dcv}){
										print OUT "$pid\t$sample\t$canty\t$site\t$maf{$canty}{$sample}{$site}{line}\t$munum\t$i\n";
												last;
										}elsif($i+100<$munum){
												if($dcv[$i+100] < $maf{$canty}{$sample}{$site}{dcv}){$i+=99;}
										}
								}
						}
				}
		}
		print "end $canty\n";$|=1;
}

close OUT;


#####################################################################################################
sub header2hash ( $ ){
		my $header = $_[0];
		my @colm = split(/\t/,$header);
		my %out = ();
		for(my $i=0; $i < scalar(@colm); $i++){
				$out{$colm[$i]}=$i;
		}
		return(%out);
}
