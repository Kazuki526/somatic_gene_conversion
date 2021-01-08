#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;


my %gene_conversion=();
open(CAND, "candidate_conversion_variants.tsv") or die  "ERROR::cannot open candidate_conversion_variants.tsv\n";
#open(CAND, "candidate_conversion_variants_test.tsv") or die  "ERROR::cannot open candidate_conversion_variants_test.tsv\n";
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

my $snum=1;
open(OUT,">hetero_germ_linked_candidate.tsv");
print OUT "$candlist_header\tgerm_variant\tdistance\tlinked_read\tlinked_state\n";
foreach my $sid (sort keys %gene_conversion){
		print "$snum: doing $sid\n";$|=1;
		$snum++;
		my $dir="bam_check/$sid";
		my ($tfile,$nfile)=("$dir/$gene_conversion{$sid}{filename}","$dir/$gene_conversion{$sid}{nfilename}");
		my @candidate = split(/,/,$gene_conversion{$sid}{site});
		my @lines = split(/\n/,$gene_conversion{$sid}{line});
		my %vcf=&read_vcf("$dir/$sid");
		if(scalar(keys %vcf) == 0){next;}
		for(my $i=0; $i < scalar(@candidate);$i++){
				my ($chr,$posi,$ref,$alt)=split(/:/,$candidate[$i]);
				foreach my $germ (keys %vcf){
						my @germ_site=split(/:/,$germ);
						if(($germ_site[0] eq $chr) && ($germ_site[1] > $posi-600) && ($germ_site[1] < $posi+600)){
								my $out = &same_allele_check($candidate[$i],$germ,$tfile);
								if($out ne ""){print OUT "$lines[$i]\t$germ\t$out\n";}
						}
				}
		}
}
close OUT;
exit;
		


#######################################################################################################
sub read_vcf ( $ ){
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
				$vcf{$site}="$norm[2]:$norm[3]:$norm[4]\t$tumor[2]:$tumor[3]:$tumor[4]";
		}
		close SNP;
		open(IND,"$vcf_name.indel.vcf");
		while(<IND>){
				if($_=~ /^#/){next;}
				chomp;
				my @line = split(/\t/,);
				if($line[8] ne "GT:GQ:DP:RD:AD:FREQ:DP4"){die "ERROR::$vcf_name.snp.vcf have not FORMAT=GT:GQ:DP:RD:AD:FREQ:DP4\n$_\n";}
				if($line[9] !~ /^0\/1:/){next;}
				my($posi,$ref,$alt)=($line[1],$line[3],$line[4]);
				if(((length($ref)!=1)&&(length($alt)!=1))||((length($ref)==1)&&(length($alt)==1))){die "ERROR::what line??\n$vcf_name.indel:\n$_\n";}
				if(length($ref)==1){ #insertion
						$ref="-";$alt=substr($alt,1);
				}else{ #deletion
						$posi++;$ref=substr($ref,1);$alt="-";
				}
				my $site="$line[0]:$posi:$ref:$alt";
				my @norm = split(/:/,$line[9]);
				my @tumor= split(/:/,$line[10]);
				$vcf{$site}="$norm[2]:$norm[3]:$norm[4]\t$tumor[2]:$tumor[3]:$tumor[4]";
		}
		close IND;
		return(%vcf);
}


#######################################################################################################
sub same_allele_check ( $ $ $ ){
		my @cand_site=split(/:/,$_[0]);
		my @germ_site=split(/:/,$_[1]);
		my $tbam=$_[2];
		my $cand_region = "$cand_site[0]:". ($cand_site[1]-500) ."-". ($cand_site[1]+500);
		my ($former,$former_posi,$later_posi);
		if($cand_site[1] < $germ_site[1]){($former,$former_posi,$later_posi)=("cand",$cand_site[1],$germ_site[1]);
		}elsif($cand_site[1] > $germ_site[1]){($former,$former_posi,$later_posi)=("germ",$germ_site[1],$cand_site[1]);
		}else{return("");}
		my %pair_read=();
		my @link_state=(0,0,0,0);  # cand,germ = (refref, refalt, altref, altalt)
		open(SAM,"samtools view $tbam $cand_region|");
		while(<SAM>){
				chomp;
				my @sam = split(/\t/,);
				if(($sam[6] ne "=") || ($sam[4]<10) || ($sam[5] eq "*") || ($sam[5] =~/H/) || ($sam[1] >= 1024)){next;}
				my $length = &cigar2length($sam[5]);
				if(($former_posi > $sam[3]) && ($former_posi < $sam[3]+$length-1)){
						if(($later_posi>$sam[3]) && ($later_posi < $sam[3]+$length-1)){
								my $cand = &ref_or_alt("this_cand",$sam[3],$sam[5],$sam[9],@cand_site);
								my $germ = &ref_or_alt("this_germ",$sam[3],$sam[5],$sam[9],@germ_site);
								if(($cand eq "") ||($germ eq "")){next;}
								my $linkage = &linkage_posi($cand,$germ);
								if($linkage eq "ERROR"){die "ERROR::ref alt error on $sam[0]. cand:$cand, germ:$germ\ncand=>@cand_site\ngerm=>@germ_site\n";}
								$link_state[$linkage]++;
						}else{
								if($sam[0] eq "HWI-ST115_0183:1:1108:12494:141620#0"){print "now $former, $former_posi, $later_posi, $sam[3]><$sam[7]\n";}
								my $pair_cigar="";
								for(my $tab=11;$tab<scalar(@sam);$tab++){
										if($sam[$tab] =~ /^MC:Z:(.+)$/){$pair_cigar=$1;}else{next;}
								}
								if($pair_cigar eq "*"){next;}
								my$pair_length = &cigar2length($pair_cigar);
								if($sam[0] eq "HWI-ST115_0183:1:1108:12494:141620#0"){print "now $former, $former_posi, $later_posi, $sam[3]><$sam[7], $pair_length\n$_\n";}
								if(($later_posi > $sam[7]) && ($later_posi < $sam[7]+$pair_length-1)){
										if($former eq "cand"){
												$pair_read{$sam[0]}=&ref_or_alt("pthis_cand",$sam[3],$sam[5],$sam[9],@cand_site);
										}else{
												$pair_read{$sam[0]}=&ref_or_alt("pthis_germ",$sam[3],$sam[5],$sam[9],@germ_site);
										}
								}
						}
				}elsif(defined $pair_read{$sam[0]}){
						if($former eq "cand"){
								my $germ=&ref_or_alt("pair_germ",$sam[3],$sam[5],$sam[9],@germ_site);
								if(($germ eq "")||($pair_read{$sam[0]} eq "")){next;}
								my $linkage = &linkage_posi($pair_read{$sam[0]},$germ);
								if($linkage eq "ERROR"){die "ERROR::ref alt error on $sam[0]. cand:$pair_read{$sam[0]}, germ:$germ\ncand=>@cand_site\ngerm=>@germ_site\n";}
								$link_state[$linkage]++;

						}else{
								my $cand=&ref_or_alt("pair_cand",$sam[3],$sam[5],$sam[9],@cand_site);
								if(($pair_read{$sam[0]} eq "")||($cand eq "")){next;}
								my $linkage = &linkage_posi($cand,$pair_read{$sam[0]});
								if($linkage eq "ERROR"){die "ERROR::ref alt error on $sam[0]. cand:$cand, germ:$pair_read{$sam[0]}\ncand=>@cand_site\ngerm=>@germ_site\n";}
								$link_state[$linkage]++;
						}
				}
		}
		close SAM;
		my $linked_read=$link_state[0]+$link_state[1]+$link_state[2]+$link_state[3];
		if($linked_read>1){
				my $out = abs($cand_site[1]-$germ_site[1]) ."\t". $linked_read ."\t". join(":",@link_state);
				return($out);
		}else{
				return("");
		}
}

sub cigar2length ( $ ){
		my $cigar = $_[0];
		my $length=0;
		while($cigar =~ /./){
				$cigar =~ s/^(\d+)([HSIDM])//;
				my($leng,$operator)=($1,$2);
				if($operator =~ /[IHS]/){next;}
				$length+=$leng;
		}
		return($length);
}
sub ref_or_alt( $  $ $ $   $ $ $ $){
		my ($state,$start,$cigar,$seq,$chr,$posi,$ref,$alt)=@_;
		if($alt eq "-"){ # deletion
				if($cigar !~ /D/){return("ref");
				}else{
						my %posi_seq = &posi2seq($start,$cigar,$seq);
						my $return="alt";
						for(my$p=$posi;$p<$posi+length($ref);$p++){
								if(!defined($posi_seq{$p})){$return="";
								}elsif($posi_seq{$p} ne "-"){$return="ref";}
						}
						return($return);
				}
		}elsif($ref eq "-"){ # insertion
				if($cigar !~ /I/){return("ref");
				}else{
						my %posi_seq = &posi2seq($start,$cigar,$seq);
						if(defined $posi_seq{"$posi-"}){
								if($posi_seq{"$posi-"} eq $alt){
										return("alt");
								}else{
										return("ref");
								}
						}else{
								return("ref");
						}
				}
		}else{ # single nucleotide mutation
				my %posi_seq = &posi2seq($start,$cigar,$seq);
				if(!defined $posi_seq{$posi}){die "ERROR::$state= $start, $cigar, $seq, $posi, $ref, $alt\n";}
				if($posi_seq{$posi} eq $ref){
						return("ref");
				}elsif($posi_seq{$posi} eq $alt){
						return("alt");
				}else{return("");}
		}
}
sub posi2seq ( $ $ $ ){
		my ($start,$cigar,$seq)=@_;
		my %posi2seq=();
		my ($posi,$seq_posi)=($start,0);
		while($cigar){
				$cigar =~ s/^(\d+)([HSIDM])//;
				my($leng,$operator)=($1,$2);
				if($operator eq "H"){
				}elsif($operator eq "S"){$seq_posi+=$leng;
				}elsif($operator eq "I"){
						my $p=$posi-1;
						$posi2seq{"$p-"}=substr($seq,$seq_posi,$leng);
						$seq_posi+=$leng;
				}elsif($operator eq "D"){
						while($leng>0){
								$posi2seq{$posi}="-";
								$leng--;$posi++;
						}
				}else{if($operator ne "M"){die "ERROR::consider operator:$operator \n";}
						while($leng>0){
								$posi2seq{$posi}=substr($seq,$seq_posi,1);
								$leng--;$posi++;$seq_posi++;
						}
				}
		}
		return(%posi2seq);
}

sub linkage_posi ( $ $ ){
		my($cand,$germ)=@_;
		if(($cand eq "ref") && ( $germ eq "ref")){return(0);
		}elsif(($cand eq "ref") && ( $germ eq "alt")){return(1);
		}elsif(($cand eq "alt") && ( $germ eq "ref")){return(2);
		}elsif(($cand eq "alt") && ( $germ eq "alt")){return(3);
		}else{return("ERROR");
		}
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

