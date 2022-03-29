library(tidyverse)
library(pipeR)
library(gridExtra)
library(RcppRoll)
library(ggsignif)
loadNamespace('cowplot')

#setwd("/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion/")
setwd('~/Dropbox/work/somatic_gene_conversion/')

#patient_list = read_tsv("patient_list.tsv")
sample_list = read_tsv("sample_list.tsv")
all_maf = read_tsv(paste0("all_pass_with_dist_position_CPE.maf.gz"))
all_maf%>>%count(patient_id)

purity_cutoff=0.7
tvaf_cutoff=0.8

################# mutations on allele count = 2 #######################
ac2_all_maf = all_maf %>>% filter(ascat_major==1,ascat_minor==1)%>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>purity_cutoff)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
  mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  #filter(tVAF>0.8)%>>%
  mutate(gene_conversion=ifelse(p<0.001 & tVAF>tvaf_cutoff,1,0))

############### number of mutation list file #############################
all_maf %>>% #filter(ascat_major==1,ascat_minor==1)%>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>0.7)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  mutate(gene_conversion=ifelse(ascat_major==1 & ascat_minor==1 & p<0.001 & tVAF>0.8,1,0))%>>%
  mutate(ac2=ifelse(ascat_major==1 & ascat_minor==1, 1,0))%>>%
  group_by(patient_id)%>>%summarise(n=n(),ac2=sum(ac2),GCn=sum(gene_conversion))%>>%
  arrange(desc(n))%>>%write_df("revise2/mut_num.tsv")


###################################################################################
#######################  focused most mutated patients  ############################
top_maf=read_tsv("revise2/all_pass_with_dcv_and_signature.maf.gz")
triplet_num=read_tsv("revise2/count_triplet_num_merged.tsv",col_names=c("patient_id","triplet","triplet_site_num"))

top_maf%>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>0.7)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  mutate(gene_conversion=ifelse(ascat_major==1 & ascat_minor==1 & p<0.001 & tVAF>0.8,1,0),
         ac2=ifelse(ascat_major==1 & ascat_minor==1,1,0))%>>%
  group_by(patient_id,triplet,alt_forsig)%>>%
  summarise(n=n(),ac2=sum(ac2),GCn=sum(gene_conversion))%>>%ungroup()%>>%
  inner_join(triplet_num)%>>%mutate(n=as.double(n))%>>%
  mutate(expect_GC=n*n/triplet_site_num)%>>%(?.)%>>%
  group_by(patient_id)%>>%
  summarise(expected_gene_conversion=sum(expect_GC),observed_gene_conversion=sum(GCn),
            SNV_num=sum(n),AC2_mut=sum(ac2))%>>%
  dplyr::arrange(-observed_gene_conversion)%>>%
  #filter(patient_id!="TCGA-FM-8000")%>>%
  write_df("revise2/expected_GC.tsv")
