library(tidyverse)
library(pipeR)
library(ggsignif)
loadNamespace('cowplot')
setwd('~/work/SomaticGeneConversion/extract_raw_maf/somatic_conversion/')
#setwd('/Volumes/HDPX-UT/SomaticGeneConversion/extract_raw_maf/somatic_conversion/')


write_df= function(x, path, delim='\t', na='NA', append=FALSE, col_names=!append, ...) {
  file = if (grepl('gz$', path)) {
    gzfile(path, ...)
  } else if (grepl('bz2$', path)) {
    bzfile(path, ...)
  } else if (grepl('xz$', path)) {
    xzfile(path, ...)
  } else {path}
  utils::write.table(x, file,
                     append=append, quote=FALSE, sep=delim, na=na,
                     row.names=FALSE, col.names=col_names)
}

sample_list = read_tsv("sample_list.tsv")

all_maf=read_tsv("all_pass_with_dist_position_CPE.maf.gz")


### for further linkage analysis candidate list
binom_p=all_maf%>>%
  inner_join(sample_list%>>%filter(purity>0.5,dcv_median95>dcv_sd95)%>>%
               dplyr::select(tumor_sample_id,screening,purity),
             by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75) %>>%
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  dplyr::select(sample_id,chr,start,p)
all_tbl = read_tsv("hetero_germ_linked_candidate.tsv")%>>%
  left_join(binom_p)

futher_focal=all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>0.5,tVAF>0.75,p>=0.999)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>%
  mutate(linked_confirm=ifelse(`alt-ref` >0 & `alt-alt` >0, "confirmed","not confirmed"))%>>%
  group_by(sample_id,chr,start,ref,alt)%>>%
  mutate(linked_confirm=ifelse(any(linked_confirm %in% "confirmed" ),"confirmed","not_confirmed"))%>>%
  ungroup()%>>%filter(linked_confirm=="not_confirmed")%>>%
  mutate(gstart=str_split(germ_variant,":",simplify = T)[,2])%>>%
  group_by(sample_id,chr,start,ref,alt,purity,mutect_dcv_posi,mutect_mut_num,tVAF)%>>%
  summarise(min=min(gstart), min_germ=germ_variant[which.min(gstart)],
            max=max(gstart), max_germ=germ_variant[which.max(gstart)])%>>%
  mutate(min_germ=ifelse(min>start,NA,min_germ),max_germ=ifelse(max<start,NA,max_germ))%>>%
  dplyr::select(sample_id,chr,start,ref,alt,purity,
                mutect_dcv_posi,mutect_mut_num,tVAF,min_germ,max_germ)

write_df(futher_focal,"~/Dropbox/work/somatic_gene_conversion/revise/further_linked_analysis_candidate.tsv")

################################################################################################################
################################################################################################################
################################################################################################################
linkage_tbl=read_tsv("~/Dropbox/work/somatic_gene_conversion/revise/further_hetero_germ_linked_candidate.tsv")
LOH_label=all_maf%>>%
  inner_join(sample_list%>>%filter(purity>0.7,dcv_median95>dcv_sd95)%>>%dplyr::select(tumor_sample_id,screening,purity),
             by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  filter(tVAF>0.8,p<0.001) %>>%
  dplyr::select(sample_id,chr,start)
linkage_tbl%>>%inner_join(LOH_label)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>% (?.%>>%count(sample_id,chr,start,ref,alt)%>>%count()) %>>%
  filter(`alt-ref` >0 & `alt-alt` >0)%>>%(?.%>>%count(sample_id,chr,start,ref,alt)%>>%count())

################################################################################################################
################################################################################################################
################################################################################################################
# list of gene covnersion mutaiton
ac2_all_maf = all_maf %>>% filter(ascat_major==1,ascat_minor==1)%>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>0.7)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  #filter(tVAF>0.8)%>>%
  mutate(gene_conversion=ifelse(p<0.001 & tVAF>0.8,1,0))
gclist=ac2_all_maf%>>%filter(gene_conversion==1)%>>%
  dplyr::select(patient_id,sample_id,chr,start,ref,alt,gene,variant_classification,impact)

linkage_tbl = read_tsv("hetero_germ_linked_candidate.tsv")%>>%
  left_join(binom_p)
have_neib_germ=linkage_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>0.7,tVAF>0.8,p<0.001)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>% count(sample_id,chr,start,ref,alt)%>>%dplyr::select(-n)%>>%
  mutate(have_neib_germ="yes")

have_evid = linkage_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>0.7,tVAF>0.8,p<0.001)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>% 
  filter(`alt-ref` >=1 & `alt-alt` >=1)%>>%count(sample_id,chr,start,ref,alt)%>>%dplyr::select(-n)%>>%
  mutate(have_evid="yes")

further_tbl=read_tsv("~/Dropbox/work/somatic_gene_conversion/revise/further_hetero_germ_linked_candidate.tsv")%>>%
  left_join(binom_p)
have_further=further_tbl %>>%
  filter(mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>0.7,tVAF>0.8,p<0.001)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>% 
  filter(`alt-ref` >=1 & `alt-alt` >=1)%>>%
  count(sample_id,chr,start,ref,alt)%>>%dplyr::select(-n)%>>%mutate(further="yes")

gclist%>>%left_join(have_neib_germ)%>>%
  left_join(have_evid)%>>%left_join(have_further)%>>%
  mutate(have_neib_germ=ifelse(!is.na(have_neib_germ),"yes","no"))%>>%
  mutate(have_evidence=ifelse(!is.na(have_evid)|!is.na(further),"yes",
                              ifelse(have_neib_germ=="yes","no","-")))%>>%
  dplyr::select(-have_evid,-further)%>>%
  arrange(desc(have_neib_germ),desc(have_evidence))%>>%
  write_df("revise/all_gene_conversion_mutation.tsv")

