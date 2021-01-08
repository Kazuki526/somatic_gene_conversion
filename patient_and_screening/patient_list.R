library(tidyverse)
library(pipeR)
library(ggsignif)
library(gridExtra)
library(purrrlyr)
setwd("/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion/")
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

maf_patient = read_tsv("../all_pass/maf_patient_list.tsv")
tcga_purity = read_tsv("/Volumes/areca42TB2/gdc/purity/by_sample/tcga_sample_purity.tsv")
#ascat_patient = read_tsv("/Volumes/areca42TB/tcga/CNA/all_patient/sample_id/patient_sample_info.tsv")
ascat_patient = read_tsv("/Volumes/areca42TB2/gdc/somatic_maf/extracted/somatic_comversion/ascat_classify/patient_sample_info.tsv")
purity = read_tsv("/Volumes/areca42TB2/gdc/purity/Dvir_purity_data.tsv") 

patient_list = maf_patient %>>%
  left_join(tcga_purity%>>%count(tumor_sample_id_full,purity)%>>%
              dplyr::rename(tumor_sample_id=tumor_sample_id_full,tcga_purity=purity)%>>%
              dplyr::select(-n))%>>%
  mutate(sample_id = str_extract(tumor_sample_id,"TCGA-..-....-..."))%>>%
  left_join(ascat_patient%>>%rename(tcga_cna_purity=tcga_purity))%>>%filter(ascat_status=="do_ascat")%>>%
  left_join(purity) %>>%
  dplyr::rename(ASCAT=ascat_purity,HE_staining=tcga_purity)%>>%
  mutate(MAX=pmax(ASCAT,HE_staining,ABSOLUTE,CPE,na.rm = T))%>>%
  filter_at(c("ASCAT","HE_staining","CPE","MAX"),all_vars(is.na(.)|.!=0))%>>%
  dplyr::select(patient_id,cancer_type,sample_id,tumor_sample_id,mutation_num,ascat_ploidy,ASCAT,HE_staining,CPE,MAX)
write_df(patient_list,"patient_list.tsv")

patient_list%>>%
  group_by(patient_id)%>>%
  filter(mutation_num==max(mutation_num))%>>%
  filter(tumor_sample_id==first(tumor_sample_id))%>>%ungroup%>>%
  summarise(n=n(),mun=sum(mutation_num),mean=mean(mutation_num),sd=sd(mutation_num))%>>%View()


############################################################################################################
# After screening
sample_list = read_tsv("sample_list.tsv")
sample_list %>>%mutate(aft_screening=ifelse(purity>0.7,ifelse(is.na(screening),1,0),0))%>>%
  group_by(cancer_type)%>>%
  summarise(n=n(),aft_screening=sum(aft_screening))%>>%
  write_df("~/Dropbox/work/somatic_gene_conversion/patient_list_bycanty.tsv")

sample_list %>>% 
  left_join(ascat_patient%>>%filter(ascat_status=="do_ascat")%>>%dplyr::select(sample_id,sample_id_full))%>>%
  dplyr::select(patient_id,cancer_type,tumor_sample_id,sample_id_full,mutation_num,purity)%>>%
  dplyr::rename(WXS_sample_id=tumor_sample_id,array_sample_id=sample_id_full)%>>%
  write_df("~/Dropbox/work/somatic_gene_conversion/sample_list.tsv")
