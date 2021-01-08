library(tidyverse)
library(pipeR)
library(gridExtra)
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
print_tbl.df = function(x,..){print(as.data.frame(x))}

purity_class = "CPE"
#patient_list = read_tsv("patient_list.tsv")
sample_list = read_tsv("sample_list.tsv")
all_maf = read_tsv(paste0("all_pass_with_dist_position_",purity_class,".maf.gz"))
all_maf%>>%count(patient_id)
driver_gene = read_tsv("~/Dropbox/cooperative/machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
  dplyr::rename(gene=`Gene Symbol`,role =`Role in Cancer`)%>>%dplyr::select(gene,role)%>>%
  filter(str_detect(role,"TSG"))

###################################### allele count == 2 ######################################
all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  #(?.%>>%count(patient_id)%>>%count())%>>%
  #(?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
  group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  #(?.%>>%count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10))%>>%
  group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
  mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
  mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter((allele_num==2 & p>=0.999)|(allele_num==1),tVAF>tvaf_cutoff)%>>%View
###### for Linkage analysis ##############
all_maf%>>%filter(ascat_minor==1,ascat_major==1)%>>%
  inner_join(sample_list%>>%filter(purity>0.5,dcv_median95>dcv_sd95)%>>%dplyr::select(tumor_sample_id,screening,purity),
             by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75) %>>%
  mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
  mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
  mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
  dplyr::select(sample_id,cancer_type,chr,start,ref,alt,gene,variant_type,consequence,impact,
                purity,screening,mutect_dcv_posi,mutect_mut_num,tVAF)%>>%
  write_df("candidate_conversion_variants.tsv")

#cutoff by binom
LOH_filtering = function(purity_cutoff,tvaf_cutoff){
  sample_list = sample_list%>>%
    mutate(screening=ifelse(purity>purity_cutoff,ifelse(dcv_median95/2 > dcv_sd95,NA,"dcv_out"),"purity_out"))
  all_maf %>>%
    filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
    inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
    mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
    mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
    #(?.%>>%count(patient_id)%>>%count())%>>%
    #(?.%>>%count())%>>%(?.%>>%count(variant_type))%>>%
    group_by(variant_type)%>>%mutate(all=n())%>>%ungroup()%>>%
    filter(allele_num<=2,allele_num>0)%>>% 
    #(?.%>>%count(genotype, variant_type,all)%>>%mutate(ratio=round(n/all*1000)/10))%>>%
    group_by(genotype,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
    mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
    mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
    mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
    mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
    filter((allele_num==2 & p>=0.999)|(allele_num==1),tVAF>tvaf_cutoff)%>>%
    count(genotype, variant_type,bef)%>>%mutate(ratio=round(n/bef*1000)/10)%>>%
    mutate(freq=n/sum(n))%>>%
    dplyr::select(genotype,variant_type,n,freq)
}
LOH_filtering(purity_cutoff = purity_cutoff,tvaf_cutoff = tvaf_cutoff)

for(p in c(0.6,0.7,0.8)){
  print(paste0("purity == ",p))
  outbl=LOH_filtering(p,0.75)%>>%dplyr::rename(`tVAF>0.75`= n, `tVAF>0.75 freq`=freq)
  for(tvaf in c(0.8,0.85)){
    outbl=left_join(outbl,LOH_filtering(p,tvaf)%>>%dplyr::rename(!!paste0("tVAF>",tvaf):=n,!!paste0("tVAF>",tvaf," freq"):=freq))
  }
  outbl
  write_df(outbl,paste0("~/Dropbox/work/somatic_gene_conversion/cutoff_table/purity",p*100,".tsv"))
}
#gene conversionとされる変異の数のみ
GC_number = function(purity_cutoff,tvaf_cutoff){
  sample_list = sample_list%>>%
    mutate(screening=ifelse(purity>purity_cutoff,ifelse(dcv_median95/2 > dcv_sd95,NA,"dcv_out"),"purity_out"))
  all_maf %>>%
    filter(ascat_major==1,ascat_minor==1)%>>% 
    filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
    inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
    mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
    mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
    mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
    mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
    filter(((allele_num==2 & p>=0.999)|(allele_num==1)) & tVAF>tvaf_cutoff)%>>%
    count()
}
tibble(purity_cutoff=c(rep(0.6,5),rep(0.7,5),rep(0.8,5)),
       tvaf_cutoff=rep(c(0.75,0.8,0.85,0.9,0.95),3))%>>%
  mutate(ratio=purrr::map2(purity_cutoff,tvaf_cutoff,~GC_number(.x,.y)))%>>%unnest(cols = c(ratio))%>>%#(?.)%>>%
  tidyr::pivot_wider(names_from = "tvaf_cutoff",values_from = "n")%>>%
  write_df("~/Dropbox/work/somatic_gene_conversion/cutoff_table/GC_number.tsv")


#LOH mutationにおけるgene conversion由来のものの割合のみでみる
GC_ratio=function(purity_cutoff,tvaf_cutoff){
  sample_list = sample_list%>>%
    mutate(screening=ifelse(purity>purity_cutoff,ifelse(dcv_median95/2 > dcv_sd95,NA,"dcv_out"),"purity_out"))
  all_maf %>>%
    filter(allele_num<=2,allele_num>0)%>>% 
    filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
    inner_join(sample_list%>>%filter(is.na(screening))%>>%dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
    mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A")))%>>%
    mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
    mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
    mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
    mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
    filter((allele_num==2 & p>=0.999)|(allele_num==1),tVAF>tvaf_cutoff)%>>%
    mutate(GC=ifelse(genotype=="AB",1,0))%>>%
    summarise(GC_ratio=paste0(sum(GC)," / ",n()," (",signif(sum(GC)/n()*100,4),"%)"))
}
tibble(purity_cutoff=c(rep(0.6,5),rep(0.7,5),rep(0.8,5)),
       tvaf_cutoff=rep(c(0.75,0.8,0.85,0.9,0.95),3))%>>%
  mutate(ratio=purrr::map2(purity_cutoff,tvaf_cutoff,~GC_ratio(.x,.y)))%>>%unnest(cols = c(ratio))%>>%#(?.)%>>%
  tidyr::pivot_wider(names_from = "tvaf_cutoff",values_from = "GC_ratio")%>>%
  write_df("~/Dropbox/work/somatic_gene_conversion/cutoff_table/GC_ration_in_LOH.tsv")








