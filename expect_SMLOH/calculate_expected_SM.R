library(tidyverse)
library(pipeR)
library(gridExtra)
library(RcppRoll)
library(ggsignif)
loadNamespace('cowplot')

#setwd("/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion/")
setwd('~/Dropbox/work/somatic_gene_conversion/')

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

count_mutations=top_maf%>>%
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
  mutate(expect_BM=ac2*n/triplet_site_num/2)%>>%(?.)%>>%
  group_by(patient_id)%>>%
  summarise(expected_biallelic_mutation=sum(expect_BM),observed_gene_conversion=sum(GCn),
            SNV_num=sum(n),AC2_mut=sum(ac2))%>>%
  dplyr::arrange(-SNV_num)
count_mutations%>>%
  write_df("revise2/expected_GC.tsv")

count_mutations%>>%mutate(order=row_number(desc(SNV_num)))%>>%
  filter(order<=50)%>%
  ggplot(aes(x=reorder(patient_id,order)))+
  geom_bar(aes(y=observed_gene_conversion,
               fill="Observed number of SMLOH,Conv"),stat="identity")+
  geom_line(aes(y=expected_biallelic_mutation,group="Expected number of biallelic mutations",
                color="Expected number of biallelic mutations"))+
  geom_point(aes(y=expected_biallelic_mutation,color="Expected number of biallelic mutations"))+
  theme_classic()+
  ylab("Number of SMs")+xlab("Most mutated 50 patients")+
  scale_fill_manual(values = "#F8766D")+
  scale_color_manual(values = "gray50")+
  theme(axis.ticks.x = element_blank(),legend.title = element_blank(),
        legend.position = c(1,1),legend.justification = c(1,1),
        legend.text = element_text(size=18),axis.text.x= element_text(size=8,angle=-45,hjust = 0),
        axis.title=element_text(size=22),axis.text=element_text(size=16,color="black"))
ggsave("~/Dropbox/work/somatic_gene_conversion/revise2/biallelic_expectation_top50.pdf",height = 7,width = 8)

count_mutations%>>%mutate(order=row_number(desc(SNV_num)))%>>%
  filter(order<=300)%>%
  ggplot(aes(x=order))+
  geom_bar(aes(y=observed_gene_conversion,
               fill="Observed number of SMLOH,Conv"),stat="identity")+
  geom_line(aes(y=expected_biallelic_mutation,group="Expected number of biallelic mutations",
                color="Expected number of biallelic mutations"))+
  #geom_point(aes(y=expected_biallelic_mutation,color="Expected number of biallelic mutations"))+
  theme_classic()+
  ylab("Number of SMs")+xlab("Most mutated 300 patients")+
  scale_fill_manual(values = "#F8766D")+
  scale_color_manual(values = "gray50")+
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300),expand = c(0.01,0))+
  scale_y_continuous(expand=c(0.01,0))+
  theme(legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=18),
        axis.title=element_text(size=22),axis.text=element_text(size=16,color="black"))
ggsave("~/Dropbox/work/somatic_gene_conversion/revise2/biallelic_expectation_top300.pdf",height = 7,width = 10)

count_mutations%>>%mutate(order=row_number(desc(SNV_num)))%>>%
  mutate(expect_cumsum=cumsum(expected_biallelic_mutation),
         observe_cumsum=cumsum(observed_gene_conversion))%>>%
  ggplot(aes(x=order))+
  geom_line(aes(y=expect_cumsum,group="Expected number of biallelic mutations",
                color="Expected number of biallelic mutations"))+
  geom_point(aes(y=expect_cumsum,color="Expected number of biallelic mutations"),size=0.5)+
  geom_line(aes(y=observe_cumsum,group="Observed number of SMLOH,Conv",
                color="Expected number of biallelic mutations"))+
  geom_point(aes(y=observe_cumsum,color="Observed number of SMLOH,Conv"),size=0.5)+
  theme_classic()+
  ylab("Number of SMs")+xlab("Patients (order by number of mutations)")+
  scale_color_manual(values = c("gray50","#F8766D"))+
  theme(axis.ticks.x = element_blank(),legend.title = element_blank(),
        legend.position = c(1,0.7),legend.justification = c(1,1),
        legend.text = element_text(size=18),axis.text.x= element_text(size=16),
        axis.title=element_text(size=22),axis.text=element_text(size=16,color="black"))
ggsave("~/Dropbox/work/somatic_gene_conversion/revise2/biallelic_expectation_cumsum.pdf",height = 6,width = 8)



######################################################################sup fig for how to calculate
sbs_rank=c("C>A","C>G","C>T","T>A","T>C","T>G")
compo_rank=c("A-A","A-C","A-G","A-T","C-A","C-C","C-G","C-T",
             "G-A","G-C","G-G","G-T","T-A","T-C","T-G","T-T")
top_maf%>>%filter(patient_id =="TCGA-AX-A2HC")%>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>0.7)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(ac2=ifelse(ascat_major==1 & ascat_minor==1,1,0))%>>%
  group_by(patient_id,triplet,alt_forsig)%>>%
  summarise(n=n(),ac2=sum(ac2))%>>%ungroup()%>>%
  left_join(triplet_num)%>>%mutate(n=as.double(n))%>>%
  mutate(expectation = ac2*n/triplet_site_num/2)%>>%
  mutate(SBS=paste0(substr(triplet,2,2),">",alt_forsig))%>>%
  mutate(compo=paste0(substr(triplet,1,1),"-",substr(triplet,3,3)))%>>%
  ggplot(aes(x=compo))+
  geom_bar(aes(y=n/2000,fill=SBS),stat="identity")+
  geom_point(aes(y=expectation,color="Expectation"),size=1)+
  facet_grid(.~SBS)+
  scale_y_continuous(limits = c(0,9),
                     sec.axis = sec_axis(~ .*2000,name="Number of mutations"))+
  scale_fill_manual(values=c(`C>A`="cyan",`C>G`="black",`C>T`="red",
                             `T>A`="grey",`T>C`="green",`T>G`="pink"))+
  scale_color_manual(values = "gray50")+
  ylab("Expected number of biallelic mutations")+
  xlab("")+
  theme_classic()+
  theme(axis.title=element_text(size=24),legend.position="none",
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(angle = 90,size=6,vjust = 0.5, hjust=1),
        strip.placement = "outside",strip.text.x = element_text(size=20))
ggsave("~/Dropbox/work/somatic_gene_conversion/revise2/example_for_calculatetion.pdf",height = 6,width = 9)
