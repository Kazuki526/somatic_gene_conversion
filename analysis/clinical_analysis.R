library(tidyverse)
library(pipeR)
loadNamespace('cowplot')
#setwd('/Volumes/HDPX-UT/SomaticGeneConversion/extract_raw_maf/somatic_conversion/')
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
print_tbl.df = function(x,..){print(as.data.frame(x))}


sample_list = read_tsv("sample_list.tsv")
all_maf = read_tsv(paste0("all_pass_with_dist_position_CPE.maf.gz"))
driver_gene = read_tsv("~/Dropbox/cooperative/cancer_driver_machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
  dplyr::rename(gene=`Gene Symbol`,role =`Role in Cancer`,tumor_type=`Tumour Types(Somatic)`)%>>%
  dplyr::select(gene,role,tumor_type)
chr_length = read_tsv("~/Dropbox/work/grch38datas/chr_arm_pq.tsv")%>>%
  group_by(chr)%>>%summarise(start=min(start),end=max(end))
purity_cutoff=0.7
tvaf_cutoff=0.8
ac2_all_maf = all_maf %>>% filter(ascat_major==1,ascat_minor==1)%>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>purity_cutoff)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(variant_type=ifelse(variant_type=="SNP","SNV","indel"))%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  #filter(tVAF>0.8)%>>%
  mutate(gene_conversion=ifelse(p<0.001 & tVAF>tvaf_cutoff,1,0))
################################# #################################
clinic_info=read_tsv("revise2/clinic_data.tsv")%>>%
  filter(!is.na(OS_time),!is.na(PFI_time))%>>%
  filter(OS_time!="#N/A",PFI_time!="#N/A")%>>%
  mutate(OS_time=as.integer(OS_time),PFI_time=as.integer(PFI_time),
         OS=as.numeric(OS),PFI=as.numeric(PFI))
clinic_GC_tbl=sample_list%>>%filter(is.na(screening),purity>purity_cutoff)%>>%
  dplyr::select(patient_id,tumor_sample_id,purity)%>>%
  inner_join(clinic_info)%>>%
  left_join(ac2_all_maf%>>%group_by(patient_id)%>>%summarise(n=n(),GCn=sum(gene_conversion)))%>>%
  mutate(n=ifelse(is.na(n),0,n),GCn=ifelse(is.na(GCn),0,GCn))

############ stage #################

clinic_GC_tbl%>>%filter(n>0)%>>%filter(str_detect(stage,"Stage"))%>>%
  mutate(stage=ifelse(str_detect(stage,"Stage III"),"Stage III",stage))%>>%
  mutate(stage=ifelse(str_detect(stage,"Stage II[ABC]"),"Stage II",stage))%>>%
  mutate(stage=ifelse(str_detect(stage,"Stage I[AB]"),"Stage I",stage))%>>%
  group_by(stage)%>>%
  summarise(GCprop=mean(GCn/n),GCsd=sd(GCn/n),n=n())%>>%
  ggplot(aes(x=stage,y=GCprop))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=GCprop-GCsd,ymax=GCprop+GCsd),width=0.2)+
  geom_text(aes(x=stage,y=-max(GCsd),label=n))+
  theme_classic()+xlab("")+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+
  theme(axis.title=element_text(size=24),axis.text.y = element_text(size = 12,color="black"),
        axis.text.x = element_text(size = 10,color="black",angle = -45,hjust = 0))
ggsave("revise2/byStage_GCrage.pdf")

clinic_GC_tbl%>>%filter(n>0)%>>%filter(str_detect(stage,"Stage"))%>>%
  mutate(stage=ifelse(str_detect(stage,"Stage III"),"Stage III",stage))%>>%
  mutate(stage=ifelse(str_detect(stage,"Stage II[ABC]"),"Stage II",stage))%>>%
  mutate(stage=ifelse(str_detect(stage,"Stage I[AB]"),"Stage I",stage))%>>%
  group_by(stage)%>>%
  summarise(GCprop=sum(GCn)/sum(n),n=n())%>>%
  ggplot(aes(x=stage,y=GCprop))+
  geom_bar(stat="identity")+
  geom_text(aes(x=stage,y=-0.001,label=n))+
  theme_classic()+xlab("")+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+
  theme(axis.title=element_text(size=24),axis.text.y = element_text(size = 12,color="black"),
        axis.text.x = element_text(size = 10,color="black",angle = -45,hjust = 0))
ggsave("revise2/byStage_GCrage_merged.pdf")


############ kaplan myer #################
library(survival)
ggplot(clinic_GC_tbl)+geom_point(aes(x=n,y=GCn))
# number of mutation >2500 の患者のみで比較してみる？上下の境は全体平均の0.0038(54vs35)
KM_group=clinic_GC_tbl%>>%filter(n>2000)%>>%mutate(GCrate=GCn/n)%>>%
  mutate(HL_GCrate=ifelse(GCrate>0.0038,"Gene conversion High","Gene conversion Low"))

pdf("revise2/OS_KM.pdf",width = 6,height = 4)
#with(data = KM_group,Surv(OS_time,OS))  
survdiff(Surv(OS_time/365.25,OS)~HL_GCrate,data=KM_group)
OS_survfit=survfit(Surv(OS_time/365.25,OS)~HL_GCrate,data=KM_group)
plot(OS_survfit,las=1,xlab="Survival Time (years)", ylab="Overall Survivial",col=2:3)
legend("bottomleft",legend=c("High gene conversion","Low gene conversion"),lty=1, col=2:3)
text(x=15,y=0.1,label="p = 0.6")
dev.off()

pdf("revise2/PFI_KM.pdf",width = 6,height = 4)
survdiff(Surv(PFI_time/365.25,PFI)~HL_GCrate,data=KM_group)
PFI_survfit=survfit(Surv(PFI_time/365.25,PFI)~HL_GCrate,data=KM_group)
plot(PFI_survfit,las=1,xlab="Survival Time (years)", ylab="Progression free interval",col=2:3)
legend("bottomleft",legend=c("High gene conversion","Low gene conversion"),lty=1, col=2:3)
text(x=10,y=0.1,label="p = 0.5")
dev.off()


######################### CMS ###############################
CMS_tbl=read_tsv("revise2/CRCSC_data/clinical_molecular_public_all.txt")
CMS_tbl%>>%filter(cms_label!="NOLBL")%>>%dplyr::rename(patient_id=sample)%>>%
  dplyr::select(patient_id,cms_label)%>>%
  inner_join(clinic_GC_tbl)%>>%group_by(cms_label)%>>%
  summarise(GCrate=sum(GCn)/sum(n),GCsd=sd(GCn/n))%>>%(?.)%>>%
  ggplot()+geom_bar(aes(x=cms_label,y=GCrate),stat = "identity")+
  theme_classic()+xlab("colorectal cancer consensus molecular subtype (CMS)")+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+
  theme(axis.title=element_text(size=24),axis.text = element_text(size = 18,color="black"))
ggsave("revise2/byCMS_GCrate_mereged.pdf")

CMS_tbl%>>%filter(cms_label!="NOLBL")%>>%dplyr::rename(patient_id=sample)%>>%
  dplyr::select(patient_id,cms_label)%>>%
  inner_join(clinic_GC_tbl)%>>%filter(n>0)%>>%
  group_by(cms_label)%>>%
  summarise(GCrate=mean(GCn/n),GCsd=sd(GCn/n))%>>%(?.)%>>%
  ggplot()+geom_bar(aes(x=cms_label,y=GCrate),stat = "identity")+
  geom_errorbar(aes(x=cms_label,ymin=GCrate-GCsd,ymax=GCrate+GCsd),width=0.2)+
  theme_classic()+xlab("colorectal cancer consensus molecular subtype (CMS)")+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+
  theme(axis.title=element_text(size=24),axis.text = element_text(size = 18,color="black"))
ggsave("revise2/byCMS_GCrate.pdf")

