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
driver_gene = read_tsv("~/Dropbox/cooperative/machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
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
  mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
  mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
  mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  #filter(tVAF>0.8)%>>%
  mutate(gene_conversion=ifelse(p>=0.999 & tVAF>tvaf_cutoff,1,0))

################################ distribution of GC rate ######################################
################################### GR hotspot search #########################################
allmut=count(ac2_all_maf)$n
allGC =sum(ac2_all_maf$gene_conversion)
GC_fisher=function(data){
  alln=data$alln
  gcn=data$gcn
  fisher.test(matrix(c(gcn, alln-gcn,allGC-gcn, allmut-alln-allGC+gcn),nrow=2),
              alternative = "g")$p.value
}
# 1Mbp window
hotspot_1mb=ac2_all_maf %>>%
  mutate(window=start %/% 10^6)%>>%
  group_by(chr,window)%>>%summarise(gcn=sum(gene_conversion),alln=n())%>>%
  ungroup()%>>%nest(data=c(gcn,alln))%>>%
  mutate(fisher_p=purrr::map(data,~GC_fisher(.)))%>>%
  unnest(cols = c(data, fisher_p))
hotspot_1mb$FDR=p.adjust(hotspot_1mb$fisher_p,"fdr")
hotspot_1mb%>>%
  mutate(Significance=ifelse(FDR<0.05,"FDR<0.05","Not significant"))%>>%
  filter(alln>50)%>>%
  mutate(gene_conversion_rate=gcn/alln,
         chr=factor(chr,c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10","chr11",
                          "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")))%>>%
  ggplot()+geom_bar(aes(x=window,y=gene_conversion_rate,fill=Significance),stat="identity")+
  scale_fill_manual(values=c("#e41a1c","#377eb8"))+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+xlab("Position (Mb)")+
  facet_wrap(~chr,scale="free_x")+theme_classic()
ggsave("~/Dropbox/work/somatic_gene_conversion/revise/window_GCrate.pdf",width=12,height = 6)
hotspot_1mb%>>%
  mutate(Window=paste0(window*1000000+1,"-",(window+1)*1000000))%>>%
write_df("~/Dropbox/work/somatic_gene_conversion/revise/hotspot_1mb_window.tsv")

## 10Mbp window
hotspot_10mb=ac2_all_maf %>>%
  mutate(window=start %/% 10^7)%>>%
  group_by(chr,window)%>>%summarise(gcn=sum(gene_conversion),alln=n())%>>%
  ungroup()%>>%nest(data=c(gcn,alln))%>>%
  mutate(fisher_p=purrr::map(data,~GC_fisher(.)))%>>%
  unnest(cols = c(data, fisher_p))
hotspot_10mb$FDR=p.adjust(hotspot_10mb$fisher_p,"fdr")
hotspot_10mb%>>%
  mutate(Significance=ifelse(fisher_p<0.05,"FDR<0.05","Not significant"))%>>%
  filter(alln>50)%>>%
  mutate(gene_conversion_rate=gcn/alln,
         chr=factor(chr,c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10","chr11",
                          "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")))%>>%
  ggplot()+geom_bar(aes(x=window,y=gene_conversion_rate,fill=Significance),stat="identity")+
  ylab("Gene conversion rate")+xlab("Position (10Mb)")+
  facet_wrap(~chr,scale="free_x")+theme_classic()
ggsave("~/Dropbox/work/somatic_gene_conversion/revise/window10_GCrate.pdf",width=12,height = 6)
write_df()

################################### mutation signature #########################################
mutsig=read_tsv("~/Dropbox/work/somatic_gene_conversion/revise/all_pass_SBS.tsv.gz")
mutsig=mutsig %>>%
  filter(genotype=="AB")%>>%
  inner_join(ac2_all_maf%>>%dplyr::select(sample_id,chr,start,ref,alt,gene_conversion))
sbs_rank=c("C>A","C>G","C>T","T>A","T>C","T>G")
compo_rank=c("A-A","A-C","A-G","A-T","C-A","C-C","C-G","C-T",
             "G-A","G-C","G-G","G-T","T-A","T-C","T-G","T-T")
mutsig_signif=mutsig %>>%
  mutate(SBS=paste0(sbs_ref,">",sbs_alt))%>>%
  mutate(SBS=factor(SBS,levels=sbs_rank))%>>%
  mutate(Components=paste0(str_sub(sbs_pattern,1,1),"-",str_sub(sbs_pattern,3,3)))%>>%
  mutate(Components=factor(Components,levels=compo_rank))%>>%
  group_by(SBS,Components)%>>%
  summarise(gcn=sum(gene_conversion),alln=n(),ratio=sum(gene_conversion)/n())%>>%
  nest(data=c(gcn,alln))%>>%
  mutate(fisher_p=purrr::map(data,~GC_fisher(.)))%>>%
  unnest(cols = c(data, fisher_p))
mutsig_signif$FDR=p.adjust(mutsig_signif$fisher_p,"fdr")
mutsig_signif%>>%
  ggplot(aes(x=Components,y=ratio,fill=SBS))+
  geom_bar(stat="identity")+
  facet_grid(.~SBS)+
  scale_y_continuous(limits = c(0,0.0085),expand = c(0,0))+
  scale_fill_manual(values=c(`C>A`="cyan",`C>G`="black",`C>T`="red",
                             `T>A`="grey",`T>C`="green",`T>G`="pink"))+
  geom_text(data=tibble(lab="*",SBS=factor("C>T",levels=sbs_rank)),
            aes(x="T-G",y=0.007,label=lab),size=12)+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+
  theme_classic()+
  theme(axis.title=element_text(size=24),legend.position = "none",
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(angle = 90,size=6,vjust = 0.5, hjust=1),
        strip.placement = "outside",strip.text.x = element_text(size=20))
ggsave("~/Dropbox/work/somatic_gene_conversion/revise/GCrate_signature.pdf",width=12,height = 4)
## TCG>T is significantly high GC rate
## SBS10b (POLE exonuclease domain mutation) is most assosiated with TCG>T
## check PLOE mutation and GC rate
POLEmut=all_maf%>>%filter(gene=="POLE",(impact=="HIGH" ))%>>%
  count(patient_id,impact)%>>%dplyr::select(-n)%>>%
  mutate(polemut="mut")

ac2_all_maf%>>%left_join(POLEmut)%>>%mutate(polemut=ifelse(is.na(polemut),"Others","Truncating"))%>>%
  group_by(polemut)%>>%summarise(GCrate=sum(gene_conversion)/n())%>>%
  mutate(polemut=factor(polemut,levels = c("Truncating","Others")))%>>%
  ggplot(aes(x=polemut,y=GCrate))+
  geom_bar(stat="identity",fill="#e41a1c")+
  scale_y_continuous(limits = c(0,0.0065),expand = c(0,0))+
  theme_classic()+xlab("POLE mutation")+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+
  theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=22),
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(size=20,color="black"))
ggsave("~/Dropbox/work/somatic_gene_conversion/revise/GCrate_POLE.pdf",width=5,height = 5)
############################ GC rate by gene and cancer type #####################################
ac2_all_maf%>>%count(consequence,impact)
# by gene
gcrate_gene=ac2_all_maf%>>%#filter(impact!="MODIFIER")%>>%
  group_by(gene)%>>%summarise(gcn=sum(gene_conversion),alln=n())%>>%
  ungroup()%>>%nest(data=c(gcn,alln))%>>%
  mutate(fisher_p=purrr::map(data,~GC_fisher(.)))%>>%
  unnest(cols = c(data, fisher_p))
gcrate_gene$FDR=p.adjust(gcrate_gene$fisher_p,"fdr")
gcrate_gene%>>%arrange(FDR)%>>%filter(FDR<0.05)%>>%
  write_df("~/Dropbox/work/somatic_gene_conversion/revise/GCrate_gene.tsv")

#by cancer type
gcrate_ct=ac2_all_maf%>>%#filter(impact!="MODIFIER")%>>%
  group_by(cancer_type)%>>%summarise(gcn=sum(gene_conversion),alln=n())%>>%
  mutate(gcrate=gcn/alln)%>>%
  ungroup()%>>%nest(data=c(gcn,alln))%>>%
  mutate(fisher_p=purrr::map(data,~GC_fisher(.)))%>>%
  unnest(cols = c(data, fisher_p))
gcrate_ct$FDR=p.adjust(gcrate_ct$fisher_p,"fdr")
gcrate_ct%>>%arrange(FDR)

gcrate_ct%>>%
  mutate(Significance=ifelse(FDR<0.05,"FDR<0.05","Not significant"))%>>%
  ggplot(aes(x=reorder(cancer_type,desc(gcrate)),y=gcrate,fill=Significance))+
  geom_bar(stat = "identity")+
  scale_y_continuous(limits = c(0,0.014),expand = c(0,0))+
  scale_fill_manual(values=c("#e41a1c","#377eb8"))+
  theme_classic()+xlab("Cancer type")+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+
  theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=22),
        axis.text.y = element_text(size=16,color="black"),
        axis.text.x = element_text(angle = 90,size=10,color="black", vjust = 0.5, hjust=1),
        legend.text = element_text(size=16),legend.title = element_text(size=18),
        strip.placement = "outside",strip.text.x = element_text(size=20),
        legend.position = c(1,1),legend.justification = c(1,1))
ggsave("~/Dropbox/work/somatic_gene_conversion/revise/GCrate_cantype.pdf",width=9,height = 5)
gcrate_ct%>>%arrange(FDR)%>>%#filter(FDR<0.05)%>>%
  write_df("~/Dropbox/work/somatic_gene_conversion/revise/GCrate_cantype.tsv")


gcrate_gene%>>%arrange(FDR)%>>%filter(FDR<0.05)%>>%
  inner_join(driver_gene)
  