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
print_tbl.df = function(x,..){print(as.data.frame(x))}



#patient_list = read_tsv("patient_list.tsv")
sample_list = read_tsv("sample_list.tsv")
all_maf = read_tsv(paste0("all_pass_with_dist_position_CPE.maf.gz"))
all_maf%>>%count(patient_id)
driver_gene = read_tsv("~/Dropbox/cooperative/machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
  dplyr::rename(gene=`Gene Symbol`,role =`Role in Cancer`)%>>%dplyr::select(gene,role)%>>%
  filter(str_detect(role,"TSG"))
chr_length = read_tsv("~/Dropbox/work/grch38datas/chr_arm_pq.tsv")%>>%
  group_by(chr)%>>%summarise(start=min(start),end=max(end))
purity_cutoff=0.7
tvaf_cutoff=0.8
#########################################################################################################
#coverageの確認
quantile(all_maf$t_depth,probs = c(0.01,0.05,0.1,0.5,1,5,10)/100)
all_maf%>>%filter(ascat_major==1,ascat_minor==1)%>>%
  inner_join(sample_list%>>%filter(purity>0.7,dcv_median95>dcv_sd95)%>>%dplyr::select(tumor_sample_id,screening,purity),
             by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(l10=ifelse(t_depth<=10,1,0),
         l15=ifelse(t_depth<=15,1,0),
         l20=ifelse(t_depth<=20,1,0))%>>%
  summarise(l10=sum(l10)/n(),l15=sum(l15)/n(),l20=sum(l20)/n())

all_maf%>>%
  filter(ascat_major==1,ascat_minor==1)%>>%
  inner_join(sample_list%>>%filter(purity>0.7,dcv_median95>dcv_sd95)%>>%dplyr::select(tumor_sample_id,screening,purity),
             by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.8) %>>%
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,p<0.001)%>>%
  mutate(l10=ifelse(t_depth<=10,1,0),
         l15=ifelse(t_depth<=15,1,0),
         l20=ifelse(t_depth<=20,1,0))%>>%
  summarise(l10=sum(l10)/n(),l15=sum(l15)/n(),l20=sum(l20)/n())

###################################### allele count == 2 ######################################
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


################################ distribution of GC rate ######################################
################# by patient gene conversion rate
by_pGC=ac2_all_maf %>>% group_by(cancer_type,sample_id,purity)%>>%
  summarise(gene_conversion=sum(gene_conversion),all=n())%>>%(?.%>>%arrange(desc(gene_conversion)))%>>%
  ggplot()+geom_point(aes(x=all,y=gene_conversion))+
  #scale_x_log10()+scale_y_log10()+
  theme_classic()+xlab("Number of somatic mutations")+
  ylab(expression(paste("Number of ",{SM["LOH,Conv"]},sep="")))+
  theme(axis.title=element_text(size=24),axis.text = element_text(size = 18,color="black"))
by_pGC
ggsave("~/Dropbox/work/somatic_gene_conversion/by_patient_num_of_mut_conv.pdf",by_pGC,height = 6,width = 6)

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


#########################################################################################################
#exact test of high indel proportion
fisher.test(matrix(c(1134589-2299, 191738-2679, 2299, 2679),nrow=2))

# indel on homopolymer vs large indel
indel_tbl=read_tsv("~/Dropbox/work/somatic_gene_conversion/revise/indel_homoplymre.tsv.gz")
LOH_label=all_maf%>>%
  inner_join(sample_list%>>%filter(purity>0.7,dcv_median95>dcv_sd95)%>>%dplyr::select(tumor_sample_id,screening,purity),
             by=c("sample_id"="tumor_sample_id"))%>>%
  filter(allele_num>0,allele_num<=2)%>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.8)%>>%
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  dplyr::select(sample_id,chr,start)
indel_class_tbl=indel_tbl %>>%
  mutate(indel_class=ifelse(indel_leng>=5,"Large indel (>5bp)","Other"))%>>%
  mutate(homopolymer=ifelse(variant_type=="DEL",ifelse(later_homply==0,0,indel_leng+later_homply),later_homply))%>>%
  mutate(indel_class=ifelse(homopolymer>=3,"Homopolymer (>3bp)",indel_class))%>>%
  mutate(variant_type=ifelse(variant_type=="INS","Insertion","Deletion"))%>>%
  group_by(genotype,variant_type,indel_class)%>>%
  mutate(N=n())%>>%ungroup()%>>%
  inner_join(LOH_label)%>>%
  count(genotype,variant_type,indel_class,N)%>>%
  mutate(LOH_ratio=n/N)
#calculate pvalue
indel_p=function(indtbl){
  focal_tbl=indtbl%>>%dplyr::select(N,n)
  #print(focal_tbl)
  return(chisq.test(focal_tbl)$p.value)
}
indel_class_tbl%>>%
  dplyr::select(genotype,variant_type,indel_class,N,n,LOH_ratio)%>>%
  nest(data=c(indel_class,N,n,LOH_ratio))%>>%
  mutate(p=purrr::map(data,~indel_p(.)))%>>%
  dplyr::select(-data)%>>%unnest(cols = p)


Conv=indel_class_tbl%>%
  filter(genotype=="AB")%>>%
  ggplot(aes(x=indel_class,y=LOH_ratio,fill=indel_class))+geom_bar(stat="identity")+
  theme_classic()+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+
  scale_y_continuous(limits = c(0,0.031),expand = c(0,0))+
  facet_wrap(.~ variant_type,strip.position = "bottom",ncol = 2)+
  #ggtitle("Gene conversion",)+
  scale_fill_brewer(palette="Set2")+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=24),
        legend.position = "top",legend.justification = "center",
        legend.title = element_blank(),legend.text = element_text(size=18),
        axis.text.x = element_blank(),axis.text.y = element_text(size=18,color="black"),
        axis.ticks.x = element_blank(),strip.placement = "outside",
        strip.background = element_blank(),strip.text.x = element_text(size=24),
        plot.title= element_text(size=24))
Conv
ggsave("~/Dropbox/work/somatic_gene_conversion/revise/indel_class_conv.pdf",Conv,width=9,height = 6)
Del=indel_class_tbl%>%
  filter(genotype=="A")%>>%
  ggplot(aes(x=indel_class,y=LOH_ratio,fill=indel_class))+geom_bar(stat="identity")+
  theme_classic()+
  ylab(expression(paste("Proportion of ",{SM["LOH,Del"]},sep="")))+
  scale_y_continuous(limits = c(0,0.5),expand = c(0,0))+
  facet_wrap(.~ variant_type,strip.position = "bottom",ncol = 2)+
  ggtitle("Deletion",)+scale_fill_brewer(palette="Set2")+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=24),
        legend.position = "top",legend.justification = "center",
        legend.title = element_blank(),legend.text = element_text(size=15),
        axis.text.x = element_blank(),axis.text.y = element_text(size=18,color="black"),
        axis.ticks.x = element_blank(),strip.placement = "outside",
        strip.background = element_blank(),strip.text.x = element_text(size=20),
        plot.title= element_text(size=24))
UPD=indel_class_tbl%>%
  filter(genotype=="AA")%>>%
  ggplot(aes(x=indel_class,y=LOH_ratio,fill=indel_class))+geom_bar(stat="identity")+
  theme_classic()+
  ylab(expression(paste("Proportion of ",{SM["LOH,UPD"]},sep="")))+
  scale_y_continuous(limits = c(0,0.4),expand = c(0,0))+
  facet_wrap(.~ variant_type,strip.position = "bottom",ncol = 2)+
  ggtitle("UPD",)+scale_fill_brewer(palette="Set2")+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=24),
        legend.position = "none",legend.justification = "center",
        legend.title = element_blank(),legend.text = element_text(size=16),
        axis.text.x = element_blank(),axis.text.y = element_text(size=18,color="black"),
        axis.ticks.x = element_blank(),strip.placement = "outside",
        strip.background = element_blank(),strip.text.x = element_text(size=20),
        plot.title= element_text(size=24))

cowplot::plot_grid(Del,UPD, nrow= 2)
#                   labels = c("A","B","C"),label_size = 40)
ggsave("~/Dropbox/work/somatic_gene_conversion/revise/indel_class.pdf",width=6,height = 10)
write_df(indel_class_tbl,"~/Dropbox/work/somatic_gene_conversion/revise/indel_class_tbl.tsv")


### recombination rate and gene conversion rate
chr_arm=read_tsv("~/Dropbox/work/grch38datas/chr_arm_pq.tsv")
recombination_rate=read_tsv("~/Dropbox/work/grch38datas/recomb-hg38/genetic_map_GRCh38.tsv")
recomb_GC=recombination_rate%>>%
  filter(chr!="chrX",chr!="chrY")%>>%
  left_join(chr_arm%>>%dplyr::rename(arm_start=start,arm_end=end))%>>%
  filter(arm_start<pos,arm_end>pos)%>>%
  mutate(arm_posi=(pos-arm_start)/(arm_end-arm_start))%>>%
  mutate(arm_posi_class=(arm_posi%/%0.01)/100)%>>%
  group_by(chr,arm,arm_posi_class)%>>%filter(pos==max(pos))%>>%ungroup()%>>%
  group_by(chr)%>>%filter(!(chr=="chr21" & arm=="p")) %>>%
  mutate(recomb_rate=(pos_cm-ifelse(is.na(lag(pos_cm,5)),0,lag(pos_cm,5)))/(pos-ifelse(is.na(lag(pos,5)),0,lag(pos,5))))%>>%
  mutate(bef_pos=ifelse(is.na(lag(pos,5)),0,lag(pos,5)))%>>%filter(arm_posi>=0.05)%>>%
  ungroup()%>>%group_by(arm,arm_posi_class)%>>%
  summarise(recomb_rate=sum(recomb_rate*(pos-bef_pos))/sum(pos-bef_pos))%>>%ungroup()%>>%
  group_by(arm)%>>%mutate(recomb_ratio=recomb_rate/mean(recomb_rate))%>>%
  left_join(ac2_all_maf %>>%dplyr::select(chr,start,gene_conversion)%>>%
              left_join(chr_arm%>>%dplyr::rename(arm_start=start,arm_end=end))%>>%
              filter(arm_start<start,arm_end>start)%>>%
              mutate(arm_posi=(start-arm_start)/(arm_end-arm_start))%>>%
              mutate(arm_posi_class=(arm_posi%/%0.01)/100)%>>%
              group_by(chr,arm,arm_posi_class)%>>%
              summarise(gene_conversion=sum(gene_conversion),num=n())%>>%
              group_by(chr,arm)%>>%mutate(gene_conversion=roll_sum(gene_conversion,n=5,align="right",fill=0),
                                          num=roll_sum(num,n=5,align="right",fill=0))%>>%
              group_by(arm,arm_posi_class)%>>%summarise(gene_conversion_rate=sum(gene_conversion)/sum(num))%>>%
              filter(!is.na(gene_conversion_rate))%>>%
              ungroup()%>>%group_by(arm)%>>%mutate(gene_conversion_ratio=gene_conversion_rate/mean(gene_conversion_rate)))%>>%
  ungroup()%>>%mutate(arm=ifelse(arm=="p","Short arm","Long arm"),arm_posi_class=arm_posi_class-0.025)%>>%
  mutate(arm=factor(arm,levels = c("Short arm","Long arm")))%>>%
  ggplot()+
  geom_line(aes(x=arm_posi_class,y=gene_conversion_ratio,color="GC"))+
  geom_line(aes(x=arm_posi_class,y=recomb_ratio,color="MR"))+
  scale_y_continuous(name="Relative frequency")+
  scale_x_continuous(limits = c(0,1))+#,breaks = c(0.00,0.25,0.50,0.75,1.00),labels = c(1.00,0.75,0.50,0.25,0.00))+
  xlab("Relative distance from centromere")+
  facet_wrap(~arm)+theme_classic()+
  scale_color_manual(name=NULL,values = c(GC="red",MR="blue"),labels=c("Gene conversion","Meiotic recombination"))+
  theme(axis.title.y=element_text(size=24),
        axis.title.x=element_text(size=24,margin = margin(t=5)),
        axis.text = element_text(size = 18,color="black"),
        panel.background =element_rect(colour = "black"), #panel.background = element_blank(),axis.line = element_line(),
        strip.text.x = element_text(size=20),strip.background = element_rect(color="white"),
        strip.placement = "inside",panel.spacing = unit(2, "lines"),
        legend.text = element_text(size=18),
        legend.position = c(0.02,0.05),legend.justification = c(0,0))
recomb_GC
ggsave("~/Dropbox/work/somatic_gene_conversion/recomb_conv_rate_in_arm.pdf",recomb_GC,height = 6,width = 12)



