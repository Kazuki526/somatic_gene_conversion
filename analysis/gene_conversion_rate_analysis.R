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

###################################### allele count == 2 ######################################
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



################# by patient gene conversion rate
by_pGC=ac2_all_maf %>>% group_by(cancer_type,sample_id,purity)%>>%
  summarise(gene_conversion=sum(gene_conversion),all=n())%>>%
  ggplot()+geom_point(aes(x=all,y=gene_conversion))+
  scale_x_log10()+scale_y_log10()+
  theme_classic()+xlab("Number of somatic mutations")+
  ylab(expression(paste("Number of ",{SM["LOH,Conv"]},sep="")))+
  theme(axis.title=element_text(size=24),axis.text = element_text(size = 18,color="black"))
ggsave("~/Dropbox/work/somatic_gene_conversion/by_patient_num_of_mut_conv.pdf",by_pGC,height = 6,width = 6)

### recombination rate and gene conversion rate
chr_arm=read_tsv("~/Dropbox/work/grch38datas/chr_arm_pq.tsv")
recombination_rate=read_tsv("~/Dropbox/work/grch38datas/recomb-hg38/genetic_map_GRCh38.tsv")
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
#########################################################################################################
# indel on homopolymer vs large indel
indel_tbl=read_tsv("~/Dropbox/work/somatic_gene_conversion/revise/indel_homoplymre.tsv.gz")
LOH_label=all_pass%>>%
  inner_join(sample_list%>>%filter(purity>0.7,dcv_median95>dcv_sd95)%>>%dplyr::select(tumor_sample_id,screening,purity),
             by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.8) %>>%
  mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
  mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
  mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
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
# A        Deletion     4.68e- 2
# A        Insertion    6.09e- 1
# AA       Deletion     2.83e-14
# AA       Insertion    7.02e- 4
# AB       Deletion     1.80e-26
# AB       Insertion    4.96e- 1

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

################### patient with TSG trunc LOH by gene conversion and other #################
CN_variance = read_tsv("by_patient_CN_variance.tsv")
TSG_truncGC= all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>purity_cutoff)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A"))) %>>%filter(genotype=="AB")%>>%
  left_join(driver_gene%>>%filter(str_detect(role,"TSG"))%>>%mutate(role="TSG")) %>>%
  mutate(role=ifelse(is.na(role),"non TSG",role))%>>%
  filter(impact=="MODERATE"|impact=="HIGH"|variant_classification=="Silent")%>>%
  mutate(variant_type=ifelse(impact=="HIGH","Truncating",ifelse(variant_type=="SNP",
                                                                ifelse(variant_classification=="Silent","Silent","Missense"),"inframe_indel"))) %>>%
  filter(variant_type=="Truncating",role=="TSG")%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
  mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
  mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter((allele_num==2 & p>=0.999)|(allele_num==1),tVAF>tvaf_cutoff)
temp=CN_variance%>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>purity_cutoff)%>>%
               dplyr::select(patient_id))%>>%
  left_join(TSG_truncGC%>>%group_by(patient_id)%>>%
              summarise(tsg_trunc_loh="yes"))%>>%
  mutate(prop_1_1= `1_1`/(`1_0`+`2_0`+`1_1`+other))

temp%>>%ggplot()+geom_boxplot(aes(x=tsg_trunc_loh,y=prop_1_1))
#permutation test
permutation=function(){
  p=0
  focal=temp%>>%filter(tsg_trunc_loh=="yes")%>>%summarise(mean(prop_1_1))%>>%{.[1,1]}
  for(n in 1:10000){
    if(mean(sample(temp$prop_1_1,80))>focal){p=p+1}
  }
  print(p/10000)
}
permutation()

#level of somatic copy number variation
p1=temp%>>%
  mutate(scnv_level=prop_1_1%/%0.1)%>>%
  mutate(scnv_level_name=paste0(scnv_level*10+10,"-",scnv_level*10,"%"))%>>%
  mutate(scnv_level_name=ifelse(scnv_level<5,"50%-",scnv_level_name))%>>%
  mutate(scnv_level_name=factor(scnv_level_name,
                                levels = c("100-90%","90-80%","80-70%","70-60%","60-50%","50%-")))%>>%
  #mutate(tsg_trunc_loh=ifelse(is.na(tsg_trunc_loh),"Others","Focal samples"))%>>%
  #mutate(tsg_trunc_loh=factor(tsg_trunc_loh,levels = c("Others","Focal samples")))%>>%
  ggplot()+geom_bar(aes(x=scnv_level_name))+ #,fill=tsg_trunc_loh))+
  theme_classic()+  
  scale_y_continuous(expand = c(0,0))+
  ylab("Number of patient samples")+
  xlab(expression(paste({italic(k)["(1,1)"]},sep="")))+
  #scale_fill_manual(breaks=c("Focal samples","Others"),values = c("#F8766D","gray30"))+
  theme(axis.ticks.x = element_blank(),
        #legend.position = "top",legend.title = element_blank(),legend.text = element_text(size=18),
        axis.title.x = element_text(family = "serif",size=32),
        axis.title.y=element_text(size=24),axis.text=element_text(size=18,color="black"))
p1
ggsave("~/Dropbox/work/somatic_gene_conversion/TSGGCvsOther_SCNlevel.pdf",p1,height = 8,width = 10)
temp%>>%mutate(tsg_trunc_loh=ifelse(is.na(tsg_trunc_loh),0,1))%>>%
  mutate(scnv_level=prop_1_1%/%0.1)%>>%
  mutate(scnv_level_name=paste0(scnv_level*10+10,"%-",scnv_level*10,"%"))%>>%
  mutate(scnv_level_name=ifelse(scnv_level<5,"50%-",scnv_level_name))%>>%
  mutate(scnv_level_name=factor(scnv_level_name,
                                levels = c("100%-90%","90%-80%","80%-70%","70%-60%","60%-50%","50%-")))%>>%
  group_by(scnv_level_name)%>>%
  summarise(tsg_trunc_GC = sum(tsg_trunc_loh),sample_n=n())

p2=temp%>>%mutate(tsg_trunc_loh=ifelse(is.na(tsg_trunc_loh),0,1))%>>%
  mutate(scnv_level=prop_1_1%/%0.1)%>>%
  mutate(scnv_level_name=paste0(scnv_level*10+10,"-",scnv_level*10,"%"))%>>%
  mutate(scnv_level_name=ifelse(scnv_level<5,"50%-",scnv_level_name))%>>%
  mutate(scnv_level_name=factor(scnv_level_name,
                                levels = c("100-90%","90-80%","80-70%","70-60%","60-50%","50%-")))%>>%
  group_by(scnv_level,scnv_level_name)%>>%
  summarise(`(1,0)  `=sum(`1_0`),`(2,0)  `=sum(`2_0`),`(1,1)  `=sum(`1_1`),`Others  `=sum(other))%>>%
  tidyr::pivot_longer(cols = c(-scnv_level_name,-scnv_level),names_to = "CN",values_to ="length")%>>%
  mutate(CN=factor(CN,levels = c("Others  ","(2,0)  ","(1,0)  ","(1,1)  ")))%>>%
  ggplot()+geom_bar(aes(x=scnv_level_name,y=length,fill=CN),position = "fill",stat="identity")+
  theme_classic()+  
  scale_y_continuous(expand = c(0,0))+
  labs(fill=expression(paste("(",{alpha[1]},",",{alpha[2]},")",sep="")))+
  xlab(expression(paste({italic(k)["(1,1)"]},sep="")))+
  ylab("Proportion")+
  scale_fill_manual(breaks=c("(1,1)  ","(1,0)  ","(2,0)  ","Others  "),values=c("#e41a1c","#8dd3c7","#ffffb3","#bebada"))+
  theme(axis.ticks.x = element_blank(),
        legend.position = "top",legend.title = element_text(size=18),
        legend.text = element_text(size=18),axis.title.x = element_text(family = "serif",size=32),
        axis.title.y=element_text(size=24),axis.text=element_text(size=18,color="black"))
p2
ggsave("~/Dropbox/work/somatic_gene_conversion/CNvsOther_SCNlevel.pdf",p2,height = 8,width = 10)
cowplot::plot_grid(p1,p2,ncol=1,labels=c("a","b"),label_size = 40)
ggsave("~/Dropbox/work/somatic_gene_conversion/SCNlevel_plots.pdf",height = 16,width = 8)

k_GCrate=temp%>>%mutate(tsg_trunc_loh=ifelse(is.na(tsg_trunc_loh),0,1))%>>%
  mutate(scnv_level=prop_1_1%/%0.1)%>>%
  mutate(scnv_level_name=paste0(scnv_level*10+10,"-",scnv_level*10,"%"))%>>%
  mutate(scnv_level_name=ifelse(scnv_level<5,"50%-",scnv_level_name))%>>%
  mutate(scnv_level_name=factor(scnv_level_name,
                                levels = c("100-90%","90-80%","80-70%","70-60%","60-50%","50%-")))%>>%
  mutate(mean_tsg_trunc_loh=mean(tsg_trunc_loh),mean_prop_1_1=mean(prop_1_1))%>>%
  group_by(scnv_level_name)%>>%
  summarise(ratio=mean(tsg_trunc_loh),
            expectation=first(mean_tsg_trunc_loh)*(mean(prop_1_1)/first(mean_prop_1_1)))%>>%(?.)%>>%
  ggplot()+
  geom_bar(aes(x=scnv_level_name,y=ratio,fill="Observation"),stat = "identity")+
  geom_line(aes(x=scnv_level_name,y=expectation,group="Expectation",color="Expectation"))+
  geom_point(aes(x=scnv_level_name,y=expectation,color="Expectation"))+
  theme_classic()+  
  scale_y_continuous(expand = c(0,0),limits = c(0,0.043))+
  scale_fill_manual(values = "#F8766D")+
  scale_color_manual(values = '#e7298a')+
  #xlab(expression(paste("Ranked by proportion of (",{alpha[1]},",",{alpha[2]},")=(1,1)",sep="")))+
  #ylab(expression(paste("Proportion corrected by (",{alpha[1]},",",{alpha[2]},")=(1,1) length",sep="")))+
  xlab(expression(paste({italic(k)["(1,1)"]},sep="")))+
  ylab(expression(paste("Proportion of patient samples with ",{STM["LOH,Conv"]}," in TSG",sep="")))+
  theme(axis.ticks.x = element_blank(),legend.title = element_blank(),
        legend.position = c(0.5,1),legend.justification = c(1,1),
        legend.text = element_text(size=20),axis.title.x = element_text(family = "serif",size=32),
        axis.title.y=element_text(size=22),axis.text=element_text(size=16,color="black"))
k_GCrate
ggsave("~/Dropbox/work/somatic_gene_conversion/CNandTSGGCvsOther_SCNlevel.pdf",k_GCrate,height = 6,width = 8)

#k(1,1) = 100-90% is significantly high proportion?  =>  X-squared test
Xsq=temp%>>%mutate(tsg_trunc_loh=ifelse(is.na(tsg_trunc_loh),0,1))%>>%
  mutate(scnv_level=prop_1_1%/%0.1)%>>%
  mutate(scnv_level_name=paste0(scnv_level*10+10,"-",scnv_level*10,"%"))%>>%
  mutate(scnv_level_name=ifelse(scnv_level<5,"50%-",scnv_level_name))%>>%
  mutate(scnv_level_name=factor(scnv_level_name,
                                levels = c("100-90%","90-80%","80-70%","70-60%","60-50%","50%-")))%>>%
  mutate(mean_tsg_trunc_loh=mean(tsg_trunc_loh),mean_prop_1_1=mean(prop_1_1))%>>%
  group_by(scnv_level_name)%>>%
  summarise(ratio=sum(tsg_trunc_loh),N=n(),
            expectation=first(mean_tsg_trunc_loh)*(mean(prop_1_1)/first(mean_prop_1_1)))
prop.test(x=Xsq$ratio,n=Xsq$N,p=Xsq$expectation)
#p=0.001292


#### gene conversion rate in TSG and nonTSG
trunc_tbl=data.frame(variant_type=factor(c("Truncating"),levels=c("Truncating","Missense","Silent","inframe_indel")),
                     role=factor(c("TSG  "),levels=c("TSG  ","non TSG")),ratio=0,g=1)
misil_tbl=data.frame(variant_type=factor(c("Missense","Silent"),levels=c("Truncating","Missense","Silent","inframe_indel")),
                     role=factor("TSG  ",levels=c("TSG  ","non TSG")),ratio=0,g=1)
TSG_all_GCrate= all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9) %>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>purity_cutoff)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A"))) %>>%filter(genotype=="AB")%>>%
  left_join(driver_gene%>>%filter(str_detect(role,"TSG"))%>>%mutate(role="TSG")) %>>%
  mutate(role=ifelse(is.na(role),"non TSG","TSG  "))%>>%
  mutate(role=factor(role,levels=c("TSG  ","non TSG")))%>>%
  filter(impact=="MODERATE"|impact=="HIGH"|variant_classification=="Silent")%>>%
  mutate(variant_type=ifelse(impact=="HIGH","Truncating",ifelse(variant_type=="SNP",
                                                                ifelse(variant_classification=="Silent","Silent","Missense"),"inframe_indel"))) %>>%
  mutate(variant_type=factor(variant_type,levels=c("Truncating","Missense","Silent","inframe_indel")))%>>%
  filter(variant_type!="inframe_indel")%>>%
  filter(allele_num<=2,allele_num>0)%>>% 
  group_by(genotype,role,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
  mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
  mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter((allele_num==2 & p>=0.999)|(allele_num==1),tVAF>tvaf_cutoff)%>>%
  count(genotype, role, variant_type,bef)%>>%mutate(ratio=n/bef)%>>%(?.)%>>%
  ggplot(aes(x=role,y=ratio,fill=role))+geom_bar(stat="identity")+facet_wrap(.~ variant_type,strip.position = "bottom")+
  geom_signif(data=trunc_tbl,aes(group=g),
              xmin=1,xmax=2,y_position=0.027,annotations="**",textsize = 6)+
  geom_signif(data=misil_tbl,aes(group=g),
              xmin=1,xmax=2,y_position=0.012,annotations="NS",textsize = 6)+
  theme_classic()+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+
  scale_y_continuous(limits = c(0,0.02999),expand = c(0,0))+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=24),
        legend.position = c(0.5,1),legend.justification = c(0.5,1),
        legend.direction = "horizontal",
        legend.title = element_blank(),legend.text = element_text(size=20),
        axis.text.x = element_blank(),axis.text.y = element_text(size=18,color="black"),
        axis.ticks.x = element_blank(),strip.placement = "outside",
        strip.background = element_blank(),strip.text.x = element_text(size=20))
TSG_all_GCrate
ggsave("~/Dropbox/work/somatic_gene_conversion/TSGvsAllgene.pdf",TSG_all_GCrate,height = 6,width = 9)



cowplot::plot_grid(cowplot::plot_grid(by_pGC,recomb_GC,nrow = 1,rel_widths = c(6,10),
                                      labels = c("a","b"),label_size = 40),
                   cowplot::plot_grid(TSG_all_GCrate,NULL,k_GCrate,nrow=1,rel_widths = c(7.5,1,7.5),
                                      labels=c("c","","d"),label_size = 40),
                   nrow = 2)
ggsave("~/Dropbox/work/somatic_gene_conversion/fig4.pdf",width=16,height = 12)
