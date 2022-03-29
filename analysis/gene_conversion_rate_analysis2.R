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




#########################################################################################
#### gene conversion rate in TSG and nonTSG
driver_gene = read_tsv("~/Dropbox/cooperative/cancer_driver_machine_learning/gene_list/CGC_v89_without_fusion.tsv")%>>%
  dplyr::rename(gene=`Gene Symbol`,role =`Role in Cancer`,tumor_type=`Tumour Types(Somatic)`)%>>%
  dplyr::select(gene,role,tumor_type)

trunc_tbl=data.frame(variant_type=factor(c("Truncating"),levels=c("Truncating","Missense","Silent","inframe_indel")),
                     role=factor(c("TSGs  "),levels=c("TSGs  ","non-TSGs")),ratio=0,g=1)
misil_tbl=data.frame(variant_type=factor(c("Missense","Silent"),levels=c("Truncating","Missense","Silent","inframe_indel")),
                     role=factor("TSGs  ",levels=c("TSGs  ","non-TSGs")),ratio=0,g=1)
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
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter((allele_num==2 & p<0.001)|(allele_num==1),tVAF>tvaf_cutoff)%>>%
  count(genotype, role, variant_type,bef)%>>%mutate(ratio=n/bef)%>>%
  mutate(low=qbinom(0.025,bef,ratio)/bef,up=qbinom(0.975,bef,ratio)/bef)%>>%(?.)%>>%
  ggplot(aes(x=role,y=ratio,fill=role))+
  geom_bar(stat="identity")+facet_wrap(.~ variant_type,strip.position = "bottom")+
  geom_errorbar(aes(ymin=low,ymax=up),width=0.2)+
  #geom_signif(data=trunc_tbl,aes(group=g),
  #            xmin=1,xmax=2,y_position=0.023,annotations="**",textsize = 6)+
  #geom_signif(data=misil_tbl,aes(group=g),
  #            xmin=1,xmax=2,y_position=0.012,annotations="NS",textsize = 6)+
  theme_classic()+
  ylab(expression(paste("Proportion of ",{SM["LOH,Conv"]},sep="")))+
  #scale_y_continuous(limits = c(0,0.0245),expand = c(0,0))+
  theme(axis.title.x = element_blank(),axis.title.y=element_text(size=24),
        legend.position = c(0.5,1),legend.justification = c(0.5,1),
        legend.direction = "horizontal",
        legend.title = element_blank(),legend.text = element_text(size=20),
        axis.text.x = element_blank(),axis.text.y = element_text(size=18,color="black"),
        axis.ticks.x = element_blank(),strip.placement = "outside",
        strip.background = element_blank(),strip.text.x = element_text(size=20))
TSG_all_GCrate
ggsave("~/Dropbox/work/somatic_gene_conversion/TSGvsAllgene.pdf",TSG_all_GCrate,height = 6,width = 8)


#pie-chart of TSG SMLOH
LOHinTSG=all_maf %>>%
  filter(mutect_dcv_posi/mutect_mut_num > 0.1, mutect_dcv_posi/mutect_mut_num < 0.9,
         allele_num<=2,allele_num>0) %>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>0.7)%>>%
               dplyr::select(tumor_sample_id,purity),by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(genotype=ifelse(ascat_major==2,"AA",ifelse(ascat_minor==1,"AB","A"))) %>>%
  left_join(driver_gene%>>%filter(str_detect(role,"TSG"))%>>%mutate(role="TSG")) %>>%
  mutate(role=ifelse(is.na(role),"non TSG","TSG  "))%>>%
  mutate(role=factor(role,levels=c("TSG  ", "non TSG")))%>>%
  filter(impact=="MODERATE"|impact=="HIGH"|variant_classification=="Silent")%>>%
  mutate(variant_type=ifelse(impact=="HIGH","Truncating",ifelse(variant_type=="SNP",
                                                                ifelse(variant_classification=="Silent","Silent","Missense"),"inframe_indel"))) %>>%
  mutate(variant_type=factor(variant_type,levels=c("Truncating","Missense","Silent","inframe_indel")))%>>%
  filter(variant_type!="inframe_indel")%>>%
  group_by(genotype,role,variant_type)%>>%mutate(bef=n())%>>%ungroup()%>>%
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter((allele_num==2 & p<0.001)|(allele_num==1),tVAF>0.8) %>%
  count(genotype, role, variant_type,bef)%>>%mutate(ratio=n/bef)%>>%(?.)
LOHinTSG%>>%
  mutate(LOH_type=ifelse(genotype=="AB","Gene conversion",
                         ifelse(genotype=="AA","UPD","Deletion")))%>>%
  mutate(LOH_type=factor(LOH_type,levels = c("Gene conversion","Deletion","UPD")))%>>%
  ggplot(aes(x=1,y=ratio,fill=variant_type))+
  coord_polar(theta="y")+
  facet_grid(role ~ LOH_type,switch = "y")+
  geom_bar(stat = "identity",color="black",position="fill")+
  theme_void()+
  scale_fill_brewer(palette="Set1")+
  theme(#axis.text=element_blank(),axis.ticks = element_blank(),
    legend.title = element_blank(),legend.text = element_text(size=14),
    strip.background = element_rect(color="black"),
    strip.text = element_text(size=18),strip.text.y = element_text(angle=-90))
ggsave("~/Dropbox/work/somatic_gene_conversion/revise/TSG_piechart.pdf",width=8,height = 5)

############################ GC rate by gene and cancer type #####################################
allmut=count(ac2_all_maf)$n
allGC =sum(ac2_all_maf$gene_conversion)
GC_fisher=function(data){
  alln=data$alln
  gcn=data$gcn
  fisher.test(matrix(c(gcn, alln-gcn,allGC-gcn, allmut-alln-allGC+gcn),nrow=2),
              alternative = "g")$p.value
}
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
  mutate(low=qbinom(0.025,alln,gcrate)/alln,up=qbinom(0.975,alln,gcrate)/alln)%>>%
  ggplot(aes(x=reorder(cancer_type,desc(gcrate)),y=gcrate,fill=Significance))+
  geom_bar(stat = "identity")+geom_errorbar(aes(ymin=low,ymax=up),width=0.2)+
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
  mutate(p=pbinom(t_depth-t_alt,t_depth,1-purity*0.5))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter((allele_num==2 & p<0.001)|(allele_num==1),tVAF>tvaf_cutoff)
temp=CN_variance%>>%
  inner_join(sample_list%>>%filter(is.na(screening),purity>purity_cutoff)%>>%
               dplyr::select(patient_id))%>>%
  left_join(TSG_truncGC%>>%group_by(patient_id)%>>%
              summarise(tsg_trunc_loh="yes"))%>>%
  mutate(prop_1_1= `1_1`/(`1_0`+`2_0`+`1_1`+other))


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
  summarise(ratio=mean(tsg_trunc_loh),n=n(),
            expectation=first(mean_tsg_trunc_loh)*(mean(prop_1_1)/first(mean_prop_1_1)))%>>%
  mutate(low=qbinom(0.025,n,ratio)/n,up=qbinom(0.975,n,ratio)/n)%>>%(?.)%>>%
  ggplot()+
  geom_bar(aes(x=scnv_level_name,y=ratio,fill="Observation"),stat = "identity")+
  geom_errorbar(aes(x=scnv_level_name,ymin=low,ymax=up),width=0.2)+
  geom_line(aes(x=scnv_level_name,y=expectation,group="Expectation",color="Expectation"))+
  geom_point(aes(x=scnv_level_name,y=expectation,color="Expectation"))+
  theme_classic()+  
  scale_y_continuous(expand = c(0,0),limits = c(0,0.043))+
  scale_fill_manual(values = "#F8766D")+
  scale_color_manual(values = '#e7298a')+
  #xlab(expression(paste("Ranked by proportion of (",{alpha[1]},",",{alpha[2]},")=(1,1)",sep="")))+
  #ylab(expression(paste("Proportion corrected by (",{alpha[1]},",",{alpha[2]},")=(1,1) length",sep="")))+
  xlab(expression(paste({italic(k)["(1,1)"]},sep="")))+
  #ylab(expression(paste("Proportion of patient samples with ",{STM["LOH,Conv"]}," in TSG",sep="")))+
  ylab("")+
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
#p=0.0009044




########################################################################
# chromothripsis
ascat=read_tsv("all_patient_ascat.tsv.gz")
ascat %>>% count(sample,chr)%>>%
  mutate(num_breakpoint=n-1)%>>%
  mutate(num_breakpoint=ifelse(num_breakpoint>25,26,num_breakpoint))%>>%
  ggplot()+
  geom_histogram(aes(x=num_breakpoint),binwidth = 1)+
  theme_classic()+
  xlim(c(-0.5,26.5))+xlab("Number of breakpoints per chromosome")+
  ylab("Number of chromosomes")+
  theme(axis.title = element_text(size=24),
        axis.text  = element_text(size=18,color="black"))
ggsave("revise/breakpoint_distribution.pdf")


chromothripsis_chr=ascat %>>%
  count(sample,chr)%>>%
  mutate(num_breakpoint=n-1,chr=paste0("chr",chr))%>>%
  dplyr::rename(patient_id=sample)


ac2_all_maf%>>%left_join(chromothripsis_chr)%>>%
  mutate(num_breakpoint=ifelse(num_breakpoint>25,26,num_breakpoint))%>>%
  group_by(num_breakpoint)%>>%
  summarise(gc_rate=sum(gene_conversion)/n(),gc=sum(gene_conversion),N=n())%>>%
  mutate(low=qbinom(0.025,N,gc_rate)/N,up=qbinom(0.975,N,gc_rate)/N)%>>%(?.)%>>%
  ggplot(aes(x=num_breakpoint,y=gc_rate))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin=low,ymax=up),width=0.2)+
  theme_classic()+
  xlim(c(-0.5,26.5))+xlab("Number of breakpoints per chromosome")+
  ylab(expression(paste("Proportion of ",{SM["LOH,UPD"]},sep="")))+
  theme(axis.title = element_text(size=24),
        axis.text  = element_text(size=18,color="black"))
ggsave("revise/breakpoint_gcrate.pdf")



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
POLEmut=all_maf%>>%
  group_by(patient_id)%>>%mutate(mutnum=n())%>>%
  filter(gene=="POLE",(impact=="HIGH"))%>>%
  #filter(gene=="POLE",start==132676598)%>>%
  count(patient_id,mutnum)%>>%dplyr::select(-n)%>>%
  mutate(polemut="mut")

ac2_all_maf%>>%left_join(POLEmut)%>>%mutate(polemut=ifelse(is.na(polemut),"Others","Truncating"))%>>%
  #group_by(patient_id)%>>%summarise(gc=sum(gene_conversion),n=n(),polemut=first(polemut))%>>%View
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