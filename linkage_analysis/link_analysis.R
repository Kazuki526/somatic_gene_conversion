library(tidyverse)
library(pipeR)
library(ggsignif)
loadNamespace('cowplot')
setwd('/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion/')


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

all_maf=read_tsv(paste0("all_pass_with_dist_position_CPE.maf.gz"))
binom_p=all_maf%>>%
  inner_join(sample_list%>>%filter(purity>0.5,dcv_median95>dcv_sd95)%>>%dplyr::select(tumor_sample_id,screening,purity),
             by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75) %>>%
  mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
  mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
  mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
  dplyr::select(sample_id,chr,start,p)


all_tbl = read_tsv("hetero_germ_linked_candidate.tsv")%>>%
  left_join(binom_p)
purity_cutoff=0.7
tvaf_cutoff=0.8


confirmed=all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>purity_cutoff,tVAF>tvaf_cutoff,p>=0.999)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>% (?.%>>%count(sample_id,chr,start,ref,alt)%>>%count()) %>>%
  filter(`alt-ref` >0 & `alt-alt` >0)%>>%(?.%>>%count(sample_id,chr,start,ref,alt)%>>%count())
  #dplyr::select(sample_id,chr,start,ref,alt,purity,germ_variant,distance,linked_read,`ref-ref`,`ref-alt`,`alt-ref`,`alt-alt`,tVAF)%>>%
  #arrange(sample_id,chr,start)%>>%
  #filter(!(`alt-ref` >0 & `alt-alt` >0))%>>%View



all_tbl%>>%inner_join(confirmed%>>%count(sample_id,chr,start,ref,alt))%>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>purity_cutoff,tVAF>tvaf_cutoff,p>=0.999)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter((`alt-ref` + `alt-alt`>5) & !(`alt-ref` >0 & `alt-alt` >0))%>>%View

hist=all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>purity_cutoff,tVAF>tvaf_cutoff,p>=0.999)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>%
  mutate(linked_convert=ifelse(`alt-ref` >0 & `alt-alt` >0, "confirmed","not confirmed"))%>>%
#  mutate(linked_convert=ifelse(`alt-ref` >4 & `alt-alt` >4, "confirmed","not confirmed"))%>>%
#  dplyr::select(-variant_type,-mutect_dcv_posi,-mutect_mut_num,-linked_read,-screening)%>>%
#  group_by(sample_id,chr,start)%>>%filter(any(linked_convert=="confirmed"))%>>%
#  write_df("~/Dropbox/work/somatic_gene_conversion/confirmed_gene_conversion.tsv")
  ggplot()+geom_histogram(aes(x=distance,fill=linked_convert),color="black",binwidth = 10) +
  theme_bw()+xlab("Distance to candidate mutation")+ylab("Count")
  

ratio=all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>purity_cutoff,tVAF>tvaf_cutoff,p>=0.999)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  #filter(`alt-ref` + `alt-alt`>5, `alt-ref` + `alt-alt`<=10)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>%
  mutate(distance_class=ifelse(distance<=200,"1-200","201-"))%>>%
  #mutate(distance_class=ifelse(distance<=100,"1-100",ifelse(distance<=200,"101-200","201-")))%>>%
  #mutate(distance_class=factor(distance_class,levels = c("1-100","101-200","201-")))%>>%
  mutate(linked_convert=ifelse(`alt-ref` >0 & `alt-alt` >0, "confirmed","not confirmed"))%>>%
  count(distance_class,linked_convert)%>>%ungroup()%>>%
  group_by(distance_class)%>>%mutate(Proportion = n/sum(n))%>>%(?.)%>>%
  filter(linked_convert=="confirmed")%>>%
  ggplot(aes(x=distance_class,y=Proportion,fill=linked_convert))+geom_bar(color="black",stat = "identity",width = 0.8)+
#  geom_signif(stat = "identity",aes(x=2,xend=3,y=1.1,yend=1.1,annotation="**"),textsize = 5)+
#  geom_signif(data=data.frame(x=c(1.05,2.05),xend=c(1.95,2.95), y=c(1.1,1.1),yend=c(1.1,1.1),ann=c("NS","**")),
#              stat="identity",aes(x=x,xend=xend, y=y,yend=yend,annotation=ann))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8),limits = c(0,0.85),expand = c(0,0))+
  scale_fill_manual(values=c("red"))+
  guides(fill=F)+
  xlab("Distance")+
  theme_classic()+
  theme(axis.ticks =element_blank(),
        axis.title = element_text(size=20),legend.text = element_text(size=16),
        legend.position = "top",legend.justification = "center",
        axis.text =element_text(size=16,color="black"),panel.grid.major.x = element_blank())
ratio
ggsave("~/Dropbox/work/somatic_gene_conversion/confirmed_distance.pdf",ratio,height = 6,width = 4)
#test 1-100 vs 101-200 : p-value = 0.1037948
fisher.test(matrix(c(76,37,533,335),nrow = 2))
#test 101-200 vs 201- : p-value = 0.009443576
fisher.test(matrix(c(191,187,79,48),nrow = 2))

#few reads (5< =<10)
all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
       purity>purity_cutoff,tVAF>tvaf_cutoff,p>=0.999)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>%
  #filter(`alt-ref` + `alt-alt`>5,`alt-ref` + `alt-alt`<=10)%>>%
  #filter(`alt-ref` + `alt-alt`>10)%>>%
  #mutate(linked_convert=ifelse((`alt-ref` >0 & `alt-alt` >0)& (`alt-ref`/(`alt-ref`+`alt-alt`) >=0.1 & `alt-alt`/(`alt-ref`+`alt-alt`) >=0.1),"confirmed","not confirmed"))%>>%
  mutate(linked_convert=ifelse(`alt-ref` >0 & `alt-alt` >0,"confirmed","not confirmed"))%>>%
  group_by(sample_id,chr,start)%>>%
  summarise(linked_convert=ifelse(any(linked_convert%in%"confirmed"),"confirmed","not confirmed"))%>>%
  ungroup()%>>%count(linked_convert)
#test 5< =<10
fisher.test(matrix(c(499,262,323,187),nrow = 2))

########################################################################################################









cowplot::plot_grid(hist+theme(legend.position ="none" ),NULL,ratio,nrow = 1,rel_widths = c(10,1,10))
ggsave("~/Dropbox/work/somatic_gene_conversion/confirmed_distance.pdf",height = 3,width = 10)

all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>purity_cutoff,tVAF>tvaf_cutoff)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>%
  mutate(linked_convert=ifelse(`alt-ref` >0 & `alt-alt` >0, "confirmed","not confirmed"))%>>%
  group_by(sample_id,chr,start,ref,alt,tVAF)%>>%
  summarise(linked_convert=ifelse(any(linked_convert=="confirmed"),"confirmed","not confirmed"))%>>%
  ungroup()%>>%(?.%>>%count(linked_convert))%>>%
  mutate(tVAF=ifelse(tVAF>1,1,tVAF))%>>%
  ggplot()+geom_histogram(aes(x=tVAF,fill=linked_convert),color="black",position = "fill",bins=25)+
  scale_y_continuous(breaks = c(0,0.5,1),limits = c(0,1),expand = c(0,0))+
  scale_x_continuous(breaks = c(0.75,0.8,0.85,0.9,0.95,1),expand=c(0,0))+
  ylab("Ratio")+theme_classic()+
  theme(axis.ticks =element_blank(),legend.title = element_blank(),
        axis.text.x =element_text(size=10),panel.grid.major.x = element_blank())
ggsave("~/Dropbox/work/somatic_gene_conversion/confirmed_tVAF.pdf",height = 3,width = 6)


############################################################################################################
#cutoffを色々変えて見た
confirmed_ratio=function(purity_cutoff,tvaf_cutoff){
  all_tbl %>>%
    filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
           purity>purity_cutoff,tVAF>tvaf_cutoff,p>=0.999)%>>%
    tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
    filter(`alt-ref` + `alt-alt`>5)%>>% 
    mutate(confirmed=ifelse(`alt-ref` >0 & `alt-alt` >0,1,0))%>>%
    group_by(sample_id,chr,start,ref,alt)%>>%summarise(confirmed=max(confirmed))%>>%ungroup()%>>%
    summarise(ratio=paste0(signif(sum(confirmed)/n()*100,4),"% ( ",sum(confirmed)," / ",n()," )"))
}
tibble(purity_cutoff=c(rep(0.5,5),rep(0.6,5),rep(0.7,5),rep(0.8,5),rep(0.9,5)),
       tvaf_cutoff=rep(c(0.75,0.8,0.85,0.9,0.95),5))%>>%
  mutate(ratio=purrr::map2(purity_cutoff,tvaf_cutoff,~confirmed_ratio(.x,.y)))%>>%unnest(cols=c(ratio))%>>%#(?.)%>>%
  tidyr::pivot_wider(names_from = "tvaf_cutoff",values_from = "ratio")%>>%
  write_df("~/Dropbox/work/somatic_gene_conversion/cutoff_table/confirmed_ration.tsv")


confirm_length=function(purity_cutoff,tvaf_cutoff){
  all_tbl %>>%
    filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
           purity>purity_cutoff,tVAF>tvaf_cutoff,p>=0.999)%>>%
    tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
    filter(`alt-ref` + `alt-alt`>5)%>>% 
    group_by(sample_id,chr,start,ref,alt)%>>%summarise(distance=max(distance))%>>%ungroup()%>>%
    summarise(length=round(mean(distance),2))
}
tibble(purity_cutoff=c(rep(0.6,3),rep(0.7,3),rep(0.8,3),tvaf_cutoff=rep(c(0.75,0.8,0.85),3)))%>>%
  mutate(ratio=purrr::map2(purity_cutoff,tvaf_cutoff,~confirm_length(.x,.y)))%>>%unnest()%>>%#(?.)%>>%
  tidyr::pivot_wider(names_from = "tvaf_cutoff",values_from = "length")%>>%
  write_df("~/Dropbox/work/somatic_gene_conversion/cutoff_table/link_analysis_mean_distance.tsv")


#########################################################################################################
#linkage anlysis
# distance and having evidence
binom_p=all_pass%>>%
  inner_join(sample_list%>>%filter(purity>0.5,dcv_median95>dcv_sd95)%>>%dplyr::select(tumor_sample_id,screening,purity),
             by=c("sample_id"="tumor_sample_id"))%>>%
  mutate(tVAF = t_alt / (t_depth * (purity*allele_num/(purity*allele_num+2*(1-purity)))))%>>%
  filter(tVAF>0.75) %>>%
  mutate(binom_purity=ifelse(purity>0.99,0.99,purity))%>>%
  mutate(binom_purity=ifelse(allele_num==1,binom_purity/(2-binom_purity),binom_purity))%>>%
  mutate(p=pbinom(t_alt,t_depth,binom_purity*0.5))%>>%
  dplyr::select(sample_id,chr,start,p)


all_tbl = read_tsv("hetero_germ_linked_candidate.tsv")%>>%
  left_join(binom_p)

confirmed=all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>0.7,tVAF>0.8,p>=0.999)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref` + `alt-alt`>5)%>>% (?.%>>%count(sample_id,chr,start,ref,alt)%>>%count()) %>>%
  #filter((`alt-ref` >(`alt-ref`+`alt-alt`)*0.1) & (`alt-alt` >(`alt-ref`+`alt-alt`)*0.1))%>>%
  filter(`alt-ref` >1 & `alt-alt` >1)%>>%
  (?.%>>%count(sample_id,chr,start,ref,alt)%>>%count())

dbin=40
confirmedornot=all_tbl %>>%
  filter(is.na(screening),mutect_dcv_posi/mutect_mut_num>0.1,mutect_dcv_posi/mutect_mut_num<0.9,
         purity>0.7,tVAF>0.8,p>=0.999)%>>%
  tidyr::separate(col=linked_state,into=c("ref-ref","ref-alt","alt-ref","alt-alt"),sep=":",convert=T)%>>%
  filter(`alt-ref`+`alt-alt`>5)%>>%
  mutate(distcl=distance%/%dbin)%>%
  mutate(distance_class=paste0(distcl*dbin+1,"-",distcl*dbin+dbin))%>>%
  mutate(linked_convert=ifelse(`alt-ref` >0 & `alt-alt` >0, "confirmed","not confirmed"))
confirm_permutation=function(data=confirmedornot){
  times=0
  focal= confirmedornot%>%
    count(distcl,distance_class,linked_convert)%>>%ungroup()%>>%
    group_by(distance_class)%>>%mutate(Proportion = n/sum(n))%>>%
    filter(linked_convert=="confirmed")%>>%
    {as.list(coef(lm(Proportion ~ distcl, data=.)))$distcl}
  for(i in 1:10000){
    if(i %% 1000 == 0){print(paste0("permutation ",i," times now"))}
    lean = confirmedornot%>%
      {mutate(.,linked_convert=sample(.$linked_convert,length(.$linked_convert)))}%>%
      count(distcl,distance_class,linked_convert)%>>%ungroup()%>>%
      group_by(distance_class)%>>%mutate(Proportion = n/sum(n))%>>%
      filter(linked_convert=="confirmed")%>>%
      {as.list(coef(lm(Proportion ~ distcl, data=.)))$distcl}
    if(lean > focal){times=times+1}
  }
  return(times/10000)
}
confirm_permutation()


confirmedornot %>>%
  count(distcl,distance_class,linked_convert)%>>%ungroup()%>>%
  group_by(distance_class)%>>%mutate(Proportion = n/sum(n))%>>%
  filter(linked_convert=="confirmed")%>>%(?.)%>>%
  ggplot(aes(x=reorder(distance_class,distcl),y=Proportion,fill=linked_convert))+
  geom_bar(color="black",stat = "identity",width = 0.8)+
  geom_smooth(mapping = aes(x = distcl+1, y = Proportion),method = "lm",se=F, color="black")+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1.05),expand = c(0,0))+
  scale_fill_manual(values=c("red"))+
  guides(fill=F)+
  xlab("Distance")+
  geom_text(aes(x=1.5,y=0.9,label="p=0.0266"),size=5)+
  theme_classic()+
  theme(axis.ticks =element_blank(),
        axis.title = element_text(size=20),legend.text = element_text(size=16),
        legend.position = "top",legend.justification = "center",
        axis.text.x =element_text(size=10,color="black", angle=-45,hjust=0),
        axis.text.y =element_text(size=16,color="black"), panel.grid.major.x = element_blank())
ggsave("~/Dropbox/work/somatic_gene_conversion/revise/linkage_distance_proportion.pdf",height = 6,width=6)

confirmedornot %>>%
  count(distcl,distance_class,linked_convert)%>>%ungroup()%>>%
  tidyr::pivot_wider(names_from = "linked_convert",values_from = "n")%>>%
  dplyr::select(-distcl)%>>%
  write_df("~/Dropbox/work/somatic_gene_conversion/revise/linked_distance_table.tsv")
