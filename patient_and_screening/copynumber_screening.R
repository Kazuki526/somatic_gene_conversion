library(tidyverse)
library(pipeR)
library(gridExtra)
loadNamespace('cowplot')
setwd("/Volumes/areca42TB2/gdc/somatic_maf/extract_raw_maf/somatic_conversion/patients_mutect/")
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

patient_list = read_tsv("../patient_list.tsv")
sample_list = patient_list%>>%
  group_by(patient_id)%>>%
  filter(mutation_num==max(mutation_num))%>>%
  filter(tumor_sample_id==first(tumor_sample_id))%>>%ungroup


# calculate Beta/Alpha (BA)
BA_status=function(file,purity_class="CPE"){
  read_tsv(file,col_types = "ccccicccciiidddd")%>>%
    filter(allele_num!=0, t_depth>=10, n_depth>=10) %>>%
    {if(purity_class=="CPE"){.%>>%mutate(CPE=ifelse(is.na(CPE),ifelse(is.na(HE_staining),ASCAT,HE_staining),CPE))}else{.}}%>>%
    mutate(purity=get(purity_class))%>>%
    filter(!is.na(purity))%>>%
    mutate(BA =t_depth/n_depth/((1-purity)+allele_num*purity/2))%>>%
    dplyr::select(sample_id,purity,BA)%>>%
    group_by(sample_id)%>>%
    inner_join(sample_list%>>%dplyr::select(tumor_sample_id)%>>%
                dplyr::rename(sample_id=tumor_sample_id))%>>%
    group_by(sample_id)%>>%
    mutate(n=n(),median=median(BA),mean=mean(BA),sd=sd(BA),rank=row_number(BA))%>>%
    filter(rank<n*0.95)%>>%
    summarise(n=first(n),median=first(median),mean=first(mean),sd=first(sd),purity=first(purity),
              median95=median(BA),mean95=mean(BA),sd95=sd(BA))
}
if(0){
BA_status("patient_ACC.maf.gz")
BA_table = tibble(file=list.files(".",pattern = "^patient_[A-Z]+.maf.gz",full.names = T))%>>%
  mutate(tbl=purrr::map(file,~BA_status(.)))%>>%
  unnest()%>>%dplyr::select(-file)
write_df(BA_table,"dcv_table_CPE.tsv")
}
BA_table=read_tsv("dcv_table_CPE.tsv")
all=BA_table %>>%filter(purity>0.5)%>>%(?.)%>>%
  mutate(`sd/median`=sd95/median95)%>>%
  ggplot()+geom_histogram(aes(x=`sd/median`),bins=100)+theme_bw()
und2=BA_table %>>%filter(purity>0.5)%>>%
  mutate(`sd/median`=sd95/median95)%>>%
  filter(`sd/median`<2)%>>%
  ggplot()+geom_histogram(aes(x=`sd/median`),bins=100)+theme_bw()
plot=cowplot::plot_grid(all,und2,nrow = 2)
plot
ggsave("~/Dropbox/innan/somatic_gene_conversion_MS/supply/DCV_SD_mean_distribution.pdf",height = 10,width = 10)

sample_list %>>%dplyr::select(-ASCAT,-HE_staining,-CPE,-MAX)%>>%
  left_join(BA_table%>>%dplyr::rename(tumor_sample_id=sample_id,mutect_n=n,BA_median95=median95,BA_sd95=sd95)%>>%
              dplyr::select(tumor_sample_id, purity, mutect_n, BA_median95, BA_sd95))%>>%
  mutate(screening=ifelse(purity>0.5,ifelse(BA_median95/2 > BA_sd95,NA,"BA_out"),"purity_out")) %>>%(?.%>>%count(screening))%>>%
  summarise(mean(mutation_num),sd(mutation_num),sum(mutation_num))%>>%
  filter(is.na(screening))%>>%summarise(mean(mutation_num),sd(mutation_num),sum(mutation_num))%>>%
  summarise(mean(mutect_n),sd(mutect_n),min(mutect_n))%>>%
  write_df("../sample_list.tsv")


####################################################################################################
# Beta/Alpha plot
pick_BA = function(file,purity_class="CPE"){
  read_tsv(file,col_types = "ccccicccciiidddd")%>>%
    filter(allele_num!=0, t_depth>=10, n_depth>=10) %>>%
    {if(purity_class=="CPE"){.%>>%mutate(CPE=ifelse(is.na(CPE),ifelse(is.na(HE_staining),ASCAT,HE_staining),CPE))}else{.}}%>>%
    mutate(purity=get(purity_class))%>>%
    filter(!is.na(purity),purity>0.5)%>>%
    mutate(BA =t_depth/n_depth/((1-purity)+allele_num*purity/2))%>>%
    dplyr::select(sample_id,purity,BA)%>>%
    group_by(sample_id)%>>%mutate(BA_correct=BA/median(BA))%>>%ungroup()
}
all_BA = tibble(file=list.files(".",pattern = "^patient_[A-Z]+.maf.gz",full.names = T))%>>%
  mutate(tbl=purrr::map(file,~pick_BA(.)))%>>%
  unnest()%>>%dplyr::select(-file)

all_BA_plot=all_BA%>>%mutate(BA_correct=ifelse(BA_correct>5,5,BA_correct))%>>%
  ggplot()+geom_histogram(aes(x=BA_correct),bins = 100)+
  theme_classic()+ylab("Count")+xlab("")

sample1="TCGA-CV-7104-01A-11D-2012-08"
sample2="TCGA-UU-A93S-01A-21D-A41F-09"
sample3="TCGA-D8-A1JF-01A-11D-A13L-09"
sample4="TCGA-A2-A259-01A-11D-A16D-09"
sample5="TCGA-D8-A27F-01A-11D-A16D-09"

sd_calculate=function(sampleN){
  all_BA%>>%filter(sample_id == sampleN)%>>%
    filter(row_number(BA) < n()*0.95)%>>%
    mutate(BA_correct=BA/median(BA))%>>%
    {round(sd(.$BA_correct),3)}
}
sd_calculate(sample4)


sample_plot=function(sampleN,fil="black"){
  all_BA%>>%filter(sample_id == sampleN)%>>%
    ggplot()+geom_histogram(aes(x=BA_correct),binwidth = 0.05,fill=fil)+
    theme_classic()+ylab("Count")+
    xlab(expression(paste({beta},"/",{alpha},sep = "")))+
    ggtitle(paste0(sampleN," (S.D.=",sd_calculate(sampleN),")"))+
    scale_y_continuous(expand = c(0,1))+
    xlim(0,6)+
    theme(axis.text = element_text(color="black"))
}
sample_plot(sample1,"#7fc97f")
sample_plot("TCGA-UU-A93S-01A-21D-A41F-09")
sample_plot("TCGA-A2-A259-01A-11D-A16D-09")

by_samples=cowplot::plot_grid(sample_plot(sample1,"#1b9e77"),
                              sample_plot(sample3,"#d95f02"),
                              sample_plot(sample2,"#7570b3"),
                              sample_plot(sample4,"#e7298a"),
                              nrow = 2,labels = c("a","c","b","d"))
by_samples


sd1=sd_calculate(sample1)
sd2=sd_calculate(sample2)
sd3=sd_calculate(sample3)
sd4=sd_calculate(sample4)
sd_plot=all_BA%>>%
  filter(purity>0.7)%>>%
  group_by(sample_id)%>>%
  filter(row_number(BA) < n()*0.95)%>>%
  mutate(BA_correct=BA/median(BA))%>>%
  summarise(sd_BA=sd(BA_correct))%>>%
  mutate(sd_BA=ifelse(sd_BA>2,2,sd_BA))%>>%
  ggplot()+geom_histogram(aes(x=sd_BA),binwidth = 0.01)+
  theme_classic()+ylab("Count")+
  xlab(expression(paste("S.D. of ",{beta},"/",{alpha},sep = "")))+
  scale_y_continuous(expand = c(0,1))+
  geom_vline(xintercept = 0.5,size=1,linetype="dashed")+
  annotate("segment",x=sd1,xend=sd1,y=600,yend=400,color="#1b9e77",size=1,arrow=arrow(angle=20,length=unit(5,"mm")))+
  annotate("segment",x=sd2,xend=sd2,y=400,yend=200,color="#d95f02",size=1,arrow=arrow(angle=20,length=unit(5,"mm")))+
  annotate("segment",x=sd3,xend=sd3,y=250,yend=50,color="#7570b3",size=1,arrow=arrow(angle=20,length=unit(5,"mm")))+
  annotate("segment",x=sd4,xend=sd4,y=250,yend=50,color="#e7298a",size=1,arrow=arrow(angle=20,length=unit(5,"mm")))+
  theme(axis.text = element_text(color="black"))
sd_plot

.plot=cowplot::plot_grid(by_samples,sd_plot,nrow=2,rel_heights = c(2,1.5),labels=c("","e"))
.plot
ggsave("~/Dropbox/innan/somatic_gene_conversion_MS2/supply/figs/BA_distributions.pdf",.plot,height = 6,width = 10)

#####################################
plot=all_BA%>>%
  filter(purity>0.7)%>>%
  group_by(sample_id)%>>%
  filter(row_number(BA) < n()*0.95)%>>%
  mutate(BA_correct=BA/median(BA))%>>%
  summarise(sd_BA=sd(BA_correct))%>>%View
#  filter(sd_BA>0.5)%>>%
  mutate(plot=purrr::map(sample_id,~sample_plot(.)))%>>%
  dplyr::select(plot)
ggsave("~/Dropbox/work/somatic_gene_conversion/large_sd_sample_BA.pdf",
       gridExtra::marrangeGrob(plot$plot,nrow = 8,ncol = 1),width =10,height = 18)  
