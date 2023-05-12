### Microbial GWA 
### 2023-3-20
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SNV/pipeline/4_R")

## 1 Preparation
### 1.1 Import
#install.packages("survminer")
rm(list = ls())
options(scipen = 200)

source("functions.R")
knitr::opts_chunk$set(echo = TRUE)
library("R.utils")
library(showtext)
showtext_auto()

#library("survival")
#library("ggplot2")
#library("survminer")
#library(gridExtra)
#library("grid")
library(reshape2)	  
#library("RColorBrewer")
#library("plyr")

###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

all_basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

mel_basic<-subset(all_basic,cancer_type=="melanoma",drop=T)
NSCLC_basic<-subset(all_basic,cancer_type=="NSCLC",drop=T)
RCC_basic<-subset(all_basic,cancer_type=="RCC",drop=T)


###################################################################################################
### 3 Results
### 3.1 Clean results
#### 3.1.1 vSV associations

## vsv
## merge result tables

melanoma_vsv_resp_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/melanoma_vsv_resp_adjAbun_res.csv",header=T)
melanoma_vsv_resp_adjAbun_overall_sig<-subset(melanoma_vsv_resp_adjAbun_overall,p<0.05,drop=T)
melanoma_vsv_resp_adjAbun_datasets<-read.csv("08.Microbial_GWAS/datasets/mel_vsv_resp_pvalue.csv",header=T)
melanoma_vsv_resp<-merge(melanoma_vsv_resp_adjAbun_overall_sig,melanoma_vsv_resp_adjAbun_datasets,by.x="X",by.y="X")
melanoma_vsv_resp_final<-left_join(melanoma_vsv_resp, vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_resp_final,"08.Microbial_GWAS/combine/melanoma_vsv_resp.csv",row.names=F)

melanoma_vsv_os_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/melanoma_vsv_os_adjAbun_res.csv",header=T)
melanoma_vsv_os_adjAbun_overall_sig<-subset(melanoma_vsv_os_adjAbun_overall,p<0.05,drop=T)
melanoma_vsv_os_adjAbun_datasets<-read.csv("08.Microbial_GWAS/datasets/mel_vsv_os_pvalue.csv",header=T)
melanoma_vsv_os<-merge(melanoma_vsv_os_adjAbun_overall_sig,melanoma_vsv_os_adjAbun_datasets,by.x="X",by.y="X")
melanoma_vsv_os_final<-left_join(melanoma_vsv_os, vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_os_final,"08.Microbial_GWAS/combine/melanoma_vsv_os.csv",row.names=F)

melanoma_vsv_pfs_12_months_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/melanoma_vsv_pfs_12_months_adjAbun_res.csv",header=T)
melanoma_vsv_pfs_12_months_adjAbun_overall_sig<-subset(melanoma_vsv_pfs_12_months_adjAbun_overall,p<0.05,drop=T)
melanoma_vsv_pfs_12_months_adjAbun_datasets<-read.csv("08.Microbial_GWAS/datasets/mel_vsv_pfs_12_months_pvalue.csv",header=T)
melanoma_vsv_pfs_12_months<-merge(melanoma_vsv_pfs_12_months_adjAbun_overall_sig,melanoma_vsv_pfs_12_months_adjAbun_datasets,by.x="X",by.y="X")
melanoma_vsv_pfs_12_months_final<-left_join(melanoma_vsv_pfs_12_months, vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_pfs_12_months_final,"08.Microbial_GWAS/combine/melanoma_vsv_pfs_12_months.csv",row.names=F)

melanoma_vsv_irAEs_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/melanoma_vsv_irAEs_adjAbun_res.csv",header=T)
melanoma_vsv_irAEs_adjAbun_overall_sig<-subset(melanoma_vsv_irAEs_adjAbun_overall,p<0.05,drop=T)
melanoma_vsv_irAEs_final<-left_join(melanoma_vsv_irAEs_adjAbun_overall_sig, vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_irAEs_final,"08.Microbial_GWAS/combine/melanoma_vsv_irAEs.csv",row.names=F)


melanoma_dsv_resp_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/melanoma_dsv_resp_adjAbun_res.csv",header=T)
melanoma_dsv_resp_adjAbun_overall_sig<-subset(melanoma_dsv_resp_adjAbun_overall,p<0.05,drop=T)
melanoma_dsv_resp_adjAbun_datasets<-read.csv("08.Microbial_GWAS/datasets/mel_dsv_resp_pvalue.csv",header=T)
melanoma_dsv_resp<-merge(melanoma_dsv_resp_adjAbun_overall_sig,melanoma_dsv_resp_adjAbun_datasets,by.x="X",by.y="X")
melanoma_dsv_resp_final<-left_join(melanoma_dsv_resp, dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_resp_final,"08.Microbial_GWAS/combine/melanoma_dsv_resp.csv",row.names=F)

melanoma_dsv_os_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/melanoma_dsv_os_adjAbun_res.csv",header=T)
melanoma_dsv_os_adjAbun_overall_sig<-subset(melanoma_dsv_os_adjAbun_overall,p<0.05,drop=T)
melanoma_dsv_os_adjAbun_datasets<-read.csv("08.Microbial_GWAS/datasets/mel_dsv_os_pvalue.csv",header=T)
melanoma_dsv_os<-merge(melanoma_dsv_os_adjAbun_overall_sig,melanoma_dsv_os_adjAbun_datasets,by.x="X",by.y="X")
melanoma_dsv_os_final<-left_join(melanoma_dsv_os, dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_os_final,"08.Microbial_GWAS/combine/melanoma_dsv_os.csv",row.names=F)

melanoma_dsv_pfs_12_months_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/melanoma_dsv_pfs_12_months_adjAbun_res.csv",header=T)
melanoma_dsv_pfs_12_months_adjAbun_overall_sig<-subset(melanoma_dsv_pfs_12_months_adjAbun_overall,p<0.05,drop=T)
melanoma_dsv_pfs_12_months_adjAbun_datasets<-read.csv("08.Microbial_GWAS/datasets/mel_dsv_pfs_12_months_pvalue.csv",header=T)
melanoma_dsv_pfs_12_months<-merge(melanoma_dsv_pfs_12_months_adjAbun_overall_sig,melanoma_dsv_pfs_12_months_adjAbun_datasets,by.x="X",by.y="X")
melanoma_dsv_pfs_12_months_final<-left_join(melanoma_dsv_pfs_12_months, dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_pfs_12_months_final,"08.Microbial_GWAS/combine/melanoma_dsv_pfs_12_months.csv",row.names=F)

melanoma_dsv_irAEs_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/melanoma_dsv_irAEs_adjAbun_res.csv",header=T)
melanoma_dsv_irAEs_adjAbun_overall_sig<-subset(melanoma_dsv_irAEs_adjAbun_overall,p<0.05,drop=T)
melanoma_dsv_irAEs_final<-left_join(melanoma_dsv_irAEs_adjAbun_overall_sig, dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_irAEs_final,"08.Microbial_GWAS/combine/melanoma_dsv_irAEs.csv",row.names=F)

#######################  NSCLC

NSCLC_vsv_resp_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/NSCLC_vsv_resp_adjAbun_res.csv",header=T)
NSCLC_vsv_resp_adjAbun_overall_sig<-subset(NSCLC_vsv_resp_adjAbun_overall,p<0.05,drop=T)
NSCLC_vsv_resp_adjAbun_datasets<-read.csv("08.Microbial_GWAS/datasets/mel_vsv_resp_pvalue.csv",header=T)
NSCLC_vsv_resp<-merge(NSCLC_vsv_resp_adjAbun_overall_sig,NSCLC_vsv_resp_adjAbun_datasets,by.x="X",by.y="X")
NSCLC_vsv_resp_final<-left_join(NSCLC_vsv_resp, vsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_vsv_resp_final,"08.Microbial_GWAS/combine/NSCLC_vsv_resp.csv",row.names=F)

NSCLC_vsv_os_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/NSCLC_vsv_os_adjAbun_res.csv",header=T)
NSCLC_vsv_os_adjAbun_overall_sig<-subset(NSCLC_vsv_os_adjAbun_overall,p<0.05,drop=T)
NSCLC_vsv_os_adjAbun_datasets<-read.csv("08.Microbial_GWAS/datasets/mel_vsv_os_pvalue.csv",header=T)
NSCLC_vsv_os<-merge(NSCLC_vsv_os_adjAbun_overall_sig,NSCLC_vsv_os_adjAbun_datasets,by.x="X",by.y="X")
NSCLC_vsv_os_final<-left_join(NSCLC_vsv_os, vsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_vsv_os_final,"08.Microbial_GWAS/combine/NSCLC_vsv_os.csv",row.names=F)


NSCLC_dsv_resp_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/NSCLC_dsv_resp_adjAbun_res.csv",header=T)
NSCLC_dsv_resp_adjAbun_overall_sig<-subset(NSCLC_dsv_resp_adjAbun_overall,p<0.05,drop=T)
NSCLC_dsv_resp_adjAbun_datasets<-read.csv("08.Microbial_GWAS/datasets/mel_dsv_resp_pvalue.csv",header=T)
NSCLC_dsv_resp<-merge(NSCLC_dsv_resp_adjAbun_overall_sig,NSCLC_dsv_resp_adjAbun_datasets,by.x="X",by.y="X")
NSCLC_dsv_resp_final<-left_join(NSCLC_dsv_resp, dsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_dsv_resp_final,"08.Microbial_GWAS/combine/NSCLC_dsv_resp.csv",row.names=F)

NSCLC_dsv_os_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/NSCLC_dsv_os_adjAbun_res.csv",header=T)
NSCLC_dsv_os_adjAbun_overall_sig<-subset(NSCLC_dsv_os_adjAbun_overall,p<0.05,drop=T)
NSCLC_dsv_os_adjAbun_datasets<-read.csv("08.Microbial_GWAS/datasets/mel_dsv_os_pvalue.csv",header=T)
NSCLC_dsv_os<-merge(NSCLC_dsv_os_adjAbun_overall_sig,NSCLC_dsv_os_adjAbun_datasets,by.x="X",by.y="X")
NSCLC_dsv_os_final<-left_join(NSCLC_dsv_os, dsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_dsv_os_final,"08.Microbial_GWAS/combine/NSCLC_dsv_os.csv",row.names=F)


#######################  RCC

RCC_vsv_resp_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/RCC_vsv_resp_adjAbun_res.csv",header=T)
RCC_vsv_resp_adjAbun_overall_sig<-subset(RCC_vsv_resp_adjAbun_overall,p<0.05,drop=T)
RCC_vsv_resp_final<-left_join(RCC_vsv_resp_adjAbun_overall_sig, vsv_info, by = c("X" = "SV_Name"))
write.csv(RCC_vsv_resp_final,"08.Microbial_GWAS/combine/RCC_vsv_resp.csv",row.names=F)

RCC_dsv_resp_adjAbun_overall<-read.csv("08.Microbial_GWAS/Rdata/RCC_dsv_resp_adjAbun_res.csv",header=T)
RCC_dsv_resp_adjAbun_overall_sig<-subset(RCC_dsv_resp_adjAbun_overall,p<0.05,drop=T)
RCC_dsv_resp_final<-left_join(RCC_dsv_resp_adjAbun_overall_sig, dsv_info, by = c("X" = "SV_Name"))
write.csv(RCC_dsv_resp_final,"08.Microbial_GWAS/combine/RCC_dsv_resp.csv",row.names=F)



#### the species count

dsv_os<-read.csv("08.Microbial_GWAS/combine/dsv_os_adjAbun.sig.anno.csv",header=T)$Taxonomy_Name
vsv_os<-read.csv("08.Microbial_GWAS/combine/vsv_os_adjAbun.sig.anno.csv",header=T)$Taxonomy_Name
os<-c(dsv_os,vsv_os)
length(unique(os))


vsv_resp<-read.csv("08.Microbial_GWAS/combine/vsv_response_adjAbun.sig.anno.csv",header=T)$Taxonomy_Name
dsv_resp<-read.csv("08.Microbial_GWAS/combine/dsv_response_adjAbun.sig.anno.csv",header=T)$Taxonomy_Name
resp<-c(vsv_resp,dsv_resp)
length(unique(resp))



