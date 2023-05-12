### Microbial GWA 
### 2022-11-24
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SNV/pipeline/4_R")

## 1 Preparation
### 1.1 Import
#install.packages("survminer")
rm(list = ls())
#options(scipen = 200)

source("functions.R")
knitr::opts_chunk$set(echo = TRUE)
#library("R.utils")
#library(showtext)
#showtext_auto()

#library("survival")
#library("ggplot2")
#library("survminer")
#library(gridExtra)
#library("grid")
library(reshape2)	  
#library("RColorBrewer")
#library("plyr")
library("metafor")


###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

all_basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

mel_basic<-subset(all_basic,cancer_type=="melanoma",drop=T)
NSCLC_basic<-subset(all_basic,cancer_type=="NSCLC",drop=T)
RCC_basic<-subset(all_basic,cancer_type=="RCC",drop=T)

all_abun<-read.table("01.cleanData/mbio_all/SV_species_abun.tsv", check.names = F)

load("01.cleanData/SV_all/dsgv.RData")
load("01.cleanData/SV_all/vsgv.RData")

vsgv<-vsgv_sub
dsgv<-dsgv_sub

vsgv$"Streptococcus vestibularis F0396:_temp"<-vsgv_sub$"Streptococcus vestibularis F0396:2_3 and 19 segments"

PRJEB22863_basic<-subset(all_basic,dataset=="PRJEB22863",drop=T)
PRJNA397906_basic<-subset(all_basic,dataset=="PRJNA397906",drop=T)
PRJNA541981_basic<-subset(all_basic,dataset=="PRJNA541981",drop=T)
PRJNA751792_basic<-subset(all_basic,dataset=="PRJNA751792",drop=T)
PRJNA762360_basic<-subset(all_basic,dataset=="PRJNA762360",drop=T)
PRJNA770295_basic<-subset(all_basic,dataset=="PRJNA770295",drop=T)
PRJEB43119_basic<-subset(all_basic,dataset=="PRJEB43119",drop=T)

PRJEB22863_NSCLC_basic<-subset(PRJEB22863_basic,cancer_type=="NSCLC",drop=T)
PRJEB22863_RCC_basic<-subset(PRJEB22863_basic,cancer_type=="RCC",drop=T)

PRJEB22863_dsgv<-dsgv[row.names(PRJEB22863_basic),]
PRJEB22863_vsgv<-vsgv[row.names(PRJEB22863_basic),]

PRJEB22863_NSCLC_dsgv<-dsgv[row.names(PRJEB22863_NSCLC_basic),]
PRJEB22863_NSCLC_vsgv<-vsgv[row.names(PRJEB22863_NSCLC_basic),]

PRJEB22863_RCC_dsgv<-dsgv[row.names(PRJEB22863_RCC_basic),]
PRJEB22863_RCC_vsgv<-vsgv[row.names(PRJEB22863_RCC_basic),]

PRJNA397906_dsgv<-dsgv[row.names(PRJNA397906_basic),]
PRJNA397906_vsgv<-vsgv[row.names(PRJNA397906_basic),]
PRJNA541981_dsgv<-dsgv[row.names(PRJNA541981_basic),]
PRJNA541981_vsgv<-vsgv[row.names(PRJNA541981_basic),]
PRJNA751792_dsgv<-dsgv[row.names(PRJNA751792_basic),]
PRJNA751792_vsgv<-vsgv[row.names(PRJNA751792_basic),]
PRJNA762360_dsgv<-dsgv[row.names(PRJNA762360_basic),]
PRJNA762360_vsgv<-vsgv[row.names(PRJNA762360_basic),]
PRJNA770295_dsgv<-dsgv[row.names(PRJNA770295_basic),]
PRJNA770295_vsgv<-vsgv[row.names(PRJNA770295_basic),]
PRJEB43119_dsgv<-dsgv[row.names(PRJEB43119_basic),]
PRJEB43119_vsgv<-vsgv[row.names(PRJEB43119_basic),]

PRJEB22863_abun<-read.table("01.cleanData/mbio_all/PRJEB22863_SV_species_abun.tsv", check.names = F) 
PRJEB22863_NSCLC_abun<-PRJEB22863_abun[row.names(PRJEB22863_NSCLC_basic),]
PRJEB22863_RCC_abun<-PRJEB22863_abun[row.names(PRJEB22863_RCC_basic),]
PRJNA397906_abun<-read.table("01.cleanData/mbio_all/PRJNA397906_SV_species_abun.tsv", check.names = F)  
PRJNA541981_abun<-read.table("01.cleanData/mbio_all/PRJNA541981_SV_species_abun.tsv", check.names = F) 
PRJNA751792_abun<-read.table("01.cleanData/mbio_all/PRJNA751792_SV_species_abun.tsv", check.names = F) 
PRJNA762360_abun<-read.table("01.cleanData/mbio_all/PRJNA762360_SV_species_abun.tsv", check.names = F) 
PRJNA770295_abun<-read.table("01.cleanData/mbio_all/PRJNA770295_SV_species_abun.tsv", check.names = F)
PRJEB43119_abun<-read.table("01.cleanData/mbio_all/PRJEB43119_SV_species_abun.tsv", check.names = F)

load("01.cleanData/SV_all/msv_pc_cum0.6.RData")
load("01.cleanData/SV_all/msv_pc_cum0.6_info.RData")

mel_msv_pc_cum0.6<-msv_pc_cum0.6[row.names(mel_basic),]
NSCLC_msv_pc_cum0.6<-msv_pc_cum0.6[row.names(NSCLC_basic),]
RCC_msv_pc_cum0.6<-msv_pc_cum0.6[row.names(RCC_basic),]

## 2 Associations between SVs and prognosis
### 2.1 Preparation
if (!dir.exists("08.Microbial_GWAS")) {dir.create("08.Microbial_GWAS")}
if (!dir.exists("08.Microbial_GWAS/RData")) {dir.create("08.Microbial_GWAS/RData")}

## Prepare covariate table
# all
all_abun$Others<-1-rowSums(all_abun)
all_abun_clr<-abundances(x=as.data.frame(na.omit(all_abun)), transform="clr") %>%as.data.frame
all_abun_clr <- all_abun[match(rownames(all_abun), rownames(all_abun_clr)),]
rownames(all_abun_clr) <- rownames(all_abun)

mel_abun_clr<-all_abun_clr[row.names(mel_basic),]
NSCLC_abun_clr<-all_abun_clr[row.names(NSCLC_basic),]
RCC_abun_clr<-all_abun_clr[row.names(RCC_basic),]

# PRJEB22863
PRJEB22863_abun$Others<-1-rowSums(PRJEB22863_abun)
PRJEB22863_abun_clr<-abundances(x=as.data.frame(na.omit(PRJEB22863_abun)), transform="clr") %>%as.data.frame
PRJEB22863_abun_clr <- PRJEB22863_abun_clr[match(rownames(PRJEB22863_abun), rownames(PRJEB22863_abun_clr)),]
rownames(PRJEB22863_abun_clr) <- rownames(PRJEB22863_abun)
PRJEB22863_covar<-cbind(PRJEB22863_basic,PRJEB22863_abun_clr)

# PRJEB22863 NSCLC
PRJEB22863_NSCLC_abun$Others<-1-rowSums(PRJEB22863_NSCLC_abun)
PRJEB22863_NSCLC_abun_clr<-abundances(x=as.data.frame(na.omit(PRJEB22863_NSCLC_abun)), transform="clr") %>%as.data.frame
PRJEB22863_NSCLC_abun_clr <- PRJEB22863_NSCLC_abun_clr[match(rownames(PRJEB22863_NSCLC_abun), rownames(PRJEB22863_NSCLC_abun_clr)),]
rownames(PRJEB22863_NSCLC_abun_clr) <- rownames(PRJEB22863_NSCLC_abun)
PRJEB22863_NSCLC_covar<-cbind(PRJEB22863_NSCLC_basic,PRJEB22863_NSCLC_abun_clr)

# PRJEB22863 RCC
PRJEB22863_RCC_abun$Others<-1-rowSums(PRJEB22863_RCC_abun)
PRJEB22863_RCC_abun_clr<-abundances(x=as.data.frame(na.omit(PRJEB22863_RCC_abun)), transform="clr") %>%as.data.frame
PRJEB22863_RCC_abun_clr <- PRJEB22863_RCC_abun_clr[match(rownames(PRJEB22863_RCC_abun), rownames(PRJEB22863_RCC_abun_clr)),]
rownames(PRJEB22863_RCC_abun_clr) <- rownames(PRJEB22863_RCC_abun)
PRJEB22863_RCC_covar<-cbind(PRJEB22863_RCC_basic,PRJEB22863_RCC_abun_clr)


# PRJNA397906
PRJNA397906_abun$Others<-1-rowSums(PRJNA397906_abun)
PRJNA397906_abun_clr<-abundances(x=as.data.frame(na.omit(PRJNA397906_abun)), transform="clr") %>%as.data.frame
PRJNA397906_abun_clr <- PRJNA397906_abun_clr[match(rownames(PRJNA397906_abun), rownames(PRJNA397906_abun_clr)),]
rownames(PRJNA397906_abun_clr) <- rownames(PRJNA397906_abun)
PRJNA397906_covar<-cbind(PRJNA397906_basic,PRJNA397906_abun_clr)

# PRJNA541981
PRJNA541981_abun$Others<-1-rowSums(PRJNA541981_abun)
PRJNA541981_abun_clr<-abundances(x=as.data.frame(na.omit(PRJNA541981_abun)), transform="clr") %>%as.data.frame
PRJNA541981_abun_clr <- PRJNA541981_abun_clr[match(rownames(PRJNA541981_abun), rownames(PRJNA541981_abun_clr)),]
rownames(PRJNA541981_abun_clr) <- rownames(PRJNA541981_abun)
PRJNA541981_covar<-cbind(PRJNA541981_basic,PRJNA541981_abun_clr)

# PRJNA751792
PRJNA751792_abun$Others<-1-rowSums(PRJNA751792_abun)
PRJNA751792_abun_clr<-abundances(x=as.data.frame(na.omit(PRJNA751792_abun)), transform="clr") %>%as.data.frame
PRJNA751792_abun_clr <- PRJNA751792_abun_clr[match(rownames(PRJNA751792_abun), rownames(PRJNA751792_abun_clr)),]
rownames(PRJNA751792_abun_clr) <- rownames(PRJNA751792_abun)
PRJNA751792_covar<-cbind(PRJNA751792_basic,PRJNA751792_abun_clr)

# PRJNA762360
PRJNA762360_abun$Others<-1-rowSums(PRJNA762360_abun)
PRJNA762360_abun_clr<-abundances(x=as.data.frame(na.omit(PRJNA762360_abun)), transform="clr") %>%as.data.frame
PRJNA762360_abun_clr <- PRJNA762360_abun_clr[match(rownames(PRJNA762360_abun), rownames(PRJNA762360_abun_clr)),]
rownames(PRJNA762360_abun_clr) <- rownames(PRJNA762360_abun)
PRJNA762360_covar<-cbind(PRJNA762360_basic,PRJNA762360_abun_clr)

# PRJNA770295
PRJNA770295_abun$Others<-1-rowSums(PRJNA770295_abun)
PRJNA770295_abun_clr<-abundances(x=as.data.frame(na.omit(PRJNA770295_abun)), transform="clr") %>%as.data.frame
PRJNA770295_abun_clr <- PRJNA770295_abun_clr[match(rownames(PRJNA770295_abun), rownames(PRJNA770295_abun_clr)),]
rownames(PRJNA770295_abun_clr) <- rownames(PRJNA770295_abun)
PRJNA770295_covar<-cbind(PRJNA770295_basic,PRJNA770295_abun_clr)

# PRJEB43119
PRJEB43119_abun$Others<-1-rowSums(PRJEB43119_abun)
PRJEB43119_abun_clr<-abundances(x=as.data.frame(na.omit(PRJEB43119_abun)), transform="clr") %>%as.data.frame
PRJEB43119_abun_clr <- PRJEB43119_abun_clr[match(rownames(PRJEB43119_abun), rownames(PRJEB43119_abun_clr)),]
rownames(PRJEB43119_abun_clr) <- rownames(PRJEB43119_abun)
PRJEB43119_covar<-cbind(PRJEB43119_basic,PRJEB43119_abun_clr)

############################################################################################################################
###  for dataset (PRJNA762360)
###  Logistic regression model 1
#Linear model with PRJNA762360 covariates: age, gender, read count and corresponding species relative abundance.
# 1 indicate the organism name column of the info file

PRJNA762360_covar <- c("read_count","age_bins","gender_num")
PRJNA762360_prog<-PRJNA762360_basic[,c("response_code","irAEs","pfs_12_months")]

PRJNA762360_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJNA762360_prog,PRJNA762360_vsgv,PRJNA762360_basic,PRJNA762360_covar,PRJNA762360_abun_clr,info,1,"response_code")
PRJNA762360_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJNA762360_prog,PRJNA762360_dsgv,PRJNA762360_basic,PRJNA762360_covar,PRJNA762360_abun_clr,info,1,"response_code")

write.csv(PRJNA762360_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA762360_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(PRJNA762360_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA762360_dsv_resp_adjAbun_res.csv",row.names=F)

PRJNA762360_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJNA762360_prog,PRJNA762360_vsgv,PRJNA762360_basic,PRJNA762360_covar,PRJNA762360_abun_clr,info,1,"irAEs")
PRJNA762360_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJNA762360_prog,PRJNA762360_dsgv,PRJNA762360_basic,PRJNA762360_covar,PRJNA762360_abun_clr,info,1,"irAEs")

write.csv(PRJNA762360_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA762360_vsv_irAEs_adjAbun_res.csv",row.names=F)
write.csv(PRJNA762360_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA762360_dsv_irAEs_adjAbun_res.csv",row.names=F)

PRJNA762360_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJNA762360_prog,PRJNA762360_vsgv,PRJNA762360_basic,PRJNA762360_covar,PRJNA762360_abun_clr,info,1,"pfs_12_months")
PRJNA762360_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJNA762360_prog,PRJNA762360_dsgv,PRJNA762360_basic,PRJNA762360_covar,PRJNA762360_abun_clr,info,1,"pfs_12_months")

write.csv(PRJNA762360_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA762360_vsv_pfs_12_months_adjAbun_res.csv",row.names=F)
write.csv(PRJNA762360_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA762360_dsv_pfs_12_months_adjAbun_res.csv",row.names=F)


### Cox regression model 1
#Survival model with PRJNA762360_covariates: age, gender, BMI, read count and corresponding species relative abundance.

PRJNA762360_surv<-PRJNA762360_basic[,c("os","os_event")]
PRJNA762360_vsv_os_cox_adjAbun_res<-cox_btw_mats_adjAbun_vsv(PRJNA762360_surv,PRJNA762360_vsgv,PRJNA762360_basic,PRJNA762360_covar,PRJNA762360_abun_clr,info,1,"os","os_event")
PRJNA762360_dsv_os_cox_adjAbun_res<-cox_btw_mats_adjAbun_dsv(PRJNA762360_surv,PRJNA762360_dsgv,PRJNA762360_basic,PRJNA762360_covar,PRJNA762360_abun_clr,info,1,"os","os_event")

write.csv(PRJNA762360_dsv_os_cox_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA762360_dsv_os_adjAbun_res.csv",row.names=F)
write.csv(PRJNA762360_vsv_os_cox_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA762360_vsv_os_adjAbun_res.csv",row.names=F)

############################################################################################################################
###  for dataset (PRJNA751792)
###  Logistic regression model 1
#Linear model with PRJNA751792 covariates: age, gender, read count and corresponding species relative abundance.
# 1 indicate the organism name column of the info file

PRJNA751792_covar <- c("read_count","age_bins","gender_num")
PRJNA751792_prog<-PRJNA751792_basic[,c("response_code","response_code")]

PRJNA751792_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJNA751792_prog,PRJNA751792_vsgv,PRJNA751792_basic,PRJNA751792_covar,PRJNA751792_abun_clr,info,1,"response_code")
PRJNA751792_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJNA751792_prog,PRJNA751792_dsgv,PRJNA751792_basic,PRJNA751792_covar,PRJNA751792_abun_clr,info,1,"response_code")

write.csv(PRJNA751792_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA751792_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(PRJNA751792_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA751792_dsv_resp_adjAbun_res.csv",row.names=F)

### Cox regression model 1
#Survival model with PRJNA751792_covariates: age, gender, BMI, read count and corresponding species relative abundance.

PRJNA751792_surv<-PRJNA751792_basic[,c("os","os_event")]
PRJNA751792_vsv_os_cox_adjAbun_res<-cox_btw_mats_adjAbun_vsv(PRJNA751792_surv,PRJNA751792_vsgv,PRJNA751792_basic,PRJNA751792_covar,PRJNA751792_abun_clr,info,1,"os","os_event")
PRJNA751792_dsv_os_cox_adjAbun_res<-cox_btw_mats_adjAbun_dsv(PRJNA751792_surv,PRJNA751792_dsgv,PRJNA751792_basic,PRJNA751792_covar,PRJNA751792_abun_clr,info,1,"os","os_event")

write.csv(PRJNA751792_dsv_os_cox_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA751792_dsv_os_adjAbun_res.csv",row.names=F)
write.csv(PRJNA751792_vsv_os_cox_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA751792_vsv_os_adjAbun_res.csv",row.names=F)

############################################################################################################################
###  for dataset (PRJNA397906)
###  Logistic regression model 1
#Linear model with PRJNA397906 covariates: age, gender, read count and corresponding species relative abundance.
# 1 indicate the organism name column of the info file

PRJNA397906_covar <- c("read_count","age_bins","gender_num")
PRJNA397906_prog<-PRJNA397906_basic[,c("response_code","irAEs")]

PRJNA397906_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJNA397906_prog,PRJNA397906_vsgv,PRJNA397906_basic,PRJNA397906_covar,PRJNA397906_abun_clr,info,1,"response_code")
PRJNA397906_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJNA397906_prog,PRJNA397906_dsgv,PRJNA397906_basic,PRJNA397906_covar,PRJNA397906_abun_clr,info,1,"response_code")

write.csv(PRJNA397906_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA397906_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(PRJNA397906_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA397906_dsv_resp_adjAbun_res.csv",row.names=F)

############################################################################################################################
###  for dataset (PRJEB22863)
###  Logistic regression model 1
#Linear model with PRJEB22863 covariates: age, gender, read count and corresponding species relative abundance.
# 1 indicate the organism name column of the info file

PRJEB22863_covar <- c("read_count","age_bins","gender_num")
PRJEB22863_prog<-PRJEB22863_basic[,c("response_code","irAEs")]

PRJEB22863_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJEB22863_prog,PRJEB22863_vsgv,PRJEB22863_basic,PRJEB22863_covar,PRJEB22863_abun_clr,info,1,"response_code")
PRJEB22863_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJEB22863_prog,PRJEB22863_dsgv,PRJEB22863_basic,PRJEB22863_covar,PRJEB22863_abun_clr,info,1,"response_code")

write.csv(PRJEB22863_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJEB22863_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(PRJEB22863_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJEB22863_dsv_resp_adjAbun_res.csv",row.names=F)

############################################################################################################################
###  for dataset (PRJEB22863) NSCLC
###  Logistic regression model 1
#Linear model with PRJEB22863 covariates: age, gender, read count and corresponding species relative abundance.
# 1 indicate the organism name column of the info file

PRJEB22863_covar <- c("read_count","age_bins","gender_num")
PRJEB22863_prog<-PRJEB22863_NSCLC_basic[,c("response_code","irAEs")]

PRJEB22863_NSCLC_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJEB22863_prog,PRJEB22863_NSCLC_vsgv,PRJEB22863_NSCLC_basic,PRJEB22863_covar,PRJEB22863_NSCLC_abun_clr,info,1,"response_code")
PRJEB22863_NSCLC_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJEB22863_prog,PRJEB22863_NSCLC_dsgv,PRJEB22863_NSCLC_basic,PRJEB22863_covar,PRJEB22863_NSCLC_abun_clr,info,1,"response_code")

write.csv(PRJEB22863_NSCLC_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJEB22863_NSCLC_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(PRJEB22863_NSCLC_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJEB22863_NSCLC_dsv_resp_adjAbun_res.csv",row.names=F)

############################################################################################################################
###  for dataset (PRJEB22863) RCC
###  Logistic regression model 1
#Linear model with PRJEB22863 covariates: age, gender, read count and corresponding species relative abundance.
# 1 indicate the organism name column of the info file

PRJEB22863_covar <- c("read_count","age_bins","gender_num")
PRJEB22863_prog<-PRJEB22863_RCC_basic[,c("response_code","irAEs")]

PRJEB22863_RCC_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJEB22863_prog,PRJEB22863_RCC_vsgv,PRJEB22863_RCC_basic,PRJEB22863_covar,PRJEB22863_RCC_abun_clr,info,1,"response_code")
PRJEB22863_RCC_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJEB22863_prog,PRJEB22863_RCC_dsgv,PRJEB22863_RCC_basic,PRJEB22863_covar,PRJEB22863_RCC_abun_clr,info,1,"response_code")

write.csv(PRJEB22863_RCC_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJEB22863_RCC_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(PRJEB22863_RCC_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJEB22863_RCC_dsv_resp_adjAbun_res.csv",row.names=F)

############################################################################################################################
###  for dataset (PRJNA770295)
###  Logistic regression model 1
#Linear model with PRJNA770295 covariates: age, gender, read count and corresponding species relative abundance.
# 1 indicate the organism name column of the info file

PRJNA770295_covar <- c("read_count","age_bins","gender_num")
PRJNA770295_prog<-PRJNA770295_basic[,c("response_code","pfs_12_months")]

PRJNA770295_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJNA770295_prog,PRJNA770295_vsgv,PRJNA770295_basic,PRJNA770295_covar,PRJNA770295_abun_clr,info,1,"response_code")
PRJNA770295_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJNA770295_prog,PRJNA770295_dsgv,PRJNA770295_basic,PRJNA770295_covar,PRJNA770295_abun_clr,info,1,"response_code")

write.csv(PRJNA770295_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA770295_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(PRJNA770295_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA770295_dsv_resp_adjAbun_res.csv",row.names=F)

PRJNA770295_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJNA770295_prog,PRJNA770295_vsgv,PRJNA770295_basic,PRJNA770295_covar,PRJNA770295_abun_clr,info,1,"pfs_12_months")
PRJNA770295_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJNA770295_prog,PRJNA770295_dsgv,PRJNA770295_basic,PRJNA770295_covar,PRJNA770295_abun_clr,info,1,"pfs_12_months")

write.csv(PRJNA770295_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA770295_vsv_pfs_12_months_adjAbun_res.csv",row.names=F)
write.csv(PRJNA770295_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA770295_dsv_pfs_12_months_adjAbun_res.csv",row.names=F)


############################################################################################################################
###  for dataset (PRJEB43119)
###  Logistic regression model 1
#Linear model with PRJEB43119 covariates: age, gender, read count and corresponding species relative abundance.
# 1 indicate the organism name column of the info file

PRJEB43119_covar <- c("read_count","age_bins","gender_num")
PRJEB43119_prog<-PRJEB43119_basic[,c("response_code","pfs_12_months")]

PRJEB43119_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJEB43119_prog,PRJEB43119_vsgv,PRJEB43119_basic,PRJEB43119_covar,PRJEB43119_abun_clr,info,1,"response_code")
PRJEB43119_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJEB43119_prog,PRJEB43119_dsgv,PRJEB43119_basic,PRJEB43119_covar,PRJEB43119_abun_clr,info,1,"response_code")

write.csv(PRJEB43119_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJEB43119_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(PRJEB43119_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJEB43119_dsv_resp_adjAbun_res.csv",row.names=F)

PRJEB43119_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJEB43119_prog,PRJEB43119_vsgv,PRJEB43119_basic,PRJEB43119_covar,PRJEB43119_abun_clr,info,1,"pfs_12_months")
PRJEB43119_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJEB43119_prog,PRJEB43119_dsgv,PRJEB43119_basic,PRJEB43119_covar,PRJEB43119_abun_clr,info,1,"pfs_12_months")

write.csv(PRJEB43119_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJEB43119_vsv_pfs_12_months_adjAbun_res.csv",row.names=F)
write.csv(PRJEB43119_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJEB43119_dsv_pfs_12_months_adjAbun_res.csv",row.names=F)


############################################################################################################################
###  for dataset (PRJNA541981)
###  Logistic regression model 1
#Linear model with PRJNA541981 covariates: age, gender, read count and corresponding species relative abundance.
# 1 indicate the organism name column of the info file

PRJNA541981_covar <- c("read_count","age_bins","gender_num")
PRJNA541981_prog<-PRJNA541981_basic[,c("pfs_12_months","response_code")]

PRJNA541981_vsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(PRJNA541981_prog,PRJNA541981_vsgv,PRJNA541981_basic,PRJNA541981_covar,PRJNA541981_abun_clr,info,1,"pfs_12_months")
PRJNA541981_dsv_lg_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(PRJNA541981_prog,PRJNA541981_dsgv,PRJNA541981_basic,PRJNA541981_covar,PRJNA541981_abun_clr,info,1,"pfs_12_months")

write.csv(PRJNA541981_vsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA541981_vsv_pfs_12_months_adjAbun_res.csv",row.names=F)
write.csv(PRJNA541981_dsv_lg_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA541981_dsv_pfs_12_months_adjAbun_res.csv",row.names=F)


### Cox regression model 1
#Survival model with PRJNA541981_covariates: age, gender, BMI, read count and corresponding species relative abundance.

PRJNA541981_surv<-PRJNA541981_basic[,c("os","os_event")]
PRJNA541981_vsv_os_cox_adjAbun_res<-cox_btw_mats_adjAbun_vsv(PRJNA541981_surv,PRJNA541981_vsgv,PRJNA541981_basic,PRJNA541981_covar,PRJNA541981_abun_clr,info,1,"os","os_event")
PRJNA541981_dsv_os_cox_adjAbun_res<-cox_btw_mats_adjAbun_dsv(PRJNA541981_surv,PRJNA541981_dsgv,PRJNA541981_basic,PRJNA541981_covar,PRJNA541981_abun_clr,info,1,"os","os_event")

write.csv(PRJNA541981_dsv_os_cox_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA541981_dsv_os_adjAbun_res.csv",row.names=F)
write.csv(PRJNA541981_vsv_os_cox_adjAbun_res, file = "08.Microbial_GWAS/datasets/PRJNA541981_vsv_os_adjAbun_res.csv",row.names=F)



###############################################################################
### dsv meta-analysis
PRJEB22863_dsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJEB22863_dsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB22863_dsv_resp_res.table)[-1]<-c(paste("PRJEB22863.",colnames(PRJEB22863_dsv_resp_res.table)[2:4],sep = ""))

PRJEB22863_NSCLC_dsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJEB22863_NSCLC_dsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB22863_NSCLC_dsv_resp_res.table)[-1]<-c(paste("PRJEB22863.",colnames(PRJEB22863_NSCLC_dsv_resp_res.table)[2:4],sep = ""))

PRJEB22863_RCC_dsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJEB22863_RCC_dsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB22863_RCC_dsv_resp_res.table)[-1]<-c(paste("PRJEB22863.",colnames(PRJEB22863_RCC_dsv_resp_res.table)[2:4],sep = ""))

PRJNA397906_dsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA397906_dsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA397906_dsv_resp_res.table)[-1]<-c(paste("PRJNA397906.",colnames(PRJNA397906_dsv_resp_res.table)[2:4],sep = ""))

PRJNA751792_dsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA751792_dsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA751792_dsv_resp_res.table)[-1]<-c(paste("PRJNA751792.",colnames(PRJNA751792_dsv_resp_res.table)[2:4],sep = ""))

PRJNA762360_dsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA762360_dsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA762360_dsv_resp_res.table)[-1]<-c(paste("PRJNA762360.",colnames(PRJNA762360_dsv_resp_res.table)[2:4],sep = ""))

PRJNA770295_dsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA770295_dsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA770295_dsv_resp_res.table)[-1]<-c(paste("PRJNA770295.",colnames(PRJNA770295_dsv_resp_res.table)[2:4],sep = ""))

PRJEB43119_dsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJEB43119_dsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB43119_dsv_resp_res.table)[-1]<-c(paste("PRJEB43119.",colnames(PRJEB43119_dsv_resp_res.table)[2:4],sep = ""))

cbind_dsv_resp_res.table<-merge(PRJEB22863_dsv_resp_res.table,PRJNA397906_dsv_resp_res.table,by.x="X",by.y="X")
cbind_dsv_resp_res.table<-merge(cbind_dsv_resp_res.table,PRJNA751792_dsv_resp_res.table,by.x="X",by.y="X")
cbind_dsv_resp_res.table<-merge(cbind_dsv_resp_res.table,PRJNA762360_dsv_resp_res.table,by.x="X",by.y="X")
cbind_dsv_resp_res.table<-merge(cbind_dsv_resp_res.table,PRJNA770295_dsv_resp_res.table,by.x="X",by.y="X")
cbind_dsv_resp_res.table<-merge(cbind_dsv_resp_res.table,PRJEB43119_dsv_resp_res.table,by.x="X",by.y="X")

write.csv(cbind_dsv_resp_res.table,"08.Microbial_GWAS/datasets/all_dsv_resp_pvalue.csv",row.names=F)

### mel
mel_dsv_resp_res.table<-merge(PRJNA397906_dsv_resp_res.table,PRJNA762360_dsv_resp_res.table,by.x="X",by.y="X")
mel_dsv_resp_res.table<-merge(mel_dsv_resp_res.table,PRJNA770295_dsv_resp_res.table,by.x="X",by.y="X")
mel_dsv_resp_res.table<-merge(mel_dsv_resp_res.table,PRJEB43119_dsv_resp_res.table,by.x="X",by.y="X")
write.csv(mel_dsv_resp_res.table,"08.Microbial_GWAS/datasets/mel_dsv_resp_pvalue.csv",row.names=F)


### nsclc
nsclc_dsv_resp_res.table<-merge(PRJNA751792_dsv_resp_res.table,PRJEB22863_NSCLC_dsv_resp_res.table,by.x="X",by.y="X")
write.csv(nsclc_dsv_resp_res.table,"08.Microbial_GWAS/datasets/nsclc_dsv_resp_pvalue.csv",row.names=F)

#####################################################
## pfs 12 months 

PRJNA762360_dsv_pfs_12_months_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA762360_dsv_pfs_12_months_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA762360_dsv_pfs_12_months_res.table)[-1]<-c(paste("PRJNA762360.",colnames(PRJNA762360_dsv_pfs_12_months_res.table)[2:4],sep = ""))

PRJNA770295_dsv_pfs_12_months_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA770295_dsv_pfs_12_months_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA770295_dsv_pfs_12_months_res.table)[-1]<-c(paste("PRJNA770295.",colnames(PRJNA770295_dsv_pfs_12_months_res.table)[2:4],sep = ""))

PRJEB43119_dsv_pfs_12_months_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJEB43119_dsv_pfs_12_months_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB43119_dsv_pfs_12_months_res.table)[-1]<-c(paste("PRJEB43119.",colnames(PRJEB43119_dsv_pfs_12_months_res.table)[2:4],sep = ""))

PRJEB541981_dsv_pfs_12_months_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA541981_dsv_pfs_12_months_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB541981_dsv_pfs_12_months_res.table)[-1]<-c(paste("PRJEB541981.",colnames(PRJEB541981_dsv_pfs_12_months_res.table)[2:4],sep = ""))

### mel
mel_dsv_pfs_12_months_res.table<-merge(PRJNA762360_dsv_pfs_12_months_res.table,PRJNA770295_dsv_pfs_12_months_res.table,by.x="X",by.y="X")
mel_dsv_pfs_12_months_res.table<-merge(mel_dsv_pfs_12_months_res.table,PRJEB43119_dsv_pfs_12_months_res.table,by.x="X",by.y="X")
mel_dsv_pfs_12_months_res.table<-merge(mel_dsv_pfs_12_months_res.table,PRJEB541981_dsv_pfs_12_months_res.table,by.x="X",by.y="X")
write.csv(mel_dsv_pfs_12_months_res.table,"08.Microbial_GWAS/datasets/mel_dsv_pfs_12_months_pvalue.csv",row.names=F)



###############################################################################################
### os 

PRJNA541981_dsv_os_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA541981_dsv_os_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA541981_dsv_os_res.table)[-1]<-c(paste("PRJNA541981.",colnames(PRJNA541981_dsv_os_res.table)[2:4],sep = ""))

PRJNA762360_dsv_os_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA762360_dsv_os_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA762360_dsv_os_res.table)[-1]<-c(paste("PRJNA762360.",colnames(PRJNA762360_dsv_os_res.table)[2:4],sep = ""))

PRJNA751792_dsv_os_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA751792_dsv_os_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA751792_dsv_os_res.table)[-1]<-c(paste("PRJNA751792.",colnames(PRJNA751792_dsv_os_res.table)[2:4],sep = ""))

cbind_dsv_os_res.table<-merge(PRJNA541981_dsv_os_res.table,PRJNA762360_dsv_os_res.table,by.x="X",by.y="X")
cbind_dsv_os_res.table<-merge(cbind_dsv_os_res.table,PRJNA751792_dsv_os_res.table,by.x="X",by.y="X")
write.csv(cbind_dsv_os_res.table,"08.Microbial_GWAS/datasets/all_dsv_os_pvalue.csv",row.names=F)

mel_dsv_os_res.table<-merge(PRJNA541981_dsv_os_res.table,PRJNA762360_dsv_os_res.table,by.x="X",by.y="X")
write.csv(mel_dsv_os_res.table,"08.Microbial_GWAS/datasets/mel_dsv_os_pvalue.csv",row.names=F)


###############################################################################
### vsv meta-analysis
PRJEB22863_vsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJEB22863_vsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB22863_vsv_resp_res.table)[-1]<-c(paste("PRJEB22863.",colnames(PRJEB22863_vsv_resp_res.table)[2:4],sep = ""))

PRJEB22863_NSCLC_vsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJEB22863_NSCLC_vsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB22863_NSCLC_vsv_resp_res.table)[-1]<-c(paste("PRJEB22863.",colnames(PRJEB22863_NSCLC_vsv_resp_res.table)[2:4],sep = ""))

PRJEB22863_RCC_vsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJEB22863_RCC_vsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB22863_RCC_vsv_resp_res.table)[-1]<-c(paste("PRJEB22863.",colnames(PRJEB22863_RCC_vsv_resp_res.table)[2:4],sep = ""))

PRJNA397906_vsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA397906_vsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA397906_vsv_resp_res.table)[-1]<-c(paste("PRJNA397906.",colnames(PRJNA397906_vsv_resp_res.table)[2:4],sep = ""))

PRJNA751792_vsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA751792_vsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA751792_vsv_resp_res.table)[-1]<-c(paste("PRJNA751792.",colnames(PRJNA751792_vsv_resp_res.table)[2:4],sep = ""))

PRJNA762360_vsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA762360_vsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA762360_vsv_resp_res.table)[-1]<-c(paste("PRJNA762360.",colnames(PRJNA762360_vsv_resp_res.table)[2:4],sep = ""))

PRJNA770295_vsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA770295_vsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA770295_vsv_resp_res.table)[-1]<-c(paste("PRJNA770295.",colnames(PRJNA770295_vsv_resp_res.table)[2:4],sep = ""))

PRJEB43119_vsv_resp_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJEB43119_vsv_resp_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB43119_vsv_resp_res.table)[-1]<-c(paste("PRJEB43119.",colnames(PRJEB43119_vsv_resp_res.table)[2:4],sep = ""))

cbind_vsv_resp_res.table<-merge(PRJEB22863_vsv_resp_res.table,PRJNA397906_vsv_resp_res.table,by.x="X",by.y="X")
cbind_vsv_resp_res.table<-merge(cbind_vsv_resp_res.table,PRJNA751792_vsv_resp_res.table,by.x="X",by.y="X")
cbind_vsv_resp_res.table<-merge(cbind_vsv_resp_res.table,PRJNA762360_vsv_resp_res.table,by.x="X",by.y="X")
cbind_vsv_resp_res.table<-merge(cbind_vsv_resp_res.table,PRJNA770295_vsv_resp_res.table,by.x="X",by.y="X")
cbind_vsv_resp_res.table<-merge(cbind_vsv_resp_res.table,PRJEB43119_vsv_resp_res.table,by.x="X",by.y="X")
write.csv(cbind_vsv_resp_res.table,"08.Microbial_GWAS/datasets/all_vsv_resp_pvalue.csv")

### mel
mel_vsv_resp_res.table<-merge(PRJNA397906_vsv_resp_res.table,PRJNA762360_vsv_resp_res.table,by.x="X",by.y="X")
mel_vsv_resp_res.table<-merge(mel_vsv_resp_res.table,PRJNA770295_vsv_resp_res.table,by.x="X",by.y="X")
mel_vsv_resp_res.table<-merge(mel_vsv_resp_res.table,PRJEB43119_vsv_resp_res.table,by.x="X",by.y="X")
write.csv(mel_vsv_resp_res.table,"08.Microbial_GWAS/datasets/mel_vsv_resp_pvalue.csv",row.names=F)

### nsclc
nsclc_vsv_resp_res.table<-merge(PRJNA751792_vsv_resp_res.table,PRJEB22863_NSCLC_vsv_resp_res.table,by.x="X",by.y="X")
write.csv(nsclc_vsv_resp_res.table,"08.Microbial_GWAS/datasets/nsclc_vsv_resp_pvalue.csv",row.names=F)

#########################  PFS 12 months 

PRJNA762360_vsv_pfs_12_months_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA762360_vsv_pfs_12_months_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA762360_vsv_pfs_12_months_res.table)[-1]<-c(paste("PRJNA762360.",colnames(PRJNA762360_vsv_pfs_12_months_res.table)[2:4],sep = ""))

PRJNA770295_vsv_pfs_12_months_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA770295_vsv_pfs_12_months_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA770295_vsv_pfs_12_months_res.table)[-1]<-c(paste("PRJNA770295.",colnames(PRJNA770295_vsv_pfs_12_months_res.table)[2:4],sep = ""))

PRJEB43119_vsv_pfs_12_months_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJEB43119_vsv_pfs_12_months_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB43119_vsv_pfs_12_months_res.table)[-1]<-c(paste("PRJEB43119.",colnames(PRJEB43119_vsv_pfs_12_months_res.table)[2:4],sep = ""))

PRJEB541981_vsv_pfs_12_months_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA541981_vsv_pfs_12_months_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJEB541981_vsv_pfs_12_months_res.table)[-1]<-c(paste("PRJEB541981.",colnames(PRJEB541981_vsv_pfs_12_months_res.table)[2:4],sep = ""))

### mel
mel_vsv_pfs_12_months_res.table<-merge(PRJNA762360_vsv_pfs_12_months_res.table,PRJNA770295_vsv_pfs_12_months_res.table,by.x="X",by.y="X")
mel_vsv_pfs_12_months_res.table<-merge(mel_vsv_pfs_12_months_res.table,PRJEB43119_vsv_pfs_12_months_res.table,by.x="X",by.y="X")
mel_vsv_pfs_12_months_res.table<-merge(mel_vsv_pfs_12_months_res.table,PRJEB541981_vsv_pfs_12_months_res.table,by.x="X",by.y="X")
write.csv(mel_vsv_pfs_12_months_res.table,"08.Microbial_GWAS/datasets/mel_vsv_pfs_12_months_pvalue.csv",row.names=F)


###############################################################################################
### os 

PRJNA541981_vsv_os_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA541981_vsv_os_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA541981_vsv_os_res.table)[-1]<-c(paste("PRJNA541981.",colnames(PRJNA541981_vsv_os_res.table)[2:4],sep = ""))

PRJNA762360_vsv_os_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA762360_vsv_os_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA762360_vsv_os_res.table)[-1]<-c(paste("PRJNA762360.",colnames(PRJNA762360_vsv_os_res.table)[2:4],sep = ""))

PRJNA751792_vsv_os_res.table<-read.csv("08.Microbial_GWAS/datasets/PRJNA751792_vsv_os_adjAbun_res.csv",header=T)[,c(1,2,4,7)]
colnames(PRJNA751792_vsv_os_res.table)[-1]<-c(paste("PRJNA751792.",colnames(PRJNA751792_vsv_os_res.table)[2:4],sep = ""))

cbind_vsv_os_res.table<-merge(PRJNA541981_vsv_os_res.table,PRJNA762360_vsv_os_res.table,by.x="X",by.y="X")
cbind_vsv_os_res.table<-merge(cbind_vsv_os_res.table,PRJNA751792_vsv_os_res.table,by.x="X",by.y="X")

write.csv(cbind_vsv_os_res.table,"08.Microbial_GWAS/datasets/all_vsv_os_pvalue.csv",row.names=F)

mel_vsv_os_res.table<-merge(PRJNA541981_vsv_os_res.table,PRJNA762360_vsv_os_res.table,by.x="X",by.y="X")
write.csv(mel_vsv_os_res.table,"08.Microbial_GWAS/datasets/mel_vsv_os_pvalue.csv",row.names=F)




