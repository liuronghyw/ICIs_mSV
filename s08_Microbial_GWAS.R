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

mel_dsgv<-dsgv[row.names(mel_basic),]
NSCLC_dsgv<-dsgv[row.names(NSCLC_basic),]
RCC_dsgv<-dsgv[row.names(RCC_basic),]

mel_vsgv<-vsgv[row.names(mel_basic),]
NSCLC_vsgv<-vsgv[row.names(NSCLC_basic),]
RCC_vsgv<-vsgv[row.names(RCC_basic),]

#load("01.cleanData/SV_all/msv_pc_cum0.6.RData")
#load("01.cleanData/SV_all/msv_pc_cum0.6_info.RData")

#mel_msv_pc_cum0.6<-msv_pc_cum0.6[row.names(mel_basic),]
#NSCLC_msv_pc_cum0.6<-msv_pc_cum0.6[row.names(NSCLC_basic),]
#RCC_msv_pc_cum0.6<-msv_pc_cum0.6[row.names(RCC_basic),]

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

mel_abun<-all_abun[row.names(mel_basic),]
mel_abun_clr<-abundances(x=as.data.frame(na.omit(mel_abun)), transform="clr") %>%as.data.frame
mel_abun_clr <- all_abun[match(rownames(mel_abun), rownames(mel_abun_clr)),]
rownames(mel_abun_clr) <- rownames(mel_abun)

NSCLC_abun<-all_abun[row.names(NSCLC_basic),]
NSCLC_abun_clr<-abundances(x=as.data.frame(na.omit(NSCLC_abun)), transform="clr") %>%as.data.frame
NSCLC_abun_clr <- all_abun[match(rownames(NSCLC_abun), rownames(NSCLC_abun_clr)),]
rownames(NSCLC_abun_clr) <- rownames(NSCLC_abun)

RCC_abun<-all_abun[row.names(RCC_basic),]
RCC_abun_clr<-abundances(x=as.data.frame(na.omit(RCC_abun)), transform="clr") %>%as.data.frame
RCC_abun_clr <- all_abun[match(rownames(RCC_abun), rownames(RCC_abun_clr)),]
rownames(RCC_abun_clr) <- rownames(RCC_abun)


### 2.2 Logistic regression model 1
# Linear model with covariates: age, gender, BMI, read count and corresponding species relative abundance.
# 9 indicate the organism name column of the info file

######################## for differnet cancer type 
## mel_basic
## NSCLC_basic
## RCC_basic

### 2.2 Logistic regression model 1
#Linear model with covariates: age, gender, BMI, read count and corresponding species relative abundance.

covar <- c("read_count","cohort","age_bins","gender_num")
mel_prog<-mel_basic[,c("response_code","irAEs","pfs_6_months","pfs_12_months")]

mel_vsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(mel_prog,mel_vsgv,mel_basic,covar,mel_abun_clr,info,1,"response_code")
mel_dsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(mel_prog,mel_dsgv,mel_basic,covar,mel_abun_clr,info,1,"response_code")

write.csv(mel_vsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(mel_dsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_dsv_resp_adjAbun_res.csv",row.names=F)

mel_vsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(mel_prog,mel_vsgv,mel_basic,covar,mel_abun_clr,info,1,"pfs_6_months")
mel_dsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(mel_prog,mel_dsgv,mel_basic,covar,mel_abun_clr,info,1,"pfs_6_months")

write.csv(mel_vsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_vsv_pfs_6_months_adjAbun_res.csv",row.names=F)
write.csv(mel_dsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_dsv_pfs_6_months_adjAbun_res.csv",row.names=F)

mel_vsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(mel_prog,mel_vsgv,mel_basic,covar,mel_abun_clr,info,1,"pfs_12_months")
mel_dsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(mel_prog,mel_dsgv,mel_basic,covar,mel_abun_clr,info,1,"pfs_12_months")

write.csv(mel_vsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_vsv_pfs_12_months_adjAbun_res.csv",row.names=F)
write.csv(mel_dsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_dsv_pfs_12_months_adjAbun_res.csv",row.names=F)

mel_vsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(mel_prog,mel_vsgv,mel_basic,covar,mel_abun_clr,info,1,"irAEs")
mel_dsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(mel_prog,mel_dsgv,mel_basic,covar,mel_abun_clr,info,1,"irAEs")

write.csv(mel_vsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_vsv_irAEs_adjAbun_res.csv",row.names=F)
write.csv(mel_dsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_dsv_irAEs_adjAbun_res.csv",row.names=F)



### 2.4 Cox regression model 1
#Survival model with covariates: age, gender, BMI, read count and corresponding species relative abundance.

mel_surv<-mel_basic[,c("os","os_event","pfs","pfs_event")]
mel_covar<-c("read_count","cohort","age_bins","gender_num")

mel_vsv_os_cox_adjAbun_res<-cox_btw_mats_adjAbun_vsv(mel_surv,mel_vsgv,mel_basic,mel_covar,mel_abun_clr,info,1,"os","os_event")
mel_vsv_pfs_cox_adjAbun_res<-cox_btw_mats_adjAbun_vsv(mel_surv,mel_vsgv,mel_basic,mel_covar,mel_abun_clr,info,1,"pfs","pfs_event")

write.csv(mel_vsv_os_cox_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_vsv_os_adjAbun_res.csv",row.names=F)
write.csv(mel_vsv_pfs_cox_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_vsv_pfs_adjAbun_res.csv",row.names=F)

mel_dsv_os_cox_adjAbun_res<-cox_btw_mats_adjAbun_dsv(mel_surv,mel_dsgv,mel_basic,covar,mel_abun_clr,info,1,"os","os_event")
mel_dsv_pfs_cox_adjAbun_res<-cox_btw_mats_adjAbun_dsv(mel_surv,mel_dsgv,mel_basic,covar,mel_abun_clr,info,1,"pfs","pfs_event")

write.csv(mel_dsv_os_cox_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_dsv_os_adjAbun_res.csv",row.names=F)
write.csv(mel_dsv_pfs_cox_adjAbun_res, file = "08.Microbial_GWAS/RData/melanoma_dsv_pfs_adjAbun_res.csv",row.names=F)


#################################################################################################################################################
###  NSCLC

### 2.2 Logistic regression model 1
#Linear model with covariates: age, gender, BMI, read count and corresponding species relative abundance.

NSCLC_covar <- c("read_count","cohort","age_bins","gender_num")
NSCLC_prog<-NSCLC_basic[,c("response_code","pfs_6_months")]

NSCLC_vsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(NSCLC_prog,NSCLC_vsgv,NSCLC_basic,NSCLC_covar,NSCLC_abun_clr,info,1,"response_code")
NSCLC_dsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(NSCLC_prog,NSCLC_dsgv,NSCLC_basic,NSCLC_covar,NSCLC_abun_clr,info,1,"response_code")

write.csv(NSCLC_vsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/NSCLC_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(NSCLC_dsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/NSCLC_dsv_resp_adjAbun_res.csv",row.names=F)

NSCLC_vsv_pfs_6_months_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(NSCLC_prog,NSCLC_vsgv,NSCLC_basic,NSCLC_covar,NSCLC_abun_clr,info,1,"response_code")
NSCLC_dsv_pfs_6_months_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(NSCLC_prog,NSCLC_dsgv,NSCLC_basic,NSCLC_covar,NSCLC_abun_clr,info,1,"response_code")

write.csv(NSCLC_vsv_pfs_6_months_adjAbun_res, file = "08.Microbial_GWAS/RData/NSCLC_vsv_pfs_6_months_adjAbun_res.csv",row.names=F)
write.csv(NSCLC_dsv_pfs_6_months_adjAbun_res, file = "08.Microbial_GWAS/RData/NSCLC_dsv_pfs_6_months_adjAbun_res.csv",row.names=F)


### 2.4 Cox regression model 1
#Survival model with covariates: age, gender, BMI, read count and corresponding species relative abundance.

NSCLC_surv<-NSCLC_basic[,c("os","os_event","pfs","pfs_event")]
NSCLC_covar<-c("read_count","cohort","age_bins","gender_num")

NSCLC_vsv_os_cox_adjAbun_res<-cox_btw_mats_adjAbun_vsv(NSCLC_surv,NSCLC_vsgv,NSCLC_basic,NSCLC_covar,NSCLC_abun_clr,info,1,"os","os_event")
write.csv(NSCLC_vsv_os_cox_adjAbun_res, file = "08.Microbial_GWAS/RData/NSCLC_vsv_os_adjAbun_res.csv",row.names=F)

NSCLC_dsv_os_cox_adjAbun_res<-cox_btw_mats_adjAbun_dsv(NSCLC_surv,NSCLC_dsgv,NSCLC_basic,covar,NSCLC_abun_clr,info,1,"os","os_event")
write.csv(NSCLC_dsv_os_cox_adjAbun_res, file = "08.Microbial_GWAS/RData/NSCLC_dsv_os_adjAbun_res.csv",row.names=F)


################################################################################################################################################
### RCC

### 2.2 Logistic regression model 1
#Linear model with covariates: age, gender, BMI, read count and corresponding species relative abundance.

RCC_covar <- c("read_count","age_bins","gender_num")
RCC_prog<-RCC_basic[,c("response_code","pfs_6_months")]

RCC_vsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(RCC_prog,RCC_vsgv,RCC_basic,RCC_covar,RCC_abun_clr,info,1,"response_code")
RCC_dsv_resp_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(RCC_prog,RCC_dsgv,RCC_basic,RCC_covar,RCC_abun_clr,info,1,"response_code")

write.csv(RCC_vsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/RCC_vsv_resp_adjAbun_res.csv",row.names=F)
write.csv(RCC_dsv_resp_adjAbun_res, file = "08.Microbial_GWAS/RData/RCC_dsv_resp_adjAbun_res.csv",row.names=F)

RCC_vsv_pfs_6_months_adjAbun_res<-logistic_btw_mats_adjAbun_vsv(RCC_prog,RCC_vsgv,RCC_basic,RCC_covar,RCC_abun_clr,info,1,"pfs_6_months")
RCC_dsv_pfs_6_months_adjAbun_res<-logistic_btw_mats_adjAbun_dsv(RCC_prog,RCC_dsgv,RCC_basic,RCC_covar,RCC_abun_clr,info,1,"pfs_6_months")

write.csv(RCC_vsv_pfs_6_months_adjAbun_res, file = "08.Microbial_GWAS/RData/RCC_vsv_pfs_6_months_adjAbun_res.csv",row.names=F)
write.csv(RCC_dsv_pfs_6_months_adjAbun_res, file = "08.Microbial_GWAS/RData/RCC_dsv_pfs_6_months_adjAbun_res.csv",row.names=F)


###################################################################################################
### 3 Results of Model 1 
### 3.1 Clean results
#### 3.1.1 vSV associations

## vsv
## merge result tables
## mel (response, os ,pfs)
## nsclc(response,os)
## rcc (response)

melanoma_vsv_resp_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_vsv_resp_adjAbun_res.csv",header=T)
melanoma_vsv_resp_adjAbun_res.sig.edge<-subset(melanoma_vsv_resp_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)
#melanoma_vsv_resp_adjAbun_res.sig.edge<-subset(melanoma_vsv_resp_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
melanoma_vsv_resp_adjAbun_res.sig.anno.edge<-left_join(melanoma_vsv_resp_adjAbun_res.sig.edge, vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_resp_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/melanoma_vsv_response_adjAbun.sig.anno.csv",row.names=F)

melanoma_vsv_os_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_vsv_os_adjAbun_res.csv",header=T)
melanoma_vsv_os_adjAbun_res.sig.edge<-subset(melanoma_vsv_os_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)
#melanoma_vsv_os_adjAbun_res.sig.edge<-subset(melanoma_vsv_os_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
melanoma_vsv_os_adjAbun_res.sig.anno.edge<-left_join(melanoma_vsv_os_adjAbun_res.sig.edge, vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_os_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/melanoma_vsv_os_adjAbun.sig.anno.csv",row.names=F)

melanoma_vsv_pfs_12_months_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_vsv_pfs_12_months_adjAbun_res.csv",header=T)
melanoma_vsv_pfs_12_months_adjAbun_res.sig.edge<-subset(melanoma_vsv_pfs_12_months_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.1,drop=T)
#melanoma_vsv_pfs_12_months_adjAbun_res.sig.edge<-subset(melanoma_vsv_pfs_12_months_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
melanoma_vsv_pfs_12_months_adjAbun_res.sig.anno.edge<-left_join(melanoma_vsv_pfs_12_months_adjAbun_res.sig.edge, vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_pfs_12_months_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/melanoma_vsv_pfs_12_months_adjAbun.sig.anno.csv",row.names=F)

melanoma_vsv_irAEs_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_vsv_irAEs_adjAbun_res.csv",header=T)
melanoma_vsv_irAEs_adjAbun_res.sig.edge<-subset(melanoma_vsv_irAEs_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)
#melanoma_vsv_irAEs_adjAbun_res.sig.edge<-subset(melanoma_vsv_irAEs_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
melanoma_vsv_irAEs_adjAbun_res.sig.anno.edge<-left_join(melanoma_vsv_irAEs_adjAbun_res.sig.edge, vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_irAEs_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/melanoma_vsv_irAEs_adjAbun.sig.anno.csv",row.names=F)


NSCLC_vsv_resp_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/NSCLC_vsv_resp_adjAbun_res.csv",header=T)
NSCLC_vsv_resp_adjAbun_res.sig.edge<-subset(NSCLC_vsv_resp_adjAbun_res.edge,as.numeric(p) < 0.1 & as.numeric(fdr.p) < 0.1,drop=T)
#NSCLC_vsv_resp_adjAbun_res.sig.edge<-subset(NSCLC_vsv_resp_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
NSCLC_vsv_resp_adjAbun_res.sig.anno.edge<-left_join(NSCLC_vsv_resp_adjAbun_res.sig.edge, vsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_vsv_resp_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/NSCLC_vsv_response_adjAbun.sig.anno.csv",row.names=F)

NSCLC_vsv_os_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/NSCLC_vsv_os_adjAbun_res.csv",header=T)
NSCLC_vsv_os_adjAbun_res.sig.edge<-subset(NSCLC_vsv_os_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)
#NSCLC_vsv_os_adjAbun_res.sig.edge<-subset(NSCLC_vsv_os_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
NSCLC_vsv_os_adjAbun_res.sig.anno.edge<-left_join(NSCLC_vsv_os_adjAbun_res.sig.edge, vsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_vsv_os_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/NSCLC_vsv_os_adjAbun.sig.anno.csv",row.names=F)

RCC_vsv_resp_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/RCC_vsv_resp_adjAbun_res.csv",header=T)
RCC_vsv_resp_adjAbun_res.sig.edge<-subset(RCC_vsv_resp_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)
#RCC_vsv_resp_adjAbun_res.sig.edge<-subset(RCC_vsv_resp_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
RCC_vsv_resp_adjAbun_res.sig.anno.edge<-left_join(RCC_vsv_resp_adjAbun_res.sig.edge, vsv_info, by = c("X" = "SV_Name"))
write.csv(RCC_vsv_resp_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/RCC_vsv_response_adjAbun.sig.anno.csv",row.names=F)


#### 3.1.2 dSV associations
## load result data

melanoma_dsv_resp_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_dsv_resp_adjAbun_res.csv",header=T)
melanoma_dsv_resp_adjAbun_res.sig.edge<-subset(melanoma_dsv_resp_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)
#melanoma_dsv_resp_adjAbun_res.sig.edge<-subset(melanoma_dsv_resp_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
melanoma_dsv_resp_adjAbun_res.sig.anno.edge<-left_join(melanoma_dsv_resp_adjAbun_res.sig.edge, dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_resp_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/melanoma_dsv_response_adjAbun.sig.anno.csv",row.names=F)

melanoma_dsv_os_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_dsv_os_adjAbun_res.csv",header=T)
melanoma_dsv_os_adjAbun_res.sig.edge<-subset(melanoma_dsv_os_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)
#melanoma_dsv_os_adjAbun_res.sig.edge<-subset(melanoma_dsv_os_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
melanoma_dsv_os_adjAbun_res.sig.anno.edge<-left_join(melanoma_dsv_os_adjAbun_res.sig.edge, dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_os_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/melanoma_dsv_os_adjAbun.sig.anno.csv",row.names=F)

melanoma_dsv_pfs_12_months_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_dsv_pfs_12_months_adjAbun_res.csv",header=T)
melanoma_dsv_pfs_12_months_adjAbun_res.sig.edge<-subset(melanoma_dsv_pfs_12_months_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.1,drop=T)
#melanoma_dsv_pfs_12_months_adjAbun_res.sig.edge<-subset(melanoma_dsv_pfs_12_months_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
melanoma_dsv_pfs_12_months_adjAbun_res.sig.anno.edge<-left_join(melanoma_dsv_pfs_12_months_adjAbun_res.sig.edge, dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_pfs_12_months_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/melanoma_dsv_pfs_12_months_adjAbun.sig.anno.csv",row.names=F)

melanoma_dsv_irAEs_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_dsv_irAEs_adjAbun_res.csv",header=T)
melanoma_dsv_irAEs_adjAbun_res.sig.edge<-subset(melanoma_dsv_irAEs_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)
#melanoma_dsv_irAEs_adjAbun_res.sig.edge<-subset(melanoma_dsv_irAEs_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
melanoma_dsv_irAEs_adjAbun_res.sig.anno.edge<-left_join(melanoma_dsv_irAEs_adjAbun_res.sig.edge, dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_irAEs_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/melanoma_dsv_irAEs_adjAbun.sig.anno.csv",row.names=F)

NSCLC_dsv_resp_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/NSCLC_dsv_resp_adjAbun_res.csv",header=T)
NSCLC_dsv_resp_adjAbun_res.sig.edge<-subset(NSCLC_dsv_resp_adjAbun_res.edge,as.numeric(p) < 0.1 & as.numeric(fdr.p) < 0.1,drop=T)
#NSCLC_dsv_resp_adjAbun_res.sig.edge<-subset(NSCLC_dsv_resp_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
NSCLC_dsv_resp_adjAbun_res.sig.anno.edge<-left_join(NSCLC_dsv_resp_adjAbun_res.sig.edge, dsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_dsv_resp_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/NSCLC_dsv_response_adjAbun.sig.anno.csv",row.names=F)

NSCLC_dsv_os_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/NSCLC_dsv_os_adjAbun_res.csv",header=T)
NSCLC_dsv_os_adjAbun_res.sig.edge<-subset(NSCLC_dsv_os_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)
#NSCLC_dsv_os_adjAbun_res.sig.edge<-subset(NSCLC_dsv_os_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
NSCLC_dsv_os_adjAbun_res.sig.anno.edge<-left_join(NSCLC_dsv_os_adjAbun_res.sig.edge, dsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_dsv_os_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/NSCLC_dsv_os_adjAbun.sig.anno.csv",row.names=F)

RCC_dsv_resp_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/RCC_dsv_resp_adjAbun_res.csv",header=T)
RCC_dsv_resp_adjAbun_res.sig.edge<-subset(RCC_dsv_resp_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)
#RCC_dsv_resp_adjAbun_res.sig.edge<-subset(RCC_dsv_resp_adjAbun_res.edge,as.numeric(p) < 0.01,drop=T)
RCC_dsv_resp_adjAbun_res.sig.anno.edge<-left_join(RCC_dsv_resp_adjAbun_res.sig.edge, dsv_info, by = c("X" = "SV_Name"))
write.csv(RCC_dsv_resp_adjAbun_res.sig.anno.edge,"08.Microbial_GWAS/RCC_dsv_response_adjAbun.sig.anno.csv",row.names=F)


#####################################
#### 3.1.3 Volcano plot
# vsv 

all_vsv_resp_adjAbun_res<-read.csv("08.Microbial_GWAS/RData/all_vsv_resp_adjAbun_res.csv",header=T)
all_vsv_resp_adjAbun_res$MetaSigAssoc<-rep('No', nrow(all_vsv_resp_adjAbun_res))
all_vsv_resp_adjAbun_res$MetaSigAssoc[all_vsv_resp_adjAbun_res$fdr.p < 0.05]<-'Yes'

p_vsv_volcano<-ggplot(all_vsv_resp_adjAbun_res,aes(log(OR), -log10(as.numeric(p)), color=MetaSigAssoc))+
  geom_point(alpha = 0.5,size = 1)+
  xlab('OR')+
  ylab('-log10(P)')+
  scale_color_manual(name   = NULL,
                     breaks = c("Yes", "No"),
                     labels = c("Associated", "Not associated"),
                     values = c("#ff4040","#4f94cd"))+
  scale_shape_discrete(name = NULL)+
  theme_bw()+
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))


################################################## dsv
all_dsv_resp_adjAbun_res<-read.csv("08.Microbial_GWAS/RData/all_dsv_resp_adjAbun_res.csv",header=T)
#all_dsv_resp_adjAbun_res<-subset(all_dsv_resp_adjAbun_res,OR<10,drop=T)

all_dsv_resp_adjAbun_res$beta<-log(all_dsv_resp_adjAbun_res$OR)

all_dsv_resp_adjAbun_res$MetaSigAssoc<-rep('No', nrow(all_dsv_resp_adjAbun_res))
all_dsv_resp_adjAbun_res$MetaSigAssoc[all_dsv_resp_adjAbun_res$fdr.p < 0.05]<-'Yes'

p_dsv_volcano<-ggplot(all_dsv_resp_adjAbun_res,aes(beta, -log10(as.numeric(p)), color=MetaSigAssoc))+
  geom_point(alpha = 0.5,size = 1)+
  xlab('Beta coefficient')+
  ylab('-log10(P)')+
  scale_color_manual(name   = NULL,
                     breaks = c("Yes", "No"),
                     labels = c("Associated", "Not associated"),
                     values = c("#ff4040","#4f94cd"))+
  scale_shape_discrete(name = NULL)+
  theme_bw()+
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA))


## plot
p_title_volcano <- ggdraw() + 
    draw_label(
      ' ', # Volcano plot
      fontface = 'bold', x = 0, hjust = 0,size =20) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))

p_sv_ba_volcano<-plot_grid(
    p_title_volcano, 
    plot_grid(p_vsv_volcano,p_dsv_volcano,
              rel_widths = c(1, 1),align = 'hv',
              labels = c("vSVs", "dSVs"),
              ncol = 2,label_size= 20,vjust = 0),
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 2)
  )

#pdf("pics/sv_adjAbun_ba_volcano.pdf", height = 4, width = 7)
#print(p_sv_ba_volcano)
#dev.off()

tiff(file = "pics/FS5_sv_adjAbun_ba_volcano.tiff", width =2500, height =1200, res =300) 
print(p_sv_ba_volcano)
dev.off()



########################################################
#####  3.3 Visualization  （将不同的表型的结果汇总起来）
####   3.3.1 Combine vSV and dSV associations

#melanoma_vsv_resp_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/melanoma_vsv_response_adjAbun.sig.anno.csv",header=T)
melanoma_vsv_os_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/melanoma_vsv_os_adjAbun.sig.anno.csv",header=T)
melanoma_vsv_pfs_12_months_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/melanoma_vsv_pfs_12_months_adjAbun.sig.anno.csv",header=T)
melanoma_vsv_irAEs_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/melanoma_vsv_irAEs_adjAbun.sig.anno.csv",header=T)

melanoma_dsv_resp_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/melanoma_dsv_response_adjAbun.sig.anno.csv",header=T)
melanoma_dsv_os_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/melanoma_dsv_os_adjAbun.sig.anno.csv",header=T)
melanoma_dsv_pfs_12_months_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/melanoma_dsv_pfs_12_months_adjAbun.sig.anno.csv",header=T)
#melanoma_dsv_irAEs_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/melanoma_dsv_irAEs_adjAbun.sig.anno.csv",header=T)

#melanoma_vsv_resp_adjAbun_res.sig.anno.edge$Pheno<-"Melanoma\nResponse"
melanoma_vsv_os_adjAbun_res.sig.anno.edge$Pheno<-"Melanoma\nOS"
melanoma_vsv_pfs_12_months_adjAbun_res.sig.anno.edge$Pheno<-"Melanoma\nPFS >=12 months"
melanoma_vsv_irAEs_adjAbun_res.sig.anno.edge$Pheno<-"Melanoma\nirAEs"

melanoma_dsv_resp_adjAbun_res.sig.anno.edge$Pheno<-"Melanoma\nResponse"
melanoma_dsv_os_adjAbun_res.sig.anno.edge$Pheno<-"Melanoma\nOS"
melanoma_dsv_pfs_12_months_adjAbun_res.sig.anno.edge$Pheno<-"Melanoma\nPFS >=12 months"

melanoma_vsv_resp_adjAbun_sub<-melanoma_vsv_resp_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]
melanoma_vsv_os_adjAbun_sub<-melanoma_vsv_os_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]
melanoma_vsv_pfs_12_months_adjAbun_sub<-melanoma_vsv_pfs_12_months_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]
melanoma_vsv_irAEs_adjAbun_sub<-melanoma_vsv_irAEs_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]

melanoma_dsv_resp_adjAbun_sub<-melanoma_dsv_resp_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]
melanoma_dsv_os_adjAbun_sub<-melanoma_dsv_os_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]
melanoma_dsv_pfs_12_months_adjAbun_sub<-melanoma_dsv_pfs_12_months_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]


NSCLC_vsv_resp_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/NSCLC_vsv_response_adjAbun.sig.anno.csv",header=T)
NSCLC_vsv_os_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/NSCLC_vsv_os_adjAbun.sig.anno.csv",header=T)

NSCLC_dsv_resp_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/NSCLC_dsv_response_adjAbun.sig.anno.csv",header=T)
NSCLC_dsv_os_adjAbun_res.sig.anno.edge<-read.csv("08.Microbial_GWAS/NSCLC_dsv_os_adjAbun.sig.anno.csv",header=T)

NSCLC_vsv_resp_adjAbun_res.sig.anno.edge$Pheno<-"NSCLC\nResponse"
NSCLC_vsv_os_adjAbun_res.sig.anno.edge$Pheno<-"NSCLC\nOS"

NSCLC_dsv_resp_adjAbun_res.sig.anno.edge$Pheno<-"NSCLC\nResponse"
NSCLC_dsv_os_adjAbun_res.sig.anno.edge$Pheno<-"NSCLC\nOS"

NSCLC_vsv_resp_adjAbun_sub<-NSCLC_vsv_resp_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]
NSCLC_vsv_os_adjAbun_sub<-NSCLC_vsv_os_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]

NSCLC_dsv_resp_adjAbun_sub<-NSCLC_dsv_resp_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]
NSCLC_dsv_resp_adjAbun_sub<-subset(NSCLC_dsv_resp_adjAbun_sub,is.na(p)==F,drop=T)
NSCLC_dsv_os_adjAbun_sub<-NSCLC_dsv_os_adjAbun_res.sig.anno.edge[,c("Pheno","X","p","fdr.p")]
NSCLC_dsv_os_adjAbun_sub<-subset(NSCLC_dsv_os_adjAbun_sub,is.na(p)==F,drop=T)

sv_prog_adjAbun_res.sig.edge<-rbind(melanoma_vsv_os_adjAbun_sub,melanoma_vsv_pfs_12_months_adjAbun_sub,
melanoma_dsv_resp_adjAbun_sub,melanoma_dsv_os_adjAbun_sub,melanoma_dsv_pfs_12_months_adjAbun_sub,
NSCLC_vsv_os_adjAbun_sub,NSCLC_dsv_resp_adjAbun_sub,NSCLC_dsv_os_adjAbun_sub,melanoma_vsv_irAEs_adjAbun_sub)

sv_prog_adjAbun_res.sig.edge<-subset(sv_prog_adjAbun_res.sig.edge,fdr.p<0.05,drop=T)
sv_prog_adjAbun_res.sig.edge<-subset(sv_prog_adjAbun_res.sig.edge,fdr.p>0.000001,drop=T)


#NSCLC_vsv_resp_adjAbun_sub,melanoma_vsv_resp_adjAbun_sub,

colnames(sv_prog_adjAbun_res.sig.edge)[2]<-"SV"

sv_prog_adjAbun.sig.pheno  <- sv_prog_adjAbun_res.sig.edge$Pheno %>%
  as.character(.) %>%
  .[!duplicated(.)]

sv_prog_adjAbun.sig.sv <- sv_prog_adjAbun_res.sig.edge$SV %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

## circos plot
sv_prog_count<-NULL

for (phe in sv_prog_adjAbun.sig.pheno) {
  #pheno <-"C4"
  sv_prog_df<-sv_prog_adjAbun_res.sig.edge[sv_prog_adjAbun_res.sig.edge$Pheno==phe,]$SV %>%
    str_replace_all(.,"\\:\\d+_\\d+.*", "") %>%
    table %>%
    as.data.frame
  colnames(sv_prog_df)<-c("Species", "Count")
  sv_prog_df<-data.frame(Pheno = rep(phe,nrow(sv_prog_df)), sv_prog_df)
  sv_prog_count<-rbind(sv_prog_count,sv_prog_df)
}

sv_prog_count$Species<-info$Short_name[match(sv_prog_count$Species, info$organism)]
sv_prog_count <- sv_prog_count[order(sv_prog_count$Count),]

sv_prog_count_species_order<-sv_prog_count %>% group_by(Species) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]
sv_prog_count_pheno_order<-sv_prog_count %>% group_by(Pheno) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]

sv_prog_count<-sv_prog_count[order(match(sv_prog_count$Pheno, sv_prog_count_pheno_order$Pheno),decreasing = F),]

sv_prog_count_species_order_str<-as.character(sv_prog_count_species_order$Species)
sv_prog_count_pheno_order_str<-as.character(sv_prog_count_pheno_order$Pheno)


#pdf("pics/F6_SV_pheno_adjAbun.circos.pdf", width = 26, height = 26)
tiff(file = "pics/F6_SV_pheno_adjAbun.circos.tiff", width =8200, height =7500, res =300) 

circos.clear()
circos.par(start.degree = 180, "clock.wise" = T) # ,xaxis.clock.wise = F

chordDiagram(sv_prog_count,annotationTrack = "grid",
             #grid.col =  c(wes_palette("Darjeeling1", length(sv_prog_count_pheno_order_str), type = "continuous"),
             #              rep('grey',length(sv_prog_count_species_order_str))),
             grid.col =  c(c("#0072B5FF","#7876B1FF","#20854EFF","#E18727FF","#EE4C97FF","#6F99ADFF"),
                           rep('grey',length(sv_prog_count_species_order_str))),
             order = c(rev(sv_prog_count_pheno_order_str),
                       sv_prog_count_species_order_str),
             big.gap = 6,
             preAllocateTracks = list(track.margin = c(0, uh(45, "mm")), 
             #                         #track.height = max(strwidth(unlist(dimnames(sv_prog_count)))))
                                       track.height=0.1)
             )

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,cex = 2.5,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

dev.off()

#mycolor6
hist(mtcars$mpg,
     breaks = 15,
     col = "#6F99ADFF",
)

#"#6F99ADFF","#E18727FF","#BC3C29FF","#7876B1FF","#0072B5FF",#20854EFF


write.csv(sv_prog_count,"08.Microbial_GWAS/sv_prog_count.csv",row.names=F)



unique(sv_prog_count$Species)

sum(sv_prog_count$Count)




### 3.4 Examples
### 3.4.1 Heatmap of certain species

melanoma_vsv_resp_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_vsv_resp_adjAbun_res.csv",header=T)
melanoma_vsv_resp_adjAbun_res.anno.edge<-left_join(melanoma_vsv_resp_adjAbun_res.edge, vsv_info, by = c("X" = "SV_Name"))

melanoma_vsv_os_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_vsv_os_adjAbun_res.csv",header=T)
melanoma_vsv_os_adjAbun_res.anno.edge<-left_join(melanoma_vsv_os_adjAbun_res.edge, vsv_info, by = c("X" = "SV_Name"))

melanoma_vsv_pfs_12_months_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_vsv_pfs_12_months_adjAbun_res.csv",header=T)
melanoma_vsv_pfs_12_months_adjAbun_res.anno.edge<-left_join(melanoma_vsv_pfs_12_months_adjAbun_res.edge, vsv_info, by = c("X" = "SV_Name"))

melanoma_vsv_resp_adjAbun_sig_spe<-subset(melanoma_vsv_resp_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)$X
melanoma_vsv_os_adjAbun_sig_spe<-subset(melanoma_vsv_os_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)$X
melanoma_vsv_pfs_12_months_adjAbun_sig_spe<-subset(melanoma_vsv_pfs_12_months_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)$X

NSCLC_vsv_resp_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/NSCLC_vsv_resp_adjAbun_res.csv",header=T)
NSCLC_vsv_resp_adjAbun_res.anno.edge<-left_join(NSCLC_vsv_resp_adjAbun_res.edge, vsv_info, by = c("X" = "SV_Name"))

NSCLC_vsv_os_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/NSCLC_vsv_os_adjAbun_res.csv",header=T)
NSCLC_vsv_os_adjAbun_res.anno.edge<-left_join(NSCLC_vsv_os_adjAbun_res.edge, vsv_info, by = c("X" = "SV_Name"))

NSCLC_vsv_resp_adjAbun_sig_spe<-subset(NSCLC_vsv_resp_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)$X
NSCLC_vsv_os_adjAbun_sig_spe<-subset(NSCLC_vsv_os_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)$X

vsv_select<-data.frame(c(melanoma_vsv_os_adjAbun_sig_spe,melanoma_vsv_resp_adjAbun_sig_spe,
melanoma_vsv_pfs_12_months_adjAbun_sig_spe,NSCLC_vsv_resp_adjAbun_sig_spe,NSCLC_vsv_os_adjAbun_sig_spe))
colnames(vsv_select)<-"X"

melanoma_dsv_resp_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_dsv_resp_adjAbun_res.csv",header=T)
melanoma_dsv_resp_adjAbun_res.anno.edge<-left_join(melanoma_dsv_resp_adjAbun_res.edge, dsv_info, by = c("X" = "SV_Name"))

melanoma_dsv_os_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_dsv_os_adjAbun_res.csv",header=T)
melanoma_dsv_os_adjAbun_res.anno.edge<-left_join(melanoma_dsv_os_adjAbun_res.edge, dsv_info, by = c("X" = "SV_Name"))

melanoma_dsv_pfs_12_months_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/melanoma_dsv_pfs_12_months_adjAbun_res.csv",header=T)
melanoma_dsv_pfs_12_months_adjAbun_res.anno.edge<-left_join(melanoma_dsv_pfs_12_months_adjAbun_res.edge, dsv_info, by = c("X" = "SV_Name"))

melanoma_dsv_resp_adjAbun_sig_spe<-subset(melanoma_dsv_resp_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)$X
melanoma_dsv_os_adjAbun_sig_spe<-subset(melanoma_dsv_os_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)$X
melanoma_dsv_pfs_12_months_adjAbun_sig_spe<-subset(melanoma_dsv_pfs_12_months_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)$X

dsv_select<-data.frame(c(melanoma_dsv_os_adjAbun_sig_spe,melanoma_dsv_resp_adjAbun_sig_spe))
colnames(dsv_select)<-"X"

NSCLC_dsv_resp_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/NSCLC_dsv_resp_adjAbun_res.csv",header=T)
NSCLC_dsv_resp_adjAbun_res.anno.edge<-left_join(NSCLC_dsv_resp_adjAbun_res.edge, dsv_info, by = c("X" = "SV_Name"))

NSCLC_dsv_os_adjAbun_res.edge<-read.csv("08.Microbial_GWAS/RData/NSCLC_dsv_os_adjAbun_res.csv",header=T)
NSCLC_dsv_os_adjAbun_res.anno.edge<-left_join(NSCLC_dsv_os_adjAbun_res.edge, dsv_info, by = c("X" = "SV_Name"))

NSCLC_dsv_resp_adjAbun_sig_spe<-subset(NSCLC_dsv_resp_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)$X
NSCLC_dsv_os_adjAbun_sig_spe<-subset(NSCLC_dsv_os_adjAbun_res.edge,as.numeric(p) < 0.05 & as.numeric(fdr.p) < 0.05,drop=T)$X

dsv_select<-data.frame(c(melanoma_dsv_os_adjAbun_sig_spe,melanoma_dsv_resp_adjAbun_sig_spe,
melanoma_dsv_pfs_12_months_adjAbun_sig_spe,NSCLC_dsv_resp_adjAbun_sig_spe,NSCLC_dsv_os_adjAbun_sig_spe))
colnames(dsv_select)<-"X"

melanoma_vsv_resp_adjAbun_select<-merge(vsv_select,melanoma_vsv_resp_adjAbun_res.anno.edge,by.x="X",by.y="X")
melanoma_vsv_os_adjAbun_select<-merge(vsv_select,melanoma_vsv_os_adjAbun_res.anno.edge,by.x="X",by.y="X")
melanoma_vsv_pfs_12_months_adjAbun_select<-merge(vsv_select,melanoma_vsv_pfs_12_months_adjAbun_res.anno.edge,by.x="X",by.y="X")

melanoma_dsv_resp_adjAbun_select<-merge(dsv_select,melanoma_dsv_resp_adjAbun_res.anno.edge,by.x="X",by.y="X")
melanoma_dsv_os_adjAbun_select<-merge(dsv_select,melanoma_dsv_os_adjAbun_res.anno.edge,by.x="X",by.y="X")
melanoma_dsv_pfs_12_months_adjAbun_select<-merge(dsv_select,melanoma_dsv_pfs_12_months_adjAbun_res.anno.edge,by.x="X",by.y="X")

melanoma_vsv_resp_adjAbun_select$Pheno<-"Melanoma.Response"
melanoma_vsv_os_adjAbun_select$Pheno<-"Melanoma.OS"
melanoma_vsv_pfs_12_months_adjAbun_select$Pheno<-"Melanoma.PFS"

melanoma_dsv_resp_adjAbun_select$Pheno<-"Melanoma.Response"
melanoma_dsv_os_adjAbun_select$Pheno<-"Melanoma.OS"
melanoma_dsv_pfs_12_months_adjAbun_select$Pheno<-"Melanoma.PFS"

melanoma_vsv_resp_adjAbun_sub<-melanoma_vsv_resp_adjAbun_select[,c("Pheno","X","p","fdr.p","bonferroni.p","OR")]
melanoma_vsv_os_adjAbun_sub<-melanoma_vsv_os_adjAbun_select[,c("Pheno","X","p","fdr.p","bonferroni.p","HR")]
melanoma_vsv_pfs_12_months_adjAbun_sub<-melanoma_vsv_pfs_12_months_adjAbun_select[,c("Pheno","X","p","fdr.p","bonferroni.p","OR")]

melanoma_dsv_resp_adjAbun_sub<-melanoma_dsv_resp_adjAbun_select[,c("Pheno","X","p","fdr.p","bonferroni.p","OR")]
melanoma_dsv_os_adjAbun_sub<-melanoma_dsv_os_adjAbun_select[,c("Pheno","X","p","fdr.p","bonferroni.p","HR")]
melanoma_dsv_pfs_12_months_adjAbun_sub<-melanoma_dsv_pfs_12_months_adjAbun_select[,c("Pheno","X","p","fdr.p","bonferroni.p","OR")]

melanoma_vsv_resp_adjAbun_sub$SV<-c(paste("vsv_",melanoma_vsv_resp_adjAbun_sub$X,sep = ""))
melanoma_vsv_os_adjAbun_sub$SV<-c(paste("vsv_",melanoma_vsv_os_adjAbun_sub$X,sep = ""))
melanoma_vsv_pfs_12_months_adjAbun_sub$SV<-c(paste("vsv_",melanoma_vsv_pfs_12_months_adjAbun_sub$X,sep = ""))

melanoma_vsv_resp_adjAbun_sub<-rename(melanoma_vsv_resp_adjAbun_sub, beta=OR)
melanoma_vsv_os_adjAbun_sub<-rename(melanoma_vsv_os_adjAbun_sub, beta=HR)
melanoma_vsv_pfs_12_months_adjAbun_sub<-rename(melanoma_vsv_pfs_12_months_adjAbun_sub, beta=OR)

melanoma_vsv_resp_adjAbun_sub$beta<-log(melanoma_vsv_resp_adjAbun_sub$beta)
melanoma_vsv_os_adjAbun_sub$beta<-log(melanoma_vsv_os_adjAbun_sub$beta)
melanoma_vsv_pfs_12_months_adjAbun_sub$beta<-log(melanoma_vsv_pfs_12_months_adjAbun_sub$beta)

melanoma_dsv_resp_adjAbun_sub$SV<-c(paste("dsv_",melanoma_dsv_resp_adjAbun_sub$X,sep = ""))
melanoma_dsv_os_adjAbun_sub$SV<-c(paste("dsv_",melanoma_dsv_os_adjAbun_sub$X,sep = ""))
melanoma_dsv_pfs_12_months_adjAbun_sub$SV<-c(paste("dsv_",melanoma_dsv_pfs_12_months_adjAbun_sub$X,sep = ""))

melanoma_dsv_resp_adjAbun_sub<-rename(melanoma_dsv_resp_adjAbun_sub, beta=OR)
melanoma_dsv_os_adjAbun_sub<-rename(melanoma_dsv_os_adjAbun_sub, beta=HR)
melanoma_dsv_pfs_12_months_adjAbun_sub<-rename(melanoma_dsv_pfs_12_months_adjAbun_sub, beta=OR)

melanoma_dsv_resp_adjAbun_sub$beta<-log(melanoma_dsv_resp_adjAbun_sub$beta)
melanoma_dsv_os_adjAbun_sub$beta<-log(melanoma_dsv_os_adjAbun_sub$beta)
melanoma_dsv_pfs_12_months_adjAbun_sub$beta<-log(melanoma_dsv_pfs_12_months_adjAbun_sub$beta)


NSCLC_vsv_resp_adjAbun_select<-merge(vsv_select,NSCLC_vsv_resp_adjAbun_res.anno.edge,by.x="X",by.y="X")
NSCLC_vsv_os_adjAbun_select<-merge(vsv_select,NSCLC_vsv_os_adjAbun_res.anno.edge,by.x="X",by.y="X")

NSCLC_dsv_resp_adjAbun_select<-merge(dsv_select,NSCLC_dsv_resp_adjAbun_res.anno.edge,by.x="X",by.y="X")
NSCLC_dsv_os_adjAbun_select<-merge(dsv_select,NSCLC_dsv_os_adjAbun_res.anno.edge,by.x="X",by.y="X")

NSCLC_vsv_resp_adjAbun_select$Pheno<-"NSCLC.Response"
NSCLC_vsv_os_adjAbun_select$Pheno<-"NSCLC.OS"

NSCLC_dsv_resp_adjAbun_select$Pheno<-"NSCLC.Response"
NSCLC_dsv_os_adjAbun_select$Pheno<-"NSCLC.OS"

NSCLC_vsv_resp_adjAbun_sub<-NSCLC_vsv_resp_adjAbun_select[,c("Pheno","X","p","fdr.p","bonferroni.p","OR")]
NSCLC_vsv_os_adjAbun_sub<-NSCLC_vsv_os_adjAbun_select[,c("Pheno","X","p","fdr.p","bonferroni.p","HR")]

NSCLC_dsv_resp_adjAbun_sub<-NSCLC_dsv_resp_adjAbun_select[,c("Pheno","X","p","fdr.p","bonferroni.p","OR")]
NSCLC_dsv_os_adjAbun_sub<-NSCLC_dsv_os_adjAbun_select[,c("Pheno","X","p","fdr.p","bonferroni.p","HR")]

NSCLC_vsv_resp_adjAbun_sub$SV<-c(paste("vsv_",NSCLC_vsv_resp_adjAbun_sub$X,sep = ""))
NSCLC_vsv_os_adjAbun_sub$SV<-c(paste("vsv_",NSCLC_vsv_os_adjAbun_sub$X,sep = ""))

NSCLC_vsv_resp_adjAbun_sub<-rename(NSCLC_vsv_resp_adjAbun_sub, beta=OR)
NSCLC_vsv_os_adjAbun_sub<-rename(NSCLC_vsv_os_adjAbun_sub, beta=HR)

NSCLC_vsv_resp_adjAbun_sub$beta<-log(NSCLC_vsv_resp_adjAbun_sub$beta)
NSCLC_vsv_os_adjAbun_sub$beta<-log(NSCLC_vsv_os_adjAbun_sub$beta)

NSCLC_dsv_resp_adjAbun_sub$SV<-c(paste("dsv_",NSCLC_dsv_resp_adjAbun_sub$X,sep = ""))
NSCLC_dsv_os_adjAbun_sub$SV<-c(paste("dsv_",NSCLC_dsv_os_adjAbun_sub$X,sep = ""))

NSCLC_dsv_resp_adjAbun_sub<-rename(NSCLC_dsv_resp_adjAbun_sub, beta=OR)
NSCLC_dsv_os_adjAbun_sub<-rename(NSCLC_dsv_os_adjAbun_sub, beta=HR)

NSCLC_dsv_resp_adjAbun_sub$beta<-log(NSCLC_dsv_resp_adjAbun_sub$beta)
NSCLC_dsv_os_adjAbun_sub$beta<-log(NSCLC_dsv_os_adjAbun_sub$beta)

sv_prog_adjAbun_res.sig.edge<-rbind(melanoma_vsv_resp_adjAbun_sub,melanoma_vsv_os_adjAbun_sub,
melanoma_vsv_pfs_12_months_adjAbun_sub,melanoma_dsv_resp_adjAbun_sub,
melanoma_dsv_os_adjAbun_sub,melanoma_dsv_pfs_12_months_adjAbun_sub,
NSCLC_vsv_resp_adjAbun_sub,NSCLC_vsv_os_adjAbun_sub,
NSCLC_dsv_resp_adjAbun_sub,NSCLC_dsv_os_adjAbun_sub)

sv_prog_adjAbun_res.sig.edge<-as.data.frame(sv_prog_adjAbun_res.sig.edge)

sv_prog_adjAbun_res.sig.edge$Sign<-NA
sv_prog_adjAbun_res.sig.edge$Sign[sv_prog_adjAbun_res.sig.edge$p < 0.05 &
  !is.na(sv_prog_adjAbun_res.sig.edge$p)] <- "☆" 

sv_prog_adjAbun_res.sig.edge$Sign[sv_prog_adjAbun_res.sig.edge$fdr.p < 0.05 &
  !is.na(sv_prog_adjAbun_res.sig.edge$fdr.p)] <- "★" 

#sv_prog_adjAbun_res.sig.edge$Sign[sv_prog_adjAbun_res.sig.edge$bonferroni.p < 0.05 &
#  !is.na(sv_prog_adjAbun_res.sig.edge$bonferroni.p)] <- "★"

sv_prog_adjAbun_res.sig.edge<-unique(sv_prog_adjAbun_res.sig.edge)

sv_prog_adjAbun.r<- sv_prog_adjAbun_res.sig.edge[,c("Pheno", "SV", "beta")] %>% spread("Pheno", "beta") %>% data.frame(row.names = "SV")
sv_prog_adjAbun.sign<- sv_prog_adjAbun_res.sig.edge[,c("Pheno", "SV", "Sign")] %>% spread("Pheno", "Sign") %>% data.frame(row.names = "SV")

sv_prog_adjAbun_res.sig.edge<-subset(sv_prog_adjAbun_res.sig.edge,is.na(p)==F,drop=T)
# write.csv(sv_prog_adjAbun_res.sig.edge,"test.csv")

##########################################################################################################
###  Akkermansia muciniphila   ## ermansia
## Akkermansia muciniphila ATCC BAA-835:

plot.pheno  <- sv_prog_adjAbun_res.sig.edge$Pheno[grep("ermansia", sv_prog_adjAbun_res.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_adjAbun_res.sig.edge$SV[grep("ermansia", sv_prog_adjAbun_res.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("ermansia", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")
#spe.sv_adjAbun.r<-spe.sv_adjAbun.r[,-1]

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("ermansia", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")
#spe.sv_adjAbun.sign[,-1]

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Akkermansia muciniphila ATCC BAA-835:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")


i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-2.0, 0, length.out=ceiling(i/2) + 1), 
              seq(2.0/i, 2.0, length.out=floor(i/2)))

sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

#pdf("08.Microbial_GWAS/F6_shahii_sv_lm_adjAbun.heatmap.pdf", width = 6, height = 10)

tiff(file = "pics/F6_Akk_sv_lm_adjAbun.heatmap.tiff", width =2800, height =1600, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#0072B5FF","#F5F5F5","#BC3C29FF"))(i), 
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("Melanoma\nResponse","Melanoma\nPFS >=12 months","NSCLC\nResponse","NSCLC\nOS"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "grey50",notecex =3.2,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
)

dev.off()


##########################################################################################################
###  C.comes

plot.pheno  <- sv_prog_adjAbun_res.sig.edge$Pheno[grep("comes", sv_prog_adjAbun_res.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_adjAbun_res.sig.edge$SV[grep("comes", sv_prog_adjAbun_res.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("comes", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("comes", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Coprococcus comes ATCC 27758:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-1.5, 0, length.out=ceiling(i/2) + 1), 
              seq(1.5/i, 1.5, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

#pdf("08.Microbial_GWAS/F6_comes_sv_lm_adjAbun.heatmap.pdf", width = 6, height = 10)

tiff(file = "pics/F6B_comes_sv_lm_adjAbun.heatmap.tiff", width =2800, height =1200, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#0072B5FF","#F5F5F5","#BC3C29FF"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("Melanoma\nResponse","Melanoma\nOS","Melanoma\nPFS >= 12 months","NSCLC\nResponse","NSCLC\nOS"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "grey50",notecex =2,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
)

dev.off()




##########################################################################################################
###  R.bromii

plot.pheno  <- sv_prog_adjAbun_res.sig.edge$Pheno[grep("bromii", sv_prog_adjAbun_res.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_adjAbun_res.sig.edge$SV[grep("bromii", sv_prog_adjAbun_res.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("bromii", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("bromii", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Ruminococcus bromii L2-63:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-2.0, 0, length.out=ceiling(i/2) + 1), 
              seq(2.0/i, 2.0, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

#pdf("08.Microbial_GWAS/F6_comes_sv_lm_adjAbun.heatmap.pdf", width = 6, height = 10)

tiff(file = "pics/F6B_Ruminococcus_bromii_sv_lm_adjAbun.heatmap.tiff", width =2800, height =1400, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#0072B5FF","#F5F5F5","#BC3C29FF"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("Melanoma\nResponse","Melanoma\nOS","Melanoma\nPFS >= 12 months","NSCLC\nResponse","NSCLC\nOS"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "grey50",notecex =2,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,14) # ("margin.Y", "margin.X")
)

dev.off()





