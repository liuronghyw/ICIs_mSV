### Survival and barplot for sites
### 2022-11-29
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SNV/pipeline/4_R")

## 1 Preparation
##rm(list = ls())

library("survival")
library("survminer")
library(gridExtra)
library("grid")
#library(reshape2)	  
#library("RColorBrewer")
library("plyr")

source("functions.R")

###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

all_basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

mel_basic<-subset(all_basic,cancer_type=="melanoma",drop=T)
NSCLC_basic<-subset(all_basic,cancer_type=="NSCLC",drop=T)
RCC_basic<-subset(all_basic,cancer_type=="RCC",drop=T)

load("01.cleanData/SV_all/dsgv.RData")
load("01.cleanData/SV_all/vsgv.RData")

vsgv<-vsgv_sub
dsgv<-dsgv_sub

#######################################################################################
### clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c("Coprococcus comes ATCC 27758:1138_1142",
"Faecalibacterium cf. prausnitzii KLE1255:1473_1475",
"Akkermansia muciniphila ATCC BAA-835:863_864",
"Alistipes shahii WAL 8301:2839_2840",
"Alistipes shahii WAL 8301:3629_3630",
"Bacteroides massiliensis B84634 = Timone 84634 = DSM 17679 = JCM 13223:2345_2347",
"Bacteroides salyersiae WAL 10018 = DSM 18765 = JCM 12988:4235_4236",
"Dialister invisus DSM 15470:1141_1142",
"Parabacteroides distasonis ATCC 8503:1574_1575",
"Parabacteroides distasonis ATCC 8503:4655_4656",
"Subdoligranulum sp. 4_3_54A2FAA:1282_1283",
"[Eubacterium] hallii DSM 3353:730_731",
"Akkermansia muciniphila ATCC BAA-835:863_864",
"Ruminococcus bromii L2-63:1764_1768",
"Ruminococcus bromii L2-63:1884_1885",
"Ruminococcus bromii L2-63:2002_2005",
"[Eubacterium] rectale DSM 17629:1116_1118;1119_1120",
"[Eubacterium] rectale DSM 17629:2820_2821;2823_2825",
"Akkermansia muciniphila ATCC BAA-835:83_84",
"Akkermansia muciniphila ATCC BAA-835:1062_1064 and 5 segments",
"Akkermansia muciniphila ATCC BAA-835:2528_2529;2532_2550",
"Bacteroides clarus YIT 12056:146_148 and 48 segments",
"Coprococcus comes ATCC 27758:1400_1402",
"Prevotella copri DSM 18205:1534_1535",
"Ruminococcus torques L2-14:1132_1137",
"Akkermansia muciniphila ATCC BAA-835:731_734",
"Bacteroides fragilis NCTC 9343:1069_1071;1071_1073")]


########################################################################################
##  response (melano)

dsgv_pic$id<-row.names(dsgv_pic)
pic<-merge(mel_basic,dsgv_pic,by.x="id",by.y="id")

#########################################################################
#####  [Eubacterium] hallii DSM 3353:730_731
pic$temp<-pic$"[Eubacterium] hallii DSM 3353:730_731"
lg_res_model <- glm(response_code~temp+age+read_count,family=binomial(logit), data=pic)
lg_res <- summary(lg_res_model)

pic_sub<-subset(pic,is.na(temp)==F,drop=T)
ratio<-table(pic_sub$response_code,pic_sub$temp)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,2)

for(i in 1:2){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Non-delection","Delection")
rownames(ratio_new)<-c("Non_RE","RE")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

table(pic_sub$temp)

tiff(file = "pics/sites/F6_dsv_mel_resp_E.hallii_730_731_barplot.tiff", width =1000, height =1800, res =300) 

###par(family = "serif")
ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
geom_col(position = 'stack', width =0.9) +ggtitle("")+
scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
labs(x = '', y = 'Response rate',title="Melanoma\nE.hallii:730_731")+
geom_text(aes(y=new_col, label=value), size=5,colour="White")+
scale_x_discrete(limits=flevels,
labels=c("Non-Delection(n=84)","Delection(n=37)"))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
strip.text = element_text(size = 15)) +
theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.title = element_blank(), 
legend.text = element_text(size = 15),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans")
)

dev.off()


########################################################################################
##  response Akkermansia muciniphila ATCC BAA-835:863_864

pic$temp<-pic$"Akkermansia muciniphila ATCC BAA-835:863_864"
lg_res_model <- glm(response_code~cohort+temp+age+read_count,family=binomial(logit), data=pic)
lg_res <- summary(lg_res_model)

pic_sub<-subset(pic,is.na(temp)==F,drop=T)
ratio<-table(pic_sub$response_code,pic_sub$temp)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,2)

for(i in 1:2){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Non-delection","Delection")
rownames(ratio_new)<-c("Non_RE","RE")

lg_res_model <- glm(response_code~cohort+temp+age+read_count,family=binomial(logit), data=pic)
lg_res <- summary(lg_res_model)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,2)

for(i in 1:2){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Non-delection","Delection")
rownames(ratio_new)<-c("Non_RE","RE")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

table(pic_sub$response_code)

tiff(file = "pics/sites/F6_dsv_resp_mel_Akk_863_864_barplot.tiff", width =1000, height =1800, res =300) 

###par(family = "serif")
ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
geom_col(position = 'stack', width =0.9) +ggtitle("")+
scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
labs(x = '', y = 'Response rate',title="Melanoma\nA.muciniphila:863_864")+
geom_text(aes(y=new_col, label=value), size=5,colour="White")+
scale_x_discrete(limits=flevels,
labels=c("Non-Delection(n=60)","Delection(n=53)"))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
strip.text = element_text(size = 15)) +
theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.title = element_blank(), 
legend.text = element_text(size = 15),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans")
)

dev.off()



########################################################################################
##  response Ruminococcus bromii L2-63:1884_1885

pic$temp<-pic$"Ruminococcus bromii L2-63:1884_1885"
lg_res_model <- glm(response_code~cohort+temp+age+read_count,family=binomial(logit), data=pic)
lg_res <- summary(lg_res_model)

pic_sub<-subset(pic,is.na(temp)==F,drop=T)
ratio<-table(pic_sub$response_code,pic_sub$temp)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,2)

for(i in 1:2){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Non-delection","Delection")
rownames(ratio_new)<-c("Non_RE","RE")

lg_res_model <- glm(response_code~cohort+temp+age+read_count,family=binomial(logit), data=pic)
lg_res <- summary(lg_res_model)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,2)

for(i in 1:2){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Non-delection","Delection")
rownames(ratio_new)<-c("Non_RE","RE")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

table(pic_sub$response_code)

tiff(file = "pics/sites/F6_dsv_mel_resp_R.bromii_1884_1885_barplot.tiff", width =1000, height =1800, res =300) 

###par(family = "serif")
ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
geom_col(position = 'stack', width =0.9) +ggtitle("")+
scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
labs(x = '', y = 'Response rate',title="Melanoma\nR.bromii:1884_1885")+
geom_text(aes(y=new_col, label=value), size=5,colour="White")+
scale_x_discrete(limits=flevels,
labels=c("Non-Delection(n=84)","Delection(n=61)"))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
strip.text = element_text(size = 15)) +
theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.title = element_blank(), 
legend.text = element_text(size = 15),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans")
)

dev.off()



########################################################################################
##  response Ruminococcus bromii L2-63:1764_1768

pic$temp<-pic$"Ruminococcus bromii L2-63:1764_1768"
lg_res_model <- glm(response_code~cohort+temp+age+read_count,family=binomial(logit), data=pic)
lg_res <- summary(lg_res_model)

pic_sub<-subset(pic,is.na(temp)==F,drop=T)
ratio<-table(pic_sub$response_code,pic_sub$temp)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,2)

for(i in 1:2){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Non-delection","Delection")
rownames(ratio_new)<-c("Non_RE","RE")

lg_res_model <- glm(response_code~cohort+temp+age+read_count,family=binomial(logit), data=pic)
lg_res <- summary(lg_res_model)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,2)

for(i in 1:2){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Non-delection","Delection")
rownames(ratio_new)<-c("Non_RE","RE")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

table(pic_sub$response_code)

tiff(file = "pics/sites/F6_dsv_mel_resp_R.bromii_1764_1768_barplot.tiff", width =1000, height =1800, res =300) 

###par(family = "serif")
ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
geom_col(position = 'stack', width =0.9) +ggtitle("")+
scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
labs(x = '', y = 'Response rate',title="Melanoma\nR.bromii:1764_1768")+
geom_text(aes(y=new_col, label=value), size=5,colour="White")+
scale_x_discrete(limits=flevels,
labels=c("Non-Delection(n=84)","Delection(n=61)"))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
strip.text = element_text(size = 15)) +
theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.title = element_blank(), 
legend.text = element_text(size = 15),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans")
)

dev.off()


#################################################################################
## os Faecalibacterium cf. prausnitzii KLE1255:1473_1475
pic_sub<-pic
pic_sub$subtype<-pic$"Faecalibacterium cf. prausnitzii KLE1255:1473_1475"
a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_mel_OS_F.prausnitzii_1473_1475.tiff", width =1400, height =1600, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = mycolor2_green_blue, 
conf.int = F, pval =T,  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(1.5,0.95),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="Melanoma\nF.prausnitzii:1473_1475",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Non-delection","Delection"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()


#################################################################################
## pfs  Akkermansia muciniphila ATCC BAA-835:863_864

pic$temp<-pic$"Akkermansia muciniphila ATCC BAA-835:863_864"
lg_res_model <- glm(pfs_12_months~cohort+temp+age+read_count,family=binomial(logit), data=pic)
lg_res <- summary(lg_res_model)

pic_sub<-subset(pic,is.na(temp)==F,drop=T)
ratio<-table(pic_sub$pfs_12_months,pic_sub$temp)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,2)

for(i in 1:2){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Non-delection","Delection")
rownames(ratio_new)<-c("PFS_less12","PFS_larger12")

lg_res_model <- glm(pfs_12_months~cohort+temp+age+read_count,family=binomial(logit), data=pic)
lg_res <- summary(lg_res_model)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,2)

for(i in 1:2){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Non-delection","Delection")
rownames(ratio_new)<-c("PFS_less12","PFS_larger12")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

table(pic_sub$pfs_12_months)

tiff(file = "pics/sites/F6_dsv_mel_pfs_AKK863_864_barplot.tiff", width =1080, height =1800, res =300) 

###par(family = "serif")
ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
geom_col(position = 'stack', width =0.9) +ggtitle("")+
scale_fill_manual(breaks=c("PFS_less12","PFS_larger12"),labels=c("PFS<12m","PFS>=12m"),values = rev(mycolor2_yellow_green )) +
labs(x = '', y = 'Response rate',title="Melanoma\nA.muciniphila:863_864")+
geom_text(aes(y=new_col, label=value), size=5,colour="White")+
scale_x_discrete(limits=flevels,
labels=c("Non-Delection(n=56)","Delection(n=52)"))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
strip.text = element_text(size = 15)) +
theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.title = element_blank(), 
legend.text = element_text(size = 15),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans")
)

dev.off()



##################################################################################################
### NSCLC Akkermansia muciniphila ATCC BAA-835:731_734

dsgv_pic$id<-row.names(dsgv_pic)
pic<-merge(NSCLC_basic,dsgv_pic,by.x="id",by.y="id")

#########################################################################
#####  
pic$temp<-pic$"Akkermansia muciniphila ATCC BAA-835:731_734"
lg_res_model <- glm(response_code~temp+age+read_count,family=binomial(logit), data=pic)
lg_res <- summary(lg_res_model)

pic_sub<-subset(pic,is.na(temp)==F,drop=T)
ratio<-table(pic_sub$response_code,pic_sub$temp)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,2)

for(i in 1:2){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Non-delection","Delection")
rownames(ratio_new)<-c("Non_RE","RE")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

table(pic_sub$temp)

tiff(file = "pics/sites/F6_dsv_NSCLC_resp_Akk_731_734_barplot.tiff", width =1100, height =1800, res =300) 

###par(family = "serif")
ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
geom_col(position = 'stack', width =0.9) +ggtitle("")+
scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
labs(x = '', y = 'Response rate',title="NSCLC\nA.muciniphila:731_734")+
geom_text(aes(y=new_col, label=value), size=5,colour="White")+
scale_x_discrete(limits=flevels,
labels=c("Non-Delection(n=53)","Delection(n=84)"))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
strip.text = element_text(size = 15)) +
theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.title = element_blank(), 
legend.text = element_text(size = 15),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans")
)

dev.off()



#################################################################################
## os Akkermansia muciniphila ATCC BAA-835:1062_1064 and 5 segments
pic_sub<-pic
pic_sub$subtype<-pic$"Akkermansia muciniphila ATCC BAA-835:1062_1064 and 5 segments"

a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_NSCLC_OS_A.muciniphila_1062_1064 and 5 segments.tiff", width =1600, height =1600, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = mycolor2_green_blue, 
conf.int = F, pval =T,  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(1.5,0.95),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="NSCLC\nA.muciniphila:\n1062_1064 and 5 segments",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Non-delection","Delection"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()


#################################################################################
## os [Eubacterium] rectale DSM 17629:1116_1118;1119_1120
pic_sub<-pic
pic_sub$subtype<-pic$"[Eubacterium] rectale DSM 17629:1116_1118;1119_1120"

a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_NSCLC_OS_E.rectale1116_1118_1119_1120.tiff", width =1600, height =1600, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = mycolor2_green_blue, 
conf.int = F, pval =T,  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(1.5,0.95),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="NSCLC\nE.rectale:\n1116_1118;1119_1120",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Non-delection","Delection"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()



##################################################################################
### vsgv picture
all_basic$id<-row.names(all_basic)

vsgv_pic<-vsgv[,c(
"Lachnospiraceae bacterium 3_1_46FAA:911_912;912_913",
"Phascolarctobacterium sp. CAG:207:890_891",
"Bacteroides xylanisolvens XB1A:3046_3048",
"Alistipes shahii WAL 8301:991_992 and 4 segments",
"Bacteroides clarus YIT 12056:3156_3159",
"Ruminococcus bromii L2-63:41_43 and 9 segments",
"Blautia wexlerae DSM 19850:3538_3541",
"Blautia wexlerae DSM 19850:3545_3547",
"Blautia wexlerae DSM 19850:3534_3537",
"Blautia wexlerae DSM 19850:3533_3534",
"Blautia wexlerae DSM 19850:1218_1219",
"Blautia wexlerae DSM 19850:4492_4495;4495_4500",
"Blautia wexlerae DSM 19850:2787_2788 and 6 segments",
"Blautia wexlerae DSM 19850:609_610;610_611",
"Blautia wexlerae DSM 19850:4425_4426")]

vsgv_pic$id<-row.names(vsgv_pic)

### mel 
pic<-merge(mel_basic,vsgv_pic,by.x="id",by.y="id")

#################################################################################
## os Lachnospiraceae bacterium 3_1_46FAA:911_912;912_913

library(RColorBrewer)
#green_colors<-brewer.pal(6,"Greens")
green_colors<-c("#A1D99B","#74C476","#31A354","#006D2C")

pic_sub<-pic
pic_sub$temp<-pic$"Lachnospiraceae bacterium 3_1_46FAA:911_912;912_913"

qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
 pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }

a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_mel_OS_vsv_Lachnospiraceae bacterium_3_1_46FAA_911_912_912_913.tiff", width =1600, height =1600, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = green_colors, 
conf.int = F, pval =T,  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(1.5,0.95),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="Melanoma:\nLachnospiraceae bacterium 3\n911_912;912_913",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Q1","Q2","Q3","Q4"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()


#################################################################################
## os Phascolarctobacterium sp. CAG:207:890_891

library(RColorBrewer)
#green_colors<-brewer.pal(6,"Greens")
green_colors<-c("#A1D99B","#74C476","#31A354","#006D2C")

pic_sub<-pic
pic_sub$temp<-pic$"Phascolarctobacterium sp. CAG:207:890_891"

pic_sub<-subset(pic_sub,is.na(os)==F,drop=T)

qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
 pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }

a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_mel_OS_vsv_Phascolarctobacterium spCAG_207:890_891.tiff", width =1600, height =1600, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = green_colors, 
conf.int = F, pval =T,  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(1.5,0.95),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="Melanoma:\nPhascolarctobacterium sp.\n890_891",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Q1","Q2","Q3","Q4"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()



#################################################################################
## pfs Bacteroides xylanisolvens XB1A:3046_3048
pic_sub<-pic
pic_sub$temp<-pic$"Bacteroides xylanisolvens XB1A:3046_3048"

qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
 pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }

lg_res_model <- glm(pfs_12_months~cohort+subtype+age+read_count,family=binomial(logit), data=pic_sub)
lg_res <- summary(lg_res_model)

ratio<-table(pic_sub$pfs_12_months,pic_sub$subtype)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,4)

for(i in 1:4){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Q1","Q2","Q3","Q4")
rownames(ratio_new)<-c("PFS_less12","PFS_larger12")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

table(pic_sub$subtype)

tiff(file = "pics/sites/F6_vsv_mel_pfs_Bx_barplot.tiff", width =1280, height =1800, res =300) 

###par(family = "serif")
ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
geom_col(position = 'stack', width =0.9) +ggtitle("")+
scale_fill_manual(breaks=c("PFS_less12","PFS_larger12"),labels=c("PFS<12m","PFS>=12m"),values = rev(mycolor2_yellow_green )) +
labs(x = '', y = 'Response rate',title="Melanoma\nBacteroides xylanisolvens:3046_3048")+
geom_text(aes(y=new_col, label=value), size=5,colour="White")+
scale_x_discrete(limits=flevels,
labels=c("Q1(n=39)","Q2(n=38)","Q3(n=38)","Q4(n=39)"))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
strip.text = element_text(size = 15)) +
theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.title = element_blank(), 
legend.text = element_text(size = 15),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans")
)

dev.off()




#################################################################################
## irAEs Blautia wexlerae DSM 19850:3533_3534
pic_sub<-pic
pic_sub<-subset(pic_sub,is.na(irAEs)==F,drop=T)
pic_sub$temp<-pic_sub$"Blautia wexlerae DSM 19850:3533_3534"

qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
 pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }

lg_res_model <- glm(irAEs~cohort+subtype+age+read_count,family=binomial(logit), data=pic_sub)
lg_res <- summary(lg_res_model)

ratio<-table(pic_sub$irAEs,pic_sub$subtype)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,4)

for(i in 1:4){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Q1","Q2","Q3","Q4")
rownames(ratio_new)<-c("Non_irAEs","irAEs")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

table(pic_sub$subtype)

tiff(file = "pics/sites/F6_vsv_mel_irAEs_Bx_barplot.tiff", width =1280, height =1800, res =300) 

###par(family = "serif")
ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
geom_col(position = 'stack', width =0.9) +ggtitle("")+
scale_fill_manual(breaks=c("Non_irAEs","irAEs"),labels=c("irAEs (No)","irAEs (Yes)"),values =c("#FFDC91FF","#EE4C97FF")) +
labs(x = '', y = 'Response rate',title="Melanoma\nBlautia wexlerae DSM 19850:3533_3534")+
geom_text(aes(y=new_col, label=value), size=5,colour="White")+
scale_x_discrete(limits=flevels,
labels=c("Q1(n=12)","Q2(n=11)","Q3(n=11)","Q4(n=11)"))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
strip.text = element_text(size = 15)) +
theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.title = element_blank(), 
legend.text = element_text(size = 15),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans")
)

dev.off()




#################################################################################
## irAEs Blautia wexlerae DSM 19850:1218_1219
pic_sub<-pic
pic_sub<-subset(pic_sub,is.na(irAEs)==F,drop=T)
pic_sub$temp<-pic_sub$"Blautia wexlerae DSM 19850:1218_1219"

qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
 pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }

lg_res_model <- glm(irAEs~cohort+subtype+age+read_count,family=binomial(logit), data=pic_sub)
lg_res <- summary(lg_res_model)

ratio<-table(pic_sub$irAEs,pic_sub$subtype)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,4)

for(i in 1:4){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Q1","Q2","Q3","Q4")
rownames(ratio_new)<-c("Non_irAEs","irAEs")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

table(pic_sub$subtype)

tiff(file = "pics/sites/F6_vsv_mel_irAEs_Bx_barplot_2.tiff", width =1280, height =1800, res =300) 

###par(family = "serif")
ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
geom_col(position = 'stack', width =0.9) +ggtitle("")+
scale_fill_manual(breaks=c("Non_irAEs","irAEs"),labels=c("irAEs (No)","irAEs (Yes)"),values =c("#FFDC91FF","#EE4C97FF")) +
labs(x = '', y = 'Response rate',title="Melanoma\nBlautia wexlerae DSM 19850:1218_1219")+
geom_text(aes(y=new_col, label=value), size=5,colour="White")+
scale_x_discrete(limits=flevels,
labels=c("Q1(n=12)","Q2(n=11)","Q3(n=11)","Q4(n=11)"))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
strip.text = element_text(size = 15)) +
theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15), legend.title = element_blank(), 
legend.text = element_text(size = 15),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans")
)

dev.off()




##################################################################################################
## NSCLC vsv os

#Alistipes shahii WAL 8301:991_992 and 4 segments
#Bacteroides clarus YIT 12056:3156_3159
#Ruminococcus bromii L2-63:41_43 and 9 segments

library(RColorBrewer)
#green_colors<-brewer.pal(6,"Greens")
green_colors<-c("#A1D99B","#74C476","#31A354","#006D2C")

pic<-merge(NSCLC_basic,vsgv_pic,by.x="id",by.y="id")
pic_sub<-pic
pic_sub$temp<-pic$"Bacteroides clarus YIT 12056:3156_3159"

qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
 pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }

a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_NSCLC_OS_vsv_B.clarus_807_808.tiff", width =1600, height =1600, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = green_colors, 
conf.int = F, pval =T,  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(1.5,0.95),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="NSCLC:\nB.clarus:\n3156_3159",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Q1","Q2","Q3","Q4"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()




library(RColorBrewer)
#green_colors<-brewer.pal(6,"Greens")
green_colors<-c("#A1D99B","#74C476","#31A354","#006D2C")

pic<-merge(NSCLC_basic,vsgv_pic,by.x="id",by.y="id")
pic_sub<-pic
pic_sub$temp<-pic$"Ruminococcus bromii L2-63:41_43 and 9 segments"

qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
 pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }

a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_NSCLC_OS_vsv_Ruminococcus bromii L2_63_41_43 and 9 segments.tiff", width =1600, height =1600, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = green_colors, 
conf.int = F, pval =T,  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(1.5,0.95),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="NSCLC:\nR.bromii:\n41_43 and 9 segments",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Q1","Q2","Q3","Q4"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()







########################################################
#### 3.2 Replication 
#### 3.2.1 vSV

### Correlation of effect size between different cancer types 
###  response 

RCC_dsv_resp<-read.csv("08.Microbial_GWAS/Rdata/RCC_dsv_resp_adjAbun_res.csv",header=T)
RCC_vsv_resp<-read.csv("08.Microbial_GWAS/Rdata/RCC_vsv_resp_adjAbun_res.csv",header=T)
RCC_resp<-rbind(RCC_dsv_resp,RCC_vsv_resp)
RCC_resp_sub<-RCC_resp[,c("X","OR")]
colnames(RCC_resp_sub)<-c("id","RCC_OR")

NSCLC_dsv_resp<-read.csv("08.Microbial_GWAS/RData/NSCLC_dsv_resp_adjAbun_res.csv",header=T)
NSCLC_vsv_resp<-read.csv("08.Microbial_GWAS/RData/NSCLC_vsv_resp_adjAbun_res.csv",header=T)
NSCLC_resp<-rbind(NSCLC_dsv_resp,NSCLC_vsv_resp)
NSCLC_resp_sub<-NSCLC_resp[,c("X","OR")]
colnames(NSCLC_resp_sub)<-c("id","NSCLC_OR")

melanoma_dsv_resp<-read.csv("08.Microbial_GWAS/RData/melanoma_dsv_resp_adjAbun_res.csv",header=T)
melanoma_vsv_resp<-read.csv("08.Microbial_GWAS/RData/melanoma_vsv_resp_adjAbun_res.csv",header=T)
melanoma_resp<-rbind(melanoma_dsv_resp,melanoma_vsv_resp)
melanoma_resp_sub<-melanoma_resp[,c("X","OR")]
colnames(melanoma_resp_sub)<-c("id","melanoma_OR")

resp<-merge(RCC_resp_sub,NSCLC_resp_sub,by.x="id",by.y="id")
resp<-merge(resp,melanoma_resp_sub,by.x="id",by.y="id")

####################################################################
## overall and melanoma 

resp<-subset(resp,melanoma_OR<30,drop=T)
resp<-subset(resp,NSCLC_OR<30,drop=T)
es_cor<-cor.test(resp$NSCLC_OR, resp$melanoma_OR,method = "spearman")

text_r<-paste("R=",round(es_cor$estimate,digits = 2),
              ", P=0.1739",
              sep = "")

p_resp_mel<-ggplot(resp) + 
  geom_point(aes(NSCLC_OR, melanoma_OR),alpha = 0.4) +
  annotate("text",x = c(5), y = c(7.5), label = c(text_r),hjust = 0,size=8)+
  geom_hline(yintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  geom_vline(xintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  xlab('OR (NSCLC)')+
  ylab('OR (melanoma)')+
  scale_color_viridis()+
 theme_bw(base_size = 17)+
  theme(legend.position = 'none',axis.text=element_text(size=20))


####################################################################
## Melanoma and RCC

resp<-subset(resp,RCC_OR<30,drop=T)
resp<-subset(resp,melanoma_OR<30,drop=T)
es_cor<-cor.test(resp$melanoma_OR, resp$RCC_OR,method = "spearman")

text_r<-paste("R=",round(es_cor$estimate,digits = 2),
              ", P=0.02",
              sep = "")

p_resp_rcc<-ggplot(resp) + 
  geom_point(aes(melanoma_OR, RCC_OR),alpha = 0.4) +
  annotate("text",x = c(5), y = c(5), label = c(text_r),hjust = 0,size=8)+
  geom_hline(yintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  geom_vline(xintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  xlab('OR (Melanoma)')+
  ylab('OR (RCC)')+
  scale_color_viridis()+
 theme_bw(base_size = 17)+
  theme(legend.position = 'none',axis.text=element_text(size=20))



####################################################################
## RCC and NSCLC

resp<-subset(resp,NSCLC_OR<30,drop=T)
resp<-subset(resp,RCC_OR<30,drop=T)
es_cor<-cor.test(resp$RCC_OR, resp$NSCLC_OR,method = "spearman")

text_r<-paste("R=",round(es_cor$estimate,digits = 2),
              ", P=0.2243",
              sep = "")

p_resp_NSCLC<-ggplot(resp) + 
  geom_point(aes(RCC_OR, NSCLC_OR),alpha = 0.4) +
  annotate("text",x = c(5), y = c(5), label = c(text_r),hjust = 0,size=8)+
  geom_hline(yintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  geom_vline(xintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  xlab('OR (RCC)')+
  ylab('OR (NSCLC)')+
  scale_color_viridis()+
 theme_bw(base_size = 17)+
  theme(legend.position = 'none',axis.text=element_text(size=20))




#######################################################################################
## os 

NSCLC_dsv_os<-read.csv("08.Microbial_GWAS/RData/NSCLC_dsv_os_adjAbun_res.csv",header=T)
NSCLC_vsv_os<-read.csv("08.Microbial_GWAS/RData/NSCLC_vsv_os_adjAbun_res.csv",header=T)
NSCLC_os<-rbind(NSCLC_dsv_os,NSCLC_vsv_os)
NSCLC_os_sub<-NSCLC_os[,c("X","HR")]
colnames(NSCLC_os_sub)<-c("id","NSCLC_HR")

melanoma_dsv_os<-read.csv("08.Microbial_GWAS/RData/melanoma_dsv_os_adjAbun_res.csv",header=T)
melanoma_vsv_os<-read.csv("08.Microbial_GWAS/RData/melanoma_vsv_os_adjAbun_res.csv",header=T)
melanoma_os<-rbind(melanoma_dsv_os,melanoma_vsv_os)
melanoma_os_sub<-melanoma_os[,c("X","HR")]
colnames(melanoma_os_sub)<-c("id","melanoma_HR")

os<-merge(NSCLC_os_sub,melanoma_os_sub,by.x="id",by.y="id")

####################################################################
## overall and melanoma 

os<-subset(os,melanoma_HR<30,drop=T)
os<-subset(os,NSCLC_HR<30,drop=T)
es_cor<-cor.test(os$NSCLC_HR, os$melanoma_HR,method = "spearman")

text_r<-paste("R=",round(es_cor$estimate,digits = 2),
              ", P=0.08",
              sep = "")

p_os_melanoma<-ggplot(os) + 
  geom_point(aes(NSCLC_HR, melanoma_HR),alpha = 0.4) +
  annotate("text",x = c(3), y = c(5), label = c(text_r),hjust = 0,size=8)+
  geom_hline(yintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  geom_vline(xintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  xlab('HR (NSCLC)')+
  ylab('HR (melanoma)')+
  scale_color_viridis()+
  theme_bw(base_size = 17)+
  theme(legend.position = 'none',axis.text=element_text(size=20))



#### 3.2.3 Combine figure
## plot

p_title_replic <- ggdraw() + 
    #draw_label(
    #   'Replication',
    #  fontface = 'bold', x = 0, hjust = 0) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))

p_sv_prog_replic<-plot_grid(
    p_title_replic, 
    plot_grid(p_resp_mel,p_resp_rcc,p_resp_NSCLC,p_os_melanoma,
              rel_widths = c(1, 1, 1,1),align = 'hv',
              labels = c("SVs-response OR", "SVs-response OR","SVs-response OR",
                         "SVs-os HR "),
              ncol = 2,label_size	=18,vjust = 0),
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 2)
  )



tiff(file = "pics/FS7_replicate.tiff", width =4000, height =4000, res =300) 
print(p_sv_prog_replic)
dev.off()


