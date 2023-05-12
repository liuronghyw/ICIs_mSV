### Subspecies 
### 2022-11-24
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SNV/pipeline/4_R")

## 1 Preparation
### 1.1 Import
##install.packages("ggExtra")

source("functions.R")
knitr::opts_chunk$set(echo = TRUE)

### 1.2 Inputs
## Read input files.

info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",
                      sep = "\t", header = T, stringsAsFactors = F)
all_basic <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")
all_dsv <- read.table("01.cleanData/SV_all/deletionStructuralVariation_all.tsv",check.names = F)
all_vsv <- read.table("01.cleanData/SV_all/variableStructuralVariation_all.tsv",check.names = F)

## 2 Clustering analysis
load("01.cleanData/SV_all/distMat/all_msv_dist_std.RData")

all_msv_wsd_res <- lapply(all_msv_dist_std, my_cluster, ps.cutoff = 0.55)
beepr::beep("mario")

if(!dir.exists("07.Subspecies")){dir.create("07.Subspecies")}
save(all_msv_wsd_res,file = "07.Subspecies/all_msv_wsd_res_ps0.55.RData")

## Optimum cluster number
load("07.Subspecies/all_msv_wsd_res_ps0.55.RData")
cluster_n <- as.data.frame(matrix(NA, nrow = nrow(info), ncol = 13))

for (i in c(1:nrow(info))) {
  #i<-1
  if(is.na(all_msv_wsd_res[[i]]$clu_n)==F){
  all_msv_wsd_res[[i]]$clu_n$mean.pred
  cluster_n[i,]<-c(info$Short_name[i], all_msv_wsd_res[[i]]$clu_n$optimalk, max(all_msv_wsd_res[[i]]$clu_n$mean.pred[-1], na.rm = T),all_msv_wsd_res[[i]]$clu_n$mean.pred)
 }
}

colnames(cluster_n)<-c("Short_name", "Optimum_cluster_n", "Maximum_prediction_strength",
                       "PS_n_1", "PS_n_2", "PS_n_3", "PS_n_4", "PS_n_5",
                       "PS_n_6", "PS_n_7", "PS_n_8", "PS_n_9", "PS_n_10")

write.table(cluster_n,"07.Subspecies/all_cluster_n.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
#write.csv(cluster_n,"07.Subspecies/all_cluster_n.csv", row.names = F)


## Get PCs
msv_pc<-as.data.frame(matrix(NA, nrow = nrow(all_basic), ncol = 5*nrow(info)))
rownames(msv_pc)<-rownames(all_basic)
colname_suf<-c("_PC1", "_PC2", "_PC3", "_PC4", "_PC5")

for (i in c(1:nrow(info))) {
  #i <- 1
  msv_pc[match(rownames(all_msv_wsd_res[[i]]$pcoa), rownames(msv_pc)),c((5*(i-1)+1):(5*i))] <- all_msv_wsd_res[[i]]$pcoa
  colnames(msv_pc)[c((5*(i-1)+1):(5*i))]<-paste(info$Short_name[i],colname_suf, sep = "")
}

write.table(msv_pc, "07.Subspecies/all_msv_pc.tsv",sep = "\t", quote = F)

## Get cluster profile
all_msv_cluster <- as.data.frame(matrix(NA, nrow = nrow(all_basic), ncol = nrow(info)))
rownames(all_msv_cluster)<-rownames(all_basic)
colnames(all_msv_cluster)<-info$Short_name

for (i in 1:nrow(info)) {
  #i <- 1
  if(is.na(all_msv_wsd_res[[i]]$clu_n)==F){
  all_msv_cluster[match(rownames(all_msv_wsd_res[[i]]$tsne_df), rownames(all_msv_cluster)),i] <- all_msv_wsd_res[[i]]$tsne_df$Cluster
 }
}

write.table(all_msv_cluster, "07.Subspecies/all_msv_cluster.tsv",sep = "\t", quote = F)


## PCoA plot panel
## cohort
pcoa_plot.list<-list()

info_sub<-info[-c(7,9,10,11,12,28,35),]
all_msv_wsd_res_sub<-all_msv_wsd_res[-c(7,9,10,11,12,28,35)]

##all_msv_dist_std[which(names(all_msv_dist_std)%in%c("msv_Acidaminococcus intestini RyC-MR95"))]<-NULL


for (i in 1:nrow(info_sub)) {
  #i<-1
  
  pcoa<-all_msv_wsd_res_sub[[i]]$pcoa
  pcoa$dataset<-all_basic$dataset[match(rownames(pcoa), rownames(all_basic))]
  p_pcoa<-ggplot(pcoa,aes(X1,X2,fill =dataset, color =dataset))+
  geom_point(size = 2,alpha = 0.5)+
  ggtitle(info_sub$Short_name[i])+
  xlab(paste("PCo1=",round(all_msv_wsd_res_sub[[i]]$pcoa_res$eig[1],digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(all_msv_wsd_res_sub[[i]]$pcoa_res$eig[2],digits = 2),"%",sep = ""))+
  scale_color_manual(name=NULL,
                      breaks = c("PRJEB22863", "PRJNA397906",
                              "PRJNA541981","PRJNA751792",
                              "PRJNA762360","PRJNA770295","PRJEB43119"),
                      labels = c("RoutyB_2018", "FrankelAE_2017","PetersBA_2020",
                               "DerosaL_2022","McCullochJA_2022","SpencerCN_2021","LeeKA_2022"),
                      values = mycolor7)+
  scale_fill_manual(name=NULL,
                      breaks = c("PRJEB22863", "PRJNA397906",
                              "PRJNA541981","PRJNA751792",
                              "PRJNA762360","PRJNA770295","PRJEB43119"),
                      labels = c("RoutyB_2018", "FrankelAE_2017","PetersBA_2020",
                                 "DerosaL_2022","McCullochJA_2022","SpencerCN_2021","LeeKA_2022"),
                      values = mycolor7)+
  theme(plot.title = element_text(size=10, face="italic"),
        plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))
  p_pcoa<-ggExtra::ggMarginal(p_pcoa, type = "histogram", groupColour = F, groupFill = TRUE,
                            xparams = list(bins = 50, alpha = 0.5,position = 'identity', color = 'white'),
                            yparams = list(bins = 50, alpha = 0.5,position = 'identity', color = 'white'))

  pcoa_plot.list[[i]]<-p_pcoa
}


#pdf("pics/Figrue7_all_pcoa_cohort.pdf",width = 18,height = 18)
#plot_grid(plotlist=pcoa_plot.list)
#dev.off()

tiff(file = "pics/F5_all_abun_sv_PCoA.tiff", width =6000, height =6000, res =300) 
plot_grid(plotlist=pcoa_plot.list)
dev.off()



## cluster
pcoa_plot.list<-list()

for (i in 1:nrow(info_sub)) {
  #i<-1
  pcoa<-all_msv_wsd_res_sub[[i]]$pcoa
  p_pcoa<-ggplot(pcoa,aes(X1,X2,fill = as.factor(all_msv_wsd_res_sub[[i]]$tsne_df$Cluster), color = all_msv_wsd_res_sub[[i]]$tsne_df$Cluster))+
    geom_point(size = 2,alpha = 0.5)+
    ggtitle(info$Short_name[i])+
    xlab(paste("PCo1=",round(all_msv_wsd_res_sub[[i]]$pcoa_res$eig[1],digits = 2),"%",sep = ""))+
    ylab(paste("PCo2=",round(all_msv_wsd_res_sub[[i]]$pcoa_res$eig[2],digits = 2),"%",sep = ""))+
    theme(plot.title = element_text(size=10, face="italic"),
          plot.subtitle = element_text(vjust = 1), 
          plot.caption = element_text(vjust = 1), 
          axis.line.x =  element_line(),
          axis.line.y = element_line(),
          legend.position = 'none',
          legend.title = element_blank(),
          legend.key = element_rect(fill = NA), 
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          panel.background = element_rect(fill = NA))+ 
    scale_color_npg()
    p_pcoa<-ggExtra::ggMarginal(p_pcoa, type = "histogram", groupColour = F, groupFill = TRUE,
                                xparams = list(bins = 50, alpha = 0.5,position = 'identity', color = 'white'),
                               yparams = list(bins = 50, alpha = 0.5,position = 'identity', color = 'white'))
  pcoa_plot.list[[i]]<-p_pcoa
}

#pdf("pics/all_pcoa.pdf",width = 18,height = 18)
#plot_grid(plotlist=pcoa_plot.list)
#dev.off()

tiff(file = "pics/F5_all_pcoa.tiff", width =6000, height =6000, res =300) 
plot_grid(plotlist=pcoa_plot.list)
dev.off()



## tsne plot panel
tsne_plot.list<-list()

for (i in c(1:nrow(info_sub))) {
  #i<-1
  tsne_df<-all_msv_wsd_res_sub[[i]]$tsne_df
  p_msv_tsne <- ggplot(tsne_df, aes(x = X, y = Y)) +
    stat_ellipse(aes(group = Cluster, fill = Cluster, color = Cluster) ,
                 type = "norm",linetype = 2, geom = "polygon", alpha = 0.05)+
    geom_point(aes(color = Cluster), alpha = 0.5, size = 0.8)+
    ggtitle(info_sub$Short_name[i])+
    scale_color_npg()+
    scale_fill_npg()+
    theme_void()+
    theme(legend.position = 'none',
          plot.title = element_text(size=5, face="italic"))
  tsne_plot.list[[i]] <- p_msv_tsne
}

#pdf("pics/all_tsne_clusters.pdf")
#plot_grid(plotlist=tsne_plot.list)
#dev.off()

tiff(file = "pics/F5_all_tsne_clusters.tiff", width =2500, height =2500, res =300) 
plot_grid(plotlist=tsne_plot.list)
dev.off()


## 3 Prognosis and within-species diversity
### 3.1 Associations between response and clusters

all_resp<-all_basic[,c("response_code","pfs_12_months")]

mel_basic<-subset(all_basic,cancer_type=="melanoma",drop=T)
NSCLC_basic<-subset(all_basic,cancer_type=="NSCLC",drop=T)
RCC_basic<-subset(all_basic,cancer_type=="RCC",drop=T)

mel_resp<-mel_basic[,c("response_code","pfs_12_months","pfs_6_months")]
NSCLC_resp<-NSCLC_basic[,c("response_code","pfs_6_months")]
RCC_resp<-RCC_basic[,c("response_code","pfs_6_months")]

for (i in 1:nrow(cluster_n)){
  cluster_n[,2][is.na(cluster_n[,2])] <- 0
} 

all_msv_cluster_sub<-all_msv_cluster[,cluster_n$Optimum_cluster_n>1]
mel_msv_cluster_sub<-all_msv_cluster_sub[rownames(mel_basic),]
NSCLC_msv_cluster_sub<-all_msv_cluster_sub[rownames(NSCLC_basic),]
RCC_msv_cluster_sub<-all_msv_cluster_sub[rownames(RCC_basic),]

#head(all_msv_cluster_sub)
#head(all_prog)
#all<-cbind(all_prog,all_msv_cluster_sub)
#write.csv(all,"all_Test.csv")

mel_cluster_resp.res <-permKW_btw_mats(mel_resp, mel_msv_cluster_sub)
save(mel_cluster_resp.res, file = "07.Subspecies/mel_cluster_response.res.RData")

NSCLC_cluster_resp.res <-permKW_btw_mats(NSCLC_resp, NSCLC_msv_cluster_sub)
save(NSCLC_cluster_resp.res, file = "07.Subspecies/NSCLC_cluster_response.res.RData")

RCC_cluster_resp.res <-permKW_btw_mats(RCC_resp, RCC_msv_cluster_sub)
save(RCC_cluster_resp.res, file = "07.Subspecies/RCC_cluster_response.res.RData")


###### 3.2 Associations between survival and clusters

mel_survival<-mel_basic[,c("os","os_event","pfs","pfs_event")]
mel_survival$id<-row.names(mel_survival)
mel_msv_cluster_sub$id<-row.names(mel_msv_cluster_sub)

cluster_cox_os(mel_survival,mel_msv_cluster_sub,"07.Subspecies/mel_cluster_os_res.csv")
cluster_cox_pfs(mel_survival,mel_msv_cluster_sub,"07.Subspecies/mel_cluster_pfs_res.csv")


NSCLC_survival<-NSCLC_basic[,c("os","os_event")]
NSCLC_survival$id<-row.names(NSCLC_survival)
NSCLC_survival<-subset(NSCLC_survival,is.na(os)==F,drop=T)

NSCLC_msv_cluster_surv_sub<-NSCLC_msv_cluster_sub[NSCLC_survival$id,]
NSCLC_msv_cluster_sub$id<-row.names(NSCLC_msv_cluster_sub)

cluster_cox_os(NSCLC_survival,NSCLC_msv_cluster_sub,"07.Subspecies/NSCLC_cluster_os_res.csv")


################################################################ draw pictures
load("07.Subspecies/mel_cluster_response.res.RData")
load("07.Subspecies/NSCLC_cluster_response.res.RData")
load("07.Subspecies/RCC_cluster_response.res.RData")

mel_cluster_res.edge<-mel_cluster_resp.res$table
write.table(mel_cluster_res.edge, "07.Subspecies/mel_cluster_resp_res.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

mel_resp_sub<-subset(mel_cluster_res.edge,Phenotype=="response_code",drop=T)
mel_resp_sub$Phenotype<-"Melanoma\nResponse"
mel_resp_sub<-mel_resp_sub[,c("Taxa","Phenotype","p")]

mel_pfs_12_sub<-subset(mel_cluster_res.edge,Phenotype=="pfs_12_months",drop=T)
mel_pfs_12_sub$Phenotype<-"Melanoma\nPFS >= 12 months"
mel_pfs_12_sub<-mel_pfs_12_sub[,c("Taxa","Phenotype","p")]

NSCLC_cluster_res.edge<-NSCLC_cluster_resp.res$table
write.table(NSCLC_cluster_res.edge, "07.Subspecies/NSCLC_cluster_resp_res.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

NSCLC_resp_sub<-subset(NSCLC_cluster_res.edge,Phenotype=="response_code",drop=T)
NSCLC_resp_sub$Phenotype<-"NSCLC\nResponse"
NSCLC_resp_sub<-NSCLC_resp_sub[,c("Taxa","Phenotype","p")]

RCC_cluster_res.edge<-RCC_cluster_resp.res$table
write.table(RCC_cluster_res.edge, "07.Subspecies/RCC_cluster_resp_res.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

RCC_resp_sub<-subset(RCC_cluster_res.edge,Phenotype=="response_code",drop=T)
RCC_resp_sub$Phenotype<-"RCC\nResponse"
RCC_resp_sub<-RCC_resp_sub[,c("Taxa","Phenotype","p")]

mel_os<-read.csv("07.Subspecies/mel_cluster_os_res.csv",header=T)
mel_os_sub<-mel_os[,c("X","p")]
colnames(mel_os_sub)<-c("Taxa","p")
mel_os_sub$Phenotype<-"Melanoma\nOS"

mel_pfs<-read.csv("07.Subspecies/mel_cluster_pfs_res.csv",header=T)
mel_pfs_sub<-mel_pfs[,c("X","p")]
colnames(mel_pfs_sub)<-c("Taxa","p")
mel_pfs_sub$Phenotype<-"Melanoma\nPFS"

NSCLC_os<-read.csv("07.Subspecies/NSCLC_cluster_os_res.csv",header=T)
NSCLC_os_sub<-NSCLC_os[,c("X","p")]
colnames(NSCLC_os_sub)<-c("Taxa","p")
NSCLC_os_sub$Phenotype<-"NSCLC\nOS"

#prog<-rbind(resp_sub,os_sub,pfs_sub)
prog<-rbind(mel_resp_sub,NSCLC_resp_sub,RCC_resp_sub,mel_os_sub,mel_pfs_12_sub,NSCLC_os_sub)

prog.sig<-prog[prog$p<0.05,]
prog.sig<-subset(prog.sig,is.na(Taxa)==F,drop=T)

#write.table(prog.sig, "07.Subspecies/all_cluster_res.sig.tsv",
#            col.names = T, row.names = F, sep = "\t", quote = F)

## circos plot
prog.sig$Phenotype<-as.character(prog.sig$Phenotype)
prog.sig$Taxa<-as.character(prog.sig$Taxa)

prog.sig.prog  <- prog.sig$Phenotype %>%
  as.character(.) %>%
  .[!duplicated(.)]

spe_prog_count<-NULL
for (pheno in prog.sig.prog) {
  #pheno <-"C4"
  spe_prog_df<-prog.sig[prog.sig$Phenotype==pheno,]$Taxa %>%
    table %>%
    as.data.frame
  colnames(spe_prog_df)<-c("Species", "Count")
  spe_prog_df<-data.frame(pheno_info = rep(pheno,nrow(spe_prog_df)), spe_prog_df)
  spe_prog_count<-rbind(spe_prog_count,spe_prog_df)
}

spe_prog_count <- spe_prog_count[order(spe_prog_count$Count),]

spe_prog_count_species_order<-spe_prog_count %>% group_by(Species) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]
spe_prog_count_prog_order<-spe_prog_count %>% group_by(pheno_info) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]

spe_prog_count<-spe_prog_count[order(match(spe_prog_count$pheno_info, spe_prog_count_prog_order$pheno_info),decreasing = T),]

spe_prog_count_species_order_str<-as.character(spe_prog_count_species_order$Species)
spe_prog_count_prog_order_str<-as.character(spe_prog_count_prog_order$pheno_info)

##pdf("pics/F7_cluster_prog.circos.pdf", width = 14, height = 14)

tiff(file = "pics/F5_cluster_prog.circos.tiff", width =2500, height =2500, res =300)

circos.clear()
circos.par(start.degree = 0)
circos.info()

chordDiagram(spe_prog_count,annotationTrack = "grid", row.col =mycolor2_blue_yellow,
             grid.col =  c(mycolor6,
             rep('grey',length(spe_prog_count_species_order_str))),
             order = c(rev(spe_prog_count_prog_order_str),rev(spe_prog_count_species_order_str)),
             big.gap = 10,
             preAllocateTracks = list(track.margin = c(0, uh(25, "mm")), 
                                      #track.height = max(strwidth(unlist(dimnames(spe_prog_count))))
                                      track.height =0.1)
                                      )

circos.track(track.index =1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

dev.off()

## red "#BC3C29FF"
## blue "#6F99ADFF"
## yellow "#E18727FF"


### 3.2 Associations between prognosis and PCs
## 将study转换成数字型

all_basic_covar <- all_basic
covar2 <- c("read_count","cohort","age_bins","gender_num")
all_prog<-all_basic[,c("response_code","pfs_12_months","os_2_year")]

all_pc_prog.res <-lm_btw_mats(all_prog, msv_pc, all_basic_covar, covar2)
beepr::beep('mario')
save(all_pc_prog.res, file = "07.Subspecies/all_pc_prog.res.RData")

load("07.Subspecies/all_pc_prog.res.RData")
all_pc_prog.res.edge<-all_pc_prog.res$table


### 4.2 Cohort difference
all_msv_cluster<-read.table("07.Subspecies/all_msv_cluster.tsv",sep = "\t")
cluster_n<-read.table("07.Subspecies/all_cluster_n.tsv",sep = "\t", header = T)

for (i in 1:nrow(cluster_n)){
  cluster_n[,2][is.na(cluster_n[,2])] <- 0
} 

all_msv_cluster_2<-all_msv_cluster[,cluster_n$Optimum_cluster_n>=2]

all_msv_cluster_2_cohort_diff<-as.data.frame(matrix(NA,nrow = ncol(all_msv_cluster_2), ncol = 2))

for (i in 1:ncol(all_msv_cluster_2)) {
  #i<-1
  chisq_res<-chisq.test(table(all_msv_cluster_2[,i], all_basic$study))
  all_msv_cluster_2_cohort_diff[i,]<-c(colnames(all_msv_cluster_2)[i], chisq_res$p.value)
}

colnames(all_msv_cluster_2_cohort_diff) <- c('Species', 'P')
all_msv_cluster_2_cohort_diff$P<-as.numeric(all_msv_cluster_2_cohort_diff$P)
all_msv_cluster_2_cohort_diff$FDR<-p.adjust(all_msv_cluster_2_cohort_diff$P, method = 'fdr')

write.table(all_msv_cluster_2_cohort_diff, "07.Subspecies/all_msv_cluster_2_cohort_diff.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

write.csv(all_msv_cluster_2_cohort_diff, "07.Subspecies/all_msv_cluster_2_cohort_diff.csv",
            row.names = F)

############## 
library("ggmosaic")

all_gender_tbl <- table(all_basic$study,all_basic$gender) 
chisq.test(all_gender_tbl) 

#pdf("02.summary/all_gender.pdf", width = 3, height = 3)
##all_msv_cluster<-read.table("07.Subspecies/all_msv_cluster.tsv",sep = "\t",header=T)

all_basic_cluster<-cbind(all_basic, all_msv_cluster)

data<-all_basic_cluster[!is.na(all_basic_cluster$Oscillibacter.sp),]
data$Oscillibacter.sp<-as.factor(data$Oscillibacter.sp)

p_Oscillibacter.sp<-ggplot(data)+
  geom_mosaic(aes(x = product(study), fill=Oscillibacter.sp))+
  ylab("Cluster")+ggtitle("Oscillibacter.sp")+
  xlab("Cohort")+
  #scale_fill_manual(values=mycolor2_green_blue) +
  theme_tufte()+
  theme(axis.ticks.length = unit(0, "cm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.y = element_text(colour = "white"), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  scale_fill_npg()

tiff(file = "pics/FS4_Oscillibacter.sp.mosaic.tiff", width =1300, height =1100, res =300)
p_Oscillibacter.sp
dev.off()



data<-all_basic_cluster[!is.na(all_basic_cluster$B.wexlerae),]
data$B.wexlerae<-as.factor(data$B.wexlerae)

p_B.wexlerae<-ggplot(data)+
  geom_mosaic(aes(x = product(study), fill=B.wexlerae))+
  ylab("Cluster")+ggtitle("B.wexlerae")+
  xlab("Cohort")+
  #scale_fill_manual(values=mycolor2_green_blue) +
  theme_tufte()+
  theme(axis.ticks.length = unit(0, "cm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.y = element_text(colour = "white"), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  scale_fill_npg()

tiff(file = "pics/FS4_B.wexlerae.mosaic.tiff", width =1300, height =1100, res =300)
p_B.wexlerae
dev.off()



data<-all_basic_cluster[!is.na(all_basic_cluster$A.putredinis),]
data$A.putredinis<-as.factor(data$A.putredinis)

p_A.putredinis<-ggplot(data)+
  geom_mosaic(aes(x = product(study), fill=A.putredinis))+
  ylab("Cluster")+ggtitle("Escherichia coli")+
  xlab("Cohort")+
  #scale_fill_manual(values=mycolor2_green_blue) +
  theme_tufte()+
  theme(axis.ticks.length = unit(0, "cm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.y = element_text(colour = "white"), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  scale_fill_npg()

tiff(file = "pics/FS4_A.putredinis.mosaic.tiff", width =1300, height =1100, res =300)
p_A.putredinis
dev.off()




##########################################################################################
####  4 Examples

### 4.1
tsne_df<-all_msv_wsd_res[[match('E.coli', info$Short_name)]]$tsne_df
tsne_df<-data.frame(tsne_df,
                    PROG=all_prog[match(rownames(tsne_df), rownames(all_prog)),match('response_code', colnames(all_prog))])

ggplot(tsne_df, aes(x = X, y = Y)) +
    #stat_ellipse(aes(group = Cluster, fill = Cluster, color = Cluster) ,
    #             type = "norm",linetype = 2, geom = "polygon", alpha = 0.05)+
    geom_point(aes(color = PROG), alpha = 0.5)+
    ggtitle(info$Short_name[i])+
    scale_color_gradientn(colours = c(wes_palette("Darjeeling1", nrow(tsne_df), type = "continuous")))+
    #scale_fill_npg()+
    theme_void()+
    theme(legend.position = 'none',
          plot.title = element_text(face="italic"))


j<-match('E.coli', info$Short_name)
pcoa<-all_msv_wsd_res[[j]]$pcoa
pcoa<-data.frame(pcoa,
                 BA=all_basic[match(rownames(pcoa), rownames(all_basic)),match("response_code", colnames(all_basic))])


ggplot(pcoa,aes(X1,X2,fill = as.factor(all_msv_wsd_res[[j]]$tsne_df$Cluster), color = BA ))+
    geom_point(alpha = 0.5)+
    ggtitle(info$Short_name[j])+
    xlab(paste("PCo1=",round(all_msv_wsd_res[[j]]$pcoa_res$eig[1],digits = 2),"%",sep = ""))+
    ylab(paste("PCo2=",round(all_msv_wsd_res[[j]]$pcoa_res$eig[2],digits = 2),"%",sep = ""))+
  scale_color_gradientn(colours = c(wes_palette("Zissou1", 100, type = "continuous")))+
    theme(plot.title = element_text(face="italic"),
          plot.subtitle = element_text(vjust = 1), 
          plot.caption = element_text(vjust = 1), 
          axis.line.x =  element_line(),
          axis.line.y = element_line(),
          legend.position = 'right',
          legend.title = element_blank(),
          legend.key = element_rect(fill = NA), 
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          panel.background = element_rect(fill = NA))

plot(pcoa$X2,pcoa$BA)



