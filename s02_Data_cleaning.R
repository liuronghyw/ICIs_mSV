### SV data processing
### 2022-11-23
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SNV/pipeline/4_R")

source("functions.R")
##install.packages("ggmosaic")
##library("plyr")
##library("dplyr")

#####################################################################
### Read SV files
dsgv<- read.delim("00.rawData/SV/Whole_dsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")
vsgv<- read.delim("00.rawData/SV/Whole_vsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")

dsgv_anno<-read.delim("00.rawData/SV/s02.dSVs_anno.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)
vsgv_anno<-read.delim("00.rawData/SV/s03.vSVs_anno.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)

### Read database files
taxa_length <- read.table("00.rawdata/database/Species_genome_size.tsv",
                        sep = "\t", header = T,check.names = F,stringsAsFactors = F)
taxonomy<- read.csv("00.rawdata/database/representatives.genomes.taxonomy.csv",
                   sep = ",", header = T,check.names = F,stringsAsFactors = F)
taxonomy[taxonomy == ""]<-"Unknown"
colnames(taxonomy)[1]<-'X'
ncbi<-read.csv("00.rawData/database/NCBI_accession.txt", sep = "\t",header = T)
tax_relationship<-read.csv("00.rawData/database/progenome1_species_relationship.tsv",sep = "\t",header = F)


############################################################################################
#### 1 Clean SV data
#### Get clean profiles
## overlap with clinical information and abdundance file 

clinical<-read.csv("00.rawData/clinical/clinical_whole_reads.csv",header=T)
rownames(clinical)<-clinical$id

SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)

clin_dsgv_inter <- intersect(rownames(clinical),rownames(dsgv))
clin_abun_inter <- intersect(clin_dsgv_inter,colnames(SV_abun_s)[-c(1,2)])

name1<-data.frame(clinical$id)
colnames(name1)<-"id"
name1$clin<-"T"

name2<-data.frame(rownames(dsgv))
colnames(name2)<-"id"
name2$dsgv<-"T"

name<-merge(name1,name2,by.x="id",by.y="id",all.x=T,all.y=T)
write.csv(name,"test_check.csv")


clinical_sub <- clinical[clin_abun_inter,]
write.table(clinical_sub, "01.cleanData/phen/Clinical_basic_overlap.tsv",sep = '\t')

dsgv<-dsgv[clin_abun_inter,]
vsgv<-vsgv[clin_abun_inter,]

# Change SV names
colnames(dsgv) <- changeSVname(colnames(dsgv))
colnames(vsgv) <- changeSVname(colnames(vsgv))

### the ratio of NAs
dsgv_non_NA_rate<-apply(dsgv[,-1], 2, myfun<-function(x){sum(is.na(x))/1013})
vsgv_non_NA_rate<-apply(vsgv[,-1], 2, myfun<-function(x){sum(is.na(x))/1013})

write.csv(dsgv_non_NA_rate,"01.cleanData/SV/dsgv_non_NA_rate.csv")
write.csv(vsgv_non_NA_rate,"01.cleanData/SV/vsgv_non_NA_rate.csv")

## Outputs
if(!dir.exists("01.cleanData")){dir.create("01.cleanData")}
if(!dir.exists("01.cleanData/SV_all")){dir.create("01.cleanData/SV_all")}


### 2 Get name conversion table
###  Name conversion
organism<-str_replace_all(colnames(dsgv),"\\:\\d+_\\d+.*","") %>%
  .[!duplicated(.)]

Short_name<- organism %>% 
  str_replace_all('\\[','') %>%
  str_replace_all('\\]', '') %>%
  str_replace_all(' cf\\.','')

Short_name[grep(' sp\\.', organism, invert = F)] <- Short_name[grep(' sp\\.', organism, invert = F)] %>%
  str_replace_all('sp\\..*','sp')

Fst_letter<-Short_name[grep(' sp\\.', organism, invert = T)] %>%
  str_replace_all(' .*','') %>%
  str_sub(start = 1,end = 1)

Spe_name<-Short_name[grep(' sp\\.', organism, invert = T)] %>%
  str_extract_all(' .*') %>%
  str_replace_all('^ ', '') %>%
  str_replace_all(' .*', '')

Short_name[grep(' sp\\.', organism, invert = T)] <-paste(Fst_letter,'.', Spe_name, sep = '')
# write.csv(Short_name,"Short_name_test.csv")

taxa_name<-data.frame(NCBI_taxonomy_id = taxonomy$X[match(organism,taxonomy$organism)],
                      organism = as.character(organism), 
                      Short_name = as.character(Short_name), stringsAsFactors = F)

taxa_name$Short_name[match('bacterium LF-3',taxa_name$organism)]<-'bacterium LF-3'
taxa_name<-left_join(taxa_name, ncbi, by = "NCBI_taxonomy_id")

if(!dir.exists("01.cleanData/SV_info")){dir.create("01.cleanData/SV_info")}
#write.table(taxa_name, "01.cleanData/SV_info/Species_name.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

#write.csv(taxa_name, "01.cleanData/SV_info/Species_name.csv",row.names=F)
## fjy<-read.csv("01.cleanData/SV_info/article_fjy.csv",header=T)
## overlap<-merge(taxa_name,fjy,by.x="NCBI_taxonomy_id",by.y="NCBI_taxonomy_id")
## write.csv(overlap,"01.cleanData/SV_info/overlap_article_fjy.csv")

### 4.3 Get SV annotation tables
# SV annotation tables
dsgv_info_anno<-data.frame(dsgv_anno,
                           SV_ID=dsgv_anno$SV_id,
                           Taxonomy_Name = taxa_name$organism[match(str_replace_all(dsgv_anno$Taxonomy_id, '\\..*', ''),
                                                                    taxa_name$NCBI_taxonomy_id)],
                           SV_Name = changeSVname(dsgv_anno$SV_id),
                           Taxonomy_ID = dsgv_anno$Taxonomy_id,
                           SV_size = calcSVSize(dsgv_anno$SV_id))[,c(9,7,6,8,3,10,4,5)]

vsgv_info_anno<-data.frame(vsgv_anno,
                           SV_ID=vsgv_anno$SV_id,
                           Taxonomy_Name = taxa_name$organism[match(str_replace_all(vsgv_anno$Taxonomy_id, '\\..*', ''),
                                                                    taxa_name$NCBI_taxonomy_id)],
                           SV_Name = changeSVname(vsgv_anno$SV_id),
                           Taxonomy_ID = vsgv_anno$Taxonomy_id,
                           SV_size = calcSVSize(vsgv_anno$SV_id))[,c(9,7,6,8,3,10,4,5)]

write.table(dsgv_info_anno, "01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t", quote = F, col.names = T, row.names = F)
write.table(vsgv_info_anno, "01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t", quote = F, col.names = T, row.names = F)

write.csv(dsgv_info_anno, "01.cleanData/SV_info/dsgv_info_anno.csv",row.names = F)
write.csv(vsgv_info_anno, "01.cleanData/SV_info/vsgv_info_anno.csv",row.names = F)


###########################################################################################
### 4.4 Get species information table
###  Get SV number per species

species_dsgv_n<-str_replace_all(colnames(dsgv),"\\:\\d+_\\d+.*","") %>%
  table(.) %>%
  as.data.frame(.)
colnames(species_dsgv_n)<-c("Species","Deletion SVs number")
species_vsgv_n<-str_replace_all(colnames(vsgv),"\\:\\d+_\\d+.*","") %>%
  table(.) %>%
  as.data.frame(.)
colnames(species_vsgv_n)<-c("Species","Variable SVs number")

species_sgv_n<-full_join(species_dsgv_n, species_vsgv_n, by = "Species")
species_sgv_n[is.na(species_sgv_n)]<-0

NCBI_taxonomy_id<-species_sgv_n$Species %>%
  match(.,taxonomy$organism) %>%
  taxonomy$X[.]
species_sgv_n<-data.frame(NCBI_taxonomy_id, species_sgv_n)

## Get sample size per species
dsgv_infor_sample_n<-str_replace_all(colnames(dsgv),"\\:\\d+_\\d+.*","") %>%
  duplicated(.) %>%
  `!`%>%
  dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)
colnames(dsgv_infor_sample_n) <- "Sample_number"
rownames(dsgv_infor_sample_n) <- rownames(dsgv_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
dsgv_infor_sample_n<-data.frame(Species = rownames(dsgv_infor_sample_n),dsgv_infor_sample_n)

Taxonomy_name <- match(dsgv_infor_sample_n$Species,taxa_name$organism) %>%
  taxa_name$Short_name[.]
sample_n<-data.frame(Short_name=Taxonomy_name, dsgv_infor_sample_n)


### output different cohort
clinical <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")
table(clinical$dataset)

PRJEB22863<-subset(clinical,dataset=="PRJEB22863",drop=T)$id
PRJNA397906<-subset(clinical,dataset=="PRJNA397906",drop=T)$id
PRJNA541981<-subset(clinical,dataset=="PRJNA541981",drop=T)$id
PRJNA751792<-subset(clinical,dataset=="PRJNA751792",drop=T)$id
PRJNA762360<-subset(clinical,dataset=="PRJNA762360",drop=T)$id
PRJNA770295<-subset(clinical,dataset=="PRJNA770295",drop=T)$id
PRJEB43119<-subset(clinical,dataset=="PRJEB43119",drop=T)$id

mel<-subset(clinical,cancer_type=="melanoma",drop=T)$id
nsclc<-subset(clinical,cancer_type=="NSCLC",drop=T)$id
rcc<-subset(clinical,cancer_type=="RCC",drop=T)$id

PRJEB22863_dsgv<- dsgv[rownames(dsgv) %in% PRJEB22863,] 
PRJNA397906_dsgv<- dsgv[rownames(dsgv) %in% PRJNA397906,] 
PRJNA541981_dsgv<- dsgv[rownames(dsgv) %in% PRJNA541981,] 
PRJNA751792_dsgv<- dsgv[rownames(dsgv) %in% PRJNA751792,] 
PRJNA762360_dsgv<- dsgv[rownames(dsgv) %in% PRJNA762360,] 
PRJNA770295_dsgv<- dsgv[rownames(dsgv) %in% PRJNA770295,] 
PRJEB43119_dsgv<- dsgv[rownames(dsgv) %in% PRJEB43119,] 

PRJEB22863_vsgv<- vsgv[rownames(vsgv) %in% PRJEB22863,] 
PRJNA397906_vsgv<- vsgv[rownames(vsgv) %in% PRJNA397906,] 
PRJNA541981_vsgv<- vsgv[rownames(vsgv) %in% PRJNA541981,] 
PRJNA751792_vsgv<- vsgv[rownames(vsgv) %in% PRJNA751792,] 
PRJNA762360_vsgv<- vsgv[rownames(vsgv) %in% PRJNA762360,] 
PRJNA770295_vsgv<- vsgv[rownames(vsgv) %in% PRJNA770295,] 
PRJEB43119_vsgv<- vsgv[rownames(vsgv) %in% PRJEB43119,] 


## PRJNA22863 sample size per species
PRJEB22863_infor_sample_n<-str_replace_all(colnames(PRJEB22863_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJEB22863_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJEB22863_infor_sample_n) <- "PRJEB22863"
rownames(PRJEB22863_infor_sample_n) <- rownames(PRJEB22863_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJEB22863_infor_sample_n<-data.frame(Species = rownames(PRJEB22863_infor_sample_n),PRJEB22863_infor_sample_n)


## PRJNA397906 sample size per species
PRJNA397906_infor_sample_n<-str_replace_all(colnames(PRJNA397906_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJNA397906_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJNA397906_infor_sample_n) <- "PRJNA397906"
rownames(PRJNA397906_infor_sample_n) <- rownames(PRJNA397906_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJNA397906_infor_sample_n<-data.frame(Species = rownames(PRJNA397906_infor_sample_n),PRJNA397906_infor_sample_n)


## PRJNA541981 sample size per species
PRJNA541981_infor_sample_n<-str_replace_all(colnames(PRJNA541981_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJNA541981_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJNA541981_infor_sample_n) <- "PRJNA541981"
rownames(PRJNA541981_infor_sample_n) <- rownames(PRJNA541981_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJNA541981_infor_sample_n<-data.frame(Species = rownames(PRJNA541981_infor_sample_n),PRJNA541981_infor_sample_n)


## PRJNA751792 sample size per species
PRJNA751792_infor_sample_n<-str_replace_all(colnames(PRJNA751792_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJNA751792_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJNA751792_infor_sample_n) <- "PRJNA751792"
rownames(PRJNA751792_infor_sample_n) <- rownames(PRJNA751792_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJNA751792_infor_sample_n<-data.frame(Species = rownames(PRJNA751792_infor_sample_n),PRJNA751792_infor_sample_n)


## PRJNA762360 sample size per species
PRJNA762360_infor_sample_n<-str_replace_all(colnames(PRJNA762360_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJNA762360_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJNA762360_infor_sample_n) <- "PRJNA762360"
rownames(PRJNA762360_infor_sample_n) <- rownames(PRJNA762360_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJNA762360_infor_sample_n<-data.frame(Species = rownames(PRJNA762360_infor_sample_n),PRJNA762360_infor_sample_n)


## PRJNA770295 sample size per species
PRJNA770295_infor_sample_n<-str_replace_all(colnames(PRJNA770295_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJNA770295_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJNA770295_infor_sample_n) <- "PRJNA770295"
rownames(PRJNA770295_infor_sample_n) <- rownames(PRJNA770295_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJNA770295_infor_sample_n<-data.frame(Species = rownames(PRJNA770295_infor_sample_n),PRJNA770295_infor_sample_n)


## PRJEB43119 sample size per species
PRJEB43119_infor_sample_n<-str_replace_all(colnames(PRJEB43119_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJEB43119_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJEB43119_infor_sample_n) <- "PRJEB43119"
rownames(PRJEB43119_infor_sample_n) <- rownames(PRJEB43119_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJEB43119_infor_sample_n<-data.frame(Species = rownames(PRJEB43119_infor_sample_n),PRJEB43119_infor_sample_n)


## merge sample size from different dataset 

info_sample<-merge(PRJEB22863_infor_sample_n,PRJNA397906_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,PRJNA541981_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,PRJNA751792_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,PRJNA762360_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,PRJNA770295_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,PRJEB43119_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,sample_n,by.x="Species",by.y="Species")

infor_sample_n <- info_sample[,c(1,9,10,seq(2,8))]

## Merge sample size and SV number information
species_sample_n<-dplyr::full_join(species_sgv_n,infor_sample_n, by = "Species")
taxa_length$Species<-str_replace_all(taxa_length$Species, '\\..*', '')
species_sample_n$NCBI_taxonomy_id<-as.character(species_sample_n$NCBI_taxonomy_id)
species_sample_n<-dplyr::left_join(species_sample_n, taxa_length, by = c("NCBI_taxonomy_id"="Species"))
species_sample_n<-data.frame(species_sample_n,
                             SVs.number = species_sample_n[,3]+species_sample_n[,4])

## Merge all information
Informative_species_information <- match(species_sample_n$NCBI_taxonomy_id, taxonomy$X)%>%
  taxonomy[.,] %>%
  cbind(.,species_sample_n)

taxa_name_short<-taxa_name[,c(2,4)]

info <- full_join(Informative_species_information[,-11],
                                             taxa_name_short,
                                             by = 'organism')[,c(1:9,13,24,11,12,23,seq(15,21),14,22)]


colnames(info)[c(1,10:23)]<-c("NCBI_taxonomy_id","Short_name","NCBI_bioproject_accession", 
"Deletion_SVs_number", "Variable_SVs_number","SVs_number","PRJEB22863_sample_number",
"PRJNA397906_sample_number","PRJNA541981_sample_number",
"PRJNA751792_sample_number","PRJNA762360_sample_number","PRJNA770295_sample_number",
"PRJEB43119_sample_number","Total_samples_number", "Length")

write.table(info, "01.cleanData/SV_info/Informative_species_information.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info, "01.cleanData/SV_info/Informative_species_information.csv", row.names = F)


##########################################################################
### 筛选species， 从物种丰度和存在的样本比率两方面  1014 个样本，在5%的样本中出现的物种

info<-subset(info,organism!="Clostridiales bacterium VE202-14" 
   & organism!="Clostridium sp. KLE 1755" & organism!="Oscillibacter sp. KLE 1728"
   & organism!="Phascolarctobacterium sp. CAG:266"
   & organism!="Lachnospiraceae bacterium 3_1_57FAA_CT1"
   & organism!="Lachnospiraceae bacterium 9_1_43BFAA"
   & organism!="Sutterella wadsworthensis 2_1_59BFAA",drop=T)

info_sub<-subset(info,Total_samples_number>51,drop=T)

#########################################################################
###  5 Clean species abundance data 

# Read species abundance files
SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)
# Basic
all_basic <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

sv_spe_taxid<-info_sub$NCBI_taxonomy_id
sv_spe_taxid[sv_spe_taxid==245018]<-649756
sv_spe_taxid[sv_spe_taxid==1118061]<-2585118

# relative abundance of species detected with SVs
#SV_abun_s[,-c(1:2)]<-apply(SV_abun_s[,-c(1:2)], 2, myfun<-function(x){x/sum(x)})

SV_abun_s <- sv_spe_taxid %>% 
  match(., tax_relationship$V2) %>% 
  tax_relationship$V1[.] %>% 
  match(., SV_abun_s$taxonomy_id) %>% 
  SV_abun_s[., -c(1:2)] %>% 
  t %>%
  as.data.frame

colnames(SV_abun_s)<-info_sub$organism
#rownames(SV_abun_s)<-str_replace_all(rownames(SV_abun_s), ".S.bracken", "")

SV_abun_s<-SV_abun_s[match(rownames(all_basic),rownames(SV_abun_s)),]
rownames(SV_abun_s)<-rownames(all_basic)

spe_abun_mean<-apply(SV_abun_s,2,mean,na.rm=T)
write.csv(spe_abun_mean,"01.cleanData/spe_abun_mean.csv")

nr<-nrow(SV_abun_s)
nc<-ncol(SV_abun_s)

zero_ratio<-c()

for(j in 1:nc){
   count<-0
   for(i in 1:nr){
     if(SV_abun_s[i,j]==0){count=count+1}
     }
   zero_ratio[j]<-count/nr
  }

zero_ratio<-data.frame(zero_ratio)
zero_ratio$organism<-colnames(SV_abun_s)
#write.csv(zero_ratio,"01.cleanData/spe_abun_zero_ratio.csv")

zero_ratio<-subset(zero_ratio,zero_ratio<0.4,drop=T)

info_final<-merge(info_sub,zero_ratio,by.x="organism",by.y="organism")[,-24]
head(info_final)

write.table(info_final, "01.cleanData/SV_info/Informative_species_information_final.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info_final, "01.cleanData/SV_info/Informative_species_information_final.csv", row.names = F)



############################### abundance final

# Read species abundance files
SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)
# Basic
all_basic <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

sv_spe_taxid<-info_final$NCBI_taxonomy_id
sv_spe_taxid[sv_spe_taxid==245018]<-649756
sv_spe_taxid[sv_spe_taxid==1118061]<-2585118

# relative abundance of species detected with SVs
#SV_abun_s[,-c(1:2)]<-apply(SV_abun_s[,-c(1:2)], 2, myfun<-function(x){x/sum(x)})

SV_abun_s <- sv_spe_taxid %>% 
  match(., tax_relationship$V2) %>% 
  tax_relationship$V1[.] %>% 
  match(., SV_abun_s$taxonomy_id) %>% 
  SV_abun_s[., -c(1:2)] %>% 
  t %>%
  as.data.frame

colnames(SV_abun_s)<-info_final$organism
#rownames(SV_abun_s)<-str_replace_all(rownames(SV_abun_s), ".S.bracken", "")

SV_abun_s<-SV_abun_s[match(rownames(all_basic),rownames(SV_abun_s)),]
rownames(SV_abun_s)<-rownames(all_basic)

# outputs
if(!dir.exists("01.cleanData/mbio_all")){dir.create("01.cleanData/mbio_all")}
write.table(SV_abun_s, "01.cleanData/mbio_all/SV_species_abun.tsv",sep = '\t')

mean(rowSums(SV_abun_s,na.rm = T),na.rm = T)  # Mean of total relative abundance of species detected with SVs: 0.825154
rowSums(SV_abun_s,na.rm = T)[rowSums(SV_abun_s,na.rm = T)!=0] %>% min # Minimum of total relative abundance of species detected with SVs: 0.43729
max(rowSums(SV_abun_s,na.rm = T),na.rm = T)   # Maximum of total relative abundance of species detected with SVs: 0.94705


###########################################################################
## 不同的数据集

PRJEB22863_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJEB22863,] 
PRJNA397906_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJNA397906,] 
PRJNA541981_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJNA541981,] 
PRJNA751792_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJNA751792,] 
PRJNA762360_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJNA762360,] 
PRJNA770295_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJNA770295,] 
PRJEB43119_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJEB43119,] 

write.table(PRJEB22863_abun, "01.cleanData/mbio_all/PRJEB22863_SV_species_abun.tsv",sep = '\t')
write.table(PRJNA397906_abun, "01.cleanData/mbio_all/PRJNA397906_SV_species_abun.tsv",sep = '\t')
write.table(PRJNA541981_abun, "01.cleanData/mbio_all/PRJNA541981_SV_species_abun.tsv",sep = '\t')
write.table(PRJNA751792_abun, "01.cleanData/mbio_all/PRJNA751792_SV_species_abun.tsv",sep = '\t')
write.table(PRJNA762360_abun, "01.cleanData/mbio_all/PRJNA762360_SV_species_abun.tsv",sep = '\t')
write.table(PRJNA770295_abun, "01.cleanData/mbio_all/PRJNA770295_SV_species_abun.tsv",sep = '\t')
write.table(PRJEB43119_abun, "01.cleanData/mbio_all/PRJEB43119_SV_species_abun.tsv",sep = '\t')

mean(rowSums(SV_abun_s,na.rm = T),na.rm = T)  # Mean of total relative abundance of species detected with SVs: 0.825154
rowSums(SV_abun_s,na.rm = T)[rowSums(SV_abun_s,na.rm = T)!=0] %>% min # Minimum of total relative abundance of species detected with SVs: 0.43729
max(rowSums(SV_abun_s,na.rm = T),na.rm = T)   # Maximum of total relative abundance of species detected with SVs: 0.94705


###################################################################################################
###  select dsgv and vsgv information

info_final<-read.delim("01.cleanData/SV_info/Informative_species_information_final.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)
organism<-info_final$organism

select_column<-c()
for(i in 1:ncol(dsgv)){
  name<-str_replace_all(colnames(dsgv[i]),"\\:\\d+_\\d+.*","")
  if(name %in% organism){select_column<-c(select_column,i)}
  }
dsgv_sub<-dsgv[,select_column]

select_column<-c()
for(i in 1:ncol(vsgv)){
  name<-str_replace_all(colnames(vsgv[i]),"\\:\\d+_\\d+.*","")
  if(name %in% organism){select_column<-c(select_column,i)}
  }
vsgv_sub<-vsgv[,select_column]

### 当variation同时出现在dsv和vsv中，去掉vsv信息

dsv_id<-colnames(dsgv_sub)
vsv_id<-colnames(vsgv_sub)
overlap_id<- intersect(dsv_id,vsv_id)
vsv_id_left<-setdiff(vsv_id,overlap_id)

vsgv_sub<-vsgv_sub[,colnames(vsgv_sub) %in% vsv_id_left]
write.csv(overlap_id,"check_overlap_id.csv")

write.table(dsgv_sub,"01.cleanData/SV_all/deletionStructuralVariation_all.tsv",sep = '\t')
write.table(vsgv_sub,"01.cleanData/SV_all/variableStructuralVariation_all.tsv",sep = '\t')
save(dsgv_sub, file = "01.cleanData/SV_all/dsgv.RData")
save(vsgv_sub, file = "01.cleanData/SV_all/vsgv.RData")


### output different cohort
clinical <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")
table(clinical$dataset)

PRJEB22863<-subset(clinical,dataset=="PRJEB22863",drop=T)$id
PRJEB22863_NSCLC<-subset(clinical,dataset=="PRJEB22863" & cancer_type=="NSCLC",drop=T)$id
PRJEB22863_RCC<-subset(clinical,dataset=="PRJEB22863" & cancer_type=="RCC",drop=T)$id

PRJNA397906<-subset(clinical,dataset=="PRJNA397906",drop=T)$id
PRJNA541981<-subset(clinical,dataset=="PRJNA541981",drop=T)$id
PRJNA751792<-subset(clinical,dataset=="PRJNA751792",drop=T)$id
PRJNA762360<-subset(clinical,dataset=="PRJNA762360",drop=T)$id
PRJNA770295<-subset(clinical,dataset=="PRJNA770295",drop=T)$id
PRJEB43119<-subset(clinical,dataset=="PRJEB43119",drop=T)$id

PRJEB22863_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJEB22863,] 
PRJEB22863_NSCLC_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJEB22863_NSCLC,] 
PRJEB22863_RCC_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJEB22863_RCC,] 
PRJNA397906_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJNA397906,] 
PRJNA541981_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJNA541981,] 
PRJNA751792_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJNA751792,] 
PRJNA762360_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJNA762360,] 
PRJNA770295_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJNA770295,] 
PRJEB43119_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJEB43119,] 

PRJEB22863_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJEB22863,] 
PRJEB22863_NSCLC_vsgv<- dsgv_sub[rownames(vsgv_sub) %in% PRJEB22863_NSCLC,] 
PRJEB22863_RCC_vsgv<- dsgv_sub[rownames(vsgv_sub) %in% PRJEB22863_RCC,] 
PRJNA397906_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJNA397906,] 
PRJNA541981_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJNA541981,] 
PRJNA751792_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJNA751792,] 
PRJNA762360_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJNA762360,] 
PRJNA770295_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJNA770295,] 
PRJEB43119_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJEB43119,] 

mel_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% mel,] 
nsclc_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% nsclc,] 
rcc_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% rcc,] 

mel_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% mel,] 
nsclc_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% nsclc,] 
rcc_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% rcc,] 

write.table(PRJEB22863_dsgv,"01.cleanData/SV_all/deletionStructuralVariation_PRJEB22863.tsv",sep = '\t')
write.table(PRJEB22863_vsgv,"01.cleanData/SV_all/variableStructuralVariation_PRJEB22863.tsv",sep = '\t')
save(PRJEB22863_dsgv, file = "01.cleanData/SV_all/dsgv_PRJEB22863.RData")
save(PRJEB22863_vsgv, file = "01.cleanData/SV_all/vsgv_PRJEB22863.RData")

write.table(PRJEB22863_NSCLC_dsgv,"01.cleanData/SV_all/deletionStructuralVariation_NSCLC_PRJEB22863.tsv",sep = '\t')
write.table(PRJEB22863_NSCLC_vsgv,"01.cleanData/SV_all/variableStructuralVariation_NSCLC_PRJEB22863.tsv",sep = '\t')
save(PRJEB22863_NSCLC_dsgv, file = "01.cleanData/SV_all/dsgv_NSCLC_PRJEB22863.RData")
save(PRJEB22863_NSCLC_vsgv, file = "01.cleanData/SV_all/vsgv_NSCLC_PRJEB22863.RData")

write.table(PRJEB22863_RCC_dsgv,"01.cleanData/SV_all/deletionStructuralVariation_RCC_PRJEB22863.tsv",sep = '\t')
write.table(PRJEB22863_RCC_vsgv,"01.cleanData/SV_all/variableStructuralVariation_RCC_PRJEB22863.tsv",sep = '\t')
save(PRJEB22863_RCC_dsgv, file = "01.cleanData/SV_all/dsgv_RCC_PRJEB22863.RData")
save(PRJEB22863_RCC_vsgv, file = "01.cleanData/SV_all/vsgv_RCC_PRJEB22863.RData")

write.table(PRJNA397906_dsgv,"01.cleanData/SV_all/deletionStructuralVariation_PRJNA397906.tsv",sep = '\t')
write.table(PRJNA397906_vsgv,"01.cleanData/SV_all/variableStructuralVariation_PRJNA397906.tsv",sep = '\t')
save(PRJNA397906_dsgv, file = "01.cleanData/SV_all/dsgv_PRJNA397906.RData")
save(PRJNA397906_vsgv, file = "01.cleanData/SV_all/vsgv_PRJNA397906.RData")

write.table(PRJNA541981_dsgv,"01.cleanData/SV_all/deletionStructuralVariation_PRJNA541981.tsv",sep = '\t')
write.table(PRJNA541981_vsgv,"01.cleanData/SV_all/variableStructuralVariation_PRJNA541981.tsv",sep = '\t')
save(PRJNA541981_dsgv, file = "01.cleanData/SV_all/dsgv_PRJNA541981.RData")
save(PRJNA541981_vsgv, file = "01.cleanData/SV_all/vsgv_PRJNA541981.RData")

write.table(PRJNA751792_dsgv,"01.cleanData/SV_all/deletionStructuralVariation_PRJNA751792.tsv",sep = '\t')
write.table(PRJNA751792_vsgv,"01.cleanData/SV_all/variableStructuralVariation_PRJNA751792.tsv",sep = '\t')
save(PRJNA751792_dsgv, file = "01.cleanData/SV_all/dsgv_PRJNA751792.RData")
save(PRJNA751792_vsgv, file = "01.cleanData/SV_all/vsgv_PRJNA751792.RData")

write.table(PRJNA762360_dsgv,"01.cleanData/SV_all/deletionStructuralVariation_PRJNA762360.tsv",sep = '\t')
write.table(PRJNA762360_vsgv,"01.cleanData/SV_all/variableStructuralVariation_PRJNA762360.tsv",sep = '\t')
save(PRJNA762360_dsgv, file = "01.cleanData/SV_all/dsgv_PRJNA762360.RData")
save(PRJNA762360_vsgv, file = "01.cleanData/SV_all/vsgv_PRJNA762360.RData")

write.table(PRJNA770295_dsgv,"01.cleanData/SV_all/deletionStructuralVariation_PRJNA770295.tsv",sep = '\t')
write.table(PRJNA770295_vsgv,"01.cleanData/SV_all/variableStructuralVariation_PRJNA770295.tsv",sep = '\t')
save(PRJNA770295_dsgv, file = "01.cleanData/SV_all/dsgv_PRJNA770295.RData")
save(PRJNA770295_vsgv, file = "01.cleanData/SV_all/vsgv_PRJNA770295.RData")

write.table(PRJEB43119_dsgv,"01.cleanData/SV_all/deletionStructuralVariation_PRJEB43119.tsv",sep = '\t')
write.table(PRJEB43119_vsgv,"01.cleanData/SV_all/variableStructuralVariation_PRJEB43119.tsv",sep = '\t')
save(PRJEB43119_dsgv, file = "01.cleanData/SV_all/dsgv_PRJEB43119.RData")
save(PRJEB43119_vsgv, file = "01.cleanData/SV_all/vsgv_PRJEB43119.RData")


################################ 4.5 Get distance matrices
#### 4.5.1 All samples

## msv (vsv+dsv) distance
sgv<-cbind(vsgv_sub, dsgv_sub)
all_shared_sv_dis<-shared_sv_dis(sgv)
save(all_shared_sv_dis, file = "01.cleanData/SV_all/all_shared_sv_dis.RData")

## SV distance matrices of all species
all_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-7
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-vsgv_sub[,grep(spe_name,colnames(vsgv_sub))]
  dsgv_i<-dsgv_sub[,grep(spe_name,colnames(dsgv_sub))]
  all_msv_i<-cbind(vsgv_i,dsgv_i) %>%
  na.omit(.)

  all_msv_i = all_msv_i[rowSums(all_msv_i[])>0,]

  all_msv_dist_i <- as.matrix(vegdist(as.data.frame(all_msv_i),method = "canberra"))
  all_msv_dist[[i]]<-all_msv_dist_i
}

names(all_msv_dist)<-paste('msv_',info_final$organism, sep = '')
all_msv_dist_std <- lapply(all_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

if(!dir.exists("01.cleanData/SV_all/distMat")){dir.create("01.cleanData/SV_all/distMat")}

#all_msv_dist[which(names(all_msv_dist)%in%c("msv_Acidaminococcus intestini RyC-MR95"))]<-NULL
#all_msv_dist_std[which(names(all_msv_dist_std)%in%c("msv_Acidaminococcus intestini RyC-MR95"))]<-NULL

save(all_msv_dist, file = "01.cleanData/SV_all/distMat/all_msv_dist.RData")
save(all_msv_dist_std, file = "01.cleanData/SV_all/distMat/all_msv_dist_std.RData")


################################# PRJEB22863 samples
## msv (vsv+dsv) distance
PRJEB22863_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-PRJEB22863_vsgv[,grep(spe_name,colnames(PRJEB22863_vsgv))]
  dsgv_i<-PRJEB22863_dsgv[,grep(spe_name,colnames(PRJEB22863_dsgv))]
  PRJEB22863_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)
  
  PRJEB22863_msv_dist_i <- as.matrix(vegdist(as.data.frame(PRJEB22863_msv_i),method = "canberra"))
  PRJEB22863_msv_dist[[i]]<-PRJEB22863_msv_dist_i
}

names(PRJEB22863_msv_dist)<-paste('msv_',info_final$organism, sep = '')
PRJEB22863_msv_dist_std <- lapply(PRJEB22863_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(PRJEB22863_msv_dist, file = "01.cleanData/SV_all/distMat/PRJEB22863_msv_dist.RData")
save(PRJEB22863_msv_dist_std, file = "01.cleanData/SV_all/distMat/PRJEB22863_msv_dist_std.RData")


################################# PRJEB22863 (NSCLC) samples
## msv (vsv+dsv) distance
PRJEB22863_NSCLC_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-PRJEB22863_NSCLC_vsgv[,grep(spe_name,colnames(PRJEB22863_NSCLC_vsgv))]
  dsgv_i<-PRJEB22863_NSCLC_dsgv[,grep(spe_name,colnames(PRJEB22863_NSCLC_dsgv))]
  PRJEB22863_NSCLC_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)

  PRJEB22863_NSCLC_msv_i = PRJEB22863_NSCLC_msv_i[rowSums(PRJEB22863_NSCLC_msv_i[])>0,]
  
  PRJEB22863_NSCLC_msv_dist_i <- as.matrix(vegdist(as.data.frame(PRJEB22863_NSCLC_msv_i),method = "canberra"))
  PRJEB22863_NSCLC_msv_dist[[i]]<-PRJEB22863_NSCLC_msv_dist_i
}

names(PRJEB22863_NSCLC_msv_dist)<-paste('msv_',info_final$organism, sep = '')
PRJEB22863_NSCLC_msv_dist_std <- lapply(PRJEB22863_NSCLC_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(PRJEB22863_NSCLC_msv_dist, file = "01.cleanData/SV_all/distMat/PRJEB22863_NSCLC_msv_dist.RData")
save(PRJEB22863_NSCLC_msv_dist_std, file = "01.cleanData/SV_all/distMat/PRJEB22863_NSCLC_msv_dist_std.RData")



################################# PRJEB22863 (RCC) samples
## msv (vsv+dsv) distance
PRJEB22863_RCC_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-PRJEB22863_RCC_vsgv[,grep(spe_name,colnames(PRJEB22863_RCC_vsgv))]
  dsgv_i<-PRJEB22863_RCC_dsgv[,grep(spe_name,colnames(PRJEB22863_RCC_dsgv))]
  PRJEB22863_RCC_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)

  PRJEB22863_RCC_msv_i = PRJEB22863_RCC_msv_i[rowSums(PRJEB22863_RCC_msv_i[])>0,]
  
  PRJEB22863_RCC_msv_dist_i <- as.matrix(vegdist(as.data.frame(PRJEB22863_RCC_msv_i),method = "canberra"))
  PRJEB22863_RCC_msv_dist[[i]]<-PRJEB22863_RCC_msv_dist_i
}

names(PRJEB22863_RCC_msv_dist)<-paste('msv_',info_final$organism, sep = '')
PRJEB22863_RCC_msv_dist_std <- lapply(PRJEB22863_RCC_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(PRJEB22863_RCC_msv_dist, file = "01.cleanData/SV_all/distMat/PRJEB22863_RCC_msv_dist.RData")
save(PRJEB22863_RCC_msv_dist_std, file = "01.cleanData/SV_all/distMat/PRJEB22863_RCC_msv_dist_std.RData")


################################# PRJNA397906 samples
## msv (vsv+dsv) distance
PRJNA397906_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-PRJNA397906_vsgv[,grep(spe_name,colnames(PRJNA397906_vsgv))]
  dsgv_i<-PRJNA397906_dsgv[,grep(spe_name,colnames(PRJNA397906_dsgv))]
  PRJNA397906_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)
  PRJNA397906_msv_i = PRJNA397906_msv_i[rowSums(PRJNA397906_msv_i[])>0,]
  
  PRJNA397906_msv_dist_i <- as.matrix(vegdist(as.data.frame(PRJNA397906_msv_i),method = "canberra"))
  PRJNA397906_msv_dist[[i]]<-PRJNA397906_msv_dist_i
}

names(PRJNA397906_msv_dist)<-paste('msv_',info_final$organism, sep = '')
PRJNA397906_msv_dist_std <- lapply(PRJNA397906_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(PRJNA397906_msv_dist, file = "01.cleanData/SV_all/distMat/PRJNA397906_msv_dist.RData")
save(PRJNA397906_msv_dist_std, file = "01.cleanData/SV_all/distMat/PRJNA397906_msv_dist_std.RData")


################################# PRJNA541981 samples
## msv (vsv+dsv) distance
PRJNA541981_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-PRJNA541981_vsgv[,grep(spe_name,colnames(PRJNA541981_vsgv))]
  dsgv_i<-PRJNA541981_dsgv[,grep(spe_name,colnames(PRJNA541981_dsgv))]
  PRJNA541981_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)

  PRJNA541981_msv_i = PRJNA541981_msv_i[rowSums(PRJNA541981_msv_i[])>0,]
  
  PRJNA541981_msv_dist_i <- as.matrix(vegdist(as.data.frame(PRJNA541981_msv_i),method = "canberra"))
  PRJNA541981_msv_dist[[i]]<-PRJNA541981_msv_dist_i
}

names(PRJNA541981_msv_dist)<-paste('msv_',info_final$organism, sep = '')
PRJNA541981_msv_dist_std <- lapply(PRJNA541981_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(PRJNA541981_msv_dist, file = "01.cleanData/SV_all/distMat/PRJNA541981_msv_dist.RData")
save(PRJNA541981_msv_dist_std, file = "01.cleanData/SV_all/distMat/PRJNA541981_msv_dist_std.RData")


################################# PRJNA751792 samples
## msv (vsv+dsv) distance
PRJNA751792_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-PRJNA751792_vsgv[,grep(spe_name,colnames(PRJNA751792_vsgv))]
  dsgv_i<-PRJNA751792_dsgv[,grep(spe_name,colnames(PRJNA751792_dsgv))]
  PRJNA751792_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)

  PRJNA751792_msv_i = PRJNA751792_msv_i[rowSums(PRJNA751792_msv_i[])>0,]
  
  PRJNA751792_msv_dist_i <- as.matrix(vegdist(as.data.frame(PRJNA751792_msv_i),method = "canberra"))
  PRJNA751792_msv_dist[[i]]<-PRJNA751792_msv_dist_i
}

names(PRJNA751792_msv_dist)<-paste('msv_',info_final$organism, sep = '')
PRJNA751792_msv_dist_std <- lapply(PRJNA751792_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(PRJNA751792_msv_dist, file = "01.cleanData/SV_all/distMat/PRJNA751792_msv_dist.RData")
save(PRJNA751792_msv_dist_std, file = "01.cleanData/SV_all/distMat/PRJNA751792_msv_dist_std.RData")



################################# PRJNA762360 samples
## msv (vsv+dsv) distance
PRJNA762360_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-PRJNA762360_vsgv[,grep(spe_name,colnames(PRJNA762360_vsgv))]
  dsgv_i<-PRJNA762360_dsgv[,grep(spe_name,colnames(PRJNA762360_dsgv))]
  PRJNA762360_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)
  
  PRJNA762360_msv_i = PRJNA762360_msv_i[rowSums(PRJNA762360_msv_i[])>0,]

  PRJNA762360_msv_dist_i <- as.matrix(vegdist(as.data.frame(PRJNA762360_msv_i),method = "canberra"))
  PRJNA762360_msv_dist[[i]]<-PRJNA762360_msv_dist_i
}

names(PRJNA762360_msv_dist)<-paste('msv_',info_final$organism, sep = '')
PRJNA762360_msv_dist_std <- lapply(PRJNA762360_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(PRJNA762360_msv_dist, file = "01.cleanData/SV_all/distMat/PRJNA762360_msv_dist.RData")
save(PRJNA762360_msv_dist_std, file = "01.cleanData/SV_all/distMat/PRJNA762360_msv_dist_std.RData")



################################# PRJNA770295 samples
## msv (vsv+dsv) distance
PRJNA770295_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-PRJNA770295_vsgv[,grep(spe_name,colnames(PRJNA770295_vsgv))]
  dsgv_i<-PRJNA770295_dsgv[,grep(spe_name,colnames(PRJNA770295_dsgv))]
  PRJNA770295_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)
  
  PRJNA770295_msv_i = PRJNA770295_msv_i[rowSums(PRJNA770295_msv_i[])>0,]

  PRJNA770295_msv_dist_i <- as.matrix(vegdist(as.data.frame(PRJNA770295_msv_i),method = "canberra"))
  PRJNA770295_msv_dist[[i]]<-PRJNA770295_msv_dist_i
}

names(PRJNA770295_msv_dist)<-paste('msv_',info_final$organism, sep = '')
PRJNA770295_msv_dist_std <- lapply(PRJNA770295_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(PRJNA770295_msv_dist, file = "01.cleanData/SV_all/distMat/PRJNA770295_msv_dist.RData")
save(PRJNA770295_msv_dist_std, file = "01.cleanData/SV_all/distMat/PRJNA770295_msv_dist_std.RData")


################################# PRJEB43119 samples
## msv (vsv+dsv) distance
PRJEB43119_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-PRJEB43119_vsgv[,grep(spe_name,colnames(PRJEB43119_vsgv))]
  dsgv_i<-PRJEB43119_dsgv[,grep(spe_name,colnames(PRJEB43119_dsgv))]
  PRJEB43119_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)

  PRJEB43119_msv_i = PRJEB43119_msv_i[rowSums(PRJEB43119_msv_i[])>0,]
  
  PRJEB43119_msv_dist_i <- as.matrix(vegdist(as.data.frame(PRJEB43119_msv_i),method = "canberra"))
  PRJEB43119_msv_dist[[i]]<-PRJEB43119_msv_dist_i
}

names(PRJEB43119_msv_dist)<-paste('msv_',info_final$organism, sep = '')
PRJEB43119_msv_dist_std <- lapply(PRJEB43119_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(PRJEB43119_msv_dist, file = "01.cleanData/SV_all/distMat/PRJEB43119_msv_dist.RData")
save(PRJEB43119_msv_dist_std, file = "01.cleanData/SV_all/distMat/PRJEB43119_msv_dist_std.RData")


################################# mel samples
## msv (vsv+dsv) distance
mel_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-mel_vsgv[,grep(spe_name,colnames(mel_vsgv))]
  dsgv_i<-mel_dsgv[,grep(spe_name,colnames(mel_dsgv))]
  mel_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)
  
  mel_msv_dist_i <- as.matrix(vegdist(as.data.frame(mel_msv_i),method = "canberra"))
  mel_msv_dist[[i]]<-mel_msv_dist_i
}

names(mel_msv_dist)<-paste('msv_',info_final$organism, sep = '')
mel_msv_dist_std <- lapply(mel_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(mel_msv_dist, file = "01.cleanData/SV_all/distMat/mel_msv_dist.RData")
save(mel_msv_dist_std, file = "01.cleanData/SV_all/distMat/mel_msv_dist_std.RData")



################################# nsclc samples
## msv (vsv+dsv) distance
nsclc_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-nsclc_vsgv[,grep(spe_name,colnames(nsclc_vsgv))]
  dsgv_i<-nsclc_dsgv[,grep(spe_name,colnames(nsclc_dsgv))]
  nsclc_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)
  
  nsclc_msv_dist_i <- as.matrix(vegdist(as.data.frame(nsclc_msv_i),method = "canberra"))
  nsclc_msv_dist[[i]]<-nsclc_msv_dist_i
}

names(nsclc_msv_dist)<-paste('msv_',info_final$organism, sep = '')
nsclc_msv_dist_std <- lapply(nsclc_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(nsclc_msv_dist, file = "01.cleanData/SV_all/distMat/nsclc_msv_dist.RData")
save(nsclc_msv_dist_std, file = "01.cleanData/SV_all/distMat/nsclc_msv_dist_std.RData")


################################# rcc samples
## msv (vsv+dsv) distance
rcc_msv_dist<-NULL

for (i in c(1:nrow(info_final))){
  #i<-16
  file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  vsgv_i<-rcc_vsgv[,grep(spe_name,colnames(rcc_vsgv))]
  dsgv_i<-rcc_dsgv[,grep(spe_name,colnames(rcc_dsgv))]
  rcc_msv_i<-cbind(vsgv_i,dsgv_i) %>%
    na.omit(.)
  
  rcc_msv_dist_i <- as.matrix(vegdist(as.data.frame(rcc_msv_i),method = "canberra"))
  rcc_msv_dist[[i]]<-rcc_msv_dist_i
}

names(rcc_msv_dist)<-paste('msv_',info_final$organism, sep = '')
rcc_msv_dist_std <- lapply(rcc_msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

save(rcc_msv_dist, file = "01.cleanData/SV_all/distMat/rcc_msv_dist.RData")
save(rcc_msv_dist_std, file = "01.cleanData/SV_all/distMat/rcc_msv_dist_std.RData")



