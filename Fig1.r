###########Fig1#######################
##Author: tzh
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(ggplot2)
library(ggsci)
library(dendextend)
library(ggraph)
library(igraph)
library(funkyheatmap)
library(dplyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(tidyverse)
library(reshape2)

dataset<-read.xlsx("ALLsample_info1.xlsx",sheet=2)

rownames(dataset)<-paste0("CCID",dataset$CCID)
dataset<-dataset[,-c(1,2)]
dataset<-dataset[,-c(3,4)]

dataset$Data_type<-factor(dataset$Data_type,levels=c("Transcriptome","Proteome","Metabolomics"))
dataset<-dataset[order(dataset$TaxonomyID,dataset$Data_type,dataset$Cancer_type,
                       dataset$Tissue),]
dataset$TaxonomyID<-as.character(dataset$TaxonomyID)
#dataset_t <- as.data.frame(t(dataset))

col_list = list(
  TaxonomyID = c("Human" = "#003C67FF", "Mouse" = "#0073C2FF"),
  Tissue = c("raM"=  "#E64B35FF","iWAT"="#00A087FF","Blood"="#DC0000FF",
             "Colon"="#3C5488FF","xM"="#F39B7FFF","pWAT"="#8491B4FF",
             "Liver"="#91D1C2FF","gM"="#4DBBD5FF","taM"="#7E6148FF",
             "eWAT"="#B09C85FF", "qM"="#374E55FF","Hypothalamus"="#DF8F44FF",
             "hM" ="#00A1D5FF",  "dM"="#B24745FF", "BAT"= "#79AF97FF",
             "Cerebellum"="#6A6599FF","Heart"="#800000FF","Hippocampus"="#80796BFF",
             "Neocortex"="#FFA319FF", "sM"="#8A9045FF","cM"="#155F83FF",
             "Gastric"="#C16622FF","Caecal"="#58593FFF","VenaCava"="#350E20FF",    
             "PortalVein"="#3D3B25FF"),
  SampleNum= colorRamp2(c(2, 226), c("#D4E1F5", "#0C2340")),
  Cancer_type = c("pancreatic cancer"="#FED439FF", #"colon cancer"="#709AE1FF",           
                  "lung cancer"="#D2AF81FF", #"gastrointestinal cancer"="#D5E4A2FF",
                  "colorectal cancer"= "#709AE1FF", "Multiple cancers"="#8A9197FF",       
                  "gastric cancer"= "#D5E4A2FF", "Breast Cancer"= "#1A9993FF",         
                  "melanoma" = "#FD8CC1FF","glioma"="#ADE2D0FF"),
  time_info = c(Yes="skyblue4", No="lightblue"),
  sex_info = c(Yes="gray50", No="gray80"),
  Data_type = c(Transcriptome="#8844AA", Proteome="#CC88CC", Metabolomics="#FFCCFF")
)

ha=HeatmapAnnotation(df=dataset,gp = gpar(col = "white",lwd=1,
                                          lineend="round",linejoin="round",
                                          linemitre=20),
                     #border=TRUE,point=point,
                     gap = unit(1, "points"),
                     annotation_name_gp = gpar(width=1),col = col_list,height=unit(4, "mm"))
zero_row_mat=matrix(nrow=0, ncol=nrow(dataset))
Hm=Heatmap(zero_row_mat, top_annotation=ha,col = col_list,km = 2,
           rect_gp = gpar(col = "white"))
draw(Hm,annotation_legend_side = "bottom")

pdf("fig1B.pdf",width=20,height = 8)
draw(Hm,annotation_legend_side = "bottom")
dev.off()


tissue<-c("BAT","iWAT","eWAT","pWAT","qM","gM","raM","taM","hM","dM","sM","cM","xM",
          "Liver","Heart","Blood","VenaCava","PortalVein","Gastric","Cerebellum",
          "Colon","Caecal","Hypothalamus","Hippocampus","Neocortex")
tissue_fullname<-c("Brown Adipose Tissue","inguinal or subcutaneous White Adipose Tissue",
                   "epididymal or visceral White Adipose Tissue","Peritumoral white adipose tissue",
                   "quadriceps Muscle","gastrocnemius c","rectus abdominus Muscle",
                   "tibialis anterior Muscle","hindlimb muscle","diaphragm Muscle","soleus Muscle",
                   "cardiac muscle","Skeletal Muscle","","","","","","","","","","","","")
tissue_group<-c("Adipose Tissue","Adipose Tissue","Adipose Tissue","Adipose Tissue",
                "Muscle Tissue","Muscle Tissue","Muscle Tissue","Muscle Tissue",
                "Muscle Tissue","Muscle Tissue","Muscle Tissue","Muscle Tissue",
                "Muscle Tissue","","","","","","","","","","","","")

datatissue<-data.frame("Tissue"=tissue,"Tissue_FullName"=tissue_fullname,"Tissue_group"=tissue_group,
                       HumanDatasetNum=0,HumanSampleNum=0,
                       MouseDatasetNum=0,MouseSampleNum=0)

for (i in c(1:dim(datatissue)[1])) {
  datatissue$HumanDatasetNum[i]<-length(which(dataset$TaxonomyID=="Human" & dataset$Tissue==datatissue$Tissue[i]))
  datatissue$HumanSampleNum[i]<-sum(dataset$SampleNum[which(dataset$TaxonomyID=="Human" & dataset$Tissue==datatissue$Tissue[i])])
  datatissue$MouseDatasetNum[i]<-length(which(dataset$TaxonomyID=="Mouse" & dataset$Tissue==datatissue$Tissue[i]))
  datatissue$MouseSampleNum[i]<-sum(dataset$SampleNum[which(dataset$TaxonomyID=="Mouse" & dataset$Tissue==datatissue$Tissue[i])])
}
rownames(datatissue)<-datatissue$Tissue
row_info <-as_tibble(datatissue %>%select(Tissue_group, Tissue))
colnames(row_info)<-c("group","id")
row_group<-tribble(
  ~group, ~level1,
  "Adipose Tissue","Adipose Tissue",
  "Muscle Tissue","Muscle Tissue",
  "","Other"
)
palettes <- tribble(
  ~palette,             ~colours,
  "p1",            grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(101),
  "p2",          grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues") %>% c("#0C2340"))(10),
  "Human", "#003C67FF",
  "Mouse", "#0073C2FF"
)
colinfo<-tribble(
    ~group,  ~id,                ~name,       ~geom,   ~palette,  ~options,
    "",       "Tissue",          "",           "text",  NA,       list(hjust = 0, width = 4),
    "",       "Tissue_FullName", "",           "text",  NA,       list(hjust = 0, width =15.5,overlay=FALSE),
    "Human",  "HumanDatasetNum", "DatasetNum", "text",   NA,       list(hjust = 0, width = 2),
    "Human",  "HumanDatasetNum", "", "bar",   "p1",     list(width = 4,overlay=FALSE),
    "Human",  "HumanSampleNum", "SampleNum", "text",   NA,       list(hjust = 0, width = 2),
    "Human",  "HumanSampleNum",  "",  "bar",   "p2",list(width = 4),
    "Mouse",  "MouseDatasetNum", "DatasetNum", "text",   NA,       list(hjust = 0, width = 2),
    "Mouse",  "MouseDatasetNum", "", "bar",   "p1",list(width = 4, legend = FALSE),
    "Mouse",  "MouseSampleNum", "SampleNum", "text",   NA,       list(hjust = 0, width = 2),
    "Mouse",  "MouseSampleNum",  "",  "bar",   "p2",list(width = 4, legend = FALSE)
)
colgroup<-tribble(
  ~group, ~level1,~palette,
  "Human","Human","Human",
  "Mouse","Mouse","Mouse"
)
p<-funky_heatmap(
  datatissue,
  row_info = row_info,
  row_groups = row_group,
  column_info = colinfo,
  column_groups = colgroup,
  palettes = palettes,
  #scale_column=FALSE,
  position_args = position_arguments(expand_xmax = 10,col_annot_angle=0,
                                     col_space = 0.3)
)
pdf("fig1C.pdf",width = 15,height=8)
p
dev.off()


humandata<-dataset[which(dataset$TaxonomyID=="Human"),]
mousedata<-dataset[which(dataset$TaxonomyID=="Mouse"),]
Cancer_type_col<-c("pancreatic cancer"="#FED439FF", #"colon cancer"="#709AE1FF",           
                   "lung cancer"="#D2AF81FF", #"gastrointestinal cancer"="#D5E4A2FF",
                   "colorectal cancer"= "#709AE1FF", "Multiple cancers"="#8A9197FF",       
                   "gastric cancer"= "#D5E4A2FF", "Breast Cancer"= "#1A9993FF",         
                   "melanoma" = "#FD8CC1FF","glioma"="#ADE2D0FF")

cancersort<-humandata %>% group_by(Cancer_type) %>% summarise("num"=sum(SampleNum),.groups = "drop") %>%
  arrange(num) %>%
  pull(Cancer_type)
humandata$Cancer_type<-factor(humandata$Cancer_type,levels=cancersort)

cbh<-ggplot(humandata)+
  geom_bar(aes(x=Cancer_type,y=SampleNum,fill=Cancer_type), stat='identity')+
  coord_polar(theta = "y")+
  theme_bw()+
  scale_fill_manual(values=Cancer_type_col)+
  theme_minimal()+
  xlab("")+ 
  ylab("")+
  theme(legend.position = "none",
        axis.text.x = element_text(color="black", face="italic"))+
  ylim(c(0,1000))
pdf("fig1D.pdf",width=5,height = 5)
cbh
dev.off()