###########Fig2########################
##Author: tzh
library(VennDiagram)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
#org.Hs.eg.db,human org.Mm.eg.db,mouse
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(rrvgo)
library(ComplexHeatmap)
library(reshape2)
library(ggsankey)
library(ggforce)
library(stringr)


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
  Cancer_type = c("pancreatic cancer"="#FED439FF",           
                  "lung cancer"="#D2AF81FF",
                  "colorectal cancer"= "#709AE1FF", "Multiple cancers"="#8A9197FF",       
                  "gastric cancer"= "#D5E4A2FF", "Breast Cancer"= "",         
                  "melanoma" = "#FD8CC1FF","glioma"=""),
  time_info = c(Yes="skyblue4", No="lightblue"),
  sex_info = c(Yes="gray50", No="gray80"),
  Data_type = c(Transcriptome="#8844AA", Proteome="#CC88CC", Metabolomics="#FFCCFF")
)
Tissue_col<-c("Muscle"=  "#E64B35FF","raM"=  "#E64B35FF","iWAT"="#00A087FF","Blood"="#DC0000FF",
  "Colon"="#3C5488FF","xM"="#F39B7FFF","pWAT"="#8491B4FF",
  "Liver"="#91D1C2FF","gM"="#4DBBD5FF","taM"="#7E6148FF",
  "eWAT"="#B09C85FF", "qM"="#374E55FF","Hypothalamus"="#DF8F44FF",
  "hM" ="#00A1D5FF",  "dM"="#B24745FF", "BAT"= "#79AF97FF",
  "Cerebellum"="#6A6599FF","Heart"="#800000FF","Hippocampus"="#80796BFF",
  "Neocortex"="#FFA319FF", "sM"="#8A9045FF","cM"="#155F83FF",
  "Gastric"="#C16622FF","Caecal"="#58593FFF","VenaCava"="#350E20FF",    
  "PortalVein"="")

  #################################Transcriptome###############################
dataset<-read.xlsx("../ALLsample_info1.xlsx",sheet=2)
dataset<-dataset[,-1]
#dataset<-dataset[-which(dataset$CCID_orig=="new"),]
geneDEG<-read.table("../MySQL_DEGs2.txt",header=T,sep='\t')

#some genes are duplicated in some CCID, duo to ID transformation
geneDEG_gene_count <- geneDEG %>%
  count(CCID, Gene, name = "count")
length(unique(geneDEG_gene_count$Gene[which(geneDEG_gene_count$count>1)]))
geneDEG_gene_count_dup<-geneDEG_gene_count[which(geneDEG_gene_count$count>1),]
write.xlsx(geneDEG_gene_count_dup,"fig2/geneDEG_gene_count.xlsx")

geneDEG_filter <- geneDEG %>%
  group_by(CCID, Gene) %>%
  filter(n() < 2) %>%
  ungroup()
geneDEG_gene_count_filter <- geneDEG_filter %>%
  count(CCID, Gene, name = "count")
length(unique(geneDEG_gene_count_filter$Gene[which(geneDEG_gene_count_filter$count>1)]))

geneDEG<-geneDEG_filter
geneDEG<-data.frame(geneDEG,"TaxonomyID"=0,"Tissue"=0,"Tissue_class"=0,"CellLine"=0,"Cancer_type"=0)
for (id in unique(geneDEG$CCID)) {
  geneDEG$TaxonomyID[which(geneDEG$CCID == id)]<-dataset$TaxonomyID[which(dataset$CCID_orig==id)]
  geneDEG$Tissue[which(geneDEG$CCID == id)]<-dataset$Tissue[which(dataset$CCID_orig==id)]
  geneDEG$Tissue_class[which(geneDEG$CCID == id)]<-dataset$Tissue_class[which(dataset$CCID_orig==id)]
  geneDEG$CellLine[which(geneDEG$CCID == id)]<-dataset$CellLine[which(dataset$CCID_orig==id)]
  geneDEG$Cancer_type[which(geneDEG$CCID == id)]<-dataset$Cancer_type[which(dataset$CCID_orig==id)]
}

geneDEG_diff<-geneDEG[which(abs(geneDEG$logFC)>=1 & geneDEG$PValue<0.05),]
table(geneDEG_diff$CCID)

######1
geneDEG_mouse<-geneDEG[geneDEG$TaxonomyID=="Mouse",]

genemouseall<-as.data.frame(geneDEG_mouse %>% group_by(Tissue) %>% 
  summarise(n_gene = n_distinct(Gene), .groups = "drop"))
genemouseall$group<-"ALLgene"
genemouseDiff<-as.data.frame(geneDEG_diff_mouse %>% group_by(Tissue) %>% 
  summarise(n_gene = n_distinct(Gene), .groups = "drop"))
genemouseDiff$group<-"Diff"
genemousenoDiff<-genemouseDiff
genemousenoDiff$n_gene<-genemouseall$n_gene-genemouseDiff$n_gene
genemousenoDiff$group<-"No Diff"
genemouse<-rbind(genemousenoDiff,genemouseDiff)

genemouseDiff$Tissue[order(genemouseDiff$n_gene,decreasing = T)]
genemouse$Tissue<-factor(genemouse$Tissue,c("taM","dM","gM","qM","sM","hM","iWAT","Liver","eWAT","BAT","Heart","Hippocampus","Cerebellum","Neocortex","Hypothalamus"))
genemouse$group<-factor(genemouse$group,levels = c("No Diff","Diff"))
pdf("fig2A.pdf",width = 6,height = 3)
ggplot(genemouse, aes(x = Tissue, y = n_gene,fill = group)) +
  geom_col(position = "stack",width=0.7) +
  scale_fill_manual(values=c("gray","#8844AA"))+
  labs(
    x = "",
    y = "Genes_num",
    title = "Distribution of DEGs in tissues"
  ) +theme_bw()+
  scale_y_continuous(breaks = seq(0,30000,5000),limits = c(0,30000))+
  theme(axis.text.x = element_text(hjust=1,vjust=1,angle=30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

######################----------consensus gene------------#########################
############filter DEG then calculate freq##########
geneDEG_diff_mouse<-geneDEG[which(abs(geneDEG$logFC)>=1 & geneDEG$PValue<0.05 & geneDEG$TaxonomyID=="Mouse"),]
tissuegene <- geneDEG_diff_mouse %>%
  mutate(direction = ifelse(logFC > 0, "up", "down")) %>% 
  group_by(Tissue, Gene) %>%
  summarise(
    n_dataset = n_distinct(CCID),                   
    n_up = sum(direction == "up"),                      
    n_down = sum(direction == "down"),   
    .groups = "drop"
  )
tissue_dataset_count <- geneDEG_diff_mouse %>%
  group_by(Tissue) %>%
  summarise(
    dataset_count = n_distinct(CCID),
    .groups = "drop"
  )
tissuegene_result <- tissuegene %>%
  left_join(tissue_dataset_count, by = "Tissue")
tissuegene_result$prop_up<-tissuegene_result$n_up/tissuegene_result$dataset_count
tissuegene_result$prop_down<-tissuegene_result$n_down/tissuegene_result$dataset_count
tissuegene_result<-tissuegene_result[order(tissuegene_result$Tissue,tissuegene_result$n_dataset,decreasing = T),]

tissuegene_consensus<-tissuegene_result[which(tissuegene_result$prop_up > 0.5| tissuegene_result$prop_down > 0.5),]
#delete the tissue with only one dataset, because the consensus gene is not meaningful
tissuegene_consensus<-tissuegene_consensus[-which(tissuegene_consensus$tissue_dataset_count==1),]

tissuegene_consensus_countup<-tissuegene_consensus %>% group_by(Tissue) %>% 
  summarise(count=sum(prop_up > 0.5, na.rm = TRUE),
            group="up") %>% ungroup()
tissuegene_consensus_countdown<-tissuegene_consensus %>% group_by(Tissue) %>% 
  summarise(count=sum(prop_down > 0.5, na.rm = TRUE),
            group="down") %>% ungroup()

tissuegene_consensus_count<-rbind(tissuegene_consensus_countup,tissuegene_consensus_countdown)
tissuegene_consensus_count$Tissue[order(tissuegene_consensus_count$count,decreasing = T)]
tissuegene_consensus_count$Tissue<-factor(tissuegene_consensus_count$Tissue,
                                          levels = c("taM","dM","gM","qM","sM","hM","iWAT","Liver",
                                                     "eWAT","BAT","Heart","Hippocampus","Cerebellum","Neocortex"))
tissuegene_consensus_count$prop<-0
for (i in c(1:dim(tissuegene_consensus_count)[1])) {
  tissue<-tissuegene_consensus_count$Tissue[i]
  all=genemouseall$n_gene[which(genemouseall$Tissue == tissue)]
  tissuegene_consensus_count$prop[i]<-tissuegene_consensus_count$count[i]/all
}

pdf("fig2B.pdf",width = 6.5,height = 3)
ggplot(tissuegene_consensus_count, aes(x=Tissue, y=ifelse(group=="up", count, -count), fill=group) )+
  geom_bar(stat="identity",width=0.7)+
  scale_fill_manual(values=c("skyblue4","#B24745FF"))+
  labs(
    x = "",
    y = "Consensus_Gene_Num"
    #title = "Distribution of DEGs in tissues"
  ) +theme_bw()+
  scale_y_continuous(labels=abs,limits = c(-500,800))+
  theme(axis.text.x = element_text(hjust=1,vjust=1,angle=30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

tissuegene_consensus$group<-"up"
tissuegene_consensus$group[which(tissuegene_consensus$prop_down > 0.5)]<-"down"
tissuegene_consensus$Tissue_class<-tissuegene_consensus$Tissue
tissuegene_consensus$Tissue_class[which(tissuegene_consensus$Tissue_class %in% c("qM","gM","raM","taM","hM","dM","sM"))]<-"Muscle"

g1<-tissuegene_consensus$Gene[which(tissuegene_consensus$Tissue_class=="Muscle" & tissuegene_consensus$group=="up")]
g2<-tissuegene_consensus$Gene[which(tissuegene_consensus$Tissue_class=="Muscle" & tissuegene_consensus$group=="down")]
g<-intersect(g1,g2) #157genes
tissuegene_consensus_bindmuscle<-tissuegene_consensus[-which(tissuegene_consensus$Tissue_class=="Muscle" & 
                                                              tissuegene_consensus$Gene %in% g),]

for (i in unique(tissuegene_consensus_bindmuscle$Tissue_class)) { #tissuegene_consensus_bindmuscle
  for (m in c("up","down")) {
    print(i)
    print(m)
    if(m=="up"){
      gene<-tissuegene_consensus_bindmuscle[which(tissuegene_consensus_bindmuscle$group=="up"),]
    }else{
      gene<-tissuegene_consensus_bindmuscle[which(tissuegene_consensus_bindmuscle$group=="down"),]
    }
    genes<-unique(gene$Gene[which(gene$Tissue_class==i)])
    ego <- enrichGO(gene          = genes,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = "SYMBOL",
                    ont           = "ALL", 
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)
    result_df <- as.data.frame(ego)
    if(length(grep("BP",result_df$ONTOLOGY))>=2){
      simMatrix <- calculateSimMatrix(result_df$ID, orgdb="org.Mm.eg.db", ont="BP", method="Rel")
      scores <- setNames(-log10(result_df$pvalue), result_df$ID)
      reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.7, orgdb="org.Mm.eg.db")
      t<-data.frame(table(reducedTerms$parentTerm))
      t<-t[order(t$Freq,decreasing = T),]
      write.csv(simMatrix,paste0("fig2/consensus_gene/bindmuscle/",i,"_",m,"_simMatrix.csv"))
      write.csv(reducedTerms,paste0("fig2/consensus_gene/bindmuscle/",i,"_",m,"_reducedTerms.csv"))
      write.csv(t,paste0("fig2/consensus_gene/bindmuscle/",i,"_",m,"_parentTerm.csv"))
    }
  }
}

rrvgo_top2<-read.xlsx("fig2/rrvgo_parentterm.xlsx",sheet=2)
rrvgo_top2$tissue<-factor(rrvgo_top2$tissue,
                                          levels = c("BAT","iWAT","eWAT","Muscle",
                                                     "Liver","Heart","Cerebellum","Hippocampus"))
rrvgo_top2$parentTerm<-factor(rrvgo_top2$parentTerm,unique(rrvgo_top2$parentTerm[order(rrvgo_top2$tissue)]))

pdf("fig2C.pdf",width = 10,height = 6)
ggplot(data=rrvgo_top2, aes(x=ifelse(group=="up", freq, -freq), y=parentTerm, fill=tissue))+
  geom_bar(stat="identity",width=0.8)+
  geom_text(aes(x=ifelse(group=="up", -40, 2),
                y=parentTerm, label=parentTerm), 
            size=2.5,
            vjust=0,hjust=0) +
  geom_vline(aes(xintercept=0),linetype = "dashed")+
  scale_fill_manual(values=Tissue_col)+
  scale_x_continuous(labels=abs)+
  theme_classic()+ 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour ='black', size =12),
        axis.ticks.x = element_line(colour ='black', linewidth = 1),
        axis.title.x = element_blank()
  )
dev.off()


################################################Protein#####################################
protein<-read.xlsx("protein_20251103.xlsx",sheet=2)
dataset_protein<-dataset[which(dataset$Data_type=="Proteome"),]
protein$Gene<-toupper(protein$Gene)

protein<-data.frame(protein,"TaxonomyID"=0,"Tissue"=0,"Tissue_class"=0,"CellLine"=0,"Cancer_type"=0)
for (id in unique(protein$CCID)) {
  protein$TaxonomyID[which(protein$CCID == id)]<-dataset_protein$TaxonomyID[which(dataset_protein$CCID_ineachomics==id)]
  protein$Tissue[which(protein$CCID == id)]<-dataset_protein$Tissue[which(dataset_protein$CCID_ineachomics==id)]
  protein$Tissue_class[which(protein$CCID == id)]<-dataset_protein$Tissue_class[which(dataset_protein$CCID_ineachomics==id)]
  protein$CellLine[which(protein$CCID == id)]<-dataset_protein$CellLine[which(dataset_protein$CCID_ineachomics==id)]
  protein$Cancer_type[which(protein$CCID == id)]<-dataset_protein$Cancer_type[which(dataset_protein$CCID_ineachomics==id)]
}
write.xlsx(protein,"fig2/protein.xlsx")

protein$PValue<-as.numeric(protein$PValue)
protein$PValue[is.na(protein$PValue)]<-0
proteinDiff<-protein[which(abs(protein$log2FC)>=1 & protein$PValue<0.05),]
table(proteinDiff$Tissue)
write.xlsx(proteinDiff,"fig2/proteinDiff.xlsx")

proteinall<-as.data.frame(protein %>% group_by(Tissue) %>% 
                              summarise(n_gene = n_distinct(Gene), .groups = "drop"))
proteinall$group<-"ALLgene"
proDiff<-as.data.frame(proteinDiff %>% group_by(Tissue) %>% 
                               summarise(n_gene = n_distinct(Gene), .groups = "drop"))
proDiff$group<-"Diff"
proDiff<-rbind(proDiff[1:4,],c("pWAT",0,"Diff"),proDiff[5:7,])
proDiff$n_gene<-as.numeric(proDiff$n_gene)
pronoDiff<-proDiff
pronoDiff$n_gene<-proteinall$n_gene-proDiff$n_gene
pronoDiff$group<-"No Diff"
pro<-rbind(pronoDiff,proDiff)

proDiff$Tissue[order(proDiff$n_gene,decreasing = T)]
pro$Tissue<-factor(pro$Tissue,levels = c("qM","Blood","xM","gM","sM","cM","Colon","pWAT"))
pro$group<-factor(pro$group,levels = c("No Diff","Diff"))
pdf("fig2D.pdf",width = 4,height = 2)
ggplot(pro, aes(x = Tissue, y = n_gene,fill = group)) +
  geom_col(position = "stack",width=0.7) +
  scale_fill_manual(values=c("gray","#CC88CC"))+
  labs(
    x = "",
    y = "proteins_num",
    title = "Distribution of Diff proteins in tissues"
  ) +theme_bw()+
  scale_y_continuous(breaks = seq(0,600,100),limits = c(0,600))+
  theme(axis.text.x = element_text(hjust=1,vjust=1,angle=30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

proteinDiff_human<-proteinDiff[which(proteinDiff$TaxonomyID=="Human"),]
proteinDiff_mouse<-proteinDiff[which(proteinDiff$TaxonomyID=="Mouse"),]
proteinDiff_mouse$Gene<-tools::toTitleCase(tolower(proteinDiff_mouse$Gene))

proteinDiff_human_stat <- proteinDiff_human %>%
  mutate(regulation = case_when(
    log2FC > 0 ~ "Up",
    log2FC < 0 ~ "Down",
    TRUE   ~ NA_character_
  )) %>%
  filter(!is.na(regulation)) %>%
  distinct(Tissue_class, Gene, regulation) %>%  
  count(Tissue_class, regulation) %>%
  mutate(labels=paste(Tissue_class, regulation,n,sep = "_"))

proteinDiff_mouse_stat <- proteinDiff_mouse %>%
  mutate(regulation = case_when(
    log2FC > 0 ~ "Up",
    log2FC < 0 ~ "Down",
    TRUE   ~ NA_character_
  )) %>%
  filter(!is.na(regulation)) %>%
  distinct(Tissue_class, Gene, regulation) %>%  
  count(Tissue_class, regulation )%>%
  mutate(labels=paste(Tissue_class, regulation,n,sep = "_"))


p<-ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+#添加颜色
  scale_fill_manual(values = c( '#374E55FF', "#DC0000FF", "#3C5488FF",
                               '#BD956A',"#E64B35FF"))+
  geom_arc_bar(data=proteinDiff_human_stat,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=n,fill=labels))
pdf("fig2E_human.pdf",width = 4.5,height = 3)
p
dev.off()

p1<-ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+#添加颜色
  scale_fill_manual(values = c( '#155F83FF',
                                '#BD956A',"#E64B35FF"))+
  geom_arc_bar(data=proteinDiff_mouse_stat,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=n,fill=labels))
pdf("fig2E_mouse.pdf",width = 5,height = 3)
p1
dev.off()

for (i in unique(proteinDiff_human$Tissue_class)) {
  for (m in c("up","down")) {
    print(i)
    print(m)
    if(m=="up"){
      gene<-proteinDiff_human[which(proteinDiff_human$log2FC >= 1),]
    }else{
      gene<-proteinDiff_human[which(proteinDiff_human$log2FC <= -1),]
    }
    genename<-unique(gene$Gene[which(gene$Tissue_class==i)])
    if(length(genename)!=0){
      gene.df <- bitr(genename,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
      geneid <- gene.df$ENTREZID   
      
      ekg <- enrichKEGG(gene = geneid,keyType = "kegg",organism= "hsa", qvalueCutoff = 0.05, pvalueCutoff=0.05)

      result_kgg <- as.data.frame(ekg)
      result_kgg <- result_kgg %>% arrange(desc(Count))
      result_kgg$Description<-factor(result_kgg$Description,levels=rev(result_kgg$Description))
      if(dim(result_kgg)[1]!=0){
        result_kgg$genename<-"0"
        for (j in c(1:dim(result_kgg)[1])) {
          name<-gene.df$SYMBOL[match(unlist(strsplit(result_kgg$geneID[j],'/')),gene.df$ENTREZID)]
          result_kgg$genename[j]<-paste(name,collapse='/')
        }
      }
      write.csv(result_kgg,paste0("fig2/protein_enrichment/HUMAN",i,"_",m,"_kegg.csv"))
    }
  }
}

proteinDiff_mouse<-proteinDiff_mouse[-grep(";",proteinDiff_mouse$Gene),]
for (i in unique(proteinDiff_mouse$Tissue_class)) {
  for (m in c("up","down")) {
    print(i)
    print(m)
    if(m=="up"){
      gene<-proteinDiff_mouse[which(proteinDiff_mouse$log2FC >= 1),]
    }else{
      gene<-proteinDiff_mouse[which(proteinDiff_mouse$log2FC <= -1),]
    }
    genename<-unique(gene$Gene[which(gene$Tissue_class==i)])
    if(length(genename)!=0){
      gene.df <- bitr(genename,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Mm.eg.db)
      geneid <- gene.df$ENTREZID   
      
      ekg <- enrichKEGG(gene = geneid,keyType = "kegg",organism= "mmu", qvalueCutoff = 0.05, pvalueCutoff=0.05)
      
      result_kgg <- as.data.frame(ekg)
      result_kgg <- result_kgg %>% arrange(desc(Count))
      result_kgg$Description<-factor(result_kgg$Description,levels=rev(result_kgg$Description))
      if(dim(result_kgg)[1]!=0){
        result_kgg$genename<-"0"
        for (j in c(1:dim(result_kgg)[1])) {
          name<-gene.df$SYMBOL[match(unlist(strsplit(result_kgg$geneID[j],'/')),gene.df$ENTREZID)]
          result_kgg$genename[j]<-paste(name,collapse='/')
        }
      }
      write.csv(result_kgg,paste0("fig2/protein_enrichment/MOUSE",i,"_",m,"_kegg.csv"))
    }
  }
}


#################################Metabolism########################################
metabolism<-read.xlsx("metabolism_diff_data_251105.xlsx",sheet=4)
dataset_metabolism<-dataset[which(dataset$Data_type=="Metabolomics"),]
metabolism<-metabolism[which(metabolism$CCID %in% dataset_metabolism$CCID_ineachomics),]

metabolism<-data.frame(metabolism,"TaxonomyID"=0,"Tissue"=0,"Tissue_class"=0,"CellLine"=0,"Cancer_type"=0,"group"=0)
for (id in unique(metabolism$CCID)) {
  metabolism$TaxonomyID[which(metabolism$CCID == id)]<-dataset_metabolism$TaxonomyID[which(dataset_metabolism$CCID_ineachomics==id)]
  metabolism$Tissue[which(metabolism$CCID == id)]<-dataset_metabolism$Tissue[which(dataset_metabolism$CCID_ineachomics==id)]
  metabolism$Tissue_class[which(metabolism$CCID == id)]<-dataset_metabolism$Tissue_class[which(dataset_metabolism$CCID_ineachomics==id)]
  metabolism$CellLine[which(metabolism$CCID == id)]<-dataset_metabolism$CellLine[which(dataset_metabolism$CCID_ineachomics==id)]
  metabolism$Cancer_type[which(metabolism$CCID == id)]<-dataset_metabolism$Cancer_type[which(dataset_metabolism$CCID_ineachomics==id)]
  metabolism$group[which(metabolism$CCID == id)]<-dataset_metabolism$group[which(dataset_metabolism$CCID_ineachomics==id)]
}
write.xlsx(metabolism,"fig2/metabolism.xlsx")
metabolismDiff<-metabolism[which(abs(metabolism$log2FC_define)>0 & metabolism$pValue<0.05),]
write.xlsx(metabolismDiff,"fig2/metabolismDiff.xlsx")

metaall<-as.data.frame(metabolism %>% group_by(Tissue) %>% 
                            summarise(n_gene = n_distinct(Full.name), .groups = "drop"))
metaall$group<-"ALLgene"
metaDiff<-as.data.frame(metabolismDiff %>% group_by(Tissue) %>% 
                         summarise(n_gene = n_distinct(Full.name), .groups = "drop"))
metaDiff$group<-"Diff"

mergemeta<-merge(metaall,metaDiff,by.x="Tissue",by.y="Tissue",all=T)
mergemeta$n_gene.y[is.na(mergemeta$n_gene.y)]<-0
mergemeta$n_gene_nodiff<-mergemeta$n_gene.x-mergemeta$n_gene.y

metaDiff<-data.frame("Tissue"=mergemeta$Tissue,"n_gene"=mergemeta$n_gene.y)
metaDiff$group<-"Diff"
metanoDiff<-data.frame("Tissue"=mergemeta$Tissue,"n_gene"=mergemeta$n_gene_nodiff)
metanoDiff$group<-"No Diff"
meta<-rbind(metanoDiff,metaDiff)

metaDiff$Tissue[order(metaDiff$n_gene,decreasing = T)]
meta$Tissue<-factor(meta$Tissue,levels = c("xM","Blood","Liver","gM","taM","Gastric","Caecal",    
                                           "PortalVein","VenaCava"))
meta$group<-factor(meta$group,levels = c("No Diff","Diff"))

pdf("fig2H.pdf",width = 4,height = 2)
ggplot(meta, aes(x = Tissue, y = n_gene,fill = group)) +
  geom_col(position = "stack",width=0.7) +
  scale_fill_manual(values=c("gray","#FFCCFF"))+
  labs(
    x = "",
    y = "metabolism_num",
    title = "Distribution of Diff_metabolism in tissues"
  ) +theme_bw()+
  scale_y_continuous(breaks = seq(0,6500,1000),limits = c(0,6500))+
  theme(axis.text.x = element_text(hjust=1,vjust=1,angle=30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


metabolismDiff_human<-metabolismDiff[which(metabolismDiff$TaxonomyID=="Human"),]
metabolismDiff_mouse<-metabolismDiff[which(metabolismDiff$TaxonomyID=="Mouse"),]

metabolismDiff_human_stat <- metabolismDiff_human %>% 
  mutate(regulation = case_when(
    log2FC_define > 0 ~ "Up",
    log2FC_define < 0 ~ "Down",
    TRUE   ~ NA_character_
  )) %>%
  filter(!is.na(regulation)) %>%
  distinct(Tissue_class, Full.name, regulation) %>%  
  count(Tissue_class, regulation) %>%
  mutate(labels=paste(Tissue_class, regulation,n,sep = "_"))

metabolismDiff_mouse_stat <- metabolismDiff_mouse %>%
  mutate(regulation = case_when(
    log2FC_define > 0 ~ "Up",
    log2FC_define < 0 ~ "Down",
    TRUE   ~ NA_character_
  )) %>%
  filter(!is.na(regulation)) %>%
  distinct(Tissue_class, Full.name, regulation) %>%  
  count(Tissue_class, regulation) %>%
  mutate(labels=paste(Tissue_class, regulation,n,sep = "_"))

p<-ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())
  xlab("")+ylab('')+
  scale_fill_manual(values = c( '#374E55FF', "#DC0000FF",
                                '#BD956A',"#E64B35FF"))+
  geom_arc_bar(data=metabolismDiff_human_stat,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=n,fill=labels))
pdf("fig2I_human.pdf",width = 4.5,height = 3)
p
dev.off()

p1<-ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+
  scale_fill_manual(values = c( '#374E55FF', "#DC0000FF",'#155F83FF',"#58593FFF","skyblue4",
                                "#80796BFF","#91D1C2FF",'#BD956A',"#E64B35FF",
                                "#3D3B25FF","#8A9197FF","#350E20FF","#8A9197FF"))+
  geom_arc_bar(data=metabolismDiff_mouse_stat,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=n,fill=labels))
pdf("fig2I_mouse.pdf",width = 5,height = 4)
p1
dev.off()

####metabolism enrichment analysis, results from metaboanalyst####