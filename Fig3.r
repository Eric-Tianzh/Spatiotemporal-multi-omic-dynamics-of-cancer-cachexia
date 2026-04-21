###############Fig3########################
##Author: tzh
library(openxlsx)
library(sva)
library(limma)
library(edgeR)
library(pheatmap)
library(ggsankey)
library(ggplot2)
library(VennDiagram) 
library(clusterProfiler)
#org.Hs.eg.db,human org.Mm.eg.db,mouse
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(dplyr)
library(rrvgo)
library(UpSetR)
library(reshape2)
library(ComplexHeatmap)
library(ggrepel)


####load data####
#### 5 datasets have sexinfo########
sampleinfo<-read.xlsx("ALLsample_info1.xlsx")
sexinfo<-unique(sampleinfo$Dir[sampleinfo$sex_related=="1"])
sexinfo<-sexinfo[!is.na(sexinfo)]
sexsample<-unique(sampleinfo$info_file[sampleinfo$sex_related=="1"])
sexsample<-sexsample[grep("GSE",sexsample)]

datasummary<-read.xlsx("ALLsample_info1.xlsx",sheet="sexinfo_summary")
datasummary<-datasummary[which(datasummary$taxo=="mouse"),]
datasummary<-datasummary %>% pivot_stages_longer(c("taxo","cancer","tissue","sex"),"num")

pos <- position_sankey(order ="ascending",v_space =0.05,direction = "forward")
pos_text <- position_sankey(order ="ascending",v_space =0.05,nudge_x =0.1,direction = "forward")
pdf("fig3A.pdf",width = 5,height = 3)
ggplot(datasummary,
       aes(x = stage, y = num, group = node,
           connector = connector,
           edge_id = edge_id)) +
  geom_sankeynode(aes(fill = node),position = pos,show.legend =F) +
  geom_sankeyedge(aes(fill = node), position = pos,show.legend =F) +
  geom_text(aes(label = node), stat ="sankeynode",color="black",
            position = pos_text, hjust =0.1, cex =3.5) +
  scale_x_discrete(expand = expansion(add = c(0.2,0.5)),
                   position ="top")+
  scale_fill_manual(values =c("#0073C2FF",
                              "#FED439FF","#709AE1FF","#709AE1FF","#D2AF81FF",
                              "#374E55FF","#7E6148FF","#4DBBD5FF",
                              "#FFCCFF","skyblue4"))+
  theme_void()+
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit ="cm"),
        axis.text.x=element_text(
          color="black",size=10))
dev.off()


c2<-read.table(paste0("../count/",sexinfo[2],"_count.txt"),header = T,sep='\t')
c3<-read.table(paste0("../count/",sexinfo[3],"_count.txt"),header = T,sep='\t')
c4<-read.table(paste0("../count/",sexinfo[4],"_count.txt"),header = T,sep='\t')
c5<-read.table(paste0("../count/",sexinfo[5],"_count.txt"),header = T,sep='\t')
genelist_m<-c2[,c("GeneName","GeneID")]
rownames(c2)<-c2$GeneID
c2<-c2[,-c(1,2)]
rownames(c3)<-c3$GeneID
c3<-c3[,-c(1,2)]
rownames(c4)<-c4$GeneID
c4<-c4[,-c(1,2)]
rownames(c5)<-c5$GeneID
c5<-c5[,-c(1,2)]

######MOUSE transcriptome##########
c2<-c2[,-c(40)]
sample_m2<-read.xlsx("fig3/sexinfo/GSE157251_series_matrix.xlsx",sheet=2)

sample_m2<-sample_m2[-grep("ACVR2B/Fc",sample_m2$agent),]
#qmKPC female
sh2<-sample_m2$samplename[which(sample_m2$Sex=="Female")]
sh2_agent<-sample_m2$agent[which(sample_m2$Sex=="Female")]
sh2_agent<-gsub("KPC_Vehicle","KPC",sh2_agent)
c2_sub1<-c2[,which(colnames(c2) %in% sh2)]

x<-c2_sub1
group<-factor(sh2_agent)
levels(group) #"Control" "KPC" FC=levels(group)[2]/levels(group)[1]
y<-DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y=normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
genemouse_1<-qlf$table
genemouse_1$genename<-genelist_m$GeneName[match(rownames(genemouse_1),genelist_m$GeneID)]
diffgenemouse_1<-genemouse_1[which((genemouse_1$logFC > 1 | genemouse_1$logFC < -1) & genemouse_1$PValue < 0.05),]
write.xlsx(genemouse_1,"fig3/gene_qm_KPC_female_cachexia_vs_control.xlsx")
write.xlsx(diffgenemouse_1,"fig3/diffgene_qm_KPC_female_cachexia_vs_control.xlsx")


#qmKPC male
sh3<-sample_m2$samplename[which(sample_m2$Sex=="Male")]
sh3_agent<-sample_m2$agent[which(sample_m2$Sex=="Male")]
sh3_agent<-gsub("KPC_Vehicle","KPC",sh3_agent)
c2_sub2<-c2[,which(colnames(c2) %in% sh3)]

x<-c2_sub2
group<-factor(sh3_agent)
levels(group)
y<-DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y=normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
genemouse_2<-qlf$table
genemouse_2$genename<-genelist_m$GeneName[match(rownames(genemouse_2),genelist_m$GeneID)]
diffgenemouse_2<-genemouse_2[which((genemouse_2$logFC > 1 | genemouse_2$logFC < -1) & genemouse_2$PValue < 0.05),]
write.xlsx(genemouse_2,"fig3/gene_qm_KPC_male_cachexia_vs_control.xlsx")
write.xlsx(diffgenemouse_2,"fig3/diffgene_qm_KPC_male_cachexia_vs_control.xlsx")


c3<-c3[,-c(89)]
sample_m3<-read.xlsx("fig3/sexinfo/GSE276018_series_matrix.xlsx",sheet=2)
sample_m3<-sample_m3[-grep("APC",sample_m3$`!Sample_title`),]
#tamcolon female
sh4<-sample_m3$samplename[which(sample_m3$Sex=="female")]
sh4_agent<-sample_m3$treatment[which(sample_m3$Sex=="female")]
c3_sub1<-c3[,which(colnames(c3) %in% sh4)]

x<-c3_sub1
group<-factor(sh4_agent)
levels(group)
y<-DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y=normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
genemouse_3<-qlf$table
genemouse_3$genename<-genelist_m$GeneName[match(rownames(genemouse_3),genelist_m$GeneID)]
diffgenemouse_3<-genemouse_3[which((genemouse_3$logFC > 1 | genemouse_3$logFC < -1) & genemouse_3$PValue < 0.05),]
write.xlsx(genemouse_3,"fig3/gene_tam_colon_female_cachexia_vs_control.xlsx")
write.xlsx(diffgenemouse_3,"fig3/diffgene_tam_colon_female_cachexia_vs_control.xlsx")

#tamcolon male
sh5<-sample_m3$samplename[which(sample_m3$Sex=="male")]
sh5_agent<-sample_m3$treatment[which(sample_m3$Sex=="male")]
c3_sub2<-c3[,which(colnames(c3) %in% sh5)]

x<-c3_sub2
group<-factor(sh5_agent)
levels(group)
y<-DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y=normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
genemouse_4<-qlf$table
genemouse_4$genename<-genelist_m$GeneName[match(rownames(genemouse_4),genelist_m$GeneID)]
diffgenemouse_4<-genemouse_4[which((genemouse_4$logFC > 1 | genemouse_4$logFC < -1) & genemouse_4$PValue < 0.05),]
write.xlsx(genemouse_4,"fig3/gene_tam_colon_male_cachexia_vs_control.xlsx")
write.xlsx(diffgenemouse_4,"fig3/diffgene_tam_colon_male_cachexia_vs_control.xlsx")


c4<-c4[,-c(41)]
sample_m4<-read.xlsx("fig3/sexinfo/GSE222317_series_matrix.xlsx",sheet=2)
sample_m4<-sample_m4[-grep("tumor mass",sample_m4$treatment),]
#gmllc female
sh6<-sample_m4$samplename[which(sample_m4$Sex=="Female")]
sh6_agent<-sample_m4$treatment[which(sample_m4$Sex=="Female")]
c4_sub1<-c4[,which(colnames(c4) %in% sh6)]

x<-c4_sub1
group<-factor(sh6_agent)
levels(group)
y<-DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y=normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
genemouse_5<-qlf$table
genemouse_5$genename<-genelist_m$GeneName[match(rownames(genemouse_5),genelist_m$GeneID)]
diffgenemouse_5<-genemouse_5[which((genemouse_5$logFC > 1 | genemouse_5$logFC < -1) & genemouse_5$PValue < 0.05),]
write.xlsx(genemouse_5,"fig3/gene_gm_llc_female_cachexia_vs_control.xlsx")
write.xlsx(diffgenemouse_5,"fig3/diffgene_gm_llc_female_cachexia_vs_control.xlsx")


c5<-c5[,-c(40)]
sample_m5<-read.xlsx("fig3/sexinfo/GSE114820_series_matrix.xlsx",sheet=2)
#gmllc male
sh7<-sample_m5$samplename[which(sample_m5$Sex=="Male")]
sh7_agent<-sample_m5$treatment[which(sample_m5$Sex=="Male")]
c5_sub1<-c5[,which(colnames(c5) %in% sh7)]

x<-c5_sub1
group<-factor(sh7_agent)
levels(group)
y<-DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y=normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
genemouse_6<-qlf$table
genemouse_6$genename<-genelist_m$GeneName[match(rownames(genemouse_6),genelist_m$GeneID)]
diffgenemouse_6<-genemouse_6[which((genemouse_6$logFC > 1 | genemouse_6$logFC < -1) & genemouse_6$PValue < 0.05),]
write.xlsx(genemouse_6,"fig3/gene_gm_llc_male_cachexia_vs_control.xlsx")
write.xlsx(diffgenemouse_6,"fig3/diffgene_gm_llc_male_cachexia_vs_control.xlsx")

diffgenemouse_female_up<-unique(diffgenemouse_female$genename[which(diffgenemouse_female$logFC>1)])
diffgenemouse_female_down<-unique(diffgenemouse_female$genename[which(diffgenemouse_female$logFC< -1)])
diffgenemouse_male_up<-unique(diffgenemouse_male$genename[which(diffgenemouse_male$logFC>1)])
diffgenemouse_male_down<-unique(diffgenemouse_male$genename[which(diffgenemouse_male$logFC< -1)])

genelistinmouse<-list("diffgenemouse_female_up"=diffgenemouse_female_up,
        "diffgenemouse_female_down"=diffgenemouse_female_down,
        "diffgenemouse_male_up"=diffgenemouse_male_up,
        "diffgenemouse_male_down"=diffgenemouse_male_down)

df<-UpSetR::fromList(genelistinmouse)
pdf("fig3C.pdf",width=6,height=4)
ComplexUpset::upset(
  df,
  intersect=colnames(df),
  min_size = 15,
  min_degree=1,
  n_intersections=40,
  sort_intersections_by="ratio",
)
dev.off()

##fig3D
listInput <- list(diffgenemouse_female_up=diffgenemouse_female_up,diffgenemouse_female_down=diffgenemouse_female_down,
                  diffgenemouse_male_up=diffgenemouse_male_up,diffgenemouse_male_down=diffgenemouse_male_down)
inter<-get.venn.partitions(listInput)
inter<-inter[order(inter$..count..,decreasing = T),]
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter<-inter[order(inter$..count..,decreasing = T),]
write.xlsx(inter,"fig3/mouse_tissue_sex_intersect_gene.xlsx")

diffmouse_mf_up<-as.character(inter$values[which(inter$..count..=="445")])
diffmouse_mf_up<-unlist(strsplit(diffmouse_mf_up,"\\|"))
gene.df <- bitr(diffmouse_mf_up,fromType="SYMBOL",toType=c("ENTREZID"), OrgDb = org.Mm.eg.db)
geneid <- gene.df$ENTREZID   
ekg <- enrichKEGG(gene = geneid,keyType = "kegg",organism= "mmu", qvalueCutoff = 1, pvalueCutoff=1)
result_kgg <- as.data.frame(ekg)
if(dim(result_kgg)[1]!=0){
  result_kgg$genename<-"0"
  for (j in c(1:dim(result_kgg)[1])) {
    name<-gene.df$SYMBOL[match(unlist(strsplit(result_kgg$geneID[j],'/')),gene.df$ENTREZID)]
    result_kgg$genename[j]<-paste(name,collapse='/')
  }
}
result_kgg_mfup<-result_kgg

diffmouse_mf_down<-as.character(inter$values[which(inter$..count..=="264")])
diffmouse_mf_down<-unlist(strsplit(diffmouse_mf_down,"\\|"))
gene.df <- bitr(diffmouse_mf_down,fromType="SYMBOL",toType=c("ENTREZID"), OrgDb = org.Mm.eg.db)
geneid <- gene.df$ENTREZID   
ekg <- enrichKEGG(gene = geneid,keyType = "kegg",organism= "mmu", qvalueCutoff = 1, pvalueCutoff=1)
result_kgg <- as.data.frame(ekg)
if(dim(result_kgg)[1]!=0){
  result_kgg$genename<-"0"
  for (j in c(1:dim(result_kgg)[1])) {
    name<-gene.df$SYMBOL[match(unlist(strsplit(result_kgg$geneID[j],'/')),gene.df$ENTREZID)]
    result_kgg$genename[j]<-paste(name,collapse='/')
  }
}
result_kgg_mfdown<-result_kgg


diffgenemouse<-rbind(diffgenemouse_1,diffgenemouse_2,diffgenemouse_3,diffgenemouse_4,diffgenemouse_5,diffgenemouse_6)
diffgenemouse_up<-diffgenemouse[which(diffgenemouse$logFC > 1),]
diffgenemouse_down<-diffgenemouse[which(diffgenemouse$logFC < -1),]

diffmouse_mf_up_data<-diffgenemouse_up[which(diffgenemouse_up$genename %in% diffmouse_mf_up),]
diffmouse_mf_up_data <- diffmouse_mf_up_data %>%
  group_by(genename) %>%
  mutate(mean_FC = mean(logFC))

diffmouse_mf_up_data_freq<-diffmouse_mf_up_data[,c("genename","Freq","mean_FC")]
diffmouse_mf_up_data_freq<-diffmouse_mf_up_data_freq[!duplicated(diffmouse_mf_up_data_freq),]
diffmouse_mf_up_data_freq<-diffmouse_mf_up_data_freq[order(diffmouse_mf_up_data_freq$Freq,diffmouse_mf_up_data_freq$mean_FC,decreasing = T),]
#upcommongene<-diffmouse_mf_up_data_freq$genename[which(diffmouse_mf_up_data_freq$Freq>3)]
diffmouse_mf_up_data_freq<-diffmouse_mf_up_data_freq[-grep("ENSMU",diffmouse_mf_up_data_freq$genename),]
diffmouse_mf_up_data_freq<-diffmouse_mf_up_data_freq[-grep("^Gm",diffmouse_mf_up_data_freq$genename),]
upcommongene<-diffmouse_mf_up_data_freq$genename[1:50]


diffmouse_mf_down_data<-diffgenemouse_down[which(diffgenemouse_down$genename %in% diffmouse_mf_down),]
diffmouse_mf_down_data <- diffmouse_mf_down_data %>%
  group_by(genename) %>%
  mutate(mean_FC = mean(logFC))

length(unique(diffmouse_mf_down_data$genename))
diffmouse_mf_down_data_freq<-diffmouse_mf_down_data[,c("genename","Freq","mean_FC")]
diffmouse_mf_down_data_freq$mean_FC<-abs(diffmouse_mf_down_data_freq$mean_FC)
diffmouse_mf_down_data_freq<-diffmouse_mf_down_data_freq[!duplicated(diffmouse_mf_down_data_freq),]
diffmouse_mf_down_data_freq<-diffmouse_mf_down_data_freq[order(diffmouse_mf_down_data_freq$Freq,diffmouse_mf_down_data_freq$mean_FC,decreasing = T),]
diffmouse_mf_down_data_freq<-diffmouse_mf_down_data_freq[-grep("ENSMU",diffmouse_mf_down_data_freq$genename),]
diffmouse_mf_down_data_freq<-diffmouse_mf_down_data_freq[-grep("^Gm",diffmouse_mf_down_data_freq$genename),]
diffmouse_mf_down_data_freq$rank<-seq(1,dim(diffmouse_mf_down_data_freq)[1],1)
downcommongene<-diffmouse_mf_down_data_freq$genename[1:50]

commongene<-c(upcommongene,downcommongene)
commongene_FC<-data.frame("Female_KPC"=genemouse_1$logFC[match(commongene,genemouse_1$genename)],
                             "Female_C26"=genemouse_3$logFC[match(commongene,genemouse_3$genename)],
                             "Female_LLC"=genemouse_5$logFC[match(commongene,genemouse_5$genename)],
                             "Male_KPC"=genemouse_2$logFC[match(commongene,genemouse_2$genename)],
                             "Male_C26"=genemouse_4$logFC[match(commongene,genemouse_4$genename)],
                             "Male_LLC"=genemouse_6$logFC[match(commongene,genemouse_6$genename)])
rownames(commongene_FC)<-commongene

commongene_pvalue<-data.frame("Female_KPC"=genemouse_1$PValue[match(commongene,genemouse_1$genename)],
                                 "Female_C26"=genemouse_3$PValue[match(commongene,genemouse_3$genename)],
                                 "Female_LLC"=genemouse_5$PValue[match(commongene,genemouse_5$genename)],
                                 "Male_KPC"=genemouse_2$PValue[match(commongene,genemouse_2$genename)],
                                 "Male_C26"=genemouse_4$PValue[match(commongene,genemouse_4$genename)],
                                 "Male_LLC"=genemouse_6$PValue[match(commongene,genemouse_6$genename)])
rownames(commongene_pvalue)<-commongene

sign_level <- function(p,f) {
  stars <- ""
  stars[p < 0.001 & abs(f) > 1] <- "***"
  stars[p >= 0.001 & p < 0.01 & abs(f) > 1] <- "**"
  stars[p >= 0.01 & p < 0.05 & abs(f) > 1] <- "*"
  return(stars)
}
star_matrix <- matrix(sign_level(commongene_pvalue,commongene_FC),nrow = nrow(commongene_pvalue))
star_matrix[is.na(star_matrix)] <- ""
pdf("fig3B.pdf",height=16,width = 5.5)
pheatmap(commongene_FC,
         display_numbers = star_matrix,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_number = 18,
         fontsize_row = 10, 
         fontsize_col = 10,
         color = colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(500),
         breaks=seq(-6,6, length.out = 501),
         legend_breaks = c(-4,-2,0,2,4,6),  
         legend_labels = c("-4","-2","0","2","4","6"),
         na_col = "white",
         border_color = NA,
         main = "Top common genes of mouse in Diff tissues and sex") 
dev.off()

###############---sex specific function enrichment analysis---###############
sampleinfo<-read.xlsx("ALLsample_info1.xlsx")
sexinfo<-unique(sampleinfo$Dir[sampleinfo$sex_related=="1"])
sexinfo<-sexinfo[!is.na(sexinfo)]
c2<-read.table(paste0("../count/",sexinfo[2],"_count.txt"),header = T,sep='\t')
c3<-read.table(paste0("../count/",sexinfo[3],"_count.txt"),header = T,sep='\t')
c4<-read.table(paste0("../count/",sexinfo[4],"_count.txt"),header = T,sep='\t')
c5<-read.table(paste0("../count/",sexinfo[5],"_count.txt"),header = T,sep='\t')
genelist_m<-c2[,c("GeneName","GeneID")]
rownames(c2)<-c2$GeneID
c2<-c2[,-c(1,2)]
rownames(c3)<-c3$GeneID
c3<-c3[,-c(1,2)]
rownames(c4)<-c4$GeneID
c4<-c4[,-c(1,2)]
rownames(c5)<-c5$GeneID
c5<-c5[,-c(1,2)]


######MOUSE transcriptome##########
##KPC
c2<-c2[,-c(40)]
sample_m2<-read.xlsx("fig3/sexinfo/GSE157251_series_matrix.xlsx",sheet=2)
sample_m2<-sample_m2[-grep("ACVR2B/Fc",sample_m2$agent),]
#qmKPC
sh2<-sample_m2$samplename[grep("KPC",sample_m2$agent)]
sh2_agent<-sample_m2$Sex[grep("KPC",sample_m2$agent)]
#sh2_agent<-gsub("KPC_Vehicle","KPC",sh2_agent)
c2_sub1<-c2[,which(colnames(c2) %in% sh2)]

x<-c2_sub1
group<-factor(sh2_agent)
levels(group) 
y<-DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y=normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
genemouse_1<-qlf$table
genemouse_1$genename<-genelist_m$GeneName[match(rownames(genemouse_1),genelist_m$GeneID)]
diffgenemouse_1<-genemouse_1[which((genemouse_1$logFC > 1 | genemouse_1$logFC < -1) & genemouse_1$PValue < 0.05),]
write.xlsx(genemouse_1,"fig3/gene_qm_KPC_female_cachexia_vs_control.xlsx")
write.xlsx(diffgenemouse_1,"fig3/diffgene_qm_KPC_female_cachexia_vs_control.xlsx")
sh3<-sample_m2$samplename[grep("Control",sample_m2$agent)]
sh3_agent<-sample_m2$Sex[grep("Control",sample_m2$agent)]
#sh3_agent<-gsub("KPC_Vehicle","KPC",sh3_agent)
c2_sub2<-c2[,which(colnames(c2) %in% sh3)]

x<-c2_sub2
group<-factor(sh3_agent)
levels(group)
y<-DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y=normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
genemouse_2<-qlf$table
genemouse_2$genename<-genelist_m$GeneName[match(rownames(genemouse_2),genelist_m$GeneID)]
diffgenemouse_2<-genemouse_2[which((genemouse_2$logFC > 1 | genemouse_2$logFC < -1) & genemouse_2$PValue < 0.05),]
write.xlsx(genemouse_2,"fig3/gene_qm_KPC_male_cachexia_vs_control.xlsx")
write.xlsx(diffgenemouse_2,"fig3/diffgene_qm_KPC_male_cachexia_vs_control.xlsx")

##c26
c3<-c3[,-c(89)]
sample_m3<-read.xlsx("fig3/sexinfo/GSE276018_series_matrix.xlsx",sheet=2)
sample_m3<-sample_m3[-grep("APC",sample_m3$`!Sample_title`),]
#tamcolon female
sh4<-sample_m3$samplename[which(sample_m3$treatment=="tumor")]
sh4_agent<-sample_m3$Sex[which(sample_m3$treatment=="tumor")]

c3_sub1<-c3[,which(colnames(c3) %in% sh4)]
x<-c3_sub1
group<-factor(sh4_agent)
levels(group)
y<-DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y=normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
genemouse_3<-qlf$table
genemouse_3$genename<-genelist_m$GeneName[match(rownames(genemouse_3),genelist_m$GeneID)]
diffgenemouse_3<-genemouse_3[which((genemouse_3$logFC > 1 | genemouse_3$logFC < -1) & genemouse_3$PValue < 0.05),]
write.xlsx(genemouse_3,"fig3/gene_tam_colon_female_cachexia_vs_control.xlsx")
write.xlsx(diffgenemouse_3,"fig3/diffgene_tam_colon_female_cachexia_vs_control.xlsx")

#tamcolon male
sh5<-sample_m3$samplename[which(sample_m3$treatment=="control")]
sh5_agent<-sample_m3$Sex[which(sample_m3$treatment=="control")]
c3_sub2<-c3[,which(colnames(c3) %in% sh5)]

x<-c3_sub2
group<-factor(sh5_agent)
levels(group)
y<-DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y=normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)
genemouse_4<-qlf$table
genemouse_4$genename<-genelist_m$GeneName[match(rownames(genemouse_4),genelist_m$GeneID)]
diffgenemouse_4<-genemouse_4[which((genemouse_4$logFC > 1 | genemouse_4$logFC < -1) & genemouse_4$PValue < 0.05),]
write.xlsx(genemouse_4,"fig3/gene_tam_colon_male_cachexia_vs_control.xlsx")
write.xlsx(diffgenemouse_4,"fig3/diffgene_tam_colon_male_cachexia_vs_control.xlsx")

##fig3E
#cancer:diffgenemouse_1,diffgenemouse_3
#control:diffgenemouse_2,diffgenemouse_4
#male vs female in each cancer and control group
cancer_up_1<-diffgenemouse_1[which(diffgenemouse_1$logFC > 1),]
cancer_up_2<-diffgenemouse_3[which(diffgenemouse_3$logFC > 1),]
cancer_up<-rbind(cancer_up_1,cancer_up_2)
cancer_down_1<-diffgenemouse_1[which(diffgenemouse_1$logFC< -1),]
cancer_down_2<-diffgenemouse_3[which(diffgenemouse_3$logFC< -1),]
cancer_down<-rbind(cancer_down_1,cancer_down_2)

control_up_1<-diffgenemouse_2[which(diffgenemouse_2$logFC > 1),]
control_up_2<-diffgenemouse_4[which(diffgenemouse_4$logFC > 1),]
control_up<-rbind(control_up_1,control_up_2)
control_down_1<-diffgenemouse_2[which(diffgenemouse_2$logFC < -1),]
control_down_2<-diffgenemouse_4[which(diffgenemouse_4$logFC < -1),]
control_down<-rbind(control_down_1,control_down_2)

#maleup
mup_cancerup<-intersect(cancer_up$genename,diffgenemouse_male_up)
mup_cancerdown<-intersect(cancer_up$genename,diffgenemouse_male_down)
#femaleup
fup_cancerup<-intersect(cancer_down$genename,diffgenemouse_female_up)
fup_cancerdown<-intersect(cancer_down$genename,diffgenemouse_female_down)

for (i in c("m","f")) {
  for (m in c("up","down")) {
    print(i)
    print(m)
    dataname<-get(paste0(i,"up_cancer",m))
    genes<-dataname
    ###KEGG
    gene.df <- bitr(genes,fromType="SYMBOL",toType=c("ENTREZID"), OrgDb = org.Mm.eg.db)
    geneid <- gene.df$ENTREZID   
    ekg <- enrichKEGG(gene = geneid,keyType = "kegg",organism= "mmu", qvalueCutoff = 1, pvalueCutoff=1)
    result_kgg <- as.data.frame(ekg)
    if(dim(result_kgg)[1]!=0){
      result_kgg$genename<-"0"
      for (j in c(1:dim(result_kgg)[1])) {
        name<-gene.df$SYMBOL[match(unlist(strsplit(result_kgg$geneID[j],'/')),gene.df$ENTREZID)]
        result_kgg$genename[j]<-paste(name,collapse='/')
      }
    }
    write.xlsx(result_kgg,paste0("fig3/mouse_enrichement/female vs male/",i,"up_cancer",m,"_kegg_dedup.xlsx"))
  }
}

cancer_m_up_kegg<-read.xlsx("fig3/mouse_enrichement/female vs male/mup_cancerup_kegg_dedup.xlsx",sheet=2)
cancer_m_up_kegg$`-log10pval`<- -log10(cancer_m_up_kegg$pvalue)
cancer_m_up_kegg$Description<-factor(cancer_m_up_kegg$Description,
                                                     levels=rev(cancer_m_up_kegg$Description))
cancer_m_up_kegg$x<-"x"
cancer_m_up_kegg$GeneRatio<-sapply(strsplit(cancer_m_up_kegg$GeneRatio, "/"), function(s) as.numeric(s[1]) / as.numeric(s[2]))

cancer_f_up_kegg<-read.xlsx("fig3/mouse_enrichement/female vs male/fup_cancerup_kegg_dedup.xlsx",sheet=2)
cancer_f_up_kegg$`-log10pval`<- -log10(cancer_f_up_kegg$pvalue)
cancer_f_up_kegg$Description<-factor(cancer_f_up_kegg$Description,
                                     levels=rev(cancer_f_up_kegg$Description))
cancer_f_up_kegg$x<-"x"
cancer_f_up_kegg$GeneRatio<-sapply(strsplit(cancer_f_up_kegg$GeneRatio, "/"), function(s) as.numeric(s[1]) / as.numeric(s[2]))

cancer_m_down_kegg<-read.xlsx("fig3/mouse_enrichement/female vs male/mup_cancerdown_kegg_dedup.xlsx",sheet=2)
cancer_m_down_kegg$`-log10pval`<- -log10(cancer_m_down_kegg$pvalue)
cancer_m_down_kegg$Description<-factor(cancer_m_down_kegg$Description,
                                     levels=rev(cancer_m_down_kegg$Description))
cancer_m_down_kegg$x<-"x"
cancer_m_down_kegg$GeneRatio<-sapply(strsplit(cancer_m_down_kegg$GeneRatio, "/"), function(s) as.numeric(s[1]) / as.numeric(s[2]))

cancer_f_down_kegg<-read.xlsx("fig3/mouse_enrichement/female vs male/fup_cancerdown_kegg_dedup.xlsx",sheet=2)
cancer_f_down_kegg$`-log10pval`<- -log10(cancer_f_down_kegg$pvalue)
cancer_f_down_kegg$Description<-factor(cancer_f_down_kegg$Description,
                                       levels=rev(cancer_f_down_kegg$Description))
cancer_f_down_kegg$x<-"x"
cancer_f_down_kegg$GeneRatio<-sapply(strsplit(cancer_f_down_kegg$GeneRatio, "/"), function(s) as.numeric(s[1]) / as.numeric(s[2]))

cancer_f_up_kegg$group<-"fup"
cancer_f_down_kegg$group<-"fdown"
cancer_f_kegg<-rbind(cancer_f_down_kegg,cancer_f_up_kegg)
cancer_f_kegg$group<-factor(cancer_f_kegg$group,c("fup","fdown"))
pdf("fig3E_female.pdf",width =6,height=3)
ggplot(data=cancer_f_kegg,aes(x = x,y=Description))+
  geom_point(aes(size=GeneRatio,color=`-log10pval`))+
  scale_color_continuous(palette =colorRampPalette(c("#FFCCFF","#8844AA"))(50))+
  #scale_x_continuous(expand = c(0,0),breaks = seq(0,6,1),limits = c(0,6))+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.line = element_line(colour ='black', linewidth =1),
        axis.text.x = element_text(colour ='black', size =12),
        #axis.ticks.x = element_line(colour ='black', linewidth = 1),
        axis.title.x = element_text(colour ='black', size =12)
  )+
  facet_grid(.~group)
dev.off()

cancer_m_up_kegg$group<-"mup"
cancer_m_down_kegg$group<-"mdown"
cancer_m_kegg<-rbind(cancer_m_down_kegg,cancer_m_up_kegg)
cancer_m_kegg$group<-factor(cancer_m_kegg$group,c("mup","mdown"))
pdf("fig3E_male.pdf",width =6,height=5)
ggplot(data=cancer_m_kegg,aes(x = x,y=Description))+
  geom_point(aes(size=GeneRatio,color=`-log10pval`))+
  scale_color_continuous(palette =colorRampPalette(c("skyblue4","#003C67FF"))(10))+
  #scale_x_continuous(expand = c(0,0),breaks = seq(0,6,1),limits = c(0,6))+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.line = element_line(colour ='black', linewidth =1),
        axis.text.x = element_text(colour ='black', size =12),
        #axis.ticks.x = element_line(colour ='black', linewidth = 1),
        axis.title.x = element_text(colour ='black', size =12)
  )+
  facet_grid(.~group)
dev.off()
