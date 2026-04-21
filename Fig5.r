###########Fig5########################
##Author: tzh


metadata_all<-read.xlsx("metabolism_diff_data_tzh_251105.xlsx",sheet=4)
metadata_all<-metadata_all[which(metadata_all$group=="cachexia vs control"),]
metadata_all<-metadata_all[which(metadata_all$TaxonomyID=="Homo sapiens"),]
metadata_all$pValue<-as.numeric(metadata_all$pValue)
metadata_all$log2FC_define<-as.numeric(metadata_all$log2FC_define)
metadata_all<-metadata_all[-grep("^C.*H.*",metadata_all$Full.name),]

metadata_all_diff<-metadata_all[which(abs(metadata_all$log2FC_define)>0 & metadata_all$pValue<0.05),]

human_blood_up<-unique(metadata_all$Full.name[which(metadata_all$log2FC_define>0 & metadata_all$pValue<0.05 & 
                                                 metadata_all$Tissue=="Blood (serum)")])
human_blood_down<-unique(metadata_all$Full.name[which(metadata_all$log2FC_define<0 & metadata_all$pValue<0.05 & 
                                                        metadata_all$Tissue=="Blood (serum)")])
human_muscle_up<-unique(metadata_all$Full.name[which(metadata_all$log2FC_define>0 & metadata_all$pValue<0.05 & 
                                                        metadata_all$Tissue=="Muscle")])
human_muscle_down<-unique(metadata_all$Full.name[which(metadata_all$log2FC_define<0 & metadata_all$pValue<0.05 & 
                                                        metadata_all$Tissue=="Muscle")])

v<-venn.diagram(x=list(human_blood_up,
                       human_blood_down,
                       human_muscle_up,
                       human_muscle_down),
                output=FALSE,filename = NULL,imagetype = 'png',
                fill=c("#B24745FF","#155F83FF","#F39B7FFF","#003C67FF"),
                alpha=0.8,scale=F,
                lwd=1,lty=1,col=c("#B24745FF","#155F83FF","#F39B7FFF","#003C67FF"),
                label.col ='black',cex=2,fontface = "bold",
                category.names = c("human_blood_up","human_blood_down",
                                   "human_muscle_up","human_muscle_down"),
                cat.pos = 0,cat.cex = 1, cat.fontface = "bold",
                cat.default.pos = "outer")

pdf('fig5A.pdf', width = 7, height = 7)
grid.draw(v)
dev.off()

listInput <- list(human_blood_up=human_blood_up,human_blood_down=human_blood_down,
                  human_muscle_up=human_muscle_up,human_muscle_down=human_muscle_down)
inter<-get.venn.partitions(listInput)
inter<-inter[order(inter$..count..,decreasing = T),]
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter<-inter[order(inter$..count..,decreasing = T),]

metadown<-intersect(human_blood_down,human_muscle_down)
metadown<-sort(metadown,decreasing = F)
metaup<-intersect(human_blood_up,human_muscle_up)
metaup<-sort(metaup,decreasing = F)
meta<-c(metaup,metadown)
metadata_all_diff_blood<-metadata_all_diff[which(metadata_all_diff$Tissue=="Blood (serum)"),]
metadata_all_diff_muscle<-metadata_all_diff[which(metadata_all_diff$Tissue=="Muscle"),]

metaFC<-data.frame("Human_blood"=metadata_all_diff_blood$log2FC_define[match(meta,metadata_all_diff_blood$Full.name)],
                   "Human_muscle"=metadata_all_diff_muscle$log2FC_define[match(meta,metadata_all_diff_muscle$Full.name)])
rownames(metaFC)<-meta

metapvalue<-data.frame("Human_blood"=metadata_all_diff_blood$pValue[match(meta,metadata_all_diff_blood$Full.name)],
                       "Human_muscle"=metadata_all_diff_muscle$pValue[match(meta,metadata_all_diff_muscle$Full.name)])
rownames(metapvalue)<-meta

sign_level <- function(p,f) {
  stars <- ""
  stars[p < 0.001 & abs(f) > 0] <- "***"
  stars[p >= 0.001 & p < 0.01 & abs(f) > 0] <- "**"
  stars[p >= 0.01 & p < 0.05 & abs(f) > 0] <- "*"
  return(stars)
}
star_matrix <- matrix(sign_level(metapvalue,metaFC),nrow = nrow(metapvalue))
star_matrix[is.na(star_matrix)] <- ""

pdf("fig5B.pdf",height=8,width = 5)
pheatmap(metaFC,
         display_numbers = star_matrix,
         cellwidth = 30, 
         cellheight = 20,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_number = 22,
         fontsize_row = 13, 
         fontsize_col = 13,
         color = colorRampPalette(c("#003C67FF", "white", "#B24745FF"))(500),
         breaks=seq(-2,2, length.out = 501)
         ) 
dev.off()

