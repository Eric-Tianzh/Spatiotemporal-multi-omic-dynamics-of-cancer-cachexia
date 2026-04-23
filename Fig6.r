###########Fig6#######################
##Author: gwy
rm(list=ls())
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(ggpubr)


DEG_all<-read_tsv("RNAseq_merge_diff_allgene.tsv")
creatine_genes<-c("Ckmt1","Ckb","Slc6a8","Gamt","Ckmt2","Gatm","Slc6a12","Slc6a11","Ckm","Slc6a7")
df<-DEG_all %>% dplyr::filter(GeneName %in% creatine_genes)

df_wide <- df %>%
  dplyr::select(GeneName, logFC, Dataset, Tissue, Time) %>%
  pivot_wider(
    names_from = GeneName,     
    values_from = logFC,       
    names_prefix = "logFC_"   
  )

heatmap_data <- df_wide %>%
  mutate(SampleID = paste(Tissue, Time, Dataset, sep = " | ")) %>%
  pivot_longer(
    cols = starts_with("logFC_"),
    names_to = "Gene",
    values_to = "logFC"
  ) %>%
  mutate(Gene = gsub("logFC_", "", Gene)) %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = Gene,
    values_from = logFC
  ) %>%
  column_to_rownames("SampleID") %>%
  as.matrix()

row_annotation <- df_wide %>%
  mutate(SampleID = paste(Tissue, Time, Dataset, sep = " | ")) %>%
  dplyr::select(SampleID, Tissue, Time, Dataset) %>%
  distinct() %>%
  column_to_rownames("SampleID")

row_annotation <- row_annotation[rownames(heatmap_data), ]

row_annotation$Time<-factor(row_annotation$Time,levels = c("Early",
                                         "Late",
                                         "7 days",
                                         "8 days",
                                         "10 days",
                                         "12 days",
                                         "14 days",
                                         "20 days",
                                         "21 days",
                                         "25 days",
                                         "28 days",
                                         "42 days",
                                         "56 days",
                                         "humane-endpoint"
))

row_annotation$Tissue<-factor(row_annotation$Tissue,levels = c("qM","dM","taM","gM","sM","BAT",
                                             "Liver","Heart",
                                             "Cerebellum","Neocortex","Hippocampus"))

row_annotation<-row_annotation %>% dplyr::arrange(Dataset,Tissue,Time)

heatmap_data<-heatmap_data[rownames(row_annotation),]
# Tissue
tissue_colors <- list(
  Tissue = c(
    "qM" = "#374E55FF",  
    "taM" = "#7E6148FF",
    "gM" = "#4DBBD5FF", 
    "sM" = "#8A9045FF",  
    "dM" = "#B24745FF",    
    "BAT" = "#79AF97FF",  
    "Cerebellum" = "#6A6599FF",  
    "Hippocampus" = "#80796BFF",
    "Neocortex" = "#FFA319FF",
    "Heart" = "#800000FF", 
    "Liver" = "#91D1C2FF" 
  ),
  Time=c("Early"="#F3E5F5",
    "Late"="#E1BEE7",
    "7 days"="#CE93D8",
    "8 days"="#BA68C8",
    "10 days"="#AB47BC",
    "12 days"="#9C27B0",
    "14 days"="#8E24AA",
    "20 days"="#7B1FA2", 
    "21 days"="#721D9E", 
    "25 days"="#6A1B9A", 
    "28 days"="#4A148C",
    "42 days"="#311B92",
    "56 days"="#1A237E", 
    "humane-endpoint"="#0D47A1"
  ),
  Dataset = c(
    "GSE114820" = "#1f77b4", 
    "GSE157251" = "#ff7f0e", 
    "GSE183613" = "#2ca02c",
    "GSE222317" = "#d62728", 
    "GSE271521" = "#e377c2", 
    "GSE276018" = "#c5b0d5"
  )
)

sig_df <- df %>%
  mutate(
    Significance = case_when(
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    SampleID = paste(Tissue, Time, Dataset, sep = " | ")
  ) %>%
  dplyr::select(SampleID, GeneName, Significance) %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = GeneName,
    values_from = Significance,
    values_fill = ""
  ) %>%
  column_to_rownames("SampleID")

sig_matrix <- as.matrix(sig_df[rownames(heatmap_data), colnames(heatmap_data)])
sig_matrix<-sig_matrix[rownames(heatmap_data),colnames(heatmap_data)]
sig_matrix_t <- t(sig_matrix)

pheatmap_result <- pheatmap(
  t(heatmap_data), 
  main = "Creatine metobolism gene (logFC with FDR Significance)",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-3, 3, length.out = 101), 

  cluster_rows = TRUE,
  cluster_cols = FALSE, 

  show_colnames = TRUE,
  show_rownames = TRUE,
  
  annotation_col = row_annotation,
  annotation_colors = tissue_colors,

  fontsize_row = 10,
  fontsize_col = 9,
  angle_col = 90,
  
  cellwidth = 18,
  cellheight = 14,
  border_color = NA,

  display_numbers = sig_matrix_t,
  number_color = "black",
  fontsize_number = 9,
  fontface_number = "bold",

  legend = TRUE,
  legend_breaks = c(-3, -1.5, 0, 1.5, 3),
  legend_labels = c("-3.0", "-1.5", "0", "1.5", "3.0")
)

pdf("Fig6A.pdf",width = 20,height = 10)
print(pheatmap_result)
dev.off()



expr_mat_new<-read_tsv("merge_6dataset_counts_final_log2TPM_ComBat_genesymbol_mean.txt")
expr_clean <- expr_mat_new %>% filter(!is.na(GeneName) & GeneName != "" & GeneName != " ")
expr_clean2<-expr_clean %>% column_to_rownames(var = "GeneName")

annotation_info<-read_tsv("merge_6dataset_phenotype_final.txt")
annotation_info<-annotation_info %>% dplyr::filter(DataSet=="GSE183613")
annotation_info$Tissue<-factor(annotation_info$Tissue,levels = c("qM","dM","taM","gM","sM","BAT" ,"Cerebellum","Hippocampus","Neocortex",
                                                                 "Heart" ,  "Liver"                                                      
))
# annotation_info$Time_Stage<-factor(annotation_info$Time_Stage,levels=c("Day_0",  "Early", "Middle", "Late")) 
annotation_info$Treatment<-factor(annotation_info$Treatment,levels = c("Control","mock injected",
                                                                       "non-cachectic melanoma cells",
                                                                       "cachectic melanoma cells",
                                                                       "LLC","C26","KPC"
))

annotation_info$Time<-factor(annotation_info$Time,
                             levels = c(
                               "14 days",      
                               "28 days",    
                               "42 days", 
                               "56 days"
                             ))

##Fig6B & FigS1
creatine_related_gene<-c("Ckmt1", "Ckb","Slc6a8","Gamt","Map4k4","Ckmt2","Gatm","Slc6a12","Slc6a11","Ckm","Slc6a7")
for(i in creatine_related_gene){
  Gamt_expr<-data.frame(Sample_ID=colnames(expr_clean2),Gene=unlist(expr_clean2[i,]))
  tmp<-inner_join(annotation_info,Gamt_expr)

  p<-tmp %>% ggplot(aes(x = Time, y = Gene, color = Treatment)) +
    stat_summary(fun = mean, geom = "line", size = 0.8, 
                 aes(group = Treatment
                 )) +
    stat_summary(fun.data = mean_se, geom = "errorbar", 
                 width = 0.2, alpha = 1) +
    stat_summary(fun = mean, geom = "point", size = 2.5) +
    facet_grid(DataSet ~ Tissue, scales = "free_y") +
    scale_color_manual(values =
                         c(
                            "cachectic melanoma cells" = "#FD8CC1FF", 
                            "non-cachectic melanoma cells" = "#FEB2D6", 
                            "Control" = "gray",    
                            "mock injected" = "gray"
                         )
                       
    ) +
    labs(x = "Time (Days)", y = paste0(i," Expression (log2 TPM)"),
         title = paste0(i," Expression Across Datasets and Tissues")) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=90,hjust = 1, size = 8),
      strip.text = element_text(size = 9),
      legend.position = "top"
    )
  dir.create("single_gene_GSE183613")
  pdf(paste0("single_gene_GSE183613/",i,"_expression_and_time.pdf"),width = 18,height = 5)
  print(p)
  dev.off()
}