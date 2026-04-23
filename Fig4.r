###########Fig4#######################
##Author: gwy
rm(list=ls())
library(stringr)
library(tidyverse)
library(ggplot2)
library(forcats)

library(limma)
library(sva)
library(DESeq2)
library(reshape2)
library(sva)

total<-read_tsv("RNAseq_merge_diff.tsv")
total$Time %>% table() %>% as.data.frame()

total$Time<-factor(total$Time,levels = c("Early",
                                         "Late",
                                         "7 days",
                                         "8 days",
                                         "10 days",
                                         "12 days",
                                         "14 days",
                                         "20 days",
                                         "25 days",
                                         "28 days",
                                         "42 days",
                                         "56 days",
                                         "humane-endpoint"
))
total$Tissue %>% table() %>% as.data.frame()


total$Tissue<-factor(total$Tissue,levels = c("qM","dM","taM","gM","sM","BAT",
                                             "Liver","Heart",
                                             "Cerebellum","Neocortex","Hippocampus"))
# ========== 2. 筛选显著差异基因（核心：统一阈值） ==========
# 阈值：FDR < 0.05 + |logFC| > 1（可按需调整，如FDR < 0.01）
deg_df <- total %>%
  filter(FDR < 0.05, abs(logFC) > 1) %>%  # 筛选显著DEG
  mutate(
    # 标记上调/下调：logFC>1=上调，logFC<-1=下调
    regulation = case_when(
      logFC > 1 ~ "up",
      logFC < -1 ~ "down"
    )
  ) %>%
  # 按GeneName+Tissue+Time去重（避免同一基因在同一组织-时间重复计数）
  distinct(GeneName, Tissue, Time, .keep_all = TRUE)

####common up genes####
up_gene_tissue_counts <- deg_df %>%
  # 只关注有明确调控方向的记录（排除可能的NA）
  filter(regulation=="up") %>% dplyr::select(GeneName,Tissue) %>% unique()
up_gene_Time_counts <- deg_df %>%
  # 只关注有明确调控方向的记录（排除可能的NA）
  filter(regulation=="up") %>% dplyr::select(GeneName,Time) %>% unique()


up_gene_by_tissue<-up_gene_tissue_counts  %>% 
  group_by(GeneName) %>% summarise(TissueNumber=n()) %>% dplyr::arrange(desc(TissueNumber))

up_gene_by_Time<-up_gene_Time_counts  %>% 
  group_by(GeneName) %>% summarise(TimeNumber=n()) %>% dplyr::arrange(desc(TimeNumber))

up_dt<-inner_join(up_gene_by_tissue,up_gene_by_Time)
up_dt$Total<-up_dt$TissueNumber+up_dt$TimeNumber

up_dt<-up_dt %>% dplyr::arrange(desc(Total))

up_dt<-up_dt %>% dplyr::filter(Total>=20)


up_final<-inner_join(up_dt,up_gene_tissue_counts)
up_final<-inner_join(up_final,up_gene_Time_counts)

gene_summary <- up_final %>%
  distinct(GeneName, Tissue, Time) %>%  # 去重：每个基因-组织/基因-时间只算1次
  group_by(GeneName) %>%
  summarise(
    TissueCount = n_distinct(Tissue),  # 该基因关联的组织数量（即原TissueNumber）
    TimeCount = n_distinct(Time),      # 该基因关联的时间数量（即原TimeNumber）
    .groups = "drop"
  )

# 第二步：重构数据，为每个组织/时间生成1行（用于堆叠）
# 组织数据（每个组织计1分）
tissue_stack <- up_final %>%
  distinct(GeneName, Tissue) %>%  # 每个基因-组织只保留1行
  group_by(GeneName) %>%
  mutate(
    value = 1,  # 每个组织计1分
    TissueNumber = n()  # 该基因的组织总数（用于验证）
  ) %>%
  ungroup()

# 时间数据（每个时间计1分）
time_stack <- up_final %>%
  distinct(GeneName, Time) %>%  # 每个基因-时间只保留1行
  group_by(GeneName) %>%
  mutate(
    value = 1,  # 每个时间计1分
    TimeNumber = n()  # 该基因的时间总数（用于验证）
  ) %>%
  ungroup()


tissue_stack$Tissue<-factor(tissue_stack$Tissue,levels = c("qM","dM","taM","gM","sM","BAT",
                                                   "Liver","Heart",
                                                   "Cerebellum","Neocortex","Hippocampus"))

time_order<-c("Early",
              "Late",
              "7 days",
              "8 days",
              "10 days",
              "12 days",
              "14 days",
              "20 days",
              "25 days",
              "28 days",
              "42 days",
              "56 days",
              "humane-endpoint"
)
time_stack$Time<-factor(time_stack$Time,levels = time_order)

# ========== 2. 绘制左侧：组织条形图（反向+堆叠计1） ==========
tissue_colors <- c(
  "qM" = "#374E55FF", "taM" = "#7E6148FF", "gM" = "#4DBBD5FF",
  "sM" = "#8A9045FF", "dM" = "#B24745FF", "BAT" = "#79AF97FF",
  "Cerebellum" = "#6A6599FF", "Hippocampus" = "#80796BFF",
  "Neocortex" = "#FFA319FF", "Heart" = "#800000FF", "Liver" = "#91D1C2FF"
)

time_colors <- c(
  "Early"="#F3E5F5",
  "Late"="#E1BEE7",
  "7 days"="#CE93D8",
  "8 days"="#BA68C8",
  "10 days"="#AB47BC",
  "12 days"="#9C27B0",
  "14 days"="#8E24AA",
  "20 days"="#7B1FA2",
  "25 days"="#6A1B9A",
  "28 days"="#4A148C",
  "42 days"="#311B92",
  "56 days"="#1A237E",
  "humane-endpoint"="#0D47A1"
)

dt<-inner_join(time_stack,tissue_stack)

dt$TimeNumber_add_TissueNumber<-dt$TimeNumber+dt$TissueNumber
gene_order<-dt %>% dplyr::select(GeneName,TimeNumber_add_TissueNumber) %>% 
  unique() %>% dplyr::arrange(desc(TimeNumber_add_TissueNumber)) %>% 
  dplyr::pull(GeneName)

tissue_stack$GeneName<-factor(tissue_stack$GeneName,levels = rev(gene_order))
write_csv(dt,"Top Common Up Genes 20260312.csv")
p_tissue <- ggplot(tissue_stack, aes(x = value, y = GeneName, fill = Tissue)) +
  geom_col(width = 0.8, position = "stack") +  # 堆叠：每个组织占1个单位
  # 核心：X轴反向 + 自定义范围（匹配组织总数）
  scale_x_reverse(
    name = "Tissue",
    limits = c(max(tissue_stack$TissueNumber)*1.1, 0),
    breaks = seq(0, max(tissue_stack$TissueNumber), 1)  # 刻度按1递增
  ) +
  # scale_fill_brewer(palette = "Set3", name = "Tissue") +
  scale_fill_manual(values = tissue_colors)+
  labs(y = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(hjust = 1),
    axis.title.x = element_text(hjust = 1),
    panel.grid = element_blank(),
    legend.position = "left",
    plot.margin = margin(10, 0, 10, 10)
  )

# ========== 3. 绘制右侧：时间条形图（堆叠计1） ==========
time_stack$GeneName<-factor(time_stack$GeneName,levels = rev(gene_order))
p_time <- ggplot(time_stack, aes(x = value, y = GeneName, fill = Time)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +  # 堆叠：每个时间占1个单位
  # geom_col(width = 0.8, position = position_stack(reverse = TRUE)) + 
  scale_x_continuous(
    name = "Time",
    limits = c(0, max(time_stack$TimeNumber)*1.1),
    breaks = seq(0, max(time_stack$TimeNumber), 1)  # 刻度按1递增
  ) +
  # scale_fill_brewer(palette = "Paired", name = "Time") +
  scale_fill_manual(values = time_colors)+
  labs(y = "") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "left",
    plot.margin = margin(10, 10, 10, 0)
  )

# ========== 4. 组合图表（共用Y轴） ==========
library(patchwork)
combined_plot <- p_tissue + p_time +
  plot_layout(ncol = 2, widths = c(1, 1.2), guides = "collect") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  plot_annotation(
    title = "Top Common Up Genes",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  )

print(combined_plot)

pdf("Fig4C_UP.pdf",width = 14,height = 7)
print(combined_plot)
dev.off()

####common down genes####
down_gene_tissue_counts <- deg_df %>%
  # 只关注有明确调控方向的记录（排除可能的NA）
  filter(regulation=="down") %>% dplyr::select(GeneName,Tissue) %>% unique()
down_gene_Time_counts <- deg_df %>%
  # 只关注有明确调控方向的记录（排除可能的NA）
  filter(regulation=="down") %>% dplyr::select(GeneName,Time) %>% unique()


down_gene_by_tissue<-down_gene_tissue_counts  %>% 
  group_by(GeneName) %>% summarise(TissueNumber=n()) %>% dplyr::arrange(desc(TissueNumber))

down_gene_by_Time<-down_gene_Time_counts  %>% 
  group_by(GeneName) %>% summarise(TimeNumber=n()) %>% dplyr::arrange(desc(TimeNumber))

down_dt<-inner_join(down_gene_by_tissue,down_gene_by_Time)
down_dt$Total<-down_dt$TissueNumber+down_dt$TimeNumber

down_dt<-down_dt %>% dplyr::arrange(desc(Total))

down_dt<-down_dt %>% dplyr::filter(Total>=15)


down_final<-inner_join(down_dt,down_gene_tissue_counts)
down_final<-inner_join(down_final,down_gene_Time_counts)

gene_summary <- down_final %>%
  distinct(GeneName, Tissue, Time) %>%  # 去重：每个基因-组织/基因-时间只算1次
  group_by(GeneName) %>%
  summarise(
    TissueCount = n_distinct(Tissue),  # 该基因关联的组织数量（即原TissueNumber）
    TimeCount = n_distinct(Time),      # 该基因关联的时间数量（即原TimeNumber）
    .groups = "drop"
  )

# 第二步：重构数据，为每个组织/时间生成1行（用于堆叠）
# 组织数据（每个组织计1分）
tissue_stack <- down_final %>%
  distinct(GeneName, Tissue) %>%  # 每个基因-组织只保留1行
  group_by(GeneName) %>%
  mutate(
    value = 1,  # 每个组织计1分
    TissueNumber = n()  # 该基因的组织总数（用于验证）
  ) %>%
  ungroup()

# 时间数据（每个时间计1分）
time_stack <- down_final %>%
  distinct(GeneName, Time) %>%  # 每个基因-时间只保留1行
  group_by(GeneName) %>%
  mutate(
    value = 1,  # 每个时间计1分
    TimeNumber = n()  # 该基因的时间总数（用于验证）
  ) %>%
  ungroup()


tissue_stack$Tissue<-factor(tissue_stack$Tissue,levels = c("qM","dM","taM","gM","sM","BAT",
                                                           "Liver","Heart",
                                                           "Cerebellum","Neocortex","Hippocampus"))


time_order<-c("Early",
              "Late",
              "7 days",
              "8 days",
              "10 days",
              "12 days",
              "14 days",
              "20 days",
              "25 days",
              "28 days",
              "42 days",
              "56 days",
              "humane-endpoint"
)
time_stack$Time<-factor(time_stack$Time,levels = time_order)


# ========== 2. 绘制左侧：组织条形图（反向+堆叠计1） ==========
tissue_colors <- c(
  "qM" = "#374E55FF", "taM" = "#7E6148FF", "gM" = "#4DBBD5FF",
  "sM" = "#8A9045FF", "dM" = "#B24745FF", "BAT" = "#79AF97FF",
  "Cerebellum" = "#6A6599FF", "Hippocampus" = "#80796BFF",
  "Neocortex" = "#FFA319FF", "Heart" = "#800000FF", "Liver" = "#91D1C2FF"
)

time_colors <- c(
  "Early"="#F3E5F5",
  "Late"="#E1BEE7",
  "7 days"="#CE93D8",
  "8 days"="#BA68C8",
  "10 days"="#AB47BC",
  "12 days"="#9C27B0",
  "14 days"="#8E24AA",
  "20 days"="#7B1FA2",
  "25 days"="#6A1B9A",
  "28 days"="#4A148C",
  "42 days"="#311B92",
  "56 days"="#1A237E",
  "humane-endpoint"="#0D47A1"
)


time_stack
tissue_stack

dt<-inner_join(time_stack,tissue_stack)

dt$TimeNumber_add_TissueNumber<-dt$TimeNumber+dt$TissueNumber
gene_order<-dt %>% dplyr::select(GeneName,TimeNumber_add_TissueNumber) %>% 
  unique() %>% dplyr::arrange(desc(TimeNumber_add_TissueNumber)) %>% 
  dplyr::pull(GeneName)
write_csv(dt,"Top Common Down Genes 20260312.csv")

tissue_stack$GeneName<-factor(tissue_stack$GeneName,levels = rev(gene_order))
p_tissue <- ggplot(tissue_stack, aes(x = value, y = GeneName, fill = Tissue)) +
  geom_col(width = 0.8, position = "stack") +  # 堆叠：每个组织占1个单位
  # 核心：X轴反向 + 自定义范围（匹配组织总数）
  scale_x_reverse(
    name = "Tissue",
    limits = c(max(tissue_stack$TissueNumber)*1.1, 0),
    breaks = seq(0, max(tissue_stack$TissueNumber), 1)  # 刻度按1递增
  ) +
  # scale_fill_brewer(palette = "Set3", name = "Tissue") +
  scale_fill_manual(values = tissue_colors)+
  labs(y = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(hjust = 1),
    axis.title.x = element_text(hjust = 1),
    panel.grid = element_blank(),
    legend.position = "left",
    plot.margin = margin(10, 0, 10, 10)
  )
time_stack$GeneName<-factor(time_stack$GeneName,levels = rev(gene_order))
# ========== 3. 绘制右侧：时间条形图（堆叠计1） ==========
p_time <- ggplot(time_stack, aes(x = value, y = GeneName, fill = Time)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +  # 堆叠：每个时间占1个单位
  # geom_col(width = 0.8, position = position_stack(reverse = TRUE)) + 
  scale_x_continuous(
    name = "Time",
    limits = c(0, max(time_stack$TimeNumber)*1.1),
    breaks = seq(0, max(time_stack$TimeNumber), 1)  # 刻度按1递增
  ) +
  # scale_fill_brewer(palette = "Paired", name = "Time") +
  scale_fill_manual(values = time_colors)+
  labs(y = "") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "left",
    plot.margin = margin(10, 10, 10, 0)
  )

# ========== 4. 组合图表（共用Y轴） ==========
library(patchwork)
combined_plot <- p_tissue + p_time +
  plot_layout(ncol = 2, widths = c(1, 1.2), guides = "collect") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  plot_annotation(
    title = "Top Common Down Genes",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  )

print(combined_plot)

pdf("Fig4C_DOWN.pdf",width = 14,height = 7)
print(combined_plot)
dev.off()

###Remove batch effect###
counts_matrix <- read.table("merge_6dataset_counts_final_log2TPM.txt")
metadata<-read_tsv("merge_6dataset_phenotype_final_time_stage.txt")
metadata<-metadata %>% column_to_rownames(var = "Sample_ID")
common_samples <- intersect(colnames(counts_matrix), rownames(metadata))
counts_matrix <- counts_matrix[, common_samples]
metadata <- metadata[common_samples, ]

metadata$Batch <- factor(metadata$DataSet)
combat_corrected <-ComBat(
  as.matrix(counts_matrix),
  metadata$Batch,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE,
  ref.batch = NULL,
  BPPARAM = bpparam("SerialParam")
)
write.table(combat_corrected,"merge_6dataset_counts_final_log2TPM_ComBat.txt")

###temporal expression pattern###
library(Mfuzz)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ClusterGVis)
library(clusterProfiler)
library(org.Mm.eg.db)  
library(enrichplot) 
library(patchwork)    


annotation_info<-read_tsv("merge_6dataset_phenotype_final_time_stage.txt")
annotation_info$Tissue<-factor(annotation_info$Tissue,levels = c("qM","dM","taM","gM","sM","BAT" ,"Cerebellum","Hippocampus","Neocortex",
                                                                 "Heart" ,  "Liver"                                                      
))
annotation_info$Time<-factor(annotation_info$Time,levels=as.data.frame(table(annotation_info$Time))$Var1) 
annotation_info$Treatment<-factor(annotation_info$Treatment,levels = c("Control","mock injected",
                                                                       "non-cachectic melanoma cells",
                                                                       "cachectic melanoma cells",
                                                                       "LLC","C26","KPC"
))

annotation_info<-annotation_info %>% dplyr::arrange(Treatment,Time_Stage,Tissue,Sex)

annotation_info$Time_Stage <-factor(annotation_info$Time_Stage ,levels = c("Day_0","Early","Middle","Late"))

dt<-annotation_info %>% dplyr::select(Sample_ID,Tissue,Treatment,Sex,Time_Stage) 
tissue_data <- dt %>%
  count(Time_Stage, Tissue, name = "Count") %>%
  mutate(Category = "Tissue", 
         Subcategory = Tissue,
         Value = Count)  # 正数

treatment_data <- dt %>%
  count(Time_Stage, Treatment, name = "Count") %>%
  mutate(Category = "Treatment", 
         Subcategory = Treatment,
         Value = -Count)  # 负数

combined_data <- bind_rows(tissue_data, treatment_data)

time_totals <- combined_data %>%
  group_by(Time_Stage, Category) %>%
  summarise(
    Total = sum(Value),  # 注意这里用Value，Treatment是负值
    Count = sum(abs(Value)),  # 实际样本数量
    .groups = "drop"
  )

tissue_colors <- c(
  "qM" = "#374E55FF", "taM" = "#7E6148FF", "gM" = "#4DBBD5FF",
  "sM" = "#8A9045FF", "dM" = "#B24745FF", "BAT" = "#79AF97FF",
  "Cerebellum" = "#6A6599FF", "Hippocampus" = "#80796BFF",
  "Neocortex" = "#FFA319FF", "Heart" = "#800000FF", "Liver" = "#91D1C2FF"
)

treatment_colors <- c(
  "LLC" = "#D2AF81FF",          
  "C26" = "#709AE1FF",       
  "KPC" = "#FED439FF",   
  "cachectic melanoma cells" = "#FD8CC1FF", 
  "non-cachectic melanoma cells" = "#FEB2D6", 
  "Control" = "#FFECF5",      
  "mock injected" = "#FFE6F2"  
)

p_mirror <- ggplot(combined_data, aes(x = Time_Stage, y = Value, fill = Subcategory)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = abs(Value), y = Value),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  geom_text(data = time_totals,
            aes(x = Time_Stage, y = Total, 
                label = paste0("Total: ", Count),
                vjust = ifelse(Category == "Tissue", -0.5, 1.5)),
            inherit.aes = FALSE,
            size = 4, fontface = "bold", color = "black") +
  geom_hline(yintercept = 0, color = "black", size = 0.8, linetype = "solid") +
  
  scale_fill_manual(
    values = c(tissue_colors, treatment_colors),
    breaks = c(names(tissue_colors), names(treatment_colors)),
    name = "Category"
  ) +
  scale_y_continuous(
    breaks = pretty(combined_data$Value),
    labels = function(x) abs(x), 
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  labs(
    title = "Sample distribution：Tissue and Treatment",
    subtitle = paste("Total:", nrow(dt), "| Top：Tissue | Bottom：Treatment"),
    x = "", 
    y = "Sample size"
  ) +ylim(-460,460)+
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "g ray40"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y =  element_blank(),#element_line(color = "gray90"),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  guides(fill = guide_legend(ncol = 1))

pdf("Fig4A.pdf",width = 6,height = 6)
print(p_mirror)
dev.off()


##===========================normalize=========================================================
expr_mat_new<-read.table("merge_6dataset_counts_final_log2TPM_ComBat.txt")
dim(expr_mat_new)
# [1] 78258   728
expr_mat_new2<-expr_mat_new[,annotation_info$Sample_ID]
dim(expr_mat_new2)
# 78258   728
##=========================== mean by time stage=========================================================
exp_long <- expr_mat_new2 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, 
               names_to = "Sample_ID", 
               values_to = "Expression")

exp_annotated <- exp_long %>%
  left_join(annotation_info %>% 
              dplyr::select(Sample_ID, Time_Stage),
            by = "Sample_ID") %>%
  filter(!is.na(Time_Stage))

stage_means <- exp_annotated %>%
  group_by(Gene, Time_Stage) %>%
  summarise(
    Mean_Expression = mean(Expression, na.rm = TRUE),
    SD_Expression = sd(Expression, na.rm = TRUE),
    N_Samples = n(),
    .groups = 'drop'
  )

stage_means_wide <- stage_means %>%
  dplyr::select(Gene, Time_Stage, Mean_Expression) %>%
  pivot_wider(names_from = Time_Stage, 
              values_from = Mean_Expression) %>%
  as.data.frame()

rownames(stage_means_wide) <- stage_means_wide$Gene
stage_means_wide <- stage_means_wide[, -1]
write.table(stage_means_wide,"Time_stage_expression_geneid.txt",sep="\t",row.names = T,col.names = T)


set.seed(123)
stage_means_wide<-read.table("Time_stage_expression_geneid.txt")
ck <- clusterData(
  obj = stage_means_wide,
  scaleData = TRUE,
  clusterMethod = "mfuzz",
  clusterNum = 8
)

dir.create("mfuzz_cluster_GO_BP_result")
saveRDS(ck,"mfuzz_cluster_GO_BP_result/mfuzz_cluster_result.rds")
for(i in 1:length(ck$cluster.list)){
  genelist<-ck$cluster.list[[i]]
  write_csv(as.data.frame(genelist),paste0("mfuzz_cluster_GO_BP_result/C",i,"_genelist.csv"))
}


for(i in 1:length(ck$cluster.list)){
  C1_go_enrich <- enrichGO(
    gene =ck$cluster.list[[i]],
    OrgDb = org.Mm.eg.db,   
    keyType = "ENSEMBL",    
    ont = "BP",            
    pAdjustMethod = "fdr",  
    qvalueCutoff = 0.05,   
    readable = TRUE         
  )
  write_tsv(C1_go_enrich@result,paste0("mfuzz_cluster_GO_BP_result/C",i,"_BP.tsv"))
}

enrich_res<-tibble()
for(i in 1:length(ck$cluster.list)){
  sub<-read_tsv(paste0("mfuzz_cluster_GO_BP_result/C",i,"_BP.tsv"))
  sub<-top_n(sub,n=5)
  sub$id<-paste0("C",i)
  enrich_res<-rbind(enrich_res,sub)
}
enrich_res<-enrich_res %>% dplyr::select("id","Description")
colnames(enrich_res)<-c("id","term")
head(enrich_res,4)
#   id                               term
# 1 C1              developmental process
# 2 C1   anatomical structure development
# 3 C1 multicellular organism development
# 4 C2                 system development


pdf('Fig4D.pdf',height = 10,width = 10)
visCluster(object = ck,
           plotType = "both",
           column_names_rot = 45,
           annoTermData = enrich_res,
           lineSide  = "left",
           show_row_dend = F)
dev.off()
