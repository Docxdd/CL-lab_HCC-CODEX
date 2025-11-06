library(data.table)
library(dplyr)
library(readxl)
library(readr)
library(tidyr)
library(pheatmap)

All_CN_TP <- read.csv('G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CT_radius50.csv', header = T)
names(All_CN_TP)

#分组展示
g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g<-g[,c('Class','LiverCirrhosis')]
colnames(g)[1]<-'imageid'

All_CN_TP <- left_join(All_CN_TP, g, by = "imageid")
Group_CN_TP <- All_CN_TP[All_CN_TP$LiverCirrhosis=='0',]
###

data<-All_CN_TP
#data<-Group_CN_TP
  
#####CN  的CT组成 heatmap#####
gc_csd_CN <-data %>% select(1:11,34)
names(gc_csd_CN)

celltype_cols <- 1:11  # 也可以使用列名
tissue_avgs <- colMeans(gc_csd_CN[, celltype_cols])

niche_clusters <- gc_csd_CN %>%
  group_by(cluster_kmeans) %>%
  summarise(across(all_of(celltype_cols), mean))

# Convert为matrix以便数值操作
niche_mat <- as.matrix(niche_clusters[,-1])  # 去除第一列的 cluster_kmeans
rownames(niche_mat) <- niche_clusters$cluster_kmeans

# 添加伪计数防止除以0
pseudo <- 1e-5
niche_mat_pseudo <- niche_mat + pseudo
tissue_avgs_pseudo <- tissue_avgs + pseudo

# 按行归一化
row_sums <- rowSums(niche_mat_pseudo)
niche_mat_norm <- sweep(niche_mat_pseudo, 1, row_sums, "/")

# 先将背景也归一化为频率（加了伪计数）
tissue_avgs_norm <- tissue_avgs_pseudo / sum(tissue_avgs_pseudo)

# 计算 log2 fold change
fc_mat <- log2(sweep(niche_mat_norm, 2, tissue_avgs_norm, "/"))

pheatmap(fc_mat,
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(c("#797BB7", "white","#E58579"))(100),
         breaks = seq(-2, 2, length.out = 101),
         main = "Cell Type Composition of CN")

