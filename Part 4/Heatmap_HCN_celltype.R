library(data.table)
library(dplyr)
library(readxl)
library(tidyr)
library(pheatmap)

####读入CN_CT_radius50#########
datact50 <- as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CT_radius50.csv"))
datact50 <- datact50[,c('cellid','imageid','cluster_kmeans')]
colnames(datact50)[3]<-'HCN'

g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g <- g[ ,c('Class','LiverCirrhosis','HepatitisB')]
colnames(g)[1] <- 'imageid'

datact50 <- datact50[datact50$imageid %in% g$imageid[g$HepatitisB==1], ]
g <- g[ ,-3]

datact50 <- left_join(datact50,g, by = 'imageid')

datacp50<-as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CP_radius50.csv"))
datacp50<-datacp50[,c(76,1:72)] #cellid,subcelltype proportion

data <- left_join(datact50, datacp50, by = 'cellid')

# 筛选出目标细胞
cell <- data %>%filter(HCN == 1)

# 提取特征列名
feature_cols <- colnames(cell)[!(colnames(cell) %in% c("cellid", "imageid", "HCN", "LiverCirrhosis"))]

# 找出包含 CD8T 和 CD4T 的列
Tcell_cols <- feature_cols[grepl("CD8T|CD4T", feature_cols)]
CD8T_cols <- feature_cols[grepl("CD8T", feature_cols)]

# 计算总 T cell
cell$total_Tcell <- rowSums(cell[, Tcell_cols], na.rm = TRUE)

# 为每个 CD8T 亚型计算比例
for (col in CD8T_cols) {
  new_col <- paste0(col)
  cell[[new_col]] <- cell[[col]] / cell$total_Tcell
}

# 获取新生成的 CD8T 比例列
rel_cols <- grep("^CD8T",names(cell), value = TRUE)

merged_wide <- cell %>%
  group_by(imageid) %>%
  summarise(
    LiverCirrhosis = first(LiverCirrhosis),
    across(all_of(rel_cols), mean, na.rm = TRUE),
    .groups = "drop"
  )

group0 <- merged_wide[merged_wide$LiverCirrhosis == 0, ]
group0 <- group0[ ,-2] #删除LiverCirrhosis列
group1 <- merged_wide[merged_wide$LiverCirrhosis == 1, ]
group1 <- group1[ ,-2] #删除LiverCirrhosis列

rname0<-group0$imageid
group0<-group0[,-1]
row.names(group0)<-rname0
group0<-t(group0)
  
rname1<-group1$imageid
group1<-group1[,-1]
row.names(group1)<-rname1
group1<-t(group1)


# 对每组分别聚类
heatmap_group0 <- pheatmap(group0, scale = "none", clustering_method = "complete", 
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean", 
                           silent = TRUE)
heatmap_group1 <- pheatmap(group1, scale = "none", clustering_method = "complete", 
                           clustering_distance_rows = "euclidean", 
                           clustering_distance_cols = "euclidean", 
                           silent = TRUE)
# 重新排列数据框，使每组样本按聚类结果排序
group0_sorted <- group0[heatmap_group0$tree_row$order, heatmap_group0$tree_col$order]
group1_sorted <- group1[heatmap_group1$tree_row$order, heatmap_group1$tree_col$order]

all <- merge(group0_sorted, group1_sorted, by = "row.names", all = TRUE)
rownames(all) <- all$Row.names
all <- all[, -1] # 去掉合并生成的行名列

all<- t(all)

# 创建分组信息
group_annotation <- data.frame(
       Sample = rownames(all),
       Group = rep(c("Group0", "Group1"),times = c(ncol(group0), ncol(group1)))
   )


# 设置行名为样本名
rownames(group_annotation) <- group_annotation$Sample

# 去掉Sample列，保留Group列
group_annotation <- group_annotation[, "Group", drop = FALSE]

# 设置分组条的颜色
annotation_colors <- list(
  Group = c("Group0" = "#193E8F", "Group1" = "#F09739") 
)

# 绘制聚类热图
pdf("heatmap of sub_CD8T in HCN1.pdf", width = 10, height = 8)  # 设置 PDF 文件的宽度和高度

pheatmap(all, 
         scale = "column", # "none", "row", "column"
         clustering_method = "complete", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         annotation_row = group_annotation, ## 注意这里使用 annotation_row 还是 annotation_col
         annotation_colors = annotation_colors,
         show_colnames = T, 
         show_rownames = F,
         color = colorRampPalette(c("#283277","white", "#fc8452"))(100), # 颜色渐变
         breaks = seq(-1, 1, length.out = 100), # 扩展颜色范围，增加色差
         legend_breaks = seq(-1, 1, 0.5), # 设置图例刻度
         cluster_cols = T, # 列聚类
         cluster_rows = F, # 行聚类
         fontsize = 8,         # 调整字体大小
         fontsize_col = 7,     # 调整列标签字体大小
         angle_col = 90,       # 旋转列标签45度
         cellheight = 1.75,    # 格子高度
         cellwidth = 8,        # 格子宽度
         treeheight_col = 5    # 列聚类树高度
)

dev.off()


