library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)
library(dplyr)


#转换ID
gene <- read_csv("markers of CD8T_Lamp1_non_cirrHCC vs cirrHCC.csv")
gene <- gene[gene$p_val_adj < 0.05,]
gene <- gene %>%
  mutate(group = ifelse(avg_log2FC > 0, "noncirrHCC", "cirrHCC"))


row.names(gene)<-gene$gene

go_results <- data.frame()
# 循环处理分组
for (group_type in c("noncirrHCC", "cirrHCC")) {
  
  # 筛选基因
  makergene <- gene[gene$group == group_type, ]
  genes <- makergene$gene
  genes <- na.omit(genes)
  
  # 获取基因对应的 Entrez ID
  entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
  entrezIDs <- as.character(entrezIDs)
  
  # GO富集分析
  go <- enrichGO(
    gene = entrezIDs,
    OrgDb = org.Hs.eg.db, 
    #pvalueCutoff = 1, 
    qvalueCutoff = 0.05,
    ont = "BP",  
    readable = TRUE
  )
  
  # 提取结果并添加标记列
  go_result <- go@result
  go_result$group <- group_type
  
  # 取前15个结果
  #go_result <- go_result[1:15, ]
  
  # 合并结果
  go_results <- rbind(go_results, go_result)
}

df<-go_results

write.csv(df,"GOBP pathway gene CD8T_Lamp1_non_cirrHCC vs cirrHCC.csv")

#excel挑选信号及免疫相关通路再进行作图
df <- read_excel("GOBP pathway gene CD8T_Lamp1_non_cirrHCC vs cirrHCC.xlsx", 
                     sheet = "pathway targeted")

ggplot(df, aes(x = zScore, y = Description, color = qvalue, size = Count)) +
  geom_point(aes(shape = group), alpha = 0.8, position = position_dodge(width = 0.6)) +  # 分组错开
  scale_color_gradient(low = "darkcyan", high = "lightcoral") +  # 颜色映射qvalue
  scale_size(range = c(1, 7)) +  # 控制气泡大小
  scale_shape_manual(values = c(16, 1)) +  # 设置不同形状
  labs(x = "zScore", y = "GO BP", color = "qvalue", size = "Count", shape = "Group") +
  theme_bw() +
  theme(
    plot.margin = margin(3, 0, 3, 0, "cm"), # 调整页边距 (上, 右, 下, 左)
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.position = "right",
    aspect.ratio = 2.5/1
  ) 
