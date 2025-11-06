library(readr)
library(ggplot2)
library(ggrepel)


data <- read.csv("markers of CD8T_Lamp1_non_cirrHCC vs cirrHCC.csv",row.names =1)
#data <- data[data$p_val_adj < 0.05,]

# 创建一个新的列 `rank`，表示基因按 avg_log2FC 排序后的索引
data$rank <- rank(-data$avg_log2FC, ties.method = "random")  

top <- data[order(-data$avg_log2FC), ][1:7, ]  # 上调基因
bottom <- data[order(data$avg_log2FC), ][1:7, ]  # 下调基因
highlighted_genes <- rbind(top, bottom)

#selected_genes <- c("APOA2")
#highlighted_genes <- data[rownames(data) %in% selected_genes, ]

# 绘制散点图
ggplot(data, aes(x = rank, y = avg_log2FC)) +
  geom_point(alpha = 0.6, color = "#ffbb78") +  # 绘制散点
  geom_text_repel(data = highlighted_genes, aes(label = rownames(highlighted_genes)), 
                  size = 3, color = "black", max.overlaps = Inf) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # 添加 y=0 参考线
  labs(x = "Gene Rank", y = "avg_log2FC", title = "Ranked Scatter Plot of DEGs") +
  theme_bw()+
  theme(aspect.ratio = 2/1)
  
