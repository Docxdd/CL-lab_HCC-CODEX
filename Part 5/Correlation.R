library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)

Tcell <- readRDS("T cell.rds")

#Tcell <- Tcell[,Tcell@meta.data[["Celltype"]] %in% c("Treg", "CD8T_Lamp1")]

# 提取metadata
metadata <- Tcell@meta.data

# 提取基因表达矩阵
expression_matrix <- GetAssayData(Tcell, slot = "data")


###计算每个样本中CD8T_Lamp1的RPSA与MKI67相关性##
# 提取x和y的表达量
x<-'ITGA4'
y<-'MKI67'
x_expression <- expression_matrix[x, ]
y_expression <- expression_matrix[y, ]

# 构建数据框
df <- data.frame(
  CellID = colnames(Tcell),  # 细胞ID
  Celltype = metadata$Celltype,  
  x = x_expression,  
  y = y_expression, 
  Individual.Name = metadata$Individual.Name,  # Individual.Name
  Cirrhosis = metadata$Cirrhosis  # Cirrhosis
)
rownames(df) <- df$CellID
df$CellID <- NULL

# 计算每个Individual.Name的x在Treg中的平均表达量
x_avg <- df %>%
  filter(Celltype == "CD8T_Lamp1") %>%
  group_by(Individual.Name) %>%
  summarise(x_avg = mean(x, na.rm = TRUE))

# 计算每个Individual.Name的y在CD8T_Lamp1中的平均表达量
y_avg <- df %>%
  filter(Celltype == "CD8T_Lamp1") %>%
  group_by(Individual.Name) %>%
  summarise(y_avg = mean(y, na.rm = TRUE))

# 合并两个结果
result_df <- merge(x_avg, y_avg, by = "Individual.Name")

# 添加Cirrhosis信息
cirrhosis_info <- df %>%
  select(Individual.Name, Cirrhosis) %>%
  distinct()

result_df <- merge(result_df, cirrhosis_info, by = "Individual.Name")

# 设置行名为Individual.Name
rownames(result_df) <- result_df$Individual.Name
result_df$Individual.Name <- NULL  # 移除Individual.Name列，因为已经设置为行名
result_df$Cirrhosis<-as.factor(result_df$Cirrhosis)

ggplot(result_df, aes(x = x_avg, y = y_avg)) +
  geom_point(size = 3, color = "black", alpha = 1) +  
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +  # 添加回归线
  stat_cor(method = "pearson", 
           label.x.npc = "left", label.y.npc = "top",
           size = 6
  ) + 
  labs(x = x, y = y)+
  theme_classic2() + 
  theme(aspect.ratio = 1/1)

