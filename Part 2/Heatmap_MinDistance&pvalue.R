library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyr)
library(reshape2)
library(tibble)

alldata<-as.data.frame(fread("All_CN_TP_0424.csv"))
cell<-alldata[,c('Class','Allsubtypes','celltype','XMin','YMin')]
cell <- cell %>%
  mutate(across(everything(), ~ gsub("_int", "", .)))

###clinic data
g <- read_excel("clinic data.xlsx")
g<-g[,c('Class','HepatitisB','LiverCirrhosis')]

merged_data <- left_join(cell, g, by = "Class")
merged_data<-merged_data[merged_data$HepatitisB==1,]#去除非乙肝样本

#rename
colnames(merged_data)[7]<-'Group' 
colnames(merged_data)[4]<-'x' 
colnames(merged_data)[5]<-'y' 

merged_data$x<-as.numeric(merged_data$x)
merged_data$y<-as.numeric(merged_data$y)
merged_data$Group <- as.factor(merged_data$Group)

celltypes <- c('CD4T', 'APC', 'B cells', 'Macrophages', 'CD8T','Treg','Endothelial cells','Fibroblasts','Lymphatic endothelial cells','Biliary tract cells')

# 初始化存储结果
results <- data.frame(Class = character(), Group = character(), Min_Distance = numeric(), Celltype1 = character(), Celltype2 = character(), stringsAsFactors = FALSE)

# 获取所有唯一的 Class 和 Group 组合
unique_class_group <- merged_data %>% select(Class, Group) %>% distinct()

# 遍历每个 Class 和 Group 组合
for (i in 1:nrow(unique_class_group)) {
  class <- unique_class_group$Class[i]
  group <- unique_class_group$Group[i]
  
  # 筛选当前 Class 和 Group 的数据
  class_group_data <- merged_data %>% filter(Class == class & Group == group)
  
  # 遍历所有 celltype 组合
  for (m in celltypes) {
    for (n in celltypes) {
      if (m == n) next  # 避免自身计算
      
      # 筛选 m 和 n 细胞的数据
      m_cells <- class_group_data %>% filter(celltype == m) %>% select(x, y)
      n_cells <- class_group_data %>% filter(celltype == n) %>% select(x, y)
      
      # 如果没有细胞，跳过
      if (nrow(m_cells) == 0 || nrow(n_cells) == 0) next
      
      # 计算最小欧几里得距离
      dist_matrix <- as.matrix(dist(rbind(m_cells, n_cells)))
      min_distance <- min(dist_matrix[1:nrow(m_cells), (nrow(m_cells) + 1):nrow(dist_matrix)], na.rm = TRUE)
      
      # 存储结果
      results <- rbind(results, data.frame(Class = class, Group = group, Min_Distance = min_distance, Celltype1 = m, Celltype2 = n, stringsAsFactors = FALSE))
    }
  }
}

# 计算每组的均值
group_means <- aggregate(
       Min_Distance ~ Group + Celltype1 + Celltype2,
       data = results,
       FUN = function(x) mean(x, na.rm = TRUE)
   )
group_means$mean_dist <- group_means$Min_Distance
  
# 生成 a 和 b 数据框
a <- group_means %>%
  filter(Group == 0) %>%
  select(Celltype1, Celltype2, mean_dist) %>%
  pivot_wider(names_from = Celltype2, values_from = mean_dist) %>%
  column_to_rownames(var = "Celltype1")

b <- group_means %>%
  filter(Group == 1) %>%
  select(Celltype1, Celltype2, mean_dist) %>%
  pivot_wider(names_from = Celltype2, values_from = mean_dist) %>%
  column_to_rownames(var = "Celltype1")

# 计算差值，生成数据框 c
c <- b - a
c <- c[celltypes, celltypes]
c <- as.matrix(c)
diag(c) <- 0

# 计算p 值，生成数据框 d
# 提取每个样本的 celltype 之间的最小距离
min_dist_data <- results %>%
  select(Class, Group, Min_Distance, Celltype1, Celltype2)

min_dist_data <- min_dist_data %>%
  mutate(Celltype_Pair = ifelse(Celltype1 < Celltype2, 
                                paste(Celltype1, Celltype2, sep = "_"), 
                                paste(Celltype2, Celltype1, sep = "_")))

# 只保留唯一的 celltype pair 组合
min_dist_data <- min_dist_data %>%
  distinct(Class, Group, Min_Distance, Celltype_Pair)

unique_pairs <- unique(min_dist_data$Celltype_Pair)

#检验计算p值
p_values_all <- data.frame(Celltype_Pair = character(), p_value = numeric(), stringsAsFactors = FALSE)

for (pair in unique_pairs) {
  # 过滤出当前 pair 的数据
  pair_data <- min_dist_data %>%
    filter(Celltype_Pair == pair)
  
  # 只有当 Group 至少有两个不同值时才检验
  if (length(unique(pair_data$Group)) > 1) {
    p_value <- wilcox.test(Min_Distance ~ Group, data = pair_data)$p.value
  } else {
    p_value <- NA  # 如果 Group 只有一个类别，p 值设为 NA
  }
  
  # 记录结果
  p_values_all <- rbind(p_values_all, data.frame(Celltype_Pair = pair, p_value = p_value))
}

# 拆分 Celltype_Pair 为 Celltype1 和 Celltype2
p_values_all <- p_values_all %>%
  separate(Celltype_Pair, into = c("Celltype1", "Celltype2"), sep = "_")

# 生成行为 Celltype1，列为 Celltype2 的矩阵
celltypes <- unique(c(p_values_all$Celltype1, p_values_all$Celltype2))
p_value_matrix <- matrix(NA, nrow = length(celltypes), ncol = length(celltypes), 
                         dimnames = list(celltypes, celltypes))

# 填充矩阵
for (i in 1:nrow(p_values_all)) {
  c1 <- p_values_all$Celltype1[i]
  c2 <- p_values_all$Celltype2[i]
  p_value_matrix[c1, c2] <- p_values_all$p_value[i]
  p_value_matrix[c2, c1] <- p_values_all$p_value[i]  # 对称填充
}
p_value_df <- as.data.frame(p_value_matrix)
p <- p_value_df[celltypes, celltypes]
p <- as.matrix(p)
diag(p) <- 1


# 将矩阵转换为长格式并显式添加行列名称
c_long <- reshape2::melt(c)
p_long <- reshape2::melt(p)

# 给转换后的数据框命名，并显式添加行列名作为变量
colnames(c_long) <- c("CellType1", "CellType2", "Distance")
colnames(p_long) <- c("CellType1", "CellType2", "PValue")

# 合并 c 和 p 数据框
df <- merge(c_long, p_long, by = c("CellType1", "CellType2"))

# 计算 log2 标准化后的距离
df <- df %>%
  mutate(log2_Distance = sign(Distance) * log2(abs(Distance) + 0.001),  # 处理 Distance
         PLabel = case_when(
           PValue < 0.001 ~ "***",
           PValue < 0.01  ~ "**",
           PValue < 0.05  ~ "*",
           TRUE ~ ""
         ))

# 绘制热图
ggplot(df, aes(x = CellType1, y = CellType2, fill = log2_Distance)) +
  geom_tile() +
  geom_text(aes(label = PLabel), size = 5, color = "black") +  # 显示星号
  scale_fill_gradient2(low = "#304d9b", mid = "white", high = "#ea514b") +  # 负值红色，正值绿色
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    aspect.ratio = 1/2.5
  ) +
  labs(title = "Heatmap of log2-transformed min distance difference",
       fill = "Signed log2(Distance)")
