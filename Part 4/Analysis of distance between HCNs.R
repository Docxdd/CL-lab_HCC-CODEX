library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggpubr)  
library(FNN)

# 读取数据
datact50<-as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CT_radius50.csv"))
cell<-datact50[,c('imageid','cluster_kmeans','XMin','YMin')]

##去掉样本量少的spot
frequency <- cell %>%
  group_by(imageid) %>%
  summarise(Frequency = n()) %>%
  ungroup()

frequency <- frequency %>%
  filter(Frequency >= 10)

cell<- semi_join(cell, frequency, by = "imageid")

g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g<-g[,c('Class','PatientID_1','HepatitisB','LiverCirrhosis')]

cell <- merge(cell, g, by.x = "imageid", by.y = "Class", all.x = TRUE)
cell<-cell[cell$HepatitisB==1,]
cell<-cell[,-6]#去除乙肝分组列
colnames(cell) <- c("spots", "neighborhood10", "X", "Y","patients","groups")

# 设定列名
tissue_col <- 'spots'  
neigh_col <- 'neighborhood10'  
patient_col <- 'patients'   
group_col <- 'groups'
X <- 'X'
Y <- 'Y'


# 计算每个spot中邻域neigh1和neigh2的平均最小距离
calculate_avg_min_distance <- function(data, n1, n2) {
  data %>%
    filter(neighborhood10 %in% c(n1, n2)) %>%
    group_by(spots) %>%
    # 确保spot中同时存在两种邻域
    filter(n_distinct(neighborhood10) == 2) %>%
    group_modify(~ {
      coords_n1 <- .x %>% filter(neighborhood10 == n1) %>% dplyr::select(X, Y) %>% as.matrix()
      coords_n2 <- .x %>% filter(neighborhood10 == n2) %>% dplyr::select(X, Y) %>% as.matrix()
      
      dist_n1_to_n2 <- get.knnx(coords_n2, coords_n1, k = 1)$nn.dist
      avg_dist_n1_to_n2 <- mean(dist_n1_to_n2)
      
      dist_n2_to_n1 <- get.knnx(coords_n1, coords_n2, k = 1)$nn.dist
      avg_dist_n2_to_n1 <- mean(dist_n2_to_n1)
      
      tibble(
        avg_min_dist_n1_to_n2 = avg_dist_n1_to_n2,
        avg_min_dist_n2_to_n1 = avg_dist_n2_to_n1,
        avg_min_dist = (avg_dist_n1_to_n2 + avg_dist_n2_to_n1) / 2,
        neighborhood_pair = paste(n1, "vs", n2)
      )
    }) %>%
    ungroup()
}

# 计算平均最小距离
result <- calculate_avg_min_distance(cell, 1, 3)

# 添加分组信息
# 获取每个spot的分组信息（确保每个spot只有一个分组）
spot_groups <- cell %>%
  dplyr::select(spots, groups) %>%
  distinct()  # 每个spot只保留一条记录

# 合并结果和分组信息
result_with_groups <- result %>%
  left_join(spot_groups, by = "spots") %>%
  filter(!is.na(groups))  # 移除没有分组信息的spot

write.csv(result_with_groups,'HCN1 & HCN3 avg_min_dist.csv')

# 绘图
ggplot(result_with_groups, aes(x = factor(groups), y = avg_min_dist)) +
  geom_boxplot(aes(fill = factor(groups)), width = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.5) +
  labs(
    title = "",
    x = "group",
    y = "avg_min_distance"
  ) +
  scale_fill_brewer(palette = "Set2", name = "group") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  stat_compare_means(
    method = "wilcox.test",  # 使用t检验
    label = "p.format", # 显示p值格式
    label.x = 1.5,      # 调整位置
    size = 5
  )

