library(data.table)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(FNN)
library(purrr)

# 读取数据
datact50<-as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CT_radius50.csv"))
cell<-datact50[,c('imageid','cluster_kmeans','XMin','YMin')]

#读取临床数据
g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g <- g[,c('Class','PatientID_1','HepatitisB','LiverCirrhosis')]

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

cells2 <- cell
# 找最近邻的邻域编号
cells2 <- cells2 %>%
  group_by_at(tissue_col) %>%
  group_modify(~ {
    knn_idx <- get.knn(data.frame(.x[[X]], .x[[Y]]), k = 1)$nn.index
    .x$neigh_neigh <- .x[[neigh_col]][knn_idx]
    .x
  }) %>%
  ungroup()%>%
  filter(!is.na(neigh_neigh)) 

# 按 group, patient, tissue, neighborhood 统计数量
counts <- cells2 %>%
  group_by(across(all_of(c(group_col, patient_col, tissue_col, neigh_col)))) %>%
  summarise(n = n(), .groups = "drop")

# 统计每个 neighborhood 邻接的其他 neigh_neigh 数量
neighs <- cells2 %>%
  group_by(across(all_of(c(group_col, patient_col, tissue_col, neigh_col))), neigh_neigh) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = neigh_neigh, values_from = n, values_fill = 0)

# 转换为 long 格式
neighs_long <- neighs %>%
  pivot_longer(-all_of(c(group_col, patient_col, tissue_col, neigh_col)), 
               names_to = "neigh_neigh", values_to = "value") %>%
  mutate(neigh_neigh = as.integer(neigh_neigh))

# 所有组合 (neigh1, neigh2)，只保留 neigh1 < neigh2
pairs <- combn(0:9, 2, simplify = FALSE)

# 定义获取 ratio 的函数
get_ratio <- function(neigh1, neigh2) {
  val1 <- neighs_long %>% 
    filter(!!sym(neigh_col) == neigh1, neigh_neigh == neigh2) %>%
    dplyr::select(all_of(c(tissue_col, group_col, patient_col)), value)
  
  val2 <- neighs_long %>% 
    filter(!!sym(neigh_col) == neigh2, neigh_neigh == neigh1) %>%
    dplyr::select(all_of(c(tissue_col, group_col, patient_col)), value)
  
  inters <- full_join(val1, val2, by = c(tissue_col, group_col, patient_col)) %>%
    rowwise() %>%
    mutate(inters = mean(c_across(starts_with("value")), na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(all_of(c(tissue_col, group_col)), inters)
  
  total <- counts %>%
    filter(!!sym(neigh_col) %in% c(neigh1, neigh2)) %>%
    group_by(across(all_of(c(group_col, patient_col, tissue_col)))) %>%
    summarise(total = sum(n), .groups = "drop") %>%
    dplyr::select(all_of(c(tissue_col, group_col)), total)
  
  ratio_df <- left_join(inters, total, by = c(tissue_col, group_col)) %>%
    mutate(ratio = inters / total) %>%
    dplyr::select(all_of(c(tissue_col, group_col)), ratio)
  
  colname <- paste0("HCN", neigh1, "_HCN", neigh2)
  names(ratio_df)[3] <- colname
  ratio_df
}

# 批量处理所有组合并合并结果
results_list <- map(pairs, ~ get_ratio(.x[1], .x[2]))
final_result <- reduce(results_list, full_join, by = c(tissue_col, group_col))

write.csv(final_result,'Mixing Score.csv')

##作图
tmp1  <- final_result[ ,c('spots', 'groups','HCN0_HCN1')] %>%
  pivot_longer(cols = -c(spots, groups), # 指定要转换为长格式的列
               names_to = "HCN", # 新列的名称，即评分的名称
               values_to = "mixing")%>%
  group_by(HCN) %>%
  filter(
    !is.na(mixing),
    mixing >= quantile(mixing, 0.15, na.rm = TRUE),#去除离散值
    mixing <= quantile(mixing, 0.85, na.rm = TRUE) #去除离散值
  ) %>%
  ungroup()

tmp1$groups <- as.factor(tmp1$groups)

p_values <- tmp1 %>%
  group_by(HCN) %>%
  summarise(p = wilcox.test(mixing ~ groups)$p.value)

# 保留三位小数
p_values$p_label <- sprintf("%.3f", p_values$p)

#可以使用每组 max(Frequency) 来作为标签的 y 坐标
label_df <- tmp1 %>%
  group_by(HCN) %>%
  summarise(y_pos = max(mixing) * 1.1) %>%
  left_join(p_values, by = "HCN")

dodge_width <- 0.8

ggplot(tmp1, aes(x = HCN, y = mixing, fill = groups)) + 
  geom_violin(position = position_dodge(width = dodge_width), alpha = 0.5, size = 0.1) +
  #geom_boxplot(position = position_dodge(width = dodge_width), width = 0.1, size = 0.1, 
               #alpha = 0, outlier.shape = NA, color = "black") +
  geom_text(data = label_df, aes(x = HCN, y = 0.006, label = p_label), 
            inherit.aes = FALSE, size = 3) +  # 手动添加 p 值标签
  scale_fill_manual(values = c("#304d9b", "#ea514b")) +
  labs(x = NULL, y = NULL) +
  facet_wrap(~ HCN, scales = "free", nrow = 2) +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.13))) +
  theme(
    plot.margin = unit(c(0, 0.5, 0, 0.5), 'cm'),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 6, color = "black"), 
    axis.text = element_text(size = 6, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 6),
    legend.position = "top",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    panel.spacing.x = unit(0.5, "lines"),
    strip.text = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    aspect.ratio = 2/1
  )
