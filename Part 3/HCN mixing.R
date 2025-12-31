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

##删除neighborhood10过少的样本
#tab <- table(cell$spots, cell$neighborhood10)
#spots_to_remove <- rownames(tab)[apply(tab, 1, function(x) any(x < 50))]
#cell <- cell[!cell$spots %in% spots_to_remove, ]

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
    select(all_of(c(tissue_col, group_col, patient_col)), value)
  
  val2 <- neighs_long %>% 
    filter(!!sym(neigh_col) == neigh2, neigh_neigh == neigh1) %>%
    select(all_of(c(tissue_col, group_col, patient_col)), value)
  
  inters <- full_join(val1, val2, by = c(tissue_col, group_col, patient_col)) %>%
    rowwise() %>%
    mutate(inters = mean(c_across(starts_with("value")), na.rm = TRUE)) %>%
    ungroup() %>%
    select(all_of(c(tissue_col, group_col)), inters)
  
  total <- counts %>%
    filter(!!sym(neigh_col) %in% c(neigh1, neigh2)) %>%
    group_by(across(all_of(c(group_col, patient_col, tissue_col)))) %>%
    summarise(total = sum(n), .groups = "drop") %>%
    select(all_of(c(tissue_col, group_col)), total)
  
  ratio_df <- left_join(inters, total, by = c(tissue_col, group_col)) %>%
    mutate(ratio = inters / total) %>%
    select(all_of(c(tissue_col, group_col)), ratio)
  
  colname <- paste0("HCN", neigh1, "_HCN", neigh2)
  names(ratio_df)[3] <- colname
  ratio_df
}

# 批量处理所有组合并合并结果
results_list <- map(pairs, ~ get_ratio(.x[1], .x[2]))
final_result <- reduce(results_list, full_join, by = c(tissue_col, group_col))

write.csv(final_result,'Mixing Score_无过滤.csv')

###批量做箱线图###
# 转换分组为因子
final_result$groups <- as.factor(final_result$groups)

# 所有要绘图的列（排除 spots 和 groups）
y_variables <- setdiff(colnames(final_result), c("spots", "groups"))

plot_list <- list()

for (y_var in y_variables) {
  
  # 如果整列都是 NA，则跳过
  if (all(is.na(final_result[[y_var]]))) next
  
  # 选出当前变量和 groups 列，去除 NA
  df <- final_result[, c("groups", y_var)] %>% drop_na()
  
  # 删除离散值
  lower <- quantile(df[[y_var]], 0.15, na.rm = TRUE)
  upper <- quantile(df[[y_var]], 0.85, na.rm = TRUE)
  df_filtered <- df %>%
    filter(.data[[y_var]] >= lower & .data[[y_var]] <= upper)
  
  # 跳过去除后数据太少的情况
  if (nrow(df_filtered) < 3) next
  
  # 绘图
  p <- ggplot(df_filtered, aes(
    x = groups,
    y = .data[[y_var]],
    fill = groups,
    color = groups
  )) +
    geom_boxplot(outlier.shape = NA,
                 position = position_dodge(width = 0.3),
                 width = 0.3, fill = "white") +
    scale_fill_manual(values = c("0" = "darkcyan", "1" = "lightcoral")) +
    scale_color_manual(values = c("0" = "darkcyan", "1" = "lightcoral")) +
    labs(x = "Cirrhosis Group", y = y_var) +
    theme_bw() +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      comparisons = list(c("0", "1"))
      #label.y = max(df_filtered[[y_var]], na.rm = TRUE) * 1.05
    ) +
    theme(aspect.ratio = 2/1)
  
  plot_list[[y_var]] <- p
}


pdf_file <- "Mixing score_boxplots_wilcox_test_无过滤_作图是删除前后0.15离群值20250624.pdf"
plots_per_row <- 4
plots_per_page <- 8
num_plots <- length(plot_list)
num_pages <- ceiling(num_plots / plots_per_page)

pdf(pdf_file, width = 11.69, height = 8.27)

for (page in seq_len(num_pages)) {
  start_idx <- (page - 1) * plots_per_page + 1
  end_idx <- min(start_idx + plots_per_page - 1, num_plots)
  
  plots_to_draw <- plot_list[start_idx:end_idx]
  num_rows <- ceiling(length(plots_to_draw) / plots_per_row)
  
  grid_plot <- ggarrange(plotlist = plots_to_draw,
                         ncol = plots_per_row,
                         nrow = num_rows)
  
  print(grid_plot)
}

dev.off()
cat("PDF 已保存到: ", pdf_file, "\n")


###生存分析
library(autoReg)
library(survival)
library(officer)
library(flextable)

g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g<-g[,c('Class','OS','Osday')]
colnames(g)[1]<-'spots'
colnames(g)[2]<-'fustat'
colnames(g)[3]<-'futime'


surdata <- left_join(final_result,g,by = 'spots')
surdata <- surdata[ ,c(3:ncol(surdata))]

col_names <- names(surdata)
cols_to_rename <- !(col_names %in% c("futime", "fustat"))
names(surdata) <- col_names

# 对这些列进行 log1p 转换
feature_cols <- setdiff(colnames(surdata), c("futime", "fustat"))
surdata[feature_cols] <- lapply(surdata[feature_cols], function(x) log(x + 0.01))

coxmod <- coxph(Surv(futime, fustat) ~ ., data = surdata)
summary(coxmod)

ft <- autoReg(coxmod,
              uni = TRUE,        # 单因素
              multi = TRUE,      # 多因素
              threshold = 0.05,     # 变量纳入多因素条件
              final = F)     # 是否逐步回归

myft(ft)

#保存为excel
library(openxlsx)
write.xlsx(ft, "cox_result.xlsx", rowNames = FALSE)

#保存为表格
ft_flex <- as_flextable(ft)
doc <- read_docx() %>%
  body_add_flextable(ft_flex)
print(doc, target = "cox_result.docx")


##森林图##
library(forestploter)
library(dplyr)
library(stringr)

ft <- ft %>%
  filter(`HR (multivariable)` != "")

df <- ft %>%
  select(name, `HR (univariable)`) %>%
  rename(label = name, hr_str = `HR (univariable)`) %>%
  filter(!is.na(hr_str), hr_str != "") %>%
  mutate(
    hr_str = str_replace_all(hr_str, "[–—]", "-"),
    HR = as.numeric(str_extract(hr_str, "^[0-9.]+")),
    lower = as.numeric(str_extract(hr_str, "(?<=\\()[0-9.]+")),
    upper = as.numeric(str_extract(hr_str, "(?<=-)[0-9.]+(?=,)")),
    p.value = as.numeric(str_extract(hr_str, "(?<=p[=|<])[0-9.eE-]+"))
  )

mean <- df$HR
lower <- df$lower
upper <- df$upper

tabletext <- cbind(
  df$label,
  sprintf("%.2f", df$HR),
  sprintf("%.2f–%.2f", df$lower, df$upper),
  sprintf("%.3f", df$p.value)
)
colnames(tabletext) <- c("Variable", "HR", "95% CI", "p-value")

df_plot <- as.data.frame(tabletext, stringsAsFactors = FALSE)
df_plot$forest <- paste(rep("  ", 15), collapse = " ")

tm <- forest_theme(
  base_size = 9,
  # 设置可信区间的外观
  ci_pch = c(15),           # 可信区间点的形状
  ci_col = c("#D97B77"),    # 可信区间的边框颜色
  ci_fill = c("#D97B77"),      # 可信区间的填充颜色
  ci_alpha = 0.8,        # 可信区间的透明度
  ci_lty = 1,            # 可信区间的线型
  ci_lwd = 1.5,          # 可信区间的线宽
  ci_Theight = 0.2,      # 设置T字在可信区间末端的高度，默认是NULL
  
  # 设置参考线的外观
  refline_lwd = 1,         # 参考线的线宽
  refline_lty = "dashed",  # 参考线的线型
  refline_col = "black",  # 参考线的颜色
  
  # 设置垂直线的外观
  vertline_lwd = 1,         # 垂直线的线宽，可以添加一条额外的垂直线，如果没有就不显示
  vertline_lty = "dashed",  # 垂直线的线型
  vertline_col = "black",  # 垂直线的颜色
  
  # 设置脚注的字体大小、字体样式和颜色
  footnote_cex = 0.6,            # 脚注字体大小
  footnote_fontface = "italic",  # 脚注字体样式
  footnote_col = "red4" ,         # 脚注文本的颜色
  
  #图例
  #legend_name = "Group",          # 图例标题
  #legend_value = c()   # 图例值
  
  #plot_margin = unit(c(8,5, 8, 5), "cm")  # 设置页边距
)

forest(
  df_plot[, c("Variable", "HR", "95% CI", "forest","p-value")],
  est = mean,
  lower = lower,
  upper = upper,
  sizes = 0.4,          # 黑框的大小
  ci_column = 4,       # 置信区间列的索引
  ref_line = 1,        # 参考线位置
  vert_line = c(1),       # 添加垂直线
  arrow_lab = c("Lower Risk", "Higher Risk"),  # 箭头标签
  #xlim = c(0.5, 3),    # x轴范围
  #ticks_at = c(0.5, 1,1.5, 2, 2.5,3),  # x轴刻度
  xlab = "Hazard Ratio (95% CI)",
  ticks_at = pretty(c(lower, upper), n = 5),
  xlim = c(min(lower, na.rm = TRUE) * 0.8,
           max(upper, na.rm = TRUE) * 1.2),
  footnote = "",  # 脚注
  col = c("#D97B77"), 
  fill = c("#D97B77"),
  #graphwidth = unit(50, "cm"),  # 设置图形的宽度
  #colgap = unit(10, "mm"),       # 增加列之间的间隔
  theme = tm
)


# wilcox 检验
pval <- wilcox.test(HCN0_HCN1 ~ groups, data = final_result)$p.value
p_label <- sprintf("%.3f", pval)

##作箱线图
y_max <- max(final_result$HCN0_HCN1, na.rm = TRUE) * 1.05

# 作图
ggviolin(final_result,
         x = "groups",
         y = "HCN0_HCN1",
         color = "groups",
         fill = "groups",  # 小提琴图填充色
         ylab = "HCN0_HCN1",
         add = "boxplot",
         add.params = list(fill = NA,  # 透明填充
                           color = "black", 
                           width = 0.1),
         bxp.errorbar.width = 0.2,
         width = 0.4,
         size = 0.5,
         notch = FALSE,
         outlier.shape = NA,
         font.label = list(size = 30),
         position = position_dodge(0.1),
         palette = c("#304d9b", "#ea514b")) +
  annotate("text", x = 1.5, y = y_max, label = p_label, size = 5) +  # 添加 p 值
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 0, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.position = "none")



tmp1  <- final_result[ ,c('spots', 'groups','HCN0_HCN1')] %>%
  pivot_longer(cols = -c(spots, groups), # 指定要转换为长格式的列
               names_to = "HCN", # 新列的名称，即评分的名称
               values_to = "mixing")%>%
  group_by(HCN) %>%
  filter(
    !is.na(mixing),
    mixing >= quantile(mixing, 0.15, na.rm = TRUE),
    mixing <= quantile(mixing, 0.85, na.rm = TRUE)
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
  geom_boxplot(position = position_dodge(width = dodge_width), width = 0.1, size = 0.1, 
               alpha = 0, outlier.shape = NA, color = "black") +
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


