library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(msigdbr)
library(AUCell)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(readr)
library(GSEABase)
library(dplyr)

genesets <- read_csv("genesets.csv")
genesets <- lapply(genesets, function(x) na.omit(x) |> as.character())

#读入数据
Tcell <- readRDS("T cell.rds")
aucdata<-Tcell[,Tcell@meta.data[["Celltype"]]%in%c('CD8T_Lamp1')]

meta <- aucdata@meta.data
sample_counts <- table(meta$Individual.Name)
valid_samples <- names(sample_counts[sample_counts >= 50]) #选取细胞数大于50的样本
aucdata <- subset(aucdata, subset = Individual.Name %in% valid_samples)

##AUCell
cells_AUC <- AUCell_run(aucdata@assays$RNA@data,genesets)

# 获取所有基因集的名称
gene_set_names <- names(genesets)

# 循环遍历每个基因集，提取评分并添加到 metadata 中
for (gene_set in gene_set_names) {
  # 提取当前基因集的 AUCell 评分
  AUCell_auc <- as.numeric(getAUC(cells_AUC)[gene_set, ])
  
  # 将评分添加到 aucdata 的 metadata 中
  aucdata[[gene_set]] <- AUCell_auc
}

# 提取元数据
meta_data <- aucdata@meta.data
meta_data$Cirrhosis <- as.character(meta_data$Cirrhosis)

# 计算每个样本的 AUC 均值
auc_mean_list <- split(meta_data, meta_data$Individual.Name) %>%
  lapply(function(df) {
    data.frame(
      Individual.Name = df$Individual.Name[1],  # 样本名称
      Cirrhosis = paste(unique(df$Cirrhosis), collapse = ";"),  
      t(colMeans(df[, 22:ncol(meta_data)], na.rm = TRUE))  # 计算 AUC 均值
    )
  })

# 合并列表为数据框
auc_mean_df <- do.call(rbind, auc_mean_list)

# 设置行名
rownames(auc_mean_df) <- auc_mean_df$Individual.Name
auc_mean_df <- auc_mean_df[, -1]  

write.csv(auc_mean_df,"Results of AUC.csv")

##绘制箱线图
#y_variables <- c("KEGG_APOPTOSIS", "KEGG_DRUG_METABOLISM_CYTOCHROME_P450","KEGG_ENDOCYTOSIS","KEGG_GLYCOLYSIS_GLUCONEOGENESIS")  
y_variables <- colnames(auc_mean_df)[2:ncol(auc_mean_df)]
auc_mean_df$Cirrhosis <- as.factor(auc_mean_df$Cirrhosis)

# 创建一个空列表，用于存储每个图
plot_list <- list()

# 循环遍历每个 y 变量
for (i in seq_along(y_variables)) {
  y_var <- y_variables[i]
  
  # 创建绘图对象
  p <- ggplot(auc_mean_df, aes(x = Cirrhosis,  
                             y = .data[[y_var]],  # 动态引用 y 变量
                             fill = Cirrhosis,  # 使用 "Cirrhosis" 列的值来填充箱线图的颜色
                             color = Cirrhosis  # 使用 "Cirrhosis" 列的值来设置颜色
  )) +
    geom_boxplot(outlier.shape = NA, #去除异常值
      position = position_dodge(width = 0.3),  # 用于分组箱线图，width 控制箱体的分组间隔
                 width = 0.3, fill = "white") +  # 设置箱线图的箱体宽度，设置箱体的填充颜色为白色，使箱体中间为空心
    scale_fill_manual(values = c("0" = "darkcyan", "1" = "lightcoral")) +  # 手动设置填充颜色
    scale_color_manual(values = c("0" = "darkcyan", "1" = "lightcoral")) +  # 手动设置边框颜色
    labs(x = "Cirrhosis Group", y = y_var) +  # 设置 x 轴和 y 轴的标签
    theme_bw() +  # 设置绘图主题
    stat_compare_means(method = "wilcox.test",  
                       label = "p.format",  # 显示具体 p 值
                       comparisons = list(c("0", "1")),  # 指定比较的分组
                       label.y = max(auc_mean_df[[y_var]]) # 将显著性标记放在 y 轴的上方
                       ) + 
    theme(aspect.ratio = 2/1)
  
  # 将绘图对象添加到列表中
  plot_list[[i]] <- p
}

# 设置输出 PDF 文件路径
pdf_file <- "CD8T_Lamp1+_AUCell_wilcox.pdf"
# 计算每页最多显示的图数
plots_per_row <- 4
plots_per_page <- 8  # 每页最多 8 个图（可以调整）

# 计算总页数
num_plots <- length(plot_list)
num_pages <- ceiling(num_plots / plots_per_page)

# 打开 PDF 设备
pdf(pdf_file, width = 11.69, height = 8.27)  # A4 横向大小 (单位: 英寸)

# 逐页绘制
for (page in seq_len(num_pages)) {
  start_idx <- (page - 1) * plots_per_page + 1
  end_idx <- min(start_idx + plots_per_page - 1, num_plots)
  
  # 提取当前页的图形
  plots_to_draw <- plot_list[start_idx:end_idx]
  
  # 计算行数
  num_rows <- ceiling(length(plots_to_draw) / plots_per_row)
  
  # 使用 ggarrange 组合绘图
  grid_plot <- ggarrange(plotlist = plots_to_draw, 
                         ncol = plots_per_row, nrow = num_rows)
  
  # 绘制当前页
  print(grid_plot)
}
# 关闭 PDF 设备
dev.off()
cat("PDF 已保存到: ", pdf_file, "\n")



