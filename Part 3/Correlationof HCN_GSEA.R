library(readr)
library(readxl)
library(GSEABase)
library(GSVA)
library(data.table)
library(dplyr)


data <- read_csv("HCN proportion.csv")

###读入数据
expr<- data.table::fread("G:/HCC CODEX/肝硬化与非肝硬化肝癌/结果/CODEX_RNAseq.csv",data.table = F)
id<-expr$ID
group<-expr$Cirrhosis
expr=expr[,-c(1:2)]
row.names(expr)=id
expr <- as.matrix(expr)
expr <-t(expr)#最终行是基因，列是样本

##准备基因集##
options(timeout = 1000) 
url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.1.Hs/c5.go.bp.v2023.1.Hs.symbols.gmt"
destfile <- "c5.go.bp.v2023.1.Hs.symbols.gmt"
download.file(url, destfile, mode = "wb")

gene_sets <- getGmt(destfile)

gobp_list <- lapply(gene_sets, geneIds)
names(gobp_list) <- sapply(gene_sets, function(x) x@setName)

#选择感兴趣通路
target_names <- c(
  "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
  "GOBP_NEGATIVE_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY"
)

# 从 gobp_list 中提取这两个
target_gobp_list <- gobp_list[target_names]

# 运行 GSEA
ssgsea_params <- ssgseaParam(expr, target_gobp_list)
gsva_result <- gsva(ssgsea_params)

gsva_result <- gsva_result %>% t() %>% as.data.frame()
gsva_result$Group<-group

g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g <- g[,c('Class','PatientID_1','HepatitisB')]

gsva_df <- as.data.frame(gsva_result)
gsva_df$PatientID_1 <- rownames(gsva_df)

gsva_df <- merge(gsva_df, g, by = "PatientID_1")
rownames(gsva_df) <- gsva_df$Class
gsva_df <- gsva_df[gsva_df$HepatitisB == 1, ]
gsva_df <- gsva_df[, !(colnames(gsva_df) %in% c("PatientID_1", "Class","HepatitisB"))]

write.csv(gsva_df,'GSVA score of target GOBP.csv')

hcn_data <- data[, c("Class","HCN1")]

merged_df <- merge(hcn_data, gsva_df, by.x = "Class", by.y = "row.names")
merged_df$Group <- as.factor(merged_df$Group)

#定义函数
plot_correlation <- function(data, color_vector, var1, var2, group = NULL, text_size = 4, method = "spearman") {
  # 加载所需的包
  library(ggplot2)
  library(ggExtra)
  library(dplyr)
  
  # 检查提供的变量列是否存在于数据框中
  if (!all(c(var1, var2) %in% colnames(data))) {
    stop("指定的变量列不在数据框中")
  }
  
  # 如果提供了分组信息，并且该分组变量在数据框中
  if (!is.null(group) && group %in% colnames(data)) {
    unique_groups <- unique(data[[group]])
    if (length(color_vector) != length(unique_groups)) {
      stop("颜色向量的长度必须与分组数量相等")
    }
    
    # 创建基础散点图对象
    base_plot <- ggplot(data, aes_string(x = var1, y = var2, color = group, fill = group)) +
      geom_point(size = 2, alpha = 0.6) +
      geom_smooth(method = "lm", aes_string(group = group), se = TRUE, color = "black", alpha = 0.3) +
      scale_color_manual(values = color_vector) +
      scale_fill_manual(values = color_vector) +
      theme_classic() +
      theme(legend.position = "top",
            panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 添加全框线
            axis.line = element_blank()  # 避免与边框重复
            )
    
    # 计算每个分组的相关系数及 p 值
    cor_data <- data %>%
      group_by(.data[[group]]) %>%
      summarise(
        Corr_Value = cor(.data[[var1]], .data[[var2]], method = method),
        P_Value = cor.test(.data[[var1]], .data[[var2]], method = method, exact = FALSE)$p.value
      ) %>%
      mutate(label = paste("Group:", .data[[group]],
                           "\n", method, "'s = ", sprintf("%.4f", Corr_Value),
                           "\np-value = ", sprintf("%.4f", P_Value)))
    
    # 添加每个分组的标签
    for (i in seq_along(unique_groups)) {
      group_label <- cor_data$label[i]
      base_plot <- base_plot +
        annotate("text", x = Inf, y = Inf,
                 label = group_label,
                 hjust = 1.1, vjust = 1.2 + (i - 1) * 1.5,
                 size = text_size, color = color_vector[i])
    }
    
    print(base_plot)
    
    # 添加边际密度图
    p_densigram_plot <- ggMarginal(base_plot, type = "densigram", margins = "both",
                                   groupColour = TRUE, groupFill = TRUE, size = 5)
    
    print(p_densigram_plot)
    
  } else {
    # 无分组
    base_plot <- ggplot(data, aes_string(x = var1, y = var2)) +
      geom_point(size = 2, alpha = 0.6, color = color_vector[2]) +
      geom_smooth(method = "lm", se = TRUE, color = color_vector[1], alpha = 0.3) +
      theme_classic() +
      theme(legend.position = "none",
            panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 添加全框线
            axis.line = element_blank()  # 避免与边框重复
            )
    
    # 计算整体相关性和 p 值
    cor_value <- cor(data[[var1]], data[[var2]], method = method)
    p_value <- cor.test(data[[var1]], data[[var2]], method = method, exact = FALSE)$p.value
    
    base_plot <- base_plot +
      annotate("text", x = Inf, y = Inf,
               label = paste("\n", method, "'s = ", sprintf("%.4f", cor_value),
                             "\np-value = ", sprintf("%.4f", p_value)),
               hjust = 1.1, vjust = 1.2, size = text_size, color = "black")
    
    print(base_plot)
    
    # 添加边际密度图
    p_densigram_plot <- ggMarginal(base_plot, type = "densigram", margins = "both",
                                   size = 5, fill = color_vector[2], color = color_vector[1])
    
    print(p_densigram_plot)
  }
}

plot_correlation(data = merged_df,
                      color_vector = c("#0073C2","#EE4C1E"),
                      var1 = "HCN1",
                      var2 = "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
                      group = 'Group',
                      method = "spearman",
                      text_size = 4)

plot_correlation(data = merged_df[merged_df$Group == 1,],
                 color_vector = c("#EE4C1E","#0073C2"),
                 var1 = "HCN1",
                 var2 = "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
                 method = "spearman",
                 text_size = 4)


