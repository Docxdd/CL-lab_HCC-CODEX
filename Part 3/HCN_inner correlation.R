library(readxl)
library(dplyr)
library(data.table)
library(readr)
library(tidyr)
library(tibble)
library(corrr)
library(ggplot2)
library(reshape2)
library(ggcorrplot)

data <- as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CT_radius50.csv"))
data<-data[,c('cellid','imageid','cluster_kmeans','phenotype')]
colnames(data)[2]<-'Class'
colnames(data)[3]<-'CN10'
data$CN10 <- paste0("HCN", data$CN10)

g <- read_excel("G:/HCC CODEX/clinic data.xlsx")

data<- left_join(data, g, by = "Class")
data<-data[data$HepatitisB==1,]#去除非乙肝样本

features <- c('HCN0','HCN1', 'HCN2','HCN3', 'HCN4', 'HCN5', 'HCN6', 'HCN7','HCN8', 'HCN9')

list_of_dataframes <- split(data, data$LiverCirrhosis)
# 移除 LiverCirrhosis列
list_of_dataframes <- lapply(list_of_dataframes, function(df) {
  df$LiverCirrhosis <- NULL
  return(df)
})


compute_proportion <- function(df, features) {
  # 每个 Class 的总 cell 数
  total_count <- df %>% group_by(Class) %>% summarise(Total = n())
  
  # 每个 Class 每个 CN 的数量
  cn_count <- df %>%
    filter(CN10 %in% features) %>%
    group_by(Class, CN10) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # 合并，计算比例
  prop_df <- left_join(cn_count, total_count, by = "Class") %>%
    mutate(Proportion = Count / Total) %>%
    select(-Count, -Total) %>%
    spread(key = CN10, value = Proportion, fill = 0) %>%
    column_to_rownames(var = "Class")
  
  return(prop_df)
}

# 分别计算
corcn0 <- compute_proportion(list_of_dataframes[["0"]], features)
corcn1 <- compute_proportion(list_of_dataframes[["1"]], features)

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


plot_correlation(data = corcn0,
                 color_vector = c("#EE4C1E","#0073C2"),
                 var1 = "HCN9",
                 var2 = "HCN1",
                 method = "spearman",
                 text_size = 4)

plot_correlation(data = corcn1,
                 color_vector = c("#EE4C1E","#0073C2"),
                 var1 = "HCN9",
                 var2 = "HCN1",
                 method = "spearman",
                 text_size = 4)



