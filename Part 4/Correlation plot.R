library(readxl)

a <- read_csv("HCN3_Treg_ki67 of T.csv")
a<-a[,c("imageid", "group",'proportion')]

b <- read_csv("HCN1_CD8T_CD107a of T.csv")
b<-b[,c("imageid", "group",'proportion')]

data <- merge(a, b, by = c("imageid", "group"))
colnames(data) <- c("imageid", "group", "HCN3_Treg_Ki67", "HCN1_CD8T_CD107a")

rname<-data$imageid
data<-data[,-1]
rownames(data)=rname

data$group <- as.factor(data$group)
data$group <- ifelse(data$group == 0, "non-cirrHCC", "cirrHCC")  
#此时data第2和3列分别是值，第1列是分组，行为样本
#删除包含0的行
data <- data %>%
  filter(across(everything(), ~ . != 0))

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
      theme_minimal() +
      theme(legend.position = "top")
    
    # 计算每个分组的相关系数及 p 值
    cor_data <- data %>%
      group_by(.data[[group]]) %>%
      summarise(
        Corr_Value = cor(.data[[var1]], .data[[var2]], method = method),
        P_Value = cor.test(.data[[var1]], .data[[var2]], method = method, exact = FALSE)$p.value
      ) %>%
      mutate(label = paste("Group:", .data[[group]],
                           "\n", method, "'s = ", format(Corr_Value, digits = 4),
                           "\np-value = ", format(P_Value, digits = 4)))
    
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
      theme_minimal() +
      theme(legend.position = "none")
    
    # 计算整体相关性和 p 值
    cor_value <- cor(data[[var1]], data[[var2]], method = method)
    p_value <- cor.test(data[[var1]], data[[var2]], method = method, exact = FALSE)$p.value
    
    base_plot <- base_plot +
      annotate("text", x = Inf, y = Inf,
               label = paste("\n", method, "'s = ", format(cor_value, digits = 4),
                             "\np-value = ", format(p_value, digits = 4)),
               hjust = 1.1, vjust = 1.2, size = text_size, color = "black")
    
    print(base_plot)
    
    # 添加边际密度图
    p_densigram_plot <- ggMarginal(base_plot, type = "densigram", margins = "both",
                                   size = 5, fill = color_vector[2], color = color_vector[1])
    print(p_densigram_plot)
  }
}


#有分组
plot_correlation(data=data,
                 color_vector=c("#EE4C1E","#0073C2"),
                 var1="HCN3_Treg_Ki67",
                 var2="HCN1_CD8T_CD107a",
                 group="group",
                 method = "spearman",
                 text_size=4)

#没有分组信息。如果没有分组，则可生成单纯的两变量相关性绘图，见下图
plot_correlation(data=data[data$group == 'non-cirrHCC',],
                 color_vector=c("#EE4C1E","#0073C2"),
                 var1="HCN3_Treg_Ki67",
                 var2="HCN1_CD8T_CD107a",
                 method = "spearman",
                 text_size=4)

plot_correlation(data=data[data$group == 'cirrHCC',],
                 color_vector=c("#EE4C1E","#0073C2"),
                 var1="HCN3_Treg_Ki67",
                 var2="HCN1_CD8T_CD107a",
                 method = "spearman",
                 text_size=4)

plot_correlation(data=data,
                 color_vector=c("#EE4C1E","#0073C2"),
                 var1="HCN3_Treg_Ki67",
                 var2="HCN1_CD8T_CD107a",
                 method = "spearman",
                 text_size=4)
