library(data.table)
library(dplyr)
library(readxl)
library(tidyr)
library(limma)
library(future.apply)
library(ggplot2) 
library(cowplot) 
library(ggrepel) 
library(cowplot)

alldata<-as.data.frame(fread("All_CN_TP_0424.csv"))
data<-alldata[,c(3,10:45)] #Class，marker expression
colnames(data) <- gsub("\\..*", "", colnames(data))
data <- data %>%
  group_by(Class) %>%
  summarize(across(everything(), mean))

g<-read_excel("clinic data.xlsx")
g<-g[,c('Class','LiverCirrhosis','HepatitisB')]

merged_data <- left_join(data, g, by = "Class")

merged_data<-merged_data[merged_data$HepatitisB==1,]#去除非乙肝样本
merged_data<-merged_data %>%select(-HepatitisB)#去除乙肝分组列

# 提取表达数据
expression<- merged_data[, c(2:37)] 
# 提取分组信息并创建一个因子变量
groups <- as.factor(merged_data$LiverCirrhosis)

plan(multisession)
# 计算每个基因的p值和logFC
results <- future_lapply(1:ncol(expression), function(i) {
  exp_data <- expression[, i]
  
  # 创建一个数据框，包含当前基因的表达值和分组信息
  comparison_data <- data.frame(exp = exp_data, group = groups)
  
  # 计算总体的中位值
  overall_median <- median(exp_data, na.rm = TRUE)
  
  # 执行Wilcoxon秩和检验并计算logFC
  if (length(unique(comparison_data$group)) > 1) {
    test_result <- wilcox.test(exp ~ group, data = comparison_data)
    p_value <- test_result$p.value
    
    # 计算logFC
    medians <- tapply(comparison_data$exp, comparison_data$group, median)
    if (length(medians) > 1) {
      logFC <- log2(medians[2] / medians[1])  # 假设我们比较的是第二组和第一组
    } else {
      logFC <- 0  # 只有一个组时，logFC无意义，设为0
    }
  } else {
    p_value <- NA  # 只有一个分组，则返回NA
    logFC <- 0     # logFC无意义
  }
  
  # 返回logFC、p值以及总体中位值
  return(c(logFC = logFC, p_value = p_value, overall_median = overall_median))
})

# 从结果列表中提取logFC、p值以及总体中位值
logFCs <- sapply(results, function(x) x["logFC.1"])
pvalues <- sapply(results, function(x) x["p_value"])
overall_medians <- sapply(results, function(x) x["overall_median"])

# 对返回的p值列表进行Benjamini-Hochberg校正
pvalues_adjusted <- p.adjust(pvalues, method = "BH")

# 将结果存入数据框
outRst <- data.frame(
  protein = colnames(expression),
  Log2FC = logFCs,
  pValues = pvalues,
  FDR = pvalues_adjusted,
  Overall_Median = overall_medians
)

write.csv(outRst,'marker difference analysis_wilcox.csv')


#气泡图
data<-read.csv('G:/HCC CODEX/marker difference analysis_wilcox.csv')
# 创建一个新的列来标识 FDR 是否小于 0.05
data$FDR_threshold <- ifelse(data$FDR < 0.05, "< 0.05", "≥ 0.05")

# 绘制气泡图
library(ggrepel)

ggplot(data, aes(x = Log2FC, y = -log10(FDR), size = abs(Log2FC))) +
  geom_point(aes(color = FDR_threshold), alpha = 0.6) +  # 根据 FDR_threshold 来区分颜色
  scale_color_manual(values = c("< 0.05" = "red", "≥ 0.05" = "blue")) +  # 自定义颜色
  theme_minimal() +
  labs(
    x = "Log2 Fold Change", 
    y = "-Log10(FDR)", 
    #title = "Bubble Plot for Differentially Expressed Proteins",
    color = "FDR"
  ) +
  theme(legend.position = "right") +  # 调整图例位置
  geom_text_repel(
    data = subset(data, FDR < 0.05),  # 筛选出 FDR 小于 0.05 
    aes(label = protein),  # 添加名字标签
    size = 3,           # 标签大小
    box.padding = 0.4,  # 标签和点之间的间距
    point.padding = 0.3 # 点和标签之间的间距
  )

