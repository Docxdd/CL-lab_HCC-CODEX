library(data.table)
library(dplyr)
library(readxl)
library(ggpubr)
library(tidyr)
library(scales)
library(reshape2)

alldata<-as.data.frame(fread("All_CN_TP_0424.csv"))
data<-alldata[,c('Class','celltype','Allsubtypes')]

##对特定细胞亚型分析用
data <- data %>%
  filter(grepl("CD8T", celltype))

#对数据分组肝硬化和非肝硬化
g <- read_excel("clinic data.xlsx")
g<-g[,c('Class','LiverCirrhosis','HepatitisB')]

merged_data <- left_join(data, g, by = "Class")

merged_data<-merged_data[merged_data$HepatitisB==1,]#去除非乙肝样本

wide_data <- reshape2::dcast(merged_data, Class + LiverCirrhosis ~ Allsubtypes, value.var = "Allsubtypes", fun.aggregate = length)

row_sums <- rowSums(wide_data[, 3:ncol(wide_data)])
wide_data[, 3:ncol(wide_data)] <- wide_data[, 3:ncol(wide_data)] / row_sums
wide_data <- na.omit(wide_data)

##此时wide_data第1列是Class，第2列是Class的分组，后序列是特征值
#wide_data<-wide_data[,-c(4,8,9)]

tmp1  <- wide_data %>%
  pivot_longer(cols = -c(Class, LiverCirrhosis), # 指定要转换为长格式的列
               names_to = "celltype", # 新列的名称，即评分的名称
               values_to = "Frequency") 
tmp1$LiverCirrhosis<-as.factor(tmp1$LiverCirrhosis)
colnames(tmp1)[2]<-'Group'


# 计算每个 Celltype 检验
p_values <- tmp1 %>%
  group_by(celltype) %>%
  summarise(p = wilcox.test(Frequency ~ Group)$p.value)

# 保留三位小数
p_values$p_label <- sprintf("%.3f", p_values$p)

#可以使用每组 max(Frequency) 来作为标签的 y 坐标
label_df <- tmp1 %>%
  group_by(celltype) %>%
  summarise(y_pos = max(Frequency) * 1.1) %>%
  left_join(p_values, by = "celltype")


dodge_width <- 0.8

ggplot(tmp1, aes(x = celltype, y = Frequency, fill = Group)) + 
  geom_violin(position = position_dodge(width = dodge_width), alpha = 0.5, size = 0.1) +
  geom_boxplot(position = position_dodge(width = dodge_width), width = 0.1, size = 0.1, 
               alpha = 0, outlier.shape = NA, color = "black") +
  geom_text(data = label_df, aes(x = celltype, y = y_pos, label = p_label), 
            inherit.aes = FALSE, size = 3) +  # 手动添加 p 值标签
  scale_fill_manual(values = c("#304d9b", "#ea514b")) +
  labs(x = NULL, y = NULL) +
  facet_wrap(~ celltype, scales = "free", nrow = 2) +
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
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )




