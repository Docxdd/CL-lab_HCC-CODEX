library(GSVA)
library(tidyverse)
library(data.table)
library(ggsci)
library(tidyr)
library(ggpubr)
library(readxl)
library(tibble)
library(magrittr)
library(reshape2)
library(ggplot2)

genesets_immune_function <- readRDS("genesets_immune_function.rds")
markers <- genesets_immune_function

###读入数据
g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g <- g[,c('PatientID_1','HepatitisB')]
colnames(g)[1] <- 'ID'
g <- g[g$HepatitisB == 1, ]

expr <- data.table::fread("CODEX_RNAseq.csv",data.table = F)
expr <- expr[expr$ID %in% g$ID, ]

id<-expr$ID
group<-expr$Cirrhosis
expr=expr[,-c(1:2)]
row.names(expr)=id
expr <- as.matrix(expr)
expr <-t(expr)#最终行是基因，列是样本

##ssgsea分析
ssgseaPar <- ssgseaParam(expr, markers) 
gsva_data <- gsva(ssgseaPar)
gsva_data<- gsva_data %>% t() %>% as.data.frame()
gsva_data$Group<-group

##作箱线图
gsva_data$Sample<-id

tmp1  <- gsva_data %>%
  pivot_longer(cols = -c(Sample, Group), # 指定要转换为长格式的列，这里排除了Sample和group列
               names_to = "GF", # 新列的名称，即评分的名称
               values_to = "Score") 
tmp1$Group<-as.factor(tmp1$Group)

# 计算每个 Group 的数量
group_counts <- tmp1 %>%
  group_by(Group) %>%
  summarise(Count = n()) %>%
  mutate(Group_Label = paste0(Group, " (n = ", Count, ")"))  # 将数量添加到 Group 名称中

# 将新的 Group 标签合并回原始数据
tmp1 <- tmp1 %>%
  left_join(group_counts, by = "Group")

# 计算每个 Celltype 检验
p_values <- tmp1 %>%
  group_by(GF) %>%
  summarise(p = wilcox.test(Score ~ Group)$p.value)

# 保留三位小数
p_values$p_label <- sprintf("%.3f", p_values$p)

#可以使用每组 max(Frequency) 来作为标签的 y 坐标
label_df <- tmp1 %>%
  group_by(GF) %>%
  summarise(y_pos = max(Score) * 1.1) %>%
  left_join(p_values, by = "GF")


dodge_width <- 0.9

#pdf("Immune function boxplot.pdf",width = 12,height = 8)
ggplot(tmp1, aes(x = GF, y = Score, fill = Group)) + 
  geom_violin(aes(fill = Group_Label),position = position_dodge(width = dodge_width), alpha = 0.8, size = 0.1) +
  geom_boxplot(aes(fill = Group_Label),position = position_dodge(width = dodge_width), width = 0.2, size = 0.1, 
               alpha = 0, outlier.shape = NA, color = "black") +
  geom_text(data = label_df, aes(x = GF, y = y_pos, label = p_label), 
            inherit.aes = FALSE, size = 3) +  # 手动添加 p 值标签
  scale_fill_manual(values = c("#304d9b", "#ea514b")) +
  labs(x = NULL, y = NULL) +
  facet_wrap(~ GF, scales = "free", nrow = 1) +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.13))) +
  theme(
    plot.margin = unit(c(2, 3, 2, 3), 'cm'),
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
