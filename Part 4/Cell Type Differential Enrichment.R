library(data.table)
library(dplyr)
library(readxl)
library(ggpubr)
library(tidyr)
library(ggplot2)

g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g<-g[,c('Class','LiverCirrhosis','HepatitisB')]
colnames(g)[1] <- 'imageid'

datact <- as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CT_radius50.csv"))

data <- left_join(datact, g, by = "imageid")

datacp <- as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CP_radius50.csv"))
datacp <- datacp[,c(76,1:72)] #cellid及CN_CP构成

data <- left_join(data, datacp, by = "cellid")
data<-data[data$HepatitisB==1,]#去除非乙肝样本


# 先选出 cluster_kmeans == a 的数据
a <- 3

datacn <- data %>%
  filter(cluster_kmeans == a)

# 获取CN组成列的平均值，按 imageid 统计
cncomposition <- datacn %>%
  group_by(imageid, LiverCirrhosis) %>%
  summarise(across(35:106, mean, na.rm = TRUE), .groups = "drop")  # 注意替换实际的列名范围

cncomposition_filtered <- cncomposition %>%
  select(imageid, LiverCirrhosis, matches("CD4T|CD8T"))

cellcomposition <- cncomposition_filtered %>%
  rowwise() %>%
  mutate(across(matches("CD4T|CD8T"), ~ .x / sum(c_across(matches("CD4T|CD8T"))))) %>%
  ungroup()

cellcomposition <- na.omit(cellcomposition) #删除缺失值行

write.csv(cellcomposition,'HCN3 T of all T composition.csv')

# wilcox 检验
cellcomposition <- cellcomposition[,c('imageid','LiverCirrhosis','CD4T_FOXP3+_KI67+')]
cellcomposition <- cellcomposition %>% filter(.data[[names(.)[3]]] != 0)

colnames(cellcomposition)[2] <- 'group'
colnames(cellcomposition)[3] <- 'proportion'

pval <- wilcox.test(proportion ~ group, data = cellcomposition)$p.value
p_label <- sprintf("%.3f", pval)

##作箱线图
y_max <- max(cellcomposition$proportion, na.rm = TRUE) * 1.05

# 作图
ggviolin(cellcomposition,
         x = "group",
         y = "proportion",
         color = "group",
         fill = "group",  # 小提琴图填充色
         ylab = "proportion",
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


