library(data.table)
library(dplyr)
library(readxl)
library(ggpubr)
library(tidyr)

alldata<-as.data.frame(fread("G:/HCC CODEX/All_CN_TP_0424.csv"))
data<-alldata[,c('Class','Allsubtypes','celltype')]
data <- data %>%
  mutate(across(everything(), ~ gsub("_int", "", .)))

###clinic data
g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g<-g[,c('Class','HepatitisB','LiverCirrhosis')]

merged_data <- left_join(data, g, by = "Class")
merged_data<-merged_data[merged_data$HepatitisB==1,]#去除非乙肝样本
merged_data<-merged_data[,-4]#去除乙肝分组列


aa<-merged_data[,-2]

proportions <- aa %>%
  group_by(Class, LiverCirrhosis, celltype) %>%
  summarise(Count = n()) %>%
  group_by(Class) %>%
  mutate(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(Proportion = Count / Total) %>%
  select(Class, LiverCirrhosis, celltype, Proportion) %>%
  spread(celltype, Proportion, fill = 0)

longaa <- melt(proportions, id.vars = c("Class", "LiverCirrhosis"), variable.name = "Feature", value.name = "Proportion")

tumor_proportion <- longaa %>%
  filter(Feature == "Tumor") %>%
  group_by(Class) %>%
  summarise(Tumor_Proportion = sum(Proportion))# 按照 Tumor 的比例降序排列

longaa <- longaa %>%
  left_join(tumor_proportion, by = "Class") %>%
  arrange(desc(Tumor_Proportion))#根据降序排列的 Tumor 比例对 Class 进行重新排序

custom_colorsa<-c('#9b5b33','#D9BC87','#FCC910','#DECCDF','#afc2d9','#4991c1','#c6307c','#A1CC44','#89558d','#79b99d','#4D4D9F')

plota <- ggplot(longaa, aes(x = reorder(Class, -Tumor_Proportion), y = Proportion, fill = Feature)) +
  geom_bar(stat = "identity", width = 1, just = 0.5) +
  geom_vline(xintercept = seq(1.5, length(unique(longaa$Class)) - 1, by = 1), color = "grey", size = 0.01) + 
  facet_wrap(~ LiverCirrhosis, scales = "free_x", ncol = 1) +#group presentation
  theme_minimal() +
  scale_fill_manual(values = custom_colorsa) +
  labs(title = "Cells Types", x = "Sample", y = "Proportion") +
  theme(
    axis.text.x = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    text = element_text(size = 8),
    strip.text = element_text(size = 8)  # 调整分面标题的大小
  ) +
  guides(fill = guide_legend(nrow = 3))+  # 设置图例行
  theme(aspect.ratio = 1/4)  # 调整图形的纵横比，1/3 可以设置为每个图的高度
  
plota

##subtype
bb<-merged_data
#bb<- subset(bb, celltype %in% c("APC", "CD4T","CD8T","B cells","Macrophages","Treg"))#immune cell
bb<- subset(bb, celltype %in% c("Endothelial cells", "Fibroblasts","Lymphatic endothelial cells","Biliary tract cells"))#stromal cell

proportions <- bb %>%
  group_by(Class, LiverCirrhosis, celltype) %>%
  summarise(Count = n()) %>%
  group_by(Class) %>%
  mutate(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(Proportion = Count / Total) %>%
  select(Class, LiverCirrhosis, celltype, Proportion) %>%
  spread(celltype, Proportion, fill = 0)

longbb <- melt(proportions, id.vars = c("Class", "LiverCirrhosis"), variable.name = "Feature", value.name = "Proportion")

range_proportion <- longbb %>%
  filter(Feature == "Fibroblasts") %>%
  group_by(Class) %>%
  summarise(range_Proportion = sum(Proportion))# 按照某种细胞的比例降序排列

feature_order <- c("Endothelial cells", "Lymphatic endothelial cells","Biliary tract cells","Fibroblasts")
longbb$Feature <- factor(longbb$Feature, levels = feature_order)

class_order <- range_proportion %>%
  arrange(desc(range_Proportion)) %>%
  pull(Class)
longbb$Class <- factor(longbb$Class, levels = class_order)

longbb <- longbb %>%
  left_join(range_proportion, by = "Class") %>%
  arrange(desc(range_Proportion))#根据降序排列的某种细胞比例对 Class 进行重新排序


custom_colorsb<-c('#9b5b33','#D9BC87','#FCC910','#DECCDF','#afc2d9','#4991c1','#c6307c','#A1CC44','#89558d','#79b99d','#435b95')

plotb <- ggplot(longbb, aes(x = reorder(Class, -range_Proportion), y = Proportion, fill = Feature)) +
  geom_bar(stat = "identity", width = 1, just = 0.5) +
  geom_vline(xintercept = seq(1.5, length(unique(longbb$Class)) - 1, by = 1), color = "white", size = 0.1) + 
  #facet_wrap(~ LiverCirrhosis, scales = "free_x", ncol = 1) +#group presentation
  theme_minimal() +
  scale_fill_manual(values = custom_colorsb) +
  labs(title = "Cells Types", x = "Sample", y = "Proportion") +
  theme(
    axis.text.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 5),
    legend.text = element_text(size = 4),
    legend.key.size = unit(0.1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    text = element_text(size = 8),
    panel.spacing = unit(0.2, "lines"),  # 控制分面之间的间距
    strip.text = element_text(size = 8)  # 调整分面标题的大小
  ) +
  guides(fill = guide_legend(ncol = 1))+  # 将图例设置为一行
  theme(aspect.ratio = 1/4)  # 调整图形的纵横比，1/3 可以设置为每个图的高度

plotb


##Tumor subtype

cc<-merged_data
cc<- subset(cc, celltype %in% c("Tumor"))#Tumor cell
cc<- cc %>%
  mutate(celltype = case_when(
    grepl("CD107a", Allsubtypes) ~ "CD107a+Tumor",
    grepl("c-Myc", Allsubtypes) ~ "c-Myc+Tumor",
    grepl("Ki67", Allsubtypes) ~ "Ki67+Tumor",
    grepl("Twist1", Allsubtypes) ~ "Twist1+Tumor",
    grepl("Vimentin", Allsubtypes) ~ "Vimentin+Tumor",
    grepl("E-Cadherin", Allsubtypes) ~ "E-Cadherin+Tumor",
    TRUE ~ "others_Tumor" # 默认值
  ))


proportions <- cc %>%
  group_by(Class, LiverCirrhosis, celltype) %>%
  summarise(Count = n()) %>%
  group_by(Class) %>%
  mutate(Total = sum(Count)) %>%
  ungroup() %>%
  mutate(Proportion = Count / Total) %>%
  select(Class, LiverCirrhosis, celltype, Proportion) %>%
  spread(celltype, Proportion, fill = 0)

longcc<- melt(proportions, id.vars = c("Class", "LiverCirrhosis"), variable.name = "Feature", value.name = "Proportion")

range_proportion <- longcc %>%
  filter(Feature == "others_Tumor") %>%
  group_by(Class) %>%
  summarise(range_Proportion = sum(Proportion))# 按照某种细胞的比例降序排列

feature_order <- c("CD107a+Tumor", "Twist1+Tumor", "Vimentin+Tumor", "E-Cadherin+Tumor","c-Myc+Tumor", "Ki67+Tumor","others_Tumor")
longcc$Feature <- factor(longcc$Feature, levels = feature_order)

class_order <- range_proportion %>%
  arrange(desc(range_Proportion)) %>%
  pull(Class)
longcc$Class <- factor(longcc$Class, levels = class_order)

longcc <- longcc %>%
  left_join(range_proportion, by = "Class") %>%
  arrange(desc(range_Proportion))#根据降序排列的某细胞比例对 Class 进行重新排序

custom_colorsc<-c('#9b5b33','#D9BC87','#FCC910','#DECCDF','#afc2d9','#79b99d','#435b95')##颜色视情况搭配

plotc <- ggplot(longcc, aes(x = reorder(Class, -range_Proportion), y = Proportion, fill = Feature)) +
  geom_bar(stat = "identity", width = 1, just = 0.5) +
  geom_vline(xintercept = seq(1.5, length(unique(longcc$Class)) - 1, by = 1), color = "white", size = 0.1) + 
  #facet_wrap(~ LiverCirrhosis, scales = "free_x", ncol = 1) +#group presentation
  theme_minimal() +
  scale_fill_manual(values = custom_colorsc) +
  labs(title = "Cells Types", x = "Sample", y = "Proportion") +
  theme(
    axis.text.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 5),
    legend.text = element_text(size = 4),
    legend.key.size = unit(0.1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    text = element_text(size = 8),
    panel.spacing = unit(0.2, "lines"),  # 控制分面之间的间距
    strip.text = element_text(size = 8)  # 调整分面标题的大小
  ) +
  guides(fill = guide_legend(ncol = 1))+  # 将图例设置为一行
  theme(aspect.ratio = 1/4)  # 调整图形的纵横比，1/3 可以设置为每个图的高度

plotc


