library(data.table)
library(dplyr)
library(readxl)

alldata<-as.data.frame(fread("G:/HCC CODEX/All_CN_TP_0424.csv"))
data<-alldata[,c('Class','Allsubtypes','celltype')]
data <- data %>%
  mutate(across(everything(), ~ gsub("_int", "", .)))

###clinic data
g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g<-g[,c('Class','HepatitisB','LiverCirrhosis')]

data <- left_join(data, g, by = "Class")
data<-data[data$HepatitisB==1,]#去除非乙肝样本
data<-data[,-4]#去除乙肝分组列

cell_type_proportion <- data %>%
  group_by(.[[3]]) %>%                  # 第 3 列celltype分组
  summarise(count = n()) %>%           # 统计每种细胞型的数量
  mutate(proportion = count / sum(count)) # 计算比例
colnames(cell_type_proportion)[1]<-'celltype'
cell_type_proportion<-cell_type_proportion[,-2]#删除count列

cell_type_proportion <- cell_type_proportion %>%
  mutate(group = case_when(
    celltype %in% c("APC", "B cells", "CD4T", "CD8T", "Macrophages", "Treg") ~ "Immune cell",
    celltype == "Tumor" ~ "Tumor cell",
    TRUE ~ "Stromal cell" # 其他类型归为 Stromal cell
  ))

# 计算 group所占percent
cell_type_proportion <- cell_type_proportion %>%
  group_by(group) %>%
  mutate(percent1 = sum(proportion, na.rm = TRUE)) %>%
  ungroup()
#调整顺序
cell_type_proportion <- cell_type_proportion %>%
  select(group, percent1, everything())

colnames(cell_type_proportion)[4]<-'percent2'

write.csv(cell_type_proportion,'sankydata.csv')


#
library(ggplot2)
devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)


# 数据转换为 ggplot2 可用格式
data_long <- cell_type_proportion %>%
  make_long(group, celltype,value = percent2)

sankey_colors<-c('Immune cell'= "#BBC258",
                 'Stromal cell'= "#9A9AF8",
                 'Tumor cell'="#85A8D3",
                 'APC'="#AAD7C8", 
                 'B cells'="#78D3AC", 
                 'Biliary tract cells'= '#Bed798',
                 'CD4T'= "#ADC6AD", 
                 'CD8T'="#C9D111", 
                 "Endothelial cells"= "#CFE2E6",
                 "Fibroblasts"= "#BAD8FB",
                 "Lymphatic endothelial cells"= "#87CEEB",
                 "Macrophages"= "#BDCA9F", 
                 "Treg"= "#C2DCBF",
                 "Tumor"="#85A8D3")

ggplot(data_long, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                      fill=factor(node),
                      label=node,
                      value = value
)
) +
  geom_sankey(flow.alpha = 0.5, node.color = "gray30", aes(fill = factor(node))) +
  geom_sankey_text(size=2,color="#000000",width=1)+
  scale_fill_manual(values = sankey_colors) +
  scale_x_discrete(position="top")+
  labs(x=NULL)+
  ggtitle("Snankey")+
  theme_sankey(base_size=10)+
  theme(
    legend.position="none",
    plot.title=element_text(hjust=0.5),
    axis.text.x.top=element_text(color="#000000")
  )
