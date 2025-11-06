library(data.table)
library(dplyr)
library(readxl)
library(ggpubr)
library(tidyr)
library(ggplot2)
library(factoextra)

alldata<-as.data.frame(fread("All_CN_TP_0424.csv"))
data<-alldata[,c('Class','Allsubtypes')]

###对数据分组肝硬化和非肝硬化
clidata<-read_excel("clinic data.xlsx")
g<-clidata[,c("Class","LiverCirrhosis","HepatitisB")] #imageid,group,HBV

merged_data <- left_join(data, g, by = "Class")

merged_data<-merged_data[merged_data$HepatitisB==1,]#去除非乙肝样本
merged_data<-merged_data %>%select(-HepatitisB)#去除乙肝分组列

wide_data <- reshape2::dcast(merged_data, Class + LiverCirrhosis ~ Allsubtypes, value.var = "Allsubtypes", fun.aggregate = length)

row_sums <- rowSums(wide_data[, 3:ncol(wide_data)])
wide_data[, 3:ncol(wide_data)] <- wide_data[, 3:ncol(wide_data)] / row_sums
colnames(wide_data)[2]<-'group'
wide_data$group<-as.factor(wide_data$group)
class<-wide_data$Class
row.names(wide_data)<-class
wide_data<-wide_data[,-1]
pcadata<- wide_data[, -1]#移除组别信息，只保留数值数据

# 进行PCA分析
pca_result <- prcomp(pcadata, scale. = TRUE)
# 可视化PCA结果
fviz_pca_ind(
  pca_result,
  geom.ind = "point",                    # 用点表示
  col.ind = wide_data$group,            # 按组上色
  palette = c("#6676AB", "#D97B77"),    # 自定义颜色
  addEllipses = TRUE,                   # 添加置信椭圆
  legend.title = "Group",
  title = "PCA - Biplot",
  pointshape = 16,                      # 统一形状为实心圆
  pointsize = 3,                        # 控制点大小
  alpha.ind = 0.9,                      # 点的透明度
  ellipse.alpha = 0.2                   # 置信椭圆透明度
) + 
  theme_minimal()

# 查看主成分负载
pca_result$rotation
#PCA解释方差
summary(pca_result)
