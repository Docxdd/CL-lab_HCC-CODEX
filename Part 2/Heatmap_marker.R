library(data.table)
library(readxl)
library(dplyr)
library(pheatmap)

alldata<-as.data.frame(fread("All_CN_TP_0424.csv"))
alldata<-alldata[,-1]
alldata<-alldata[,c("cellid","Allsubtypes","celltype","Class",colnames(alldata)[9:44])]##cellid,subcelltype,type,imageid,protein intensity
#colnames(alldata) <- sub("\\..*", "", colnames(alldata))
colnames(alldata) <- gsub("Nucleus Intensity", "", colnames(alldata))

clidata<-read_excel("clinic data.xlsx")
clidata<-clidata[,c("Class","LiverCirrhosis","HepatitisB")]#imageid,group,HBV

data3<- merge(alldata, clidata, by= "Class", all.x = TRUE)
data4<-data3[data3$HepatitisB==1,]
#下行代码决定是否分组
#data4<-data4[data4$LiverCirrhosis==1,]

data5<-data4[,c(3:40)]
colnames(data5)[1]<-'subtype'

data6 <- data5 %>%
  group_by(subtype, celltype) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')

subtype<-data6$subtype
data6<-data6[,-1]
rownames(data6)<-subtype

#重新行名排序
immune_cellsname<-c("APC","B cells","CD4T","CD8T","Macrophages","Treg")
immune_cells<-data6 %>%
  filter(celltype %in% immune_cellsname) %>%
  arrange(celltype)
imname<-row.names(immune_cells)
immune_cells$group<-'immune_cell'
row.names(immune_cells)<-imname


stromal_cellsname<-c("Fibroblasts","Biliary tract cells","Endothelial cells", "Lymphatic endothelial cells")
stromal_cells<-data6 %>%
  filter(celltype %in% stromal_cellsname) %>%
  arrange(celltype)
strname<-row.names(stromal_cells)
stromal_cells$group<-'stromal_cell'
rownames(stromal_cells)<-strname

tumor_cellsname<-c("Tumor")
tumor_cells<-data6 %>%
  filter(celltype %in% tumor_cellsname) %>%
  arrange(celltype)
tumname<-row.names(tumor_cells)
tumor_cells$group<-'tumor_cell'
rownames(tumor_cells)<-tumname

##############各细胞型分别展示###############
data7<-immune_cells
subtype <- rownames(data7)
celltype <- data7$celltype
group <- data7$group

# 分离行名（subtype）、类型（type）信息
subtype <- rownames(data7)
celltype <- data7$celltype

# 提取蛋白表达值
expression_data <- data7[, !(names(data7) %in% c("celltype", "group"))]
row.names(expression_data)<-row.names(data7)
expression_data<-t(expression_data)##视排版而定

#为方便排列，immune调动列顺序
expression_data <- expression_data[, c(3, 1:2, 33:45,4:32)]

#expression_data <- scale(expression_data) # 默认对列标准化
#expression_data<- t(scale(t(expression_data))) # 对行标准化

# 创建注释数据框
annotation_data <- data.frame(
  celltype = factor(celltype),
  group = factor(group)
)

rownames(annotation_data) <- subtype

#定义颜色方案
annotation_colors <- list(
  celltype = c("B cells" = "#EDB11A","APC" = "#D88851","Macrophages" = "#B31E22","Treg" = "#00A9BD", "CD4T" =  "#8197C6","CD8T" = "#334C8a"
               #"Biliary tract cells" = "#7CC6A2", "Endothelial cells" = "#6495ED", "Fibroblasts" = "#ff9200","Lymphatic endothelial cells"="#c94733"
               #"Tumor"="#501d8a"
               )
  #group = c("immune_cell" = "#DC143C", "stromal_cell" = "#191970", "tumor_cell" = "#808080") 
)


# 保存热图到文件并调整画布大小
save_pheatmap <- function(expression_data, annotation_data, annotation_colors, filename) {
  pdf(filename, width = 10, height = 10) # 调整图形尺寸，宽度和高度可以根据需要调整
  pheatmap(expression_data,
           cluster_rows = T,
           cluster_cols = FALSE,
           annotation_col = annotation_data,
           annotation_colors = annotation_colors,
           color = colorRampPalette(c("#304d9b", "white", "#ea514b"))(100), # "#8A233F"
           show_rownames = TRUE,
           show_colnames = TRUE,
           scale = "row", #'none', 'row' or 'column'
           fontsize = 8,         # 调整字体大小
           fontsize_col = 5,     # 调整列标签字体大小
           fontsize_row = 5,     # 调整行标签字体大小
           angle_col = 45,       # 旋转列标签45度
           cellheight = 5.5, #格子高度
           cellwidth = 8.5, #格子宽度
           treeheight_row = 10,
           border_color = "black",  # 边框颜色
           border = 0.001,           # 边框宽度
           legend = TRUE
           #legend_breaks = c(-2, 0, 2),
           #legend_labels = c("-2", "0", "2")
  )
  dev.off()
}

# 调用保存函数
save_pheatmap(expression_data, annotation_data, annotation_colors, "heatmap of marker in immune cells.pdf")



############下面是合并在一起展示###################

data7<-rbind(stromal_cells,tumor_cells)


# 分离行名（subtype）、类型（type）和组（group）信息
subtype <- rownames(data7)
celltype <- data7$celltype
group <- data7$group

# 提取蛋白表达值
expression_data <- data7[, !(names(data7) %in% c("celltype", "group"))]
row.names(expression_data)<-row.names(data7)
expression_data<-t(expression_data)##视排版而定

# 创建注释数据框
annotation_data <- data.frame(
  celltype = factor(celltype),
  group = factor(group)
)

rownames(annotation_data) <- subtype

#定义颜色方案
annotation_colors <- list(
  celltype = c(#"APC" = "#FFC0CB", "B cells" = "#FFF0F5", "CD4T" = "#DB7093","CD8T" = "#C71585","Macrophages" = "#DA70D6","Treg" = "#D8BFD8",
               "Biliary tract cells" = "#456990", "Endothelial cells" = "#6495ED", "Fibroblasts" = "#87CEFA","Lymphatic endothelial cells"="#B0E0E6",
               "Tumor"="#808080"),  
  group = c(#"immune_cell" = "#DC143C", 
            "stromal_cell" = "#191970", 
            "tumor_cell" = "#808080") 
)


# 保存热图到文件并调整画布大小
save_pheatmap <- function(expression_data, annotation_data, annotation_colors, filename) {
  pdf(filename, width = 10, height = 10) # 调整图形尺寸，宽度和高度可以根据需要调整
  pheatmap(expression_data,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           annotation_col = annotation_data,
           annotation_colors = annotation_colors,
           color = colorRampPalette(c("#304d9b", "white", "#ea514b"))(100),
           show_rownames = TRUE,
           show_colnames = TRUE,
           scale = "row",
           fontsize = 8,         # 调整字体大小
           fontsize_col = 3.5,     # 调整列标签字体大小
           fontsize_row = 3.5,     # 调整行标签字体大小
           angle_col = 45,       # 旋转列标签45度
           cellheight = 6, #格子高度
           cellwidth = 5, #格子宽度
           treeheight_row = 10,
           border_color = "black",  # 边框颜色
           border = 0.01,           # 边框宽度
           legend = TRUE
  )
  dev.off()
}

# 调用保存函数
save_pheatmap(expression_data, annotation_data, annotation_colors, "heatmap of marker in stromal and tumor cells.pdf")


##cell type
adata<-data5[,-1]
adata<- data5 %>%
      group_by(celltype) %>%  # 按 celltype 分组
      summarise(across(everything(), mean, na.rm = TRUE))  # 对每列（除分组列外）计算均值

rn<-adata$celltype
adata<-adata[,-1]
row.names(adata)<-rn

adata<-t(adata)

pdf("heatmap of cell marker.pdf", width = 12, height = 8) 
pheatmap(adata,
         color = colorRampPalette(c("#304d9b", "white", "#ea514b"))(100), # "#8A233F"
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "row", #'none', 'row' or 'column'
         fontsize = 8,         # 调整字体大小
         fontsize_col = 8,     # 调整列标签字体大小
         fontsize_row = 8,     # 调整行标签字体大小
         angle_col = 90,       # 旋转列标签45度
         cellheight = 10, #格子高度
         cellwidth = 8.5, #格子宽度
         border_color = "black",  # 边框颜色
         border = 0.2,           # 边框宽度
         legend = TRUE,
         cluster_rows = F,
         cluster_cols = F
         )
dev.off()
