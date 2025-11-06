library(deldir)
library(ggplot2)
library(sf)
library(tibble)
library(data.table)
library(sp)
library(readxl)
library(dplyr)
library(readr)

alldata <- as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CT_radius50.csv"))
data <- alldata[,c('cellid','imageid','XMin','YMin','phenotype','cluster_kmeans')]

data<-data[,-1]
colnames(data)[1]<-'Class'
colnames(data)[2]<-'X'
colnames(data)[3]<-'Y'
colnames(data)[5]<-'neighborhood10'
data[, 5] <- paste("HCN", data[, 5], sep="")
#data <- data %>%
  #mutate(across(everything(), ~ gsub("_int", "", .)))##如果用亚型，把里边一些无用的符号去掉

###对数据分组肝硬化和非肝硬化
g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g<-g[,c('Class','LiverCirrhosis','HepatitisB')]

merged_data <- left_join(data, g, by = "Class")
merged_data$X<-as.numeric(merged_data$X)
merged_data$Y<-as.numeric(merged_data$Y)

merged_data<-merged_data[merged_data$HepatitisB==1,]#去除非乙肝样本
data <- subset(merged_data, select = -HepatitisB)#去除乙肝分组列

table(data[data$LiverCirrhosis == 0,]$Class)


# 将 neighborhood10 映射到颜色
colors_map <- c(HCN0 = "#d2ebc8", HCN1 = "#304d9b",  HCN2 = "#d55640", HCN3 = "#6cb8d2",
                HCN4 = "#aecde1", HCN5 = "#9b5b33", HCN6 = "#ee934e", 
                HCN7 = "#b383b9", HCN8 = "#d1352b", HCN9 = "#3c77af")
neighborhood_color <- data.frame(
  neighborhood10 = names(colors_map),
  neighborhood_color = colors_map,
  stringsAsFactors = FALSE
)

##
spot<-data[data$Class=='TMA4_4_reg025',]
spot<-spot[,c('X','Y','neighborhood10')]

# 获取 Voronoi 图
tesselation <- deldir(spot$X, spot$Y)
tiles <- tile.list(tesselation)

spot <- left_join(spot, neighborhood_color, by = "neighborhood10")
neighborhood_colors <-spot$neighborhood_color


layout(matrix(c(1, 1), nrow = 1, byrow = TRUE), widths = c(1, 1))
par(mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))  # oma增加了额外的外边界空间
# 准备绘带有标记的图
pdf("1TMA4_4_reg025_symbol.pdf", width = 10, height = 10) 
plot.new()
plot.window(xlim = range(spot$X), ylim = range(spot$Y))
par(bty = 'n', xaxt = 'n', yaxt = 'n', ann = FALSE)  

# 遍历每个 Voronoi 区域
for (i in seq_along(tiles)) {
  tile <- tiles[[i]]
  
  # 确保多边形闭合
  if (!all(tile$x[1] == tile$x[length(tile$x)] && tile$y[1] == tile$y[length(tile$y)])) {
    tile$x <- c(tile$x, tile$x[1])
    tile$y <- c(tile$y, tile$y[1])
  }
  
  # 创建 sf 对象以计算面积
  sf_poly <- st_polygon(list(cbind(tile$x, tile$y)))
  poly_area <- st_area(sf_poly)
  
  # 仅对面积小于阈值的多边形进行显示和标记
  if (as.numeric(poly_area) < 60000) {  # 面积小于阈值的多边形进行显示
    nearest_point_index <- which.min(sqrt((spot$X - mean(tile$x))^2 + (spot$Y - mean(tile$y))^2))
    poly_color <- neighborhood_colors[nearest_point_index]
    poly_label <- spot$neighborhood10[nearest_point_index]
    
    # 绘制多边形
    polygon(tile$x, tile$y, col = poly_color, 
            border = 'white',#多边形边界颜色,可设置NA
            lwd = 0.01 #多边形边界宽度
            )
    
    # 计算多边形质心坐标
    centroid <- c(mean(tile$x, na.rm = TRUE), mean(tile$y, na.rm = TRUE))
    
    # 检查质心坐标是否有 NA
    if (!any(is.na(centroid))) {
      # 在多边形内标记 `neighborhood10`
      text(centroid[1], centroid[2], labels = poly_label, cex = 0.05, col = "black", font = 2)  # 调整 cex 使文本可见
    }
  }
}

dev.off()



layout(matrix(c(1, 1), nrow = 1, byrow = TRUE), widths = c(1, 1))
par(mar = c(1, 1, 1, 1), oma = c(0, 0, 0, 0))  # oma增加了额外的外边界空间
# 准备绘无标记的图
pdf("1TMA4_4_reg025.pdf", width = 10, height = 10) 
plot.new()
plot.window(xlim = range(spot$X), ylim = range(spot$Y))
par(bty = 'n', xaxt = 'n', yaxt = 'n', ann = FALSE)  

# 遍历每个 Voronoi 区域
for (i in seq_along(tiles)) {
  tile <- tiles[[i]]
  
  # 确保多边形闭合
  if (!all(tile$x[1] == tile$x[length(tile$x)] && tile$y[1] == tile$y[length(tile$y)])) {
    tile$x <- c(tile$x, tile$x[1])
    tile$y <- c(tile$y, tile$y[1])
  }
  
  # 创建 sf 对象以计算面积
  sf_poly <- st_polygon(list(cbind(tile$x, tile$y)))
  poly_area <- st_area(sf_poly)
  
  # 仅对面积小于阈值的多边形进行显示和标记
  if (as.numeric(poly_area) < 60000) {  # 面积小于阈值的多边形进行显示
    nearest_point_index <- which.min(sqrt((spot$X - mean(tile$x))^2 + (spot$Y - mean(tile$y))^2))
    poly_color <- neighborhood_colors[nearest_point_index]
    poly_label <- spot$neighborhood10[nearest_point_index]
    
    # 绘制多边形
    polygon(tile$x, tile$y, col = poly_color,
            border = NA,#多边形边界颜色,可设置NA
            lwd = 0.01 #多边形边界宽度
            )
    
  }
}

dev.off()



#lengend

pdf("lengend.pdf", width = 5, height = 3) 
plot.new()

legend("bottom",  # 使用计算出的位置
       legend = neighborhood_color$neighborhood10,  # 显示每个 neighborhood10 的标签
       fill = neighborhood_color$neighborhood_color,  # 对应的颜色
       #title = "CN",  # 图例标题
       cex = 0.5,  # 控制图例文本的大小
       pt.cex = 0.6,  # 控制图例点的大小
       #x.intersp = 0.2,  # 控制图例间距
       text.width = 0.4,  # 控制图例文本宽度
       #title.adj = 0.3,  # 控制图例标题对齐
       bty = 'n', # 没有边框
       horiz=F
         ) 
dev.off()

ggplot(spot, aes(x = X, y = Y, color = neighborhood10)) +
  geom_point(alpha = 0.8, size = 0.6) +
  scale_color_manual(values = colors_map) +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "right",
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

