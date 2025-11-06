library(Seurat)
library(ggplot2)


Tcell <- readRDS("T cell.rds")

CD4T <- Tcell[,Tcell@meta.data[["Celltype"]]%in%c('CD4T')]

CD4T<- NormalizeData(CD4T, normalization.method = "LogNormalize", scale.factor = 10000)
CD4T <- FindVariableFeatures(CD4T, selection.method = "vst", nfeatures = 2000) 
CD4T<- ScaleData(CD4T) 
CD4T<- RunPCA(CD4T, features = VariableFeatures(object =CD4T))

DimPlot(CD4T, reduction = "pca") + NoLegend()
ElbowPlot(CD4T)

CD4T<- FindNeighbors(CD4T, dims = 1:10)
CD4T<- FindClusters(CD4T, resolution =0.5)

CD4T<- RunUMAP(CD4T, dims = 1:10)
DimPlot(CD4T, reduction = "umap",label = TRUE,group.by ='seurat_clusters')


DotPlot(CD4T, 
        features = c('CD3D','CD3E','CD3G',#Tcell
                     'CD4',#CD4T
                     #'CCR7',#naive 
                     'FOXP3','CTLA4'#Treg
                     #'IFNG','CXCL13',#Th1
                     #'IL23R', #Th17
                     #'CXCR6', 'TCF7','CD6', #Memory
                     #'CCR3','CCR4' #Th2
                     
           
        ),
        cluster.idents = T,
        dot.scale = 7,
        cols = c("white","#EC3232")
        #group.by = "Cirrhosis"
)+
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1), # x-axis text
        axis.text.y = element_text(size = 6), # y-axis text
        axis.title = element_text(size = 14)) # axis titles

#注释
CD4T$Celltype <- NA

CD4T$Celltype[CD4T@active.ident %in% c(1,2,3,4,7,9,10,11,12,13)] <- "CD4T_Th"
CD4T$Celltype[CD4T@active.ident %in% c(0,5,6,8)] <- "Treg"

CD4T_Th<- CD4T[,CD4T@meta.data[["Celltype"]]%in%c('CD4T_Th')]
Treg<- CD4T[,CD4T@meta.data[["Celltype"]]%in%c('Treg')]

saveRDS(CD4T,"CD4T.rds")
saveRDS(CD4T_Th,"CD4T_Th.rds")
saveRDS(Treg,"Treg.rds")


DimPlot(CD4T, reduction = "umap",
        label = TRUE,
        label.size = 3,
        group.by ='Celltype',
        #repel = TRUE,
        raster=T,
        cols = c("CD4T_Th" = "#C12E2C", "Treg" = "#3F77B4")
)+
  theme(text = element_text(size = 12))


FeaturePlot(CD4T, 
            features = "SPP1", 
            cols = c("grey", "red"),
            raster=T
)


#DEG of Treg
Idents(Treg) <- Treg$Cirrhosis
markers <- FindMarkers(Treg, ident.1 = "0", ident.2 = "1", min.pct = 0.1,test.use = "wilcox")
write.csv(markers,'markers of Treg_noncirrHCC vs cirrHCC.csv')

data <- markers
data$gene <- rownames(data)

# 创建一个新的列来标识 FDR 是否小于 0.05
data$FDR_threshold <- ifelse(data$p_val_adj < 0.05, "< 0.05", "≥ 0.05")

# 绘制气泡图
library(ggrepel)

ggplot(data, aes(x = avg_log2FC, y = -log10(p_val_adj), size = abs(avg_log2FC))) +
  geom_point(aes(color = FDR_threshold), alpha = 0.6) +  # 根据 FDR_threshold 来区分颜色
  scale_color_manual(values = c("< 0.05" = "#ea514b", "≥ 0.05" = "#304d9b")) +  # 自定义颜色
  theme_minimal() +
  labs(
    x = "log2FC", 
    y = "-Log10FDR", 
    #title = "Bubble Plot for Differentially Expressed Proteins",
    color = "p_val_adj"
  ) +
  theme(legend.position = "right") +  # 调整图例位置
  geom_text_repel(
    data = subset(data, gene == c("SPP1")), 
    aes(label = gene),
    size = 3,
    box.padding = 0.3,  # 标签和点之间的间距
    point.padding = 0.3, # 点和标签之间的间距
    segment.color = "black",            # 设置连接线颜色
    segment.size = 1,                  # 设置连接线粗细
    max.overlaps = Inf,                  # 避免被滤掉
    nudge_x = 0.5,
    nudge_y = 1
  ) + 
  theme(aspect.ratio = 1.2/1)


##Ranked scatter plot
# 创建一个新的列 `rank`，表示基因按 avg_log2FC 排序后的索引
data$rank <- rank(-data$avg_log2FC, ties.method = "random")  

top <- data[order(-data$avg_log2FC), ][1:5, ]  # 上调基因
bottom <- data[order(data$avg_log2FC), ][1:5, ]  # 下调基因
highlighted_genes <- rbind(top, bottom)

selected_genes <- c("SPP1")
selected_genes <- data[rownames(data) %in% selected_genes, ]

highlighted_genes <- rbind(top, bottom,selected_genes)
# 绘制散点图
ggplot(data, aes(x = rank, y = avg_log2FC)) +
  geom_point(alpha = 0.8, color = "#ffbb78",size = 1) +  # 绘制散点
  geom_text_repel(data = highlighted_genes, aes(label = rownames(highlighted_genes)), 
                  size = 3, color = "black", max.overlaps = Inf
                  ) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # 添加 y=0 参考线
  labs(x = "Gene Rank", y = "avg_log2FC", title = "Ranked Scatter Plot of DEGs") +
  theme_bw()+
  theme(aspect.ratio = 1.2/1)


