options(timeout = 300)

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph", version = "2.0.3", repos = "http://cran.us.r-project.org")
}

if (!requireNamespace("monocle", quietly = TRUE)) {
  devtools::install_github("cole-trapnell-lab/monocle2")
}

library(igraph)
library(monocle)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)


data<- readRDS("T cell.rds")
data=subset(data,Celltype %in% c('CD8T_Effector','CD8T_Lamp1','CD8T_Ki67','CD8T_Naive','CD8T_Exhausted'))

set.seed(123)
data <- data[,sample(c(1:ncol(data)),10000)]

# 提取基因表达矩阵
expression_matrix=data@assays[["RNA"]]@counts

cell_metadata=data.frame(cell = colnames(expression_matrix),
                         NewSample.ID = data@meta.data[["Individual.Name"]],
                         Celltype = data@meta.data[["Celltype"]],
                         group = data@meta.data[["Cirrhosis"]],
                         cluster = data@active.ident)

rownames(cell_metadata) = cell_metadata$cell

gene_metadata = data.frame(gene = rownames(expression_matrix))
rownames(gene_metadata) = gene_metadata$gene
gene_metadata$gene_short_name = gene_metadata$gene

pd <- new('AnnotatedDataFrame', data = cell_metadata) 
fd <- new('AnnotatedDataFrame', data = gene_metadata)

cds <- newCellDataSet(expression_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

#释放内存
rm(cell_metadata,expression_matrix,gene_metadata,fd,pd,data)
gc()

options(future.globals.maxSize = 16000 * 1024^2)
#估计细胞的大小因子和分散性
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

##使用monocle选择的高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.05 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
topgenes <- head(disp.genes[order(disp_table$dispersion_empirical[match(disp.genes, disp_table$gene_id)], decreasing = TRUE)],100)

#设置基因排序过滤器
cds <- setOrderingFilter(cds, disp.genes)

#绘制基因排序
plot_ordering_genes(cds)

#降维处理（降到2个主成分）
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree') 

#对细胞进行排序（拟时序分析）
cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "State")
#指定方向, 重新运行 orderCells() 并指定 root
root_state <- which.max(tapply(pData(cds)$Pseudotime, pData(cds)$State, mean))
cds <- orderCells(cds, root_state = 4) #root_state可以指定
cds@phenoData@data[["group"]]<-as.factor(cds@phenoData@data[["group"]])

saveRDS(cds,'CD8T_Psedotime.rds')

#定义颜色
mycol = c("#D83D94", '#1f77b4','#2ca02c', '#d62728', '#9467bd', '#e377c2', '#17becf', '#aec7e8', '#ffbb78', '#f7b6d2')

###绘制基于拟时序的轨迹图###
plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE,raster=T)

###绘制组别细胞轨迹图###
plot_cell_trajectory(cds, color_by = "group",size = 0.2) + 
  #scale_color_manual(values = mycol) +
  theme(legend.key = element_rect(fill = NA),
        legend.title = element_text(color = 'black', size = 14, face = 2)) +
  guides(color = guide_legend(override.aes = list(size = 6)))+
  facet_wrap("~group", nrow = 5, ncol = 4)


##自定义基因##
targene<-c('CD86','HLA-DRA','CTLA4','IRF4',  #immunotherepy by CTLA4
           'AKT1','AKT2', 'CBLB','CD247','CD28',#TCR signal pathway
           'FAS','CD44','NKG7','KLRB1','CX3CR1', 'STAT1','STAT3', #Activation/Effector function
           'GZMA','GZMB','PRF1','IFNG','TNF','SERPINB9','CTSA','CTSB', #Cytotoxicity
           'PDCD1','HAVCR2','LAG3','TIGIT','TOX', #Exhaustion
           'SLC2A3','SLC16A1','PFKFB3','ALDOA'  #Glycolysis
)

##作热图和核密度图
mapgene<-topgenes
heatmap_obj<-plot_pseudotime_heatmap(cds[mapgene,],
                                     num_clusters =1,
                                     cores = 1,
                                     show_rownames = T,
                                     return_heatmap = T
)

heatmap_obj

###绘制密度图###
df_density <- data.frame(
  pseudotime = cds@phenoData@data[["Pseudotime"]], 
  group = pData(cds)$Celltype  # 细胞分组
)
# 绘制密度曲线
ggplot(df_density, aes(x = pseudotime, fill = group)) +
  geom_density(alpha = 0.5) +  # 透明度
  #scale_fill_manual(values = c("#5A9BD4","#E15759")) + 
  theme_minimal()+
  facet_wrap(vars(group), nrow = 6
             #scales = "free_y"
  )

##分组展示基因随时间变化曲线
# 提取 基因 的表达数据
gene_id <- "CTLA4"
expression_data <- data.frame(
  Pseudotime = pData(cds)$Pseudotime,
  Expression = exprs(cds)[gene_id, ],
  Group = pData(cds)$group
)

# 计算每个 group 在伪时间上的平滑曲线
ggplot(expression_data, aes(x = Pseudotime, y = Expression, color = Group)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +  # 每组一条平滑曲线
  scale_color_manual(values = c( "#1f77b4","#d62728")) + 
  theme_bw() +
  labs(title = paste("Expression of", gene_id, "over Pseudotime"),
       x = "Pseudotime",
       y = "Expression") +
  theme(legend.title = element_blank())

