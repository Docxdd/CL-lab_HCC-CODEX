library(readr)
library(readxl)
library(WGCNA)
library(limma)
library(tidyverse)

options(timeout = 10000)

rs <- read_csv("G:/HCC CODEX/CODEX_RNAseq.csv")
rs <- rs[ ,-2]
rs = as.matrix(rs)
rownames(rs) = rs[,1]

exp = rs[,2:ncol(rs)]
exp <- t(exp)
dimnames=list(rownames(exp),colnames(exp)) #【构建一个列表】
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames) #【从上一步的列表中构建一个新的矩阵】
data=avereps(data)
data2=as.data.frame(data) #【将新的矩阵构建为数据框】

datExpr0 <- t(data2)
gsg = goodSamplesGenes(datExpr0, verbose = 3)#【检查缺缺失值】
gsg$allOK
#【如果有缺失值就在下方删除】
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###过滤表达量【将表达量过低的基因删除】#####
meanFPKM=0.1  #【设置过滤值】
n=nrow(datExpr0)
datExpr0 <- rbind(datExpr0, mean = colMeans(datExpr0))#【读取表达量值】
datExpr0 <- datExpr0[1:(nrow(datExpr0)-1), datExpr0[nrow(datExpr0), ] > meanFPKM]#【将表达量值与过滤值进行比较】
dimnames=list(rownames(datExpr0), colnames(datExpr0))#【保存符合条件的值】
data=matrix(as.numeric(as.matrix(datExpr0)), nrow=nrow(datExpr0), dimnames=dimnames)
datExpr0=avereps(data)

###导出过滤后的基因表达矩阵
#filtered=t(datExpr0)#【将矩阵进行行列转换，转换为基因名在上，样品名在左】
filtered=datExpr0
filtered=data.frame(rownames(filtered),filtered)
names(filtered)[1]="sample"
#【将矩阵导出】
write.csv(filtered, file="rt_filter.csv")

###样品分群聚类
sampleTree = hclust(dist(datExpr0), method = "average")#【将每个相似的样品进行聚类】
pdf(file = "1_samplenet.pdf", width = 10, height = 10)#【绘制聚类图】
par(cex = 0.5)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
dev.off()

hlevel=30000#【根据聚类图，查看是否有离群值，并设置删除高度】

###【删除离群值】
pdf(file = "1_samplenet_cut.pdf", width = 10, height = 10)
par(cex = 0.5)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
abline(h = hlevel, col = "red")#【h】是设置删除高度
dev.off()
###【保留非离群样本】
clust = cutreeStatic(sampleTree, cutHeight = hlevel, minSize = 10)
table(clust)
keepSamples = (clust==1)#【提取非离群数据】
datExpr0 = datExpr0[keepSamples, ] #【构建非离群矩阵】

###载入临床性状数据
cli <- read_excel("G:/HCC CODEX/肝硬化与非肝硬化肝癌/结果/clinic data.xlsx")
cli <- cli[cli$HepatitisB == 1 , -1]
cli <- cli[!duplicated(cli$PatientID_1), ]
rn <- cli$PatientID_1
rownames(cli) <- rn

common_samples <- intersect(rownames(datExpr0), rownames(cli))
datExpr0 <- datExpr0[common_samples, ]
cli <- cli[common_samples, ]
rn <- cli$PatientID_1

cli <- cli[ ,c("Gender","Age","Differentiation","TumorSize_grade", "AFP","LiverCirrhosis")]
rownames(cli) <- rn

cli$TumorSize_grade <- ifelse(cli$TumorSize_grade == "≥5cm", 2,
                                 ifelse(cli$TumorSize_grade == "＜5cm", 1, NA))

cli$AFP <- ifelse(cli$AFP == "Negative", 0,
                  ifelse(cli$AFP == "Positive", 1, NA))

datTraits <- cli

collectGarbage()


###【绘制非离群样本聚类图】
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="2_sample_map.pdf",width=10,height=10)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",cex.lab = 1,
                    cex.axis = 1, cex.main = 1)
dev.off()

enableWGCNAThreads() #多线程工作
powers = c(1:13) #【无标度拓扑拟合指数】
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="3_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
###【无标度拓扑拟合指数散点图】
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #【若没有数值到达0.9则可以降低标准，但最低不低于0.8】
###平均连通性散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

sft #【查看最佳软阈】
softPower =sft$powerEstimate #【最佳软阈】
adjacency = adjacency(datExpr0, power = softPower)

###【进行邻接矩阵转换】
###【转换为TOM矩阵】
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


###【基因构建网络】
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="4_gene.pdf",width=12,height=9)#【绘制图形】
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

###【构建模块识别】
minModuleSize = 100 #【设置模块基因数目，每个模块最大基因数量为5000】
#【构建基因模块】
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 1, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)#【查看模块标识】
dynamicColors = labels2colors(dynamicMods) #【设置模块颜色】
table(dynamicColors) #【查看模块颜色】
#【绘制模块图形】
pdf(file="5_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

###【将相似的模块进行网络构建】
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)#【提取相似模块】
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="6_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

#【删除离群模块】
MEDissThres = 0.25  #【设置剪切高度】

#【绘制模块的网络图】
pdf(file="6_module_cut.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

###【绘制模块与性状数据热图】
nGenes = ncol(datExpr0)#【提取模块基因】
nSamples = nrow(datExpr0)#【提取模块样本】
moduleTraitCor = cor(MEs, datTraits, use = "p")#【计算模块相关性】
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)#【计算模块P值】
pdf(file="8_Module.pdf",width=10,height=8)#【绘制性状数据热图】
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

###【导出基因所在的模块】
moduleColors=dynamicColors
probes = colnames(datExpr0)#【读取基因名】
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
geneOrder =order(geneInfo0$moduleColor)#【读取模块】
geneInfo = geneInfo0[geneOrder, ]
#【保存所有模块基因】
write.table(geneInfo, file = "all_genes.txt",sep="\t",row.names=F,quote=F)

###【输出每个模块的基因】
#【使用循环对每个模块的基因进行导出】
for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

###【计算MM和GS】
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

###【批量输出性状和模块散点图】
picDir="module_trait"#【创建文件夹】
dir.create(picDir)
for (trait in traitNames){#【对每个性状和模块进行循环绘制图片】
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      pdfFile=paste("9_", trait, "_", module,".pdf",sep="")
      outPdf=paste(picDir,pdfFile,sep="\\")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}


### 提取目标基因
module = "turquoise"
inModule = (moduleColors == module)
modProbes = probes[inModule]  # turquoise模块的基因名

# 提取对应的 MM 和 GS
MM_col = paste0("MM", module)
GS_col = "GS.LiverCirrhosis"
gene_MM = geneModuleMembership[inModule, MM_col]
gene_GS = geneTraitSignificance[inModule, GS_col]

# 合并为表格
hubGeneTable = data.frame(Gene = modProbes,
                          MM = gene_MM,
                          GS = gene_GS)

# 筛选高MM + 高GS
hubGeneTable_filtered = hubGeneTable[abs(hubGeneTable$MM) > 0.75 & abs(hubGeneTable$GS) > 0.3, ]

# 查看前几行
head(hubGeneTable_filtered)

# 保存结果
write.csv(hubGeneTable, file = "Cirrhosis_turquoise_HubGenes.csv")


####Hubgene 富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggprism)
library(tidyverse)


kegg <- enrichKEGG(gene = bitr(hubGeneTable_filtered$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)
#dotplot(kegg)

# 使用setReadable将KEGG富集结果中的基因ID转换为基因符号
kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
KEGG <- as.data.frame(kegg_readable)

# 筛选 TOP 基因通路
# 对每个分类选择top，根据p.adjust排序，相同p.adjust值取富集到基因最多的通路
use_pathway <- top_n(KEGG,5, -p.adjust) %>%
      group_by(p.adjust) %>%
      top_n(1, wt = Count) %>%
      mutate(ONTOLOGY ='KEGG')%>%
  ungroup() %>%
  dplyr::arrange(qvalue) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

# 绘制富集通路图
plot_enrichment <- function() {
  p <- use_pathway %>%
    ggplot(aes(x = -log10(qvalue), y = reorder(Description, -log10(qvalue)), fill = Description)) +
    geom_col(width = 0.6, alpha = 0.8) +
    geom_text(
      aes(label = Description),
      x = 0.05, hjust = 0, size = 6
    ) +
    geom_text(
      aes(label = geneID),
      x = 0.1, hjust = 0, vjust = 2.5, size = 4, 
      show.legend = FALSE
    ) +
    geom_point(
      aes(x = -0.6, size = Count),
      shape = 21
    ) +
    geom_text(
      aes(x = -0.6, label = Count),
      size = 3
    ) +
    scale_size_continuous(name = 'Count', range = c(5, 8)) +
    labs(x = "-log10(qvalue)", y = NULL, fill = "Pathway") +
    scale_fill_manual(values = setNames(RColorBrewer::brewer.pal(n = nrow(use_pathway), name = "Set3"), use_pathway$Description)) +
    #scale_color_manual(values = setNames(RColorBrewer::brewer.pal(n = nrow(use_pathway), name = "Set3"), use_pathway$Description)) +
    scale_y_discrete(expand = expansion(mult = c(0.1, 0.1))) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),    # 去掉y轴文字
      axis.ticks.y = element_blank(),   # 去掉y轴刻度
      axis.title.y = element_blank(),   # 去掉y轴标题
      axis.title.x = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.position = "none"          # 不显示图例
    )
  return(p)
}

# 生成图形
enrichment_plot <- plot_enrichment()
enrichment_plot
