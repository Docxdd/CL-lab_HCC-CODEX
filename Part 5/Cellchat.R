library(CellChat)
library(patchwork)
library(Seurat)
library(ggpubr)

#1.读入数据并初步处理
Tcell <- readRDS("T cell.rds")
#Tcell$Cirrhosis <- as.factor(Tcell$Cirrhosis)
CD8T_Lamp1 <- Tcell[,Tcell@meta.data[["Celltype"]]%in%c('CD8T_Lamp1')]
rm(Tcell)

Treg <- readRDS("Treg.rds")
CD4T_Th <- readRDS("CD4T_Th.rds")
DC <- readRDS("J:/张泽民肝癌数据/01 HCC DC.rds")
Macro <- readRDS("J:/张泽民肝癌数据/01 HCC Macro.rds")

seurat <- merge(CD8T_Lamp1, y = c(Treg, CD4T_Th,DC, Macro))
Idents(seurat) <- seurat$Celltype
rm(CD8T_Lamp1,Treg, CD4T_Th,DC,Macro)

noncirrHCC_seurat <- seurat[,seurat@meta.data[["Cirrhosis"]] %in% c('0')]
cirrHCC_seurat <- seurat[,seurat@meta.data[["Cirrhosis"]] %in% c('1')]
seurat<-noncirrHCC_seurat

#2.创建 CellChat 对象
cellchat <- createCellChat(object = seurat, group.by = "Celltype")
rm(seurat)

#3.导入配受体数据库（人或者鼠）
CellChatDB <- CellChatDB.human
#CellChatDB <- CellChatDB.mouse
#直接使用CellChatDB全库进行细胞通讯分析：
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
#选择数据库中特定子集进行细胞通讯分析:
#CellChatDB.use <- subsetDB(CellChatDB, search = 'Secreted Signaling')#可选择Secreted Signaling、ECM-Receptor或Cell-Cell Contact
#将数据库添加到CellChat对象中(DB)
cellchat@DB <- CellChatDB.use

#4.cellchat
#将信号基因的表达矩阵子集化，节省计算成本
cellchat <- subsetData(cellchat) 

library(future)
set.seed(123)
#future::plan('multisession', workers = 2)
options(future.globals.maxSize = 16000 * 1024^2)# 增加允许的最大对象大小到1000 MiB（可以根据需要调整）

#鉴定与每个细胞亚群相关的过表达信号基因：
##基于表达该基因的细胞比例、差异倍数和p值判定。
cellchat <- identifyOverExpressedGenes(cellchat,
                                       only.pos = TRUE, #是否仅返回positive markers
                                       thresh.pc = 0, #细胞比例阈值
                                       thresh.fc = 0, #差异倍数
                                       thresh.p = 0.05) #P-Value

#识别过表达基因配体-受体互作：
cellchat <- identifyOverExpressedInteractions(cellchat)

#将基因表达数据映射到PPI网络(可跳过)：
#cellchat <- projectData(cellchat, PPI.human) #返回结果:cellchat@data.project

#计算细胞通讯概率：
cellchat <- computeCommunProb(cellchat, raw.use = TRUE) #返回结果:cellchat@options$parameter
#默认使用原始表达数据(cellchat@data.Signaling)，若想使用上一步PPI矫正数据,设置raw.use = FALSE
cellchat <- filterCommunication(cellchat, min.cells = 10)

#信号通路水平的细胞通讯表
cellchat <- computeCommunProbPathway(cellchat)

#计算细胞对间通讯的数量和概率强度
cellchat <- aggregateNet(cellchat)

# 计算中心分数
cellchat <- netAnalysis_computeCentrality(cellchat,
                                          slot.name = "netP")

saveRDS(cellchat,"CD8T_Lamp1_noncirrHCC_cellchat.rds")

#信号通路可视化
netVisual_bubble(cellchat,
                 sources.use = c("CD4T_Th","DC","Macrophage","Treg"),
                 targets.use = c("CD8T_Lamp1"),
                 #signaling = c("MHC-I","MHC-II"), 
                 remove.isolate = FALSE,
                 sort.by.source = T,
                 sort.by.target = F,
                 group = c("0","1"),
                 font.size =6)+ 
  theme(axis.text.x = element_text(angle = 90),
        aspect.ratio = 1/4  #高宽比
        ) +
  coord_flip()

#合并分析#
noncirrHCC_cellchat_r <- readRDS("CD8T_Lamp1_noncirrHCC_cellchat.rds")
cirrHCC_cellchat_r <- readRDS("CD8T_Lamp1_cirrHCC_cellchat.rds")

object.list <- list(noncirrHCC = noncirrHCC_cellchat_r,cirrHCC = cirrHCC_cellchat_r)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(cellchat,"merged_cellchat.rds")

###比较不同条件下的细胞间通信###
#比较细胞通信的总交互数量和交互强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) #总交互数量
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") # 交互强度
gg1 + gg2

#计算不同细胞类型间的交互变化,红色在第2组中增加的信号传递, 蓝色在第2组中减少的信号传递
#热图展示
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

#比较每个信号通路的整体信息流
#这部分通过比较每个信号通路的信息流来识别保守和情境特定的信号通路，信息流定义为所有细胞组对之间的通讯概率的总和。可以绘制堆叠柱状图，以展示每个信号通路在不同条件下的信息流。
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2


###识别上调和下调的信号配体-受体对###
netVisual_bubble(cellchat, 
                 sources.use = c('Treg','CD4T_Th','DC','Macrophage'), 
                 targets.use = c('CD8T_Lamp1'), 
                 comparison = c(1, 2),
                 signaling = c('MHC-I'),
                 angle.x = 90, 
                 grid.on = T, 
                 line.on = T) 





