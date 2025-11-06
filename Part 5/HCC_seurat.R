library(dplyr)
library(tibble)
library(readxl)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)

#import data
Tcell <- readRDS("I:/01 HCC T cell.rds")
cli <- read_excel("I:/cli.xlsx")

#merge cirrhosis group
Tcell@meta.data <- Tcell@meta.data %>%
  left_join(cli %>% dplyr::select(Patient, Cirrhosis), by = c("Individual.Name" = "Patient")) %>%
  column_to_rownames(var = "barcode")
Tcell@meta.data$Cirrhosis <- as.factor(Tcell@meta.data$Cirrhosis)

#DimPlot
DimPlot(Tcell, reduction = "umap",label = TRUE,label.size = 3,
        group.by ='seurat_clusters',
        raster = T
        )

#DotPlot
DotPlot(Tcell, 
        features = c('CD3D','CD3E','CD3G',#Tcell
                      'CD8A', 'CD8B',#CD8T
                      'CD4',#CD4T
                      'CCR7',#naive T
                      'FOXP3',#Treg
                      'LAMP1',
                      'TNF','IFNG', #Effetor
                      'MKI67', #Proliferation
                      'HIF1A', #Hypoxia
                      'PDCD1','CTLA4','LAG3' #Exhausted
                       ),
        cluster.idents = T,
        dot.scale = 7,
        cols = c("white","#EC3232")
        #group.by = "Cirrhosis"
        )+
    theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1), # x-axis text
        axis.text.y = element_text(size = 6), # y-axis text
        axis.title = element_text(size = 14)) # axis titles

# Identify the celltype
Tcell$Celltype <- NA

Tcell$Celltype[Tcell@active.ident %in% c(0,1)] <- "CD3_Naive"
Tcell$Celltype[Tcell@active.ident %in% c(4,8,9,10,14,26,27)] <- "CD4T"
Tcell$Celltype[Tcell@active.ident %in% c(2,11,16,17,20,25,31)] <- "CD8T_Effector"
Tcell$Celltype[Tcell@active.ident %in% c(3,5,18,21,24)] <- "CD8T_Lamp1"
Tcell$Celltype[Tcell@active.ident %in% c(6)] <- "CD8T_Ki67"
Tcell$Celltype[Tcell@active.ident %in% c(7,23)] <- "CD8T_Naive"
Tcell$Celltype[Tcell@active.ident %in% c(12,13,15,19,22,29)] <- "CD8T_Exhausted"
Tcell$Celltype[Tcell@active.ident %in% c(28,30,32)] <- "T cell_Others"

DimPlot(Tcell, reduction = "umap",
        label = TRUE,
        label.size = 3,
        group.by ='Celltype',
        repel = TRUE,
        raster=T,
        cols = c("CD3_Naive" = "#D83D94", "CD4T" = "#3F77B4", "CD8T_Effector" = "#EE811D",
                 "CD8T_Lamp1" = "#52968E", "CD8T_Ki67" = "#C12E2C", "CD8T_Naive" = "#906ABC",
                 "CD8T_Exhausted" = "#BBBC2D","T cell_Others" = "#7F7F7F")
        )+
  theme(text = element_text(size = 8))

saveRDS(Tcell,"T cell.rds")

#DEG of CD8T_Lamp1
targetcell<-Tcell[,Tcell@meta.data[["Celltype"]]%in%c('CD8T_Lamp1')]
Idents(targetcell) <- targetcell$Cirrhosis
markers <- FindMarkers(targetcell, ident.1 = "0", ident.2 = "1", min.pct = 0.15,test.use = "wilcox")
write.csv(markers,'markers of CD8T_Lamp1_non_cirrHCC vs cirrHCC.csv')



