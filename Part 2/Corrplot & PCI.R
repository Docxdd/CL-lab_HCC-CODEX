library(data.table)
library(tidyverse) 
library(tidyr)  
library(ggpubr)   
library(ggthemes)   
library(pheatmap)   
library(corrplot)    
library(psych)   
library(reshape)  
library(cowplot)  
library(gplots)   
library(stringr)

cells_type <- as.data.frame(fread("All_CN_TP_0424.csv"))

cells_type <- cells_type %>% dplyr::select(Class, celltype)
colnames(cells_type)[2] <- 'celltype'

g <- read_excel("G:/HCC CODEX/clinic data.xlsx")

data0 <- cells_type[cells_type$Class %in% g[g$HepatitisB == '1'& g$LiverCirrhosis == '0',]$Class,]
data1 <- cells_type[cells_type$Class %in% g[g$HepatitisB == '1'& g$LiverCirrhosis == '1',]$Class,]

dat <- data1
dat <- as.data.frame(with(dat, table(Class, celltype)))

cells_subtype_Freq <- spread(dat, celltype, Freq)
dat_subtype_freq <- cells_subtype_Freq[,-1]
row.names(dat_subtype_freq) <- cells_subtype_Freq$Class

dat_subtype_percent <- dat_subtype_freq / rowSums(dat_subtype_freq) * 100
row.names(dat_subtype_percent)<- cells_subtype_Freq$Class

MM<-corr.test(dat_subtype_percent, method="pearson")

temp <-  corrplot(MM$r,method = "square",order = 'hclust',
                  p.mat = MM$p, sig.level = 0.05, insig = "label_sig")

#Spatially Single-cell Protein Landscape Reveals the Distinct Tumor Microenvironment Architecture and Communications in Hepatocellular Carcinoma提供代码得到PCI文件spatialanalysis_output_celltype_class.csv和spatialanalysis_output_subtype_class.csv
data <- read_csv("G:/HCC CODEX/spatialanalysis_output_celltype_class.csv")
colnames(data)[2] <- 'celltype'

datan <- data[data$Class %in% g[g$HepatitisB == '1'& g$LiverCirrhosis == '0',]$Class,]
datac <- data[data$Class %in% g[g$HepatitisB == '1'& g$LiverCirrhosis == '1',]$Class,]

functionX <- function(dat){
  dat <- dplyr::rename(dat, subtype=celltype)  # 重命名 "Allsubtypes" 列为 "subtype"
  dat0 <- as.matrix(dat[c(-1,-2)])   # 只保留数值列
  print(sum(dat0==Inf))              # 检查是否有 Inf
  dat0[dat0==Inf] <- 100000          # 把 Inf 替换成 100000
  dat1 <- dat0/1000                  # 缩小数据规模（除以 1000）
  dat2 <- round(dat1,2)              # 四舍五入到小数点后 2 位
  dat3 <- as.data.frame(dat2)        # 转换回数据框
  dat4 <- cbind(dat[2], dat3)        # 重新合并 subtype 列
  dat5 <- aggregate(. ~ subtype, dat4, sum)  # 按 subtype 汇总
  return(dat5)
}


data <- functionX(datac) #注意改组

data2 <- aggregate(. ~ subtype, data, sum)
row.names(data2) <- data2$subtype
data2 <- data2[,-1]

data2 <- data2[order(rownames(data2)),order(colnames(data2))]

dat <- as.matrix(data2)
dat_filter <- dat
#dat_filter <- dat_filter[rowSums(dat_filter)>100,colSums(dat_filter)>100]
rowSums(dat_filter)

#按likelihood ratios 
dat_new<-matrix(nrow = length(rownames(dat_filter)), ncol=length(colnames(dat_filter)), dimnames = list(rownames(dat_filter),colnames(dat_filter)))

dat_sum <- sum(dat_filter)
for (i in 1:length(rownames(dat_filter))) {
  for (j in 1:length(colnames(dat_filter))) 
  {
    a0 <- dat_sum/sum(dat_filter[i,])/sum(dat_filter[,j])
    a<- dat_filter[i,j]*a0
    b<- log(a+0.0001)
    dat_new[i,j]<-b
  }
}

###########
mat <- t(round(dat_new,2))

bk <- c(seq(-6,-0.1,by=0.01),seq(0,6,by=0.01))
pheatmap(mat,
         cellwidth = 15,
         cellheight = 15,
         border=FALSE,
         color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),colorRampPalette(color=c("white","red"))(length(bk)/2)),
         scale = "none",
         cluster_rows = F,
         show_rownames = T,
         cluster_cols = F,
         legend_breaks=seq(-6,6,2),
         breaks=bk)

#### plot ####
temp1 = MM$r[rownames(temp$corr),colnames(temp$corr)]   
melt_corr = reshape2::melt(temp1)
melt_corr$Var1 <- factor(melt_corr$Var1)
melt_corr$Var2 <- factor(melt_corr$Var2)

mat1 = mat[rownames(temp$corr),colnames(temp$corr)]
melted.dat_notTransposed <- reshape2::melt(data.matrix(mat1))
melted.dat_notTransposed$Var1 <- factor(melted.dat_notTransposed$Var1)
melted.dat_notTransposed$Var2 <- factor(melted.dat_notTransposed$Var2)


p1 <-ggplot(melt_corr, aes(y = Var1,
                           x = Var2, fill=value)) +        
  geom_tile() +         
  scale_fill_gradient2(low = "#6676AB", high = "#D97B77", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation")+
  labs(x="",y="")+
  theme(legend.position = c(0.6,-0.15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

p1

melted.dat_notTransposed$color = ifelse(melted.dat_notTransposed$value >0, "#455655","#a1d8b1")
p2 <- ggplot(melted.dat_notTransposed,aes(y = Var1,x = Var2))+
  geom_point(aes(size = melted.dat_notTransposed$value), color = melted.dat_notTransposed$color)+   
   #scale_color_gradient("#455655" = "#455655",
                        #"#a1d8b1" = "#a1d8b1",
                       #name = "Likelihood\nRatio")+
  scale_size(range=c(-4,4),name = "Likelihood\nRatio")+
  scale_color_manual(values = c("#455655" = "#455655", "#a1d8b1" = "#a1d8b1"), 
                     name = "Sign", labels = c("#455655" = "> 0", "#a1d8b1" = "<= 0")) + 
  labs(x="",y="")+
  theme(legend.position = c(0.9,-0.15),
        legend.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"), 
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
p2

p3 = ggdraw() +
  draw_plot(p1, 0, 0, 1, 1) +
  draw_plot(p2, 0, 0, 1, 1)+ggtitle('Overlay')

p3
