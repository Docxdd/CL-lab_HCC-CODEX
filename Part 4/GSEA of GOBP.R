library(readr)
library(readxl)
library(dplyr)
library(GSEABase)
library(ggplot2)
library(dplyr)
library(ggpubr)

###读入数据
expr<- data.table::fread("G:/HCC CODEX/肝硬化与非肝硬化肝癌/结果/CODEX_RNAseq.csv",data.table = F)
id<-expr$ID
group<-expr$Cirrhosis
expr=expr[,-c(1:2)]
row.names(expr)=id
expr <- as.matrix(expr)
expr <-t(expr)#最终行是基因，列是样本

##准备基因集##
options(timeout = 1000) 
url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.1.Hs/c5.go.bp.v2023.1.Hs.symbols.gmt"
destfile <- "c5.go.bp.v2023.1.Hs.symbols.gmt"
download.file(url, destfile, mode = "wb")

gene_sets <- getGmt(destfile)

gobp_list <- lapply(gene_sets, geneIds)
names(gobp_list) <- sapply(gene_sets, function(x) x@setName)

pathways <- c(
  "GOBP_LIPID_GLYCOSYLATION",
  "GOBP_T_CELL_APOPTOTIC_PROCESS",
  "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
  "GOBP_POSITIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY"
)

target_gobp_list <- gobp_list[pathways]

# 运行 ssGSEA
ssgsea_params <- ssgseaParam(expr, target_gobp_list)
gsva_result <- gsva(ssgsea_params)

gsva_result <- gsva_result %>% t() %>% as.data.frame()
gsva_result$Group<-group

g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g <- g[,c('Class','PatientID_1','HepatitisB')]

gsva_df <- as.data.frame(gsva_result)
gsva_df$PatientID_1 <- rownames(gsva_df)

gsva_df <- merge(gsva_df, g, by = "PatientID_1")
rownames(gsva_df) <- gsva_df$Class
gsva_df <- gsva_df[gsva_df$HepatitisB == 1, ]
gsva_df <- gsva_df[, !(colnames(gsva_df) %in% c("PatientID_1", "Class","HepatitisB"))]

gsva_df <- gsva_df %>%
  mutate(Group = factor(Group, levels = c(0, 1), labels = c("Non-cirrhosis", "Cirrhosis")))

# 整理数据
plot_df <- gsva_df %>%
  dplyr::select(Group, all_of(pathways)) %>%
  tidyr::pivot_longer(cols = -Group, names_to = "Pathway", values_to = "Score")

# 绘图
ggplot(plot_df, aes(x = Group, y = Score, fill = Group)) +
  geom_violin( alpha = 0.6, width = 0.6)+
  geom_boxplot(outlier.shape = NA, alpha = 0, width = 0.3) +
  #geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_wrap(~ Pathway, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 1) +
  scale_fill_manual(values = c("Cirrhosis" = "#D73027", "Non-cirrhosis" = "#4575B4")) +
  theme_bw(base_size = 8) +
  theme(
    strip.text = element_text(size = 4),
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  ylab("GSVA score")


