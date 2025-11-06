library(data.table)
library(dplyr)
library(readxl)
library(ggpubr)
library(tidyr)
library(tibble)


datact50<-as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CT_radius50.csv"))
data<-datact50[,c('imageid','cluster_kmeans')]
data$cluster_kmeans <- paste0("CN", data$cluster_kmeans)
colnames(data)[1]<-'Class'

###对数据分组肝硬化和非肝硬化
cli <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g<-cli[,c('Class','LiverCirrhosis','HepatitisB')]

merged_data <- left_join(data, g, by = "Class")

merged_data<-merged_data[merged_data$HepatitisB==1,]#去除非乙肝样本
merged_data[, 2] <-paste0("H", merged_data[, 2])
#merged_data[, 3] <- ifelse(merged_data[, 3] == 0, "non-cirrHCC", "cirrHCC")

wide_data <- reshape2::dcast(merged_data, Class + LiverCirrhosis ~ cluster_kmeans, value.var = "cluster_kmeans", fun.aggregate = length)

row_sums <- rowSums(wide_data[, 3:ncol(wide_data)])
wide_data[, 3:ncol(wide_data)] <- wide_data[, 3:ncol(wide_data)] / row_sums
wide_data <- wide_data[ ,-2]

write.csv(wide_data,'HCN proportion.csv')

cli <- cli[cli$HepatitisB==1,]

cli_ordered <- cli %>%
  mutate(
    # 转换性别：1=Male，2=Female
    Gender = ifelse(Gender == 1, "Male", "Female"),
    
    # 转换肝硬化：0=No，1=Yes
    LiverCirrhosis = ifelse(LiverCirrhosis == 1, "Yes", "No"),
    
    # 年龄：<60 / ≥60
    Age = ifelse(Age < 60, "<60", "≥60"),
    
    # Stage：原始是 "ⅠB", "ⅢA" 等，按包含的罗马数字归类
    Stage = ifelse(grepl("Ⅰ|Ⅱ", Stage), "Ⅰ&Ⅱ", "Ⅲ&Ⅳ"),
    
    # 分化程度：用数字映射为文字 (可自定义)
    Differentiation = case_when(
      Differentiation %in% c(1,2) ~ "Ⅰ&Ⅱ",
      Differentiation  %in% c(3, 4) ~ "Ⅲ&Ⅳ",
     TRUE ~ as.character(Differentiation)
    ),
    
    # VascularTumorEmboli: 0=No，1=Yes
    VascularTumorEmboli = ifelse(VascularTumorEmboli == 1, "Yes", "No"),
    
    # ExtrahepaticMetastasis: 0=No，1=Yes
    ExtrahepaticMetastasis = ifelse(ExtrahepaticMetastasis == 1, "Yes", "No"),
    
    # 其他列保持字符
    T = as.character(T),
    TumorSize_grade = as.character(TumorSize_grade),
    AFP = as.character(AFP)
  ) %>%
  select(
    Class,
    LiverCirrhosis,
    Gender,
    Age,
    T,
    Stage,
    Differentiation,
    TumorSize_grade,
    AFP,
    VascularTumorEmboli,
    ExtrahepaticMetastasis
  )


mat <- wide_data
rownames(mat) <- mat$Class
mat <- mat[, !(colnames(mat) %in% c("Class"))]
mat <- as.matrix(mat)
mat <- t(mat)  # 行为基因，列为样本

# ==== 生成颜色映射函数 ====
make_named_colors <- function(vec, palette_func) {
  vals <- unique(na.omit(vec))  # 排除NA用于生成颜色
  cols <- palette_func(length(vals))
  color_vec <- setNames(cols, vals)
  color_vec["NA"] <- "white"    # 为 NA 添加白色
  return(color_vec)
}

# ==== 定义颜色列表，全部为命名向量 ====
anno_colors <- list(
  LiverCirrhosis = c("Yes" = "#E41A1C", "No" = "#377EB8", "NA" = "white"),
  Gender = c("Male" = "#4DAF4A", "Female" = "#984EA3", "NA" = "white"),
  Age = c("<60" = "#FF7F00", "≥60" = "#FFFF33", "NA" = "white"),
  Stage = c("Ⅰ&Ⅱ" = "#A6CEE3", "Ⅲ&Ⅳ" = "#1F78B4", "NA" = "white"),
  T = make_named_colors(cli_ordered$T, colorRampPalette(brewer.pal(8, "Set1"))),
  Differentiation = make_named_colors(cli_ordered$Differentiation, colorRampPalette(brewer.pal(8, "Set3"))),
  TumorSize_grade = make_named_colors(cli_ordered$TumorSize_grade, colorRampPalette(brewer.pal(8, "Paired"))),
  AFP = make_named_colors(cli_ordered$AFP, colorRampPalette(brewer.pal(8, "Pastel1"))),
  VascularTumorEmboli = c("Yes" = "#FB8072", "No" = "#80B1D3", "NA" = "white"),
  ExtrahepaticMetastasis = c("Yes" = "#B3DE69", "No" = "#FCCDE5", "NA" = "white")
)

# ==== 确保cli_ordered顺序和mat一致 ====
cli_ordered <- cli_ordered %>%
  filter(Class %in% colnames(mat)) %>%
  arrange(match(Class, colnames(mat)))

# ==== 创建HeatmapAnnotation ====
col_anno <- HeatmapAnnotation(
  Gender = as.character(cli_ordered$Gender),
  Age = as.character(cli_ordered$Age),
  T = as.character(cli_ordered$T),
  Stage = as.character(cli_ordered$Stage),
  Differentiation = as.character(cli_ordered$Differentiation),
  TumorSize_grade = as.character(cli_ordered$TumorSize_grade),
  AFP = as.character(cli_ordered$AFP),
  VascularTumorEmboli = as.character(cli_ordered$VascularTumorEmboli),
  ExtrahepaticMetastasis = as.character(cli_ordered$ExtrahepaticMetastasis),
  LiverCirrhosis = as.character(cli_ordered$LiverCirrhosis),
  col = anno_colors,
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_height = unit(rep(0.1, 10), "mm")  
)

# ==== 绘制热图 ====
Heatmap(
  mat,
  name = "Proportion",
  col = circlize::colorRamp2(c(0, max(mat, na.rm = TRUE)), c("white", "red")),
  top_annotation = col_anno,
  cluster_columns = T,
  cluster_rows = F,
  column_split = factor(cli_ordered$LiverCirrhosis, levels = c("No", "Yes")),  # No在左 Yes在右
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_title = NULL,
  heatmap_width = unit(10, "cm"),
  heatmap_height = unit(10, "cm"),
  row_gap = unit(0.1, "mm"),
  column_gap = unit(0.1, "mm")
)
