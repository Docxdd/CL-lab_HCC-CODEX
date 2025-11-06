library(Seurat)
library(nichenetr)
library(tidyverse)
library(ggplot2)


Tcell <- readRDS("T cell.rds")
Tcell$Cirrhosis <- as.factor(Tcell$Cirrhosis)
CD8T_Lamp1 <- Tcell[,Tcell@meta.data[["Celltype"]]%in%c('CD8T_Lamp1')]
rm(Tcell)

Treg <- readRDS("Treg.rds")
CD4T_Th <- readRDS("CD4T_Th.rds")
DC <- readRDS("01 HCC DC.rds")
Macro <- readRDS("01 HCC Macro.rds")

seu_list <- list(CD4T_Th = CD4T_Th, CD8T_Lamp1 = CD8T_Lamp1, Treg = Treg, 
                DC = DC, Macro = Macro)
rm(CD4T_Th, CD8T_Lamp1, Treg, DC, Macro)

set.seed(123)
# 遍历 Seurat 对象列表，随机抽取特定比例的细胞
for (name in names(seu_list)) {
  selected_cells <- sample(Cells(seu_list[[name]]), 
                           size = length(Cells(seu_list[[name]])) * 0.5)
  seu_list[[name]] <- subset(seu_list[[name]], cells = selected_cells)
}

seurat <- merge(x = seu_list[[1]], y = seu_list[-1])
rm(seu_list)

Idents(seurat) <- seurat$Celltype

#获取 NicheNet 需要的参考数据
weighted_networks <- readRDS("weighted_networks.rds")
lr_network <- readRDS("lr_network.rds")
ligand_target_matrix <- readRDS("ligand_target_matrix.rds")

lr_network <- lr_network %>% distinct(from, to)

#Define a set of potential ligands for both the sender-agnostic and sender-focused approach
receiver = "CD8T_Lamp1"
expressed_genes_receiver <- get_expressed_genes(receiver, seurat, pct = 0.1)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- c("Treg")
# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat, 0.1)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)

#组间差异表达基因
#DEG <- DEG <- read_csv("markers of CD8T_Lamp1_non_cirrHCC vs cirrHCC.csv")
#geneset_oi <-DEG[DEG$p_val_adj < 0.05 & abs(DEG$avg_log2FC)>= 0.25,]$gene

#组间差异表达基因
condition_oi <-  "1"
condition_reference <- "0"

seurat_obj_receiver <- subset(seurat, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "Cirrhosis",
                                  min.pct = 0.1) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.1) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define the background genes
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
length(geneset_oi)

#预测配体活性
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
#查看预测的配体活性，显示前n个活性最高的配体。
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(50, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

best_upstream_ligands <- ligand_activities %>% top_n(50, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

#visualize the ligand activity measure (AUPR) of these top-ranked ligands
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()
          )

#推断高排名配体的靶基因和受体,哪些基因和受体可能是 最优先配体（top-ranked ligands） 的下游靶点?
#Active target gene inference
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

#we identify which receptors have the highest interaction potential with the top-ranked ligands
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                    y_name = "Ligands", x_name = "Receptors",  
                    color = "mediumvioletred",
                    legend_title = "Prior interaction potential"
                    ))

#发送细胞（Sender cells）是否实际表达这些优先配体?
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6, angle = 90),  # 图例文本大小
        legend.title = element_text(size = 6, face = "bold"),  # 图例标题大小 + 加粗
        legend.key.height = unit(0.4, "cm"),  # 控制图例色条高度
        legend.key.width = unit(0.6, "cm"),  # 控制图例色条宽度
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6,  angle = 90, hjust = 0)
        )

p_ligand_target

# Receptor plot
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor

best_upstream_ligands_all %in% rownames(seurat) %>% table()

# Dotplot of sender-focused approach
p_dotplot <- DotPlot(subset(seurat, Celltype %in% sender_celltypes),
                     features = rev(best_upstream_ligands), 
                     cols = "RdYlBu",
                     dot.scale = 4) + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot


celltype_order <- levels(Idents(seurat)) 
# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = seurat,
  condition_colname = "Cirrhosis",
  condition_oi = "1",
  condition_reference = "0",
  celltype_col = "Celltype",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", 
                                        mid_color = "white",
                                        mid = 0,#median(vis_ligand_lfc), 
                                        high_color = "red",
                                        legend_title = "LFC")

p_lfc



figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none",
                        axis.ticks = element_blank(),
                        axis.title.y = element_blank(),
                        axis.title.x = element_text(size = 6),
                        axis.text.y = element_text(size = 6),
                        axis.text.x = element_text(size = 6,  angle = 90, hjust = 0)),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 6),
                    axis.text.y = element_text(size = 6),
                    axis.text.x = element_text(size = 6,  angle = 90, hjust = 0)
                    ) +
    ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none",
                axis.title.y = element_blank(),
                axis.title.x = element_text(size = 6),
                axis.text.y = element_text(size = 6),
                axis.text.x = element_text(size = 6,  angle = 90, hjust = 0)
                ),
  p_ligand_receptor + theme(legend.position = "none",
                          axis.title.y = element_blank(),
                          axis.title.x = element_text(size = 6),
                          axis.text.y = element_text(size = 6),
                          axis.text.x = element_text(size = 6,  angle = 90, hjust = 0)
                          ),
  align = "hv",
  nrow = 1,
  rel_widths = c(3, 4, 4, 20)
  )

figures_without_legend

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_receptor)),
  nrow = 1,
  align = "h", rel_widths = c(1, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

all_data <- list(
  ligand_activities = ligand_activities,
  vis_ligand_aupr = vis_ligand_aupr,
  vis_ligand_target = vis_ligand_target,
  vis_ligand_receptor_network = vis_ligand_receptor_network,
  seurat = seurat,
  sender_celltypes = sender_celltypes,
  best_upstream_ligands = best_upstream_ligands,
  celltype_order = celltype_order,
  vis_ligand_lfc = vis_ligand_lfc,
  combined_plot
)

# 保存整个列表为一个 .rds 文件
saveRDS(all_data, file = "nichenet_results_Treg.rds")




