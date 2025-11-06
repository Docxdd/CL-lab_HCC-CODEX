library(data.table)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(CCA)     
library(psych)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(igraph)


cells <- as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CT_radius50.csv"))

g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g<-g[,c('Class','LiverCirrhosis','HepatitisB')]
colnames(g)[1] <- 'imageid'
colnames(g)[2] <- 'groups'

cells <- left_join(cells, g, by = "imageid")
cells <- cells[cells$HepatitisB == 1, ]

to_group <- cells %>%  dplyr::select(imageid, groups) %>%  distinct() %>%  deframe()

group0 <- names(to_group[to_group == 0])
group1 <- names(to_group[to_group == 1])

#选择要分析的邻域和功能细胞子集
cns <- c(1,3,4,5,8,9)
subsets <- c(#"APC", 
             "CD4T", 
             "CD8T", 
             #"Macrophages", 
             #"B cells", 
             "Treg"
             #"Endothelial cells",
             #"Fibroblasts",
             #"Lymphatic endothelial cells",
             #"Biliary tract cells"
             )

#计算 log(1e-3 + neighborhood-specific cell type frequency)
nsctf <- cells %>%
  group_by(imageid, cluster_kmeans) %>%
  summarise(across(all_of(subsets), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(all_of(subsets), ~log(.x + 1e-3))) %>%
  unite("imageid_neigh", imageid, cluster_kmeans, remove = FALSE) %>%
  column_to_rownames("imageid_neigh")


# 注意 nsctf 是 log(1e-3 + 平均频率) 的数据
# 需要根据 CN 取出相应行
stats_group <- list()
for (cn_i in cns) {
  for (cn_j in cns) {
    if (cn_i < cn_j) {
      
      # 提取患者在这两个邻域的功能子集数据并合并
      x_df <- nsctf %>%
        filter(cluster_kmeans == cn_i & imageid %in% group0) %>%
        dplyr::select(all_of(c(subsets, "imageid")))
      
      y_df <- nsctf %>% 
        filter(cluster_kmeans == cn_j & imageid %in% group0) %>%
        dplyr::select(all_of(subsets), imageid)
      
      combined <- inner_join(x_df, y_df, by = "imageid", suffix = c("_x", "_y")) %>% 
        drop_na()
      
      x <- as.matrix(combined %>% dplyr::select(ends_with("_x")))
      y <- as.matrix(combined %>% dplyr::select(ends_with("_y")))
      
      if (nrow(x) > 2) {
        cca_res <- cc(x, y)
        
        # 第一对 canonical component
        x_scores <- x %*% cca_res$xcoef[, 1]
        y_scores <- y %*% cca_res$ycoef[, 1]
        
        true_corr <- cor(x_scores, y_scores)
        perm_corrs <- numeric(1000)
        
        for (i in 1:1000) {
          idx <- sample(nrow(x))
          x_perm <- x[idx, ]
          cca_perm <- cc(x_perm, y)
          x_perm_scores <- x_perm %*% cca_perm$xcoef[, 1]
          y_perm_scores <- y %*% cca_perm$ycoef[, 1]
          perm_corrs[i] <- cor(x_perm_scores, y_perm_scores)
        }
        
        stats_group[[paste0(cn_i, "_", cn_j)]] <- list(true = true_corr, perm = perm_corrs)
        cat(sprintf("Finished CN pair (%d,%d): true corr = %.4f\n", cn_i, cn_j, true_corr))
      }
    }
  }
}



# 构建边数据框
edges_df <- tibble::enframe(stats_group, name = "pair", value = "val") %>%
  separate(pair, into = c("cn1", "cn2"), sep = "_", convert = TRUE) %>%
  mutate(
    obs = map_dbl(val, ~ .x[[1]]),
    perms = map(val, ~ .x[[2]]),
    p = map2_dbl(obs, perms, ~ mean(.x > .y))
  ) %>%
  filter(p > 0.9) %>%
  dplyr::select(from = cn1, to = cn2, weight = obs)

edges_df$from <- as.character(edges_df$from)
edges_df$to <- as.character(edges_df$to)


# 构建图
g <- graph_from_data_frame(edges_df, directed = FALSE)

# 确保节点名字是字符
V(g)$name <- as.character(V(g)$name)

layout_df <- layout_with_fr(g) %>% 
  as.data.frame()
colnames(layout_df) <- c("x", "y")
layout_df$name <- V(g)$name

ggraph(g, layout = layout_df) +
  geom_edge_link(aes(edge_alpha = 5 * weight^5, edge_width = 5 * weight^5), color = "black") +
  geom_node_point(aes(color = name), size = 10) +
  geom_node_text(aes(label = name), size = 6, vjust = 0.5, hjust = 0.5) +
  scale_edge_alpha(range = c(0.1, 1)) +
  scale_edge_width(range = c(0.1, 1)) +
  theme_void() +
  ggtitle("Canonical Correlation Network (Group 1)")



