library(data.table)
library(dplyr)
library(readxl)
library(ggpubr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)


g <- read_excel("G:/HCC CODEX/clinic data.xlsx")
g<-g[,c('Class','LiverCirrhosis','HepatitisB','OS','Osday')]
colnames(g)[1] <- 'imageid'
colnames(g)[4] <- 'fustat'
colnames(g)[5] <- 'futime'

datact <- as.data.frame(fread("G:/HCC CODEX/nonsigCN windows/adata_nonsig_purity50_CN_CT_radius50.csv"))

data <- left_join(datact, g, by = "imageid")

data$CN10 <- paste0("HCN", data$cluster_kmeans)

data <-data[data$HepatitisB == 1, ]

# 按 imageid 和 CN10 统计每个 CN 的数量及比例
shannon_data <- data %>%
  group_by(imageid, CN10) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(imageid) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# 计算 Shannon 指数
shannon_index <- shannon_data %>%
  group_by(imageid) %>%
  summarise(
    Shannon = -sum(Proportion * log(Proportion)),
    .groups = "drop"
  )

# 合并分组信息
shannon_index <- left_join(shannon_index, data %>% distinct(imageid, LiverCirrhosis), by = "imageid")
# 合并生存信息
shannon_surv <- left_join(shannon_index, g[, c("imageid", "futime", "fustat")], by = "imageid")

cutoff <- quantile(shannon_surv$Shannon, probs = 0.7)
#根据 cutoff 分组
shannon_surv <- shannon_surv %>%
  mutate(
    group = ifelse(Shannon > cutoff, "High", "Low"),
    group = factor(group, levels = c("Low", "High"))
  )
table(shannon_surv$group) 

#构建生存对象
surv_obj <- Surv(time = shannon_surv$futime, event = shannon_surv$fustat)
#拟合生存曲线
fit <- survfit(surv_obj ~ group, data = shannon_surv)

ggsurvplot(
  fit,
  data = shannon_surv,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = c("#377eb8", "#e41a1c"),
  legend.title = "Shannon Group",
  legend.labs = c("Low", "High"),
  surv.median.line = NULL,
  ggtheme = theme_classic(base_size = 10) + 
    theme(
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.8),
      axis.line.y = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.3, "cm"),
      axis.title = element_text(size = 16, face = "bold", color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.position = c(0.85, 0.85),
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
      legend.key.size = unit(1.2, "lines"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold")
    )
)



# 可视化 Shannon 指数按分组分布
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

ggplot(shannon_index, aes(x = as.factor(LiverCirrhosis), y = Shannon, color = as.factor(LiverCirrhosis))) +
  geom_beeswarm(alpha = 0.9, size = 1.5, cex = 3) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format", 
                     label.y = max(shannon_index$Shannon) * 1.05, 
                     size = 5) +
  scale_color_manual(
    values = c("0" = "#377eb8", "1" = "#e41a1c"),
    labels = c("0" = "Non-cirrhosis", "1" = "Cirrhosis")
  ) +
  labs(
    x = NULL,
    y = "Shannon index",
    color = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.line = element_line(size = 0.8, color = "black"),
    axis.ticks = element_line(size = 0.8, color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    axis.text = element_text(color = "black", size = 14),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 14),
    plot.margin = margin(10, 10, 10, 10),
    panel.grid = element_blank()
  ) +
  coord_cartesian(ylim = c(0, max(shannon_index$Shannon) * 1.1))


