library(readr)
library(tidyr)
library(ggplot2)

dfPurity <- read.csv("dfPurity.csv", row.names = 1)

# totalnumber 是 dfPurity 的第一行
totalnumber <- dfPurity[1, ]
# 进行行除法（每行除以 totalnumber）
dfPurity_percent <- sweep(dfPurity, 2, as.numeric(totalnumber), "/")

dfPurity_percent$Threshold <- rownames(dfPurity_percent)
df_long <- pivot_longer(dfPurity_percent, cols = -Threshold, names_to = "Group", values_to = "Proportion")

# 将 Threshold 转为有序的因子（若你想按 90, 80, ... 排列）
df_long$Threshold <- factor(df_long$Threshold, levels = unique(dfPurity_percent$Threshold))

ggplot(df_long, aes(x = Threshold, y = Proportion, color = Group, group = Group)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Trend of Each Group under Different Purity Thresholds",
       x = "Purity Threshold",
       y = "Proportion")
