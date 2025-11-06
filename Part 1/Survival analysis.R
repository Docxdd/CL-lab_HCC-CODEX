library(readxl)
library(survival)
library(survminer)
library(autoReg)
library(dplyr)


data <- read_excel("zsdata.xlsx")

colnames(data) <- c("patientID","Gender","Age","OStime","OS","RFStime","RFS","Cirrhosis","HBV","TumorSize","Pathology","MVI","PT","INR","AFP","Hb","WBC","PLT","BCLC")

# 绘制生存曲线
survfit <- survfit(Surv(RFStime,RFS == 1) ~ Cirrhosis,  # 创建生存对象 
               data = data)

median_time <- summary(survfit)$table[, "median"]

ggsurvplot(
  survfit, 
  data = data,
  palette = c("blue", "red"),  # 设置生存曲线的颜色
  #risk.table = TRUE,
  #table.pval = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  #xlim = c(0, 2000),
  #break.x.by = 500,
  ggtheme = theme_bw(),
  conf.int.style = "step",
  linetype = c(1, 1)  # 修改线条类型
)


##UniCox and MultuCox
data <- data %>%
  mutate(
    # 将 Maxcm 分组分类
    TumorSize = case_when(
      TumorSize < 5 ~ 1,
      TumorSize >= 5 ~ 2,
      TRUE ~ NA_real_
    ),
    TumorSize = factor(TumorSize, levels = 1:2, labels = c("1", "2")),
    
    # pathology归类
    Pathology = case_when(
      Pathology %in% c("1", "2") ~ 1,
      Pathology %in% c("3", "4") ~ 2,
      TRUE ~ NA_real_
    ),
    Pathology = factor(Pathology, levels = c(1, 2), labels = c("1", "2")),
  
    # 其他分类变量转因子
    Cirrhosis = factor(Cirrhosis, levels = c("0", "1")),
    MVI = factor(MVI),
    BCLC = factor(BCLC)
  )

fit <- coxph(Surv(RFStime,RFS) ~ Gender + Age + Cirrhosis + Pathology + MVI + BCLC + AFP + TumorSize + PT + INR +  WBC + PLT, 
                data = data)
result <- autoReg(fit, uni = TRUE, 
                  threshold = 1,  
                  multi = TRUE, 
                  final = F)

myft(result)

write.csv(result, 'ZS_hospital_RFS_Cox_result.csv', row.names = FALSE, fileEncoding = "GBK")
