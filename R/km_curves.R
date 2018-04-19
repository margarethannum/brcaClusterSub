getwd()
setwd("/Users/TinyDragon/Dropbox/Margie/Documents/Coursework/Practicum/data/")

load("BRCA_survdat.Rdata")
load("tcga.iCluster.new.Rdata")
load("km.metabric.assign.Rdata")

## Load survival package
library(survival)
## Load lung data
data(lung)

## Show first 6 rows
head(lung)


## Check data
head(lung)

km.as.one <- survfit(SurvObj ~ 1, data = lung, conf.type = "log-log")
km.by.sex <- survfit(SurvObj ~ sex, data = lung, conf.type = "log-log")

plot(km.as.one)
plot(km.by.sex)

#Trying to get survival data
survTCGA <- data.frame(ID=names(res2$class), Class = res2$class)
survdat.df <- as.data.frame(survdat)
survdat.df$ID <- rownames(survdat.df)
survTCGA <- merge(survTCGA, survdat.df, by="ID")
   
survTCGA$SurvObj <- with(survTCGA, Surv(OS_OS, vital_status == 1))

km.as.one <- survfit(SurvObj ~ 1, data = survTCGA, conf.type = "log-log")
plot(km.as.one)
km.by.class <- survfit(SurvObj ~ Class, data = survTCGA, conf.type = "log-log")
plot(km.by.class)


#now try to make prettier
install.packages("survminer")
library(survminer)

plot10 <- ggsurvplot(km.by.class, data = survTCGA, 
           risk.table = F,       # show risk table.
           pval = TRUE,             # show p-value of log-rank test.
           conf.int = F,         # show confidence intervals for 
           # point estimaes of survival curves.
           xscale = c("d_y"),
           xlim = c(0,3652.5),        # present narrower X axis, but not affect
           # survival estimates.
           xlab = "Time in Years",
           legend = c("right"),
           break.time.by = 730.5,     # break X axis in time intervals by 500.
           ggtheme = theme_minimal(), # customize plot and risk table with a theme.
           risk.table.y.text.col = T, # colour risk table text annotations.
           risk.table.y.text = FALSE # show bars instead of names in text annotations
           # in legend of risk table)
)
plot10 <- plot10 + labs(
  title    = "Survival curves for 10 METABRIC subgroups in TCGA Breast Cancer Sample",   
  subtitle = "Predicted assignments to 10 METABRIC groups based on copy number and gene expression"                    
)
png("km10.4.png", width = 10, height = 4, units = 'in', res = 350)
plot10
dev.off()
#ggsave("km1.png", width = 6, height = 4, units = "in", dpi = 350)


#18 cluster KM curve
load("icluster.9.27.17.RData")

survTCGA.18 <- data.frame(ID=rownames(tcga.list3[[1]]), Class = fit.cn18$clusters)
survdat.df <- as.data.frame(survdat)
survdat.df$ID <- rownames(survdat.df)
survTCGA.18 <- merge(survTCGA.18, survdat.df, by="ID")

survTCGA.18$SurvObj <- with(survTCGA.18, Surv(OS_OS, vital_status == 1))

km.as.one.18 <- survfit(SurvObj ~ 1, data = survTCGA.18, conf.type = "log-log")
plot(km.as.one.18)
km.by.class.18 <- survfit(SurvObj ~ Class, data = survTCGA.18, conf.type = "log-log")
plot(km.by.class.18)

plot18 <- ggsurvplot(km.by.class.18, data = survTCGA.18, 
           risk.table = F,       # show risk table.
           pval = TRUE,             # show p-value of log-rank test.
           conf.int = F,         # show confidence intervals for 
           # point estimaes of survival curves.
           xscale = c("d_y"),
           xlim = c(0,3652.5),        # present narrower X axis, but not affect
           # survival estimates.
           xlab = "Time in Years",
           legend = c("right"),
           #theme = c(legend.text=element_text(size=14)),
           break.time.by = 730.5,     # break X axis in time intervals by 500.
           ggtheme = theme_minimal(), # customize plot and risk table with a theme.
           risk.table.y.text.col = T, # colour risk table text annotations.
           risk.table.y.text = FALSE # show bars instead of names in text annotations
           # in legend of risk table)
)
plot18 <- plot18 + labs(
  title    = "Survival curves for 18 novel subgroups in TCGA Breast Cancer Sample",   
  subtitle = "Novel subgroups discovered by training new classifiers on copy number, gene expression, and DNA methylation data"                    
)
plot18.2 <- plot18 + theme(legend.text=element_text(size=14))

png("km.18.4.png", width = 10, height = 5, units = 'in', res = 350)
plot18
dev.off()
#ggsave("km1.png", width = 6, height = 4, units = "in", dpi = 350)


#KM curve for regular PAM50 classification
load("BRCA_clin.supp.Rdata")
library(plyr)
clin <- rename(clin, c("Complete.TCGA.ID" = "ID"))
survTCGA.PAM <- merge(survdat.df, clin, by="ID")

survTCGA.PAM$SurvObj <- with(survTCGA.PAM, Surv(OS_OS, vital_status == 1))
km.by.class.PAM <- survfit(SurvObj ~ PAM50.mRNA, data = survTCGA.PAM, conf.type = "log-log")
plot(km.by.class.PAM)

png("km.PAM.png", width = 10, height = 12, units = 'in', res = 350)
ggsurvplot(km.by.class.PAM, data = survTCGA.PAM, 
           risk.table = F,       # show risk table.
           pval = TRUE,             # show p-value of log-rank test.
           conf.int = F,         # show confidence intervals for 
           # point estimaes of survival curves.
           xlim = c(0,8000),        # present narrower X axis, but not affect
           # survival estimates.
           break.time.by = 730.5,     # break X axis in time intervals by 500.
           ggtheme = theme_minimal(), # customize plot and risk table with a theme.
           risk.table.y.text.col = T, # colour risk table text annotations.
           risk.table.y.text = FALSE # show bars instead of names in text annotations
           # in legend of risk table)
) 

dev.off()

pamplot <- ggsurvplot(km.by.class.PAM, data = survTCGA.PAM, 
                      risk.table = F,       # show risk table.
                      pval = TRUE,             # show p-value of log-rank test.
                      conf.int = F,         # show confidence intervals for 
                      # point estimaes of survival curves.
                      xscale = c("d_y"),
                      xlim = c(0,3652.5),        # present narrower X axis, but not affect
                      # survival estimates.
                      xlab = "Time in Years",
                      legend = c("right"),
                      break.time.by = 730.5,     # break X axis in time intervals by 500.
                      ggtheme = theme_minimal(), # customize plot and risk table with a theme.
                      risk.table.y.text.col = T, # colour risk table text annotations.
                      risk.table.y.text = FALSE, 
                     # show bars instead of names in text annotations
                      # in legend of risk table)
                      legend.labs = c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-Like")
                     )
pamplot <- pamplot + labs(
  title    = "Survival curves for PAM50 groups in TCGA Breast Cancer Sample",   
  subtitle = "Intrinsic groups widely clinically used, based on ER, PR, and HER2 markers"
)

png("km.PAM4.png", width = 10, height = 4, units = 'in', res = 350)
pamplot
dev.off()
