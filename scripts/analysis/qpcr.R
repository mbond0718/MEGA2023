## SUPPLEMENTAL FIGURE 1
## qPCR analysis of ITGB3 and KLF1

## Load libraries -----------------------------------------
library(ggplot2)
library(readxl)
library(pcr)
library(cowplot)
library(ggbreak)

## Read in data -------------------------------------------
data <- read_xlsx("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/rna/qpcr.xlsx", sheet = 1)
data <- as.data.frame(data[,1:3])
data$rep <- rep(c(rep(1, 4), rep(2, 4)), 9)

## Separate targets
itgb3 <- data[data$`Target Name`=="ITGB3",]
klf1 <- data[data$`Target Name`=="KLF1",]
gapdh <- data[data$`Target Name`=="GAPDH",]

## Analyze
itgb3 <- cbind(itgb3$`Sample Name`, itgb3$CT, gapdh$CT) |> 
  as.data.frame() |> 
  na.omit() |> 
  setNames(c("Sample", "ITGB3", "GAPDH"))
samples <- itgb3$Sample |> as.character()

data <- pcr_analyze(itgb3[,2:3], 
                    group_var = samples, 
                    reference_gene = 'GAPDH', 
                    reference_group = '0')
data <- cbind(data[,1:2], 
      c(log2(data[1,5]/data[1,5]), log2(data[2,5]/data[1,5]), log2(data[3,5]/data[1,5])), 
      c(log2(data[1,7]/data[1,7]), log2(data[2,7]/data[1,7]), log2(data[3,7]/data[1,7])), 
      c(log2(data[1,8]/data[1,8]), log2(data[2,8]/data[1,8]), log2(data[3,8]/data[1,8])))
colnames(data) <- c("group", "normalized", "lfc", "lower", "upper")

ggplot(data = data, aes(x = group, y = lfc, fill = group)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + 
  scale_fill_manual(values = rep("#53A9B2", 3)) + 
  theme_classic() + theme(legend.position = "none")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup1/itgb3.pdf")

klf1 <- cbind(klf1$`Sample Name`, klf1$CT, gapdh$CT) |> 
  as.data.frame() |> 
  na.omit() |> 
  setNames(c("Sample", "KLF1", "GAPDH"))
samples <- klf1$Sample |> as.character()

data <- pcr_analyze(klf1[,2:3], 
                    group_var = samples, 
                    reference_gene = 'GAPDH', 
                    reference_group = '0')
data <- cbind(data[,1:2], 
              c(log2(data[1,5]/data[1,5]), log2(data[2,5]/data[1,5]), log2(data[3,5]/data[1,5])), 
              c(log2(data[1,7]/data[1,7]), log2(data[2,7]/data[1,7]), log2(data[3,7]/data[1,7])), 
              c(log2(data[1,8]/data[1,8]), log2(data[2,8]/data[1,8]), log2(data[3,8]/data[1,8])))
colnames(data) <- c("group", "normalized", "lfc", "lower", "upper")

ggplot(data = data, aes(x = group, y = lfc, fill = group)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + 
  scale_fill_manual(values = rep("#E0524D", 3)) + 
  theme_classic() + theme(legend.position = "none")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup1/klf1.pdf")


