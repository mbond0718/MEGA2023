## Compare feature ranks in random forest to linear model

library(ggplot2)
library(ggrepel)
library(reshape2)

randomforest <- data.frame(rank = 1:24, 
                           feature = c("C_Anchor_Max","A_Int_Max","K_Int_Max","K_Int_Sum","T_Int_Sum","A_Int_Sum","C_Anchor_Sum",
                                       "K_Anchor_Max","J_Int_Max","K_Anchor_Sum","J_Anchor_Sum","R_Int_Sum","C_Int_Max","A_Anchor_Sum",
                                       "C_Int_Sum","R_Int_Max","A_Anchor_Max","R_Anchor_Max","J_Int_Sum","R_Anchor_Sum","J_Anchor_Max",
                                       "T_Anchor_Sum", "T_Int_Max", "T_Anchor_Max"))

linearmodel <- data.frame(rank = 1:24, 
                           feature = c("C_Anchor_Max","A_Int_Max","K_Int_Max","J_Int_Max","R_Anchor_Max","A_Anchor_Sum","C_Anchor_Sum",
                                       "R_Anchor_Sum","A_Anchor_Max","K_Anchor_Max","K_Anchor_Sum","A_Int_Sum","K_Int_Sum","R_Int_Sum",
                                       "J_Anchor_Max","C_Int_Max","J_Int_Sum","J_Anchor_Sum","R_Int_Max","C_Int_Sum","T_Int_Max",
                                       "T_Anchor_Max","T_Int_Sum","T_Anchor_Sum"))

merged <- merge(randomforest, linearmodel, by = "feature")
colnames(merged) <- c("feature", "rfRank", "lmRank")
merged <- merged[order(merged$rfRank),]
merged$average <- (merged$rfRank+merged$lmRank)/2

#Scatter plot
ggplot(merged, aes(x = rfRank, y = lmRank)) + 
  geom_point() +
  geom_text_repel(label = merged$feature) + 
  geom_abline(intercept = 0, slope = 1, lty = 3) + 
  theme_classic()

#Line plot
long <- melt(merged, "feature")
long$variable <- factor(long$variable, levels = c("rfRank", "average", "lmRank"))
long$labelA <- ""
long$labelB <- ""
long[long$variable=="rfRank",]$labelA <- long[long$variable=="rfRank",]$feature
long[long$variable=="lmRank",]$labelB <- long[long$variable=="lmRank",]$feature

ggplot(long, aes(x = variable, y = value, group = feature)) + 
  geom_line(alpha = 0.5) + 
  geom_text(label = long$labelA, nudge_x = -0.35) + 
  geom_text(label = long$labelB, nudge_x = 0.35) + 
  theme_classic()


## Alternative version--plot each of the models on top of each other
lm <- read.table("~/Phanstiel Lab Dropbox/Shared Folder/Projects/MEGA/MEGA_paper/GR Resubmission 02/allFeatureCoefficients.txt", 
                 header=T)
lm$absCoef <- abs(lm$coefficient)
lm <- lm[order(lm$absCoef, decreasing = T),]
lm$feature <- factor(lm$feature, levels = lm$feature)

## Absolute coefficient plot: 
lmPlot <- 
  ggplot(lm, aes(x = feature, y = absCoef, fill = model)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = "firebrick") + 
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + 
  ylab("Absolute(Coefficient)") + xlab("")
  

#Key plots
pal <- brewer.pal(8, "YlGnBu")
cols <- c("#8ed368", pal[4:8], pal[c(6,8)], "#71797E", "#D3D3D3")

key <- (strsplit(as.character(lm$feature), "_"))
key <- 
  cbind(lapply(key, `[[`, 1), 
        lapply(key, `[[`, 2), 
        lapply(key, `[[`, 3)) |> 
  as.data.frame()
colnames(key) <- c("element", "position", "measure")
key <- apply(key, 2, as.character)
key <- melt(key)
colnames(key) <- c("x", "variable", "value")
key$variable <- factor(key$variable, levels = c("measure", "position", "element"))
#Replace values
key[key$variable=="element" & key$value=="A",]$value <- "ATAC"
key[key$variable=="element" & key$value=="C",]$value <- "CTCF"
key[key$variable=="element" & key$value=="J",]$value <- "Jun"
key[key$variable=="element" & key$value=="K",]$value <- "K27"
key[key$variable=="element" & key$value=="R",]$value <- "Rad"
key[key$variable=="element"& key$value=="T",]$value <- "RNA"
key$value <- str_replace_all(key$value, "Anchor", "Anc")
key$value <- str_replace_all(key$value, "MaxBW", "Max")
key$value <- str_replace_all(key$value, "BWMax", "Max")
key$value <- factor(key$value, levels = (c("ATAC", "CTCF", "Jun", "K27", "Rad", "RNA", "Anc", "Int", "Sum", "Max")))

lm$x <- unique(key$x)
key <- merge(key, lm)

##Generate key plot
keyplot <- 
  ggplot(key, aes(x = feature, y = variable, fill = value)) + 
  geom_tile(color = "white", alpha = 0.85) + coord_fixed() + geom_text(aes(label=value)) + 
  scale_fill_manual(values = cols) + 
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
                          axis.text.y=element_blank(), axis.ticks.y = element_blank(), legend.position = "") + 
  xlab("") + ylab("")

## Random forest plot: (eventually read in actual data)
randomforest <- data.frame(incMSC = c(40,36,30,26,25.5,25,25,24.6,24.2,24.1,24,23,
                                    22.6,22,20,19.8,19.6,18,18,16,14,14,10,6), 
                           feature = c("C_Anchor_Max","A_Int_Max","K_Int_Max","K_Int_Sum","T_Int_Sum","A_Int_Sum","C_Anchor_Sum",
                                       "K_Anchor_Max","J_Int_Max","K_Anchor_Sum","J_Anchor_Sum","R_Int_Sum","C_Int_Max","A_Anchor_Sum",
                                       "C_Int_Sum","R_Int_Max","A_Anchor_Max","R_Anchor_Max","J_Int_Sum","R_Anchor_Sum","J_Anchor_Max",
                                       "T_Anchor_Sum", "T_Int_Max", "T_Anchor_Max"))
randomforest$feature <- factor(randomforest$feature, levels = lm$feature)
randomforest$mod <- "rf"

## Random forest feature importance actual data: 
randomforest <- read.delim("~/Phanstiel Lab Dropbox/Shared Folder/Projects/MEGA/MEGA_paper/GR Resubmission 02/varimp.val.csv", header = TRUE, sep = ",")
colnames(randomforest) <- c("feature", "IncMSE", "IncNodePurity")
randomforest$feature <- factor(randomforest$feature, levels = lm$feature)
randomforest$mod <- "rf"

## Percent incMSC plot: 
rfPlot <- 
  ggplot(randomforest, aes(x = feature, y = IncMSE, fill = mod)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = "firebrick") + 
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + 
  ylab("%IncMSC") + xlab("") +
  scale_y_reverse()

## Combine everything together: 
library(cowplot)
plot_grid(plotlist = list(lmPlot, keyplot, rfPlot), 
          align = "v", 
          ncol = 1)
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/compareFeatureImportance.pdf", last_plot(), 
       width = 10, height = 10)
