## FIGURE 4 LINEAR MODEL 
## This script prepares data for linear modeling, performs LASSO feature selection, and tests on the final subset of the dataset

## Load libraries -----------------------------------------
library(caret)
library(dplyr)
library(stringr)
library(knitr)
library(corrgram)
library(leaps)
library(reshape2)
library(ggpmisc)

## Read in data -----------------------------------------
## Loop classes: 
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
gained <- loops[loops$class == "gained",] |> rownames()
lost <- loops[loops$class == "lost",] |> rownames()
static <- loops[loops$class == "static",] |> rownames()
gainedMatch <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/gainedMatch.rds") |> as.data.frame() |> 
  rownames()
lostMatch <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/lostMatch.rds") |> as.data.frame() |> 
  rownames()
matched <- c(gainedMatch, lostMatch)
diff72 <- read.csv("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diff72Loops.txt", 
                     header=F)
diff6 <- read.csv("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/diff6Loops.txt", 
                  header=F)

## Overlap matrix 
data <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/model/lmData.rds")

## Pull just differential and matched loops from data frame to include in model: 
data <- data[rownames(data) %in% c(gained, lost, gainedMatch, lostMatch),]

#Separate 0 and 72
data0 <- data[,grepl("_0",colnames(data))]
data72 <- data[,grepl("_72",colnames(data))]
cols <- colnames(data0)
cols <- unlist(strsplit(cols, "_0"))
colnames(data0) <- cols
colnames(data72) <- cols

## Add pseudocount
pseudo <- 1000
data0 <- data0+pseudo
data72 <- data72+pseudo

#get the difference between 0h and 72h, store as delta: 
delta <- log2(data72/data0)

#remove loop counts
delta <- delta[,2:ncol(delta)]
#scale all features including y
delta <- data.frame(scale(delta))
#add back in DeSeq2 fold changes: 
deseqLFC <- read.delim("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/hicSummary.txt", 
                       header = T)
deseqLFC <- deseqLFC[,c(1,4)]
rownames(deseqLFC) <- deseqLFC$Row.names
deseqLFC <- deseqLFC[rownames(deseqLFC) %in% rownames(delta),]
#add actual deseqLFC onto delta: 
delta <- cbind(deseqLFC$log2FoldChange_4320, delta)
colnames(delta)[1] <- "L_counts"

#Add column classifying loops as gained, lost, or static
delta$class <- ""
delta[rownames(delta) %in% gained,]$class <- "gained"
delta[rownames(delta) %in% lost,]$class <- "lost"
delta[rownames(delta) %in% static,]$class <- "static"
delta <- delta[!delta$class=="",]

#Split into training and testing datasets: 
set.seed(458)


model <- delta
model$class <- factor(model$class, levels = c("static", "gained", "lost"))

train <- model[rownames(model) %in% c(sample(gained, length(gained)*0.75), 
                                      sample(lost, length(lost)*0.75), 
                                      sample(matched, length(matched)*0.75)),]
test <- model[!rownames(model) %in% rownames(train),]

#1 to 1 comparisons
allFeatures <- list()
features <- colnames(model)[2:(ncol(model)-1)]
for (i in 1:length(features)){
  print(i)
  feature <- features[i]
  eq <- paste0("L_counts", " ~ 0 + ", feature)
  tmp <- lm(formula = eq, data = model[,c(1, i+1)])
  cor <- cor(model[,1], model[,i+1])
  allFeatures[[i]] <- c(summary(tmp)$r.squared, cor)
}
allFeatures <- do.call(rbind, allFeatures) |> 
  as.data.frame()
colnames(allFeatures) <- c("R2", "Cor")
allFeatures$feature <- features
allFeatures <- allFeatures[order(allFeatures$R2, decreasing = TRUE),]
allFeatures$feature <- factor(allFeatures$feature, levels = allFeatures$feature)
#NEW
allFeatures$sign <- sign(allFeatures$Cor)
allFeatures$R2sign <- allFeatures$R2*allFeatures$sign
allFeatures$model <- "all"

library(scales)
# ggplot(allFeatures, aes(x = feature, y = R2sign, fill = R2sign)) + 
#   geom_bar(stat = "identity") + scale_fill_gradientn(colours = c(brewer.pal(9, "RdBu")[9],"white",brewer.pal(9, "RdBu")[1]), values = rescale(c(-.003,0,0.2)), guide = "colorbar", limits=c(-.003,0.2)) + 
#   theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") + 
#   ylab("R2 * sign of association")

ggplot(allFeatures, aes(x = feature, y = R2sign, fill = sign)) + 
  geom_bar(stat = "identity") + scale_fill_gradientn(colours = c(brewer.pal(9, "RdBu")[9],"white",brewer.pal(9, "RdBu")[1]), values = rescale(c(-1,0,1)), guide = "colorbar", limits=c(-1,1)) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") + 
  ylab("R2 * sign of association") + ylim(-0.005,0.42)
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/allR2signNew.pdf")


#Key plots
pal <- brewer.pal(8, "YlGnBu")
cols <- c("#8ed368", pal[4:8], pal[c(6,8)], "#71797E", "#D3D3D3")
#cols <- c(pal[c(2:9)], "#71797E", "#D3D3D3")

key <- (strsplit(as.character(allFeatures$feature), "_"))
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

allFeatures$x <- unique(key$x)
key <- merge(key, allFeatures)

##Generate key plot
ggplot(key, aes(x = feature, y = variable, fill = value)) + 
  geom_tile(color = "white", alpha = 0.85) + coord_fixed() + geom_text(aes(label=value)) + 
  scale_fill_manual(values = cols) + 
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
                          axis.text.y=element_blank(), axis.ticks.y = element_blank(), legend.position = "") + 
  xlab("") + ylab("")
#ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/allKey.pdf")
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/allKeyNew.pdf")

## LASSO ##
## Set control
#ctl <- trainControl(method="cv", number = 10)
ctl <- trainControl(method="cv", number=10, selectionFunction = "oneSE")
modelR2 <- list()

## Anchor only model 
ancFeatures <- c("L_counts", features[grep("Anchor", features)], "class")
ancTrain <- train[colnames(train) %in% ancFeatures]
ancMod <- train(L_counts ~ 0 + ., data=ancTrain[,1:(ncol(ancTrain)-1)], 
                method="glmnet", 
                trControl=ctl,
                trace=0)
lambda <- (ancMod$bestTune$lambda)
idx <- which(ancMod$finalModel$lambda < lambda)[1] #first model where lambda < smallest lambda 
table(ancMod$finalModel$beta[,idx] != 0)
ancMod <- ancMod$finalModel$beta[,idx]
ancMod <- ancMod[order(ancMod)]
ancMod <- data.frame(ancMod)
ancMod$feature <- rownames(ancMod)
ancMod$model <- "all"
colnames(ancMod) <- c("coefficient", "feature", "model")

ancFeatures <- c("L_counts", ancFeatures[ancFeatures %in% ancMod[abs(ancMod$coefficient)>0,]$feature], "class")
ancTest <- test[colnames(test) %in% ancFeatures]
ancMod <- train(L_counts ~ 0 + ., data=ancTest[,1:(ncol(ancTest)-1)], 
                method="glmnet", 
                trControl=ctl,
                trace=0)
predict <- predict(ancMod)
actual <- ancTest$L_counts
predicted <- data.frame(predict, actual)
predicted$class <- test$class
r2 <- (cor(predicted$predict, predicted$actual, method = "pearson"))^2 |> round(3)
modelR2[[1]] <- r2
#Plot
ggplot(predicted, aes(x = predict, y = actual, color = class)) +
  geom_point(alpha = 0.75, size = 1.5) + 
  scale_color_manual(values = c("gray50", "#227B7F", "#931C1E")) + 
  geom_smooth(inherit.aes = FALSE, aes(x = predict, y = actual), method = "lm", color = "gray50", se = FALSE) + 
  annotate("text", x = -1.5, y = 1.5, label = paste0("R2 = ", r2)) +
  theme_classic() +  theme(legend.position = "none") + xlab("Predicted Loop LFC") + ylab("Actual Loop LFC") + 
  xlim(-2,2.5) + ylim(-2,2.5) + 
  geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3) + 
  geom_abline(slope = 1, alpha = 0.35, lty = 2)
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/anchorPredict.pdf", width = 5, height = 5)


## Interior only model 
intFeatures <- c("L_counts", features[grep("Int", features)], "class")
intTrain <- train[colnames(train) %in% intFeatures]
intMod <- train(L_counts ~ 0 + ., data=intTrain[,1:(ncol(intTrain)-1)], 
                method="glmnet", 
                trControl=ctl,
                trace=0)
lambda <- (intMod$bestTune$lambda)
idx <- which(intMod$finalModel$lambda < lambda)[1] #first model where lambda < smallest lambda 
table(intMod$finalModel$beta[,idx] != 0)
intMod <- intMod$finalModel$beta[,idx]
intMod <- intMod[order(intMod)]
intMod <- data.frame(intMod)
intMod$feature <- rownames(intMod)
intMod$model <- "all"
colnames(intMod) <- c("coefficient", "feature", "model")

intFeatures <- c("L_counts", intFeatures[intFeatures %in% intMod[abs(intMod$coefficient)>0,]$feature], "class")
intTest <- test[colnames(test) %in% intFeatures]
intMod <- train(L_counts ~ 0 + ., data=intTest[,1:(ncol(intTest)-1)], 
                method="glmnet", 
                trControl=ctl,
                trace=0)
predict <- predict(intMod)
actual <- intTest$L_counts
predicted <- data.frame(predict, actual)
predicted$class <- test$class
r2 <- (cor(predicted$predict, predicted$actual, method = "pearson"))^2 |> round(3)
modelR2[[2]] <- r2
#Plot
ggplot(predicted, aes(x = predict, y = actual, color = class)) +
  geom_point(alpha = 0.75, size = 1.5) + 
  scale_color_manual(values = c("gray50", "#227B7F", "#931C1E")) + 
  geom_smooth(inherit.aes = FALSE, aes(x = predict, y = actual), method = "lm", color = "gray50", se = FALSE) + 
  annotate("text", x = -1.5, y = 1.5, label = paste0("R2 = ", r2)) +
  theme_classic() +  theme(legend.position = "none") + xlab("Predicted Loop LFC") + ylab("Actual Loop LFC") + 
  xlim(-2,2.5) + ylim(-2,2.5) + 
  geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3) + 
  geom_abline(slope = 1, alpha = 0.35, lty = 2)
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/interiorPredict.pdf", width = 5, height = 5)

## Train on all features:
mod <- train(L_counts ~ 0 + ., data=train[,1:(ncol(train)-1)], 
             method="glmnet", 
             trControl=ctl,
             trace=0)
lambda <- (mod$bestTune$lambda)
idx <- which(mod$finalModel$lambda < lambda)[1] #first model where lambda < smallest lambda 
table(mod$finalModel$beta[,idx] != 0)
mod <- mod$finalModel$beta[,idx]
mod <- mod[order(mod)]
mod <- data.frame(mod)
mod$feature <- rownames(mod)
mod$model <- "all"
colnames(mod) <- c("coefficient", "feature", "model")

#Prepare to plot
mod <- mod[order(mod$coefficient, decreasing = TRUE),]
mod$feature <- factor(mod$feature, levels = allFeatures$feature)
#Save for comparison to RF: 
write.table(mod, file = "~/Phanstiel Lab Dropbox/Shared Folder/Projects/MEGA/MEGA_paper/GR Resubmission 02/allFeatureCoefficients.txt", row.names = F, quote = F, sep = "\t")
ggplot(mod, aes(x = feature, y = model, fill = coefficient)) + 
  geom_tile(color = "white") + coord_fixed() + scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), limits = c(-0.105, 0.105)) +  geom_text(aes(label = round(coefficient, digits = 3))) + 
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
                          axis.text.y=element_blank(), axis.ticks.y = element_blank()) + xlab("") + ylab("")
#ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/allModel.pdf")
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/allModelNew.pdf")

## Use the predict function: 
features <- rownames(mod) #only run this line once

allFeat <- c("L_counts", as.character(mod[abs(mod$coefficient)>0,]$feature), "class")
aMod <- test[colnames(test) %in% allFeat]
mod <- train(L_counts ~ 0 + ., data=aMod[,1:(ncol(aMod)-1)], 
             method="glmnet", 
             trControl=ctl,
             trace=0)
predict <- predict(mod)
actual <- aMod$L_counts

predicted <- data.frame(predict, actual)
r2 <- (cor(predicted$predict, predicted$actual, method = "pearson"))^2 |> round(3)
modelR2[[3]] <- r2

data <- data.frame(predict, actual)
data$class <- aMod$class
colnames(data) <- c("predict", "actual", "class")

ggplot(data, aes(x = predict, y = actual, color = class)) +
  geom_point(alpha = 0.75, size = 1.5) + 
  scale_color_manual(values = c("gray50", "#227B7F", "#931C1E")) + 
  geom_smooth(inherit.aes = FALSE, aes(x = predict, y = actual), method = "lm", color = "gray50", se = FALSE) + 
  annotate("text", x = -1.5, y = 1.5, label = paste0("R2 = ", r2)) +
  theme_classic() +  theme(legend.position = "none") + xlab("Predicted Loop LFC") + ylab("Actual Loop LFC") + 
  xlim(-2,2.5) + ylim(-2,2.5) + 
  geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3) + 
  geom_abline(slope = 1, alpha = 0.35, lty = 2)
#ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/allPredict.pdf")
#ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/allPredictNew.pdf", width = 5, height = 5)

#version with just coloring points that are differential at 72h 
exclude <- diff6$V1[!diff6$V1 %in% diff72$V1]
ggplot(data[!rownames(data) %in% exclude,], aes(x = predict, y = actual, color = class)) +
  geom_point(alpha = 0.75, size = 1.5) + 
  scale_color_manual(values = c("gray50", "#227B7F", "#931C1E")) + 
  geom_smooth(inherit.aes = FALSE, aes(x = predict, y = actual), method = "lm", color = "gray50", se = FALSE) + 
  annotate("text", x = -1.5, y = 1.5, label = paste0("R2 = ", r2)) +
  theme_classic() +  theme(legend.position = "none") + xlab("Predicted Loop LFC") + ylab("Actual Loop LFC") + 
  xlim(-2,2.5) + ylim(-2,2.5) + 
  geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3) + 
  geom_abline(slope = 1, alpha = 0.35, lty = 2)
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/allPredictNewExclude6h.pdf", width = 5, height = 5)

data <- 
  data.frame(unlist(modelR2), 
           c("Anchor", "Interior", "All")) |> 
  setNames(c("R2", "Model"))
data$Model <- factor(data$Model, levels = c("Anchor", "Interior", "All"))

pal <- brewer.pal(8, "YlGnBu")
ggplot(data, aes(x = Model, y = R2, fill = Model)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c(pal[8], pal[6], pal[3])) + 
  theme_classic() + xlab("") + ylim(-0.005,0.42)
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/anchorIntModels.pdf")


# Sign correlations: 
predicted$sign <- sign(predicted$predict) * sign(predicted$actual)
predicted[rownames(predicted) %in% static,]$sign |> table()
predicted[rownames(predicted) %in% gained,]$sign |> table()
predicted[rownames(predicted) %in% lost,]$sign |> table()

data <- rbind(nrow(predicted[rownames(predicted) %in% gained & predicted$sign>0,])/nrow(predicted[rownames(predicted) %in% gained,]),
              nrow(predicted[rownames(predicted) %in% static & predicted$sign>0,])/nrow(predicted[rownames(predicted) %in% static,]), 
              nrow(predicted[rownames(predicted) %in% lost & predicted$sign>0,])/nrow(predicted[rownames(predicted) %in% lost,]))
data <- as.data.frame(data)
data$model <- c("gained", "static", "lost")
colnames(data) <- c("percentage", "model")
data$model <- factor(data$model, levels = c("static", "gained", "lost"))
ggplot(data, aes(x = model, y = percentage, fill = model)) + 
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("gray50", "#227B7F", "#931C1E")) + 
  theme_classic() + theme(legend.position = "none") + 
  ylab("Percentage of loops whose direction was accurately predicted")
#ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/allDirection.pdf")









