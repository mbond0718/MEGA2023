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
gainedMatch <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/gainedMatch.rds") |> as.data.frame()
lostMatch <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/lostMatch.rds") |> as.data.frame()

## Overlap matrix 
data <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/model/lmData.rds")

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

#scale all features including y
delta <- data.frame(scale(delta))

#Add column classifying loops as gained, lost, or static
delta$class <- ""
delta[rownames(delta) %in% gained,]$class <- "gained"
delta[rownames(delta) %in% lost,]$class <- "lost"
delta[rownames(delta) %in% static,]$class <- "static"
delta <- delta[!delta$class=="",]


## Generate models with 2x static loops matched for distance and contact
set.seed(123)
gModel <- rbind(delta[rownames(delta) %in% gained,], 
                delta[rownames(delta) %in% rownames(gainedMatch),])
table(gModel$class)
gModel$class <- factor(gModel$class, levels = c("static", "gained"))

lModel <- rbind(delta[rownames(delta) %in% lost,], 
                delta[rownames(delta) %in% rownames(lostMatch),])
table(lModel$class)
lModel$class <- factor(lModel$class, levels = c("static", "lost"))
#Separate into 75% training and 25% testing datasets
gModelTrain <- gModel[rownames(gModel) %in% sample(rownames(gModel), nrow(gModel)*0.75),]
gModelTest <- gModel[!rownames(gModel) %in% rownames(gModelTrain),]

lModelTrain <- lModel[rownames(lModel) %in% sample(rownames(lModel), nrow(lModel)*0.75),]
lModelTest <- lModel[!rownames(lModel) %in% rownames(lModelTrain),]

## 1:1 R2 values from L_counts vs each feature
## For each feature individually, generate R2 values
slope <- function(x, y){
  mean_x <- mean(x)
  mean_y <- mean(y)
  nom <- sum((x - mean_x)*(y-mean_y))
  denom <- sum((x - mean_x)^2)
  m <- nom / denom
  return(m)
}

#Gained Loops
gainedFeatures <- list()
features <- colnames(gModel)[2:(ncol(gModel)-1)]
for (i in 1:length(features)){
  print(i)
  feature <- features[i]
  eq <- paste0("L_counts", " ~ 0 + ", feature)
  tmp <- lm(formula = eq, data = gModel[,c(1, i+1)])
  sl <- slope(unlist(gModel[colnames(gModel) %in% feature]), gModel$L_counts)
  gainedFeatures[[i]] <- c(summary(tmp)$r.squared, sl)
}
gainedFeatures <- do.call(rbind, gainedFeatures) |> 
  as.data.frame()
colnames(gainedFeatures) <- c("R2", "Slope")
gainedFeatures$feature <- features
gainedFeatures <- gainedFeatures[order(gainedFeatures$R2, decreasing = TRUE),]
gainedFeatures$feature <- factor(gainedFeatures$feature, levels = gainedFeatures$feature)
## Plot
ggplot(gainedFeatures, aes(x = feature, y = R2)) + 
  geom_bar(stat = "identity", fill = "#227B7F") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/gainedR2.pdf")
ggplot(gainedFeatures, aes(x = feature, y = Slope)) + 
  geom_bar(stat = "identity", fill = "#227B7F", alpha = 0.75) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/gainedSlope.pdf")


#Lost Loops
lostFeatures <- list()
features <- colnames(lModel)[2:(ncol(lModel)-1)]
for (i in 1:length(features)){
  print(i)
  feature <- features[i]
  eq <- paste0("L_counts", " ~ 0 + ", feature)
  tmp <- lm(formula = eq, data = lModel[,c(1, i+1)])
  sl <- slope(unlist(lModel[colnames(lModel) %in% feature]), lModel$L_counts)
  lostFeatures[[i]] <- c(summary(tmp)$r.squared, sl)
}
lostFeatures <- do.call(rbind, lostFeatures) |> 
  as.data.frame()
colnames(lostFeatures) <- c("R2", "Slope")
lostFeatures$feature <- features
lostFeatures <- lostFeatures[order(lostFeatures$R2, decreasing = TRUE),]
lostFeatures$feature <- factor(lostFeatures$feature, levels = lostFeatures$feature)
## Plot
ggplot(lostFeatures, aes(x = feature, y = R2)) + 
  geom_bar(stat = "identity", fill = "#931C1E") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/lostR2.pdf")
ggplot(lostFeatures, aes(x = feature, y = Slope)) + 
  geom_bar(stat = "identity", fill = "#931C1E", alpha = 0.75) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/lostSlope.pdf")


## Add information about each feature to dataframes to generate legend 
#set color palette
pal <- brewer.pal(8, "YlGnBu")
cols <- c("#8ed368", pal[4:8], pal[c(6,8)], "#71797E", "#D3D3D3")
#cols <- c(pal[c(2:9)], "#71797E", "#D3D3D3")

gainedKey <- (strsplit(as.character(gainedFeatures$feature), "_"))
gainedKey <- 
  cbind(lapply(gainedKey, `[[`, 1), 
        lapply(gainedKey, `[[`, 2), 
        lapply(gainedKey, `[[`, 3)) |> 
  as.data.frame()
colnames(gainedKey) <- c("element", "position", "measure")
gainedKey <- apply(gainedKey, 2, as.character)
gainedKey <- melt(gainedKey)
colnames(gainedKey) <- c("x", "variable", "value")
gainedKey$variable <- factor(gainedKey$variable, levels = c("measure", "position", "element"))
#Replace values
gainedKey[gainedKey$variable=="element" & gainedKey$value=="A",]$value <- "ATAC"
gainedKey[gainedKey$variable=="element" & gainedKey$value=="C",]$value <- "CTCF"
gainedKey[gainedKey$variable=="element" & gainedKey$value=="J",]$value <- "Jun"
gainedKey[gainedKey$variable=="element" & gainedKey$value=="K",]$value <- "K27"
gainedKey[gainedKey$variable=="element" & gainedKey$value=="R",]$value <- "Rad"
gainedKey[gainedKey$variable=="element"& gainedKey$value=="T",]$value <- "RNA"
gainedKey$value <- str_replace_all(gainedKey$value, "Anchor", "Anc")
gainedKey$value <- str_replace_all(gainedKey$value, "MaxBW", "Max")
gainedKey$value <- str_replace_all(gainedKey$value, "BWMax", "Max")
gainedKey$value <- factor(gainedKey$value, levels = (c("ATAC", "CTCF", "Jun", "K27", "Rad", "RNA", "Anc", "Int", "Sum", "Max")))

gainedFeatures$x <- unique(gainedKey$x)
gainedKey <- merge(gainedKey, gainedFeatures)

##Generate key plot
ggplot(gainedKey, aes(x = feature, y = variable, fill = value)) + 
  geom_tile(color = "white", alpha = 0.85) + coord_fixed() + geom_text(aes(label=value)) + 
  scale_fill_manual(values = cols) + 
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
                          axis.text.y=element_blank(), axis.ticks.y = element_blank(), legend.position = "") + 
  xlab("") + ylab("")
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/gainedKey.pdf")

lostKey <- (strsplit(as.character(lostFeatures$feature), "_"))
lostKey <- 
  cbind(lapply(lostKey, `[[`, 1), 
        lapply(lostKey, `[[`, 2), 
        lapply(lostKey, `[[`, 3)) |> 
  as.data.frame()
colnames(lostKey) <- c("element", "position", "measure")
lostKey <- apply(lostKey, 2, as.character)
lostKey <- melt(lostKey)
colnames(lostKey) <- c("x", "variable", "value")
lostKey$variable <- factor(lostKey$variable, levels = c("measure", "position", "element"))
#Replace values
lostKey[lostKey$variable=="element" & lostKey$value=="A",]$value <- "ATAC"
lostKey[lostKey$variable=="element" & lostKey$value=="C",]$value <- "CTCF"
lostKey[lostKey$variable=="element" & lostKey$value=="J",]$value <- "Jun"
lostKey[lostKey$variable=="element" & lostKey$value=="K",]$value <- "K27"
lostKey[lostKey$variable=="element" & lostKey$value=="R",]$value <- "Rad"
lostKey[lostKey$variable=="element"& lostKey$value=="T",]$value <- "RNA"
lostKey$value <- str_replace_all(lostKey$value, "Anchor", "Anc")
lostKey$value <- str_replace_all(lostKey$value, "MaxBW", "Max")
lostKey$value <- str_replace_all(lostKey$value, "BWMax", "Max")
lostKey$value <- factor(lostKey$value, levels = (c("ATAC", "CTCF", "Jun", "K27", "Rad", "RNA", "Anc", "Int", "Sum", "Max")))

lostFeatures$x <- unique(lostKey$x)
lostKey <- merge(lostKey, lostFeatures)

##Generate key plot
ggplot(lostKey, aes(x = feature, y = variable, fill = value)) + 
  geom_tile(color = "white", alpha = 0.85) + coord_fixed() + geom_text(aes(label=value)) + 
  scale_fill_manual(values = cols) + 
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
                          axis.text.y=element_blank(), axis.ticks.y = element_blank(), legend.position = "") + 
  xlab("") + ylab("")
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/lostKey.pdf")

##### LASSO
## Set control
ctl <- trainControl(method="cv", number=10)

## Train on gained data:
g <- train(L_counts ~ 0 + ., data=gModelTrain[,1:(ncol(gModelTrain)-1)], 
           method="glmnet", 
           trControl=ctl,
           trace=0)
lambda <- (g$bestTune$lambda)
idx <- which(g$finalModel$lambda < lambda)[1] #first model where lambda < smallest lambda 
table(g$finalModel$beta[,idx] != 0)
g <- g$finalModel$beta[,idx]
g <- g[order(g)]
g <- data.frame(g)
g$feature <- rownames(g)
g$model <- "gained"
colnames(g) <- c("coefficient", "feature", "model")
#Prepare to plot
g <- g[order(g$coefficient, decreasing = TRUE),]
g$feature <- factor(g$feature, levels = gainedFeatures$feature)
ggplot(g, aes(x = feature, y = model, fill = coefficient)) + 
  geom_tile(color = "white") + coord_fixed() + scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), limits = c(-0.3, 0.3)) +  geom_text(aes(label = round(coefficient, digits = 3))) + 
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
                          axis.text.y=element_blank(), axis.ticks.y = element_blank()) + xlab("") + ylab("")
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/gainedModel.pdf")

## Train on lost data with the following features: 
l <- train(L_counts ~ 0 + ., data=lModelTrain[,1:(ncol(lModelTrain)-1)], 
           method="glmnet", 
           trControl=ctl,
           trace=0)
lambda <- (l$bestTune$lambda)
idx <- which(l$finalModel$lambda < lambda)[1] #first model where lambda < smallest lambda 
table(l$finalModel$beta[,idx] != 0)
l <- l$finalModel$beta[,idx]
l <- l[order(l)]
l <- data.frame(l)
l$feature <- rownames(l)
l$model <- "lost"
colnames(l) <- c("coefficient", "feature", "model")
#Prepare to plot
l <- l[order(l$coefficient, decreasing = TRUE),]
l$feature <- factor(l$feature, levels = lostFeatures$feature)
ggplot(l, aes(x = feature, y = model, fill = coefficient)) + 
  geom_tile(color = "white") + coord_fixed() + scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), limits = c(-0.35, 0.35)) +  geom_text(aes(label = round(coefficient, digits = 3))) + 
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
                          axis.text.y=element_blank(), axis.ticks.y = element_blank()) + xlab("") + ylab("") 
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/lostModel.pdf")

## Use the predict function: 
# Gained 
features <- rownames(g)
gainedFeatures <- c("L_counts", as.character(g[abs(g$coefficient)>0.05,]$feature), "class")
# features <- features[features %in% c("L_counts", "C_Anchor_Max", "A_Anchor_Max", "J_Anchor_Max", "K_Anchor_Max", "R_Anchor_Max", "T_Anchor_Max",
#                          "C_Int_Sum", "A_Int_Max", "J_Int_Max", "K_Int_Max", "Rad_Int_Sum", "T_Int_Max")]
gMod <- gModelTest[colnames(gModelTest) %in% gainedFeatures]
g <- train(L_counts ~ 0 + ., data=gMod[,1:(ncol(gMod)-1)], 
           method="glmnet", 
           trControl=ctl,
           trace=0)
predict <- predict(g)
actual <- gMod$L_counts


predicted <- data.frame(predict, actual)
predicted$sign <- sign(predicted$predict) * sign(predicted$actual)
predicted[rownames(predicted) %in% static,]$sign |> table()
predicted[rownames(predicted) %in% gained,]$sign |> table()
r2 <- (cor(predicted$predict, predicted$actual, method = "pearson"))^2 |> round(3)

data <- rbind(nrow(predicted[rownames(predicted) %in% gained & predicted$sign>0,])/nrow(predicted[rownames(predicted) %in% gained,]),
           nrow(predicted[rownames(predicted) %in% static & predicted$sign>0,])/nrow(predicted[rownames(predicted) %in% static,]))
data <- as.data.frame(data)
data$model <- c("gained", "static")
colnames(data) <- c("percentage", "model")
data$model <- factor(data$model, levels = c("static", "gained"))
ggplot(data, aes(x = model, y = percentage, fill = model)) + 
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("gray50", "#227B7F")) + 
  theme_classic() + theme(legend.position = "none") + 
  ylab("Percentage of loops whose direction was accurately predicted")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/gainedDirection.pdf")

data <- data.frame(predict, actual)
colnames(data) <- c("predict", "actual")
data$class <- gMod$class
ggplot(data, aes(x = actual, y = predict)) + 
  geom_point(col = "#227B7F", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "#227B7F", se = FALSE) +
  annotate("text", x = -3, y = 3, label = paste0("R2 = ", r2)) +
  theme_classic() + 
  xlab("Predicted Loop FC") + ylab("Actual Loop FC") + 
  geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3)
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/gainedPredict.pdf")

# Lost 
features <- rownames(l)
LostFeatures <- c("L_counts", as.character(l[abs(l$coefficient)>0.05,]$feature), "class")
# features <- features[features %in% c("L_counts", "C_Anchor_Max", "A_Anchor_Max", "J_Anchor_Sum", "K_Anchor_Max", "R_Anchor_Max", "T_Anchor_Max",
#                                      "C_Int_Sum", "A_Int_Max", "J_Int_Max", "K_Int_Max", "Rad_Int_Max", "T_Int_Max")]
lMod <- lModelTest[colnames(lModelTest) %in% LostFeatures]
l <- train(L_counts ~ 0 + ., data=lMod[,1:(ncol(lMod)-1)], 
           method="glmnet", 
           trControl=ctl,
           trace=0)
predict <- predict(l)
actual <- lMod$L_counts

predicted <- data.frame(predict, actual)
predicted$sign <- sign(predicted$predict) * sign(predicted$actual)
predicted[rownames(predicted) %in% static,]$sign |> table()
predicted[rownames(predicted) %in% lost,]$sign |> table()
r2 <- (cor(predicted$predict, predicted$actual, method = "pearson"))^2 |> round(3)

data <- rbind(nrow(predicted[rownames(predicted) %in% lost & predicted$sign>0,])/nrow(predicted[rownames(predicted) %in% lost,]),
              nrow(predicted[rownames(predicted) %in% static & predicted$sign>0,])/nrow(predicted[rownames(predicted) %in% static,]))
data <- as.data.frame(data)
data$model <- c("lost", "static")
colnames(data) <- c("percentage", "model")
data$model <- factor(data$model, levels = c("static", "lost"))
ggplot(data, aes(x = model, y = percentage, fill = model)) + 
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("gray50", "#931C1E")) + 
  theme_classic() + theme(legend.position = "none") + 
  ylab("Percentage of loops whose direction was accurately predicted")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/lostDirection.pdf")

data <- data.frame(predict, actual)
data$class <- lMod$class
colnames(data) <- c("predict", "actual")
ggplot(data, aes(x = actual, y = predict)) + 
  geom_point(col = "#931C1E", alpha = 0.75) + 
  geom_smooth(method = "lm", color = "#931C1E", se = FALSE) +
  annotate("text", x = -3, y = 3, label = paste0("R2 = ", r2)) +
  theme_classic() + 
  xlab("Predicted Loop FC") + ylab("Actual Loop FC") + 
  geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3)
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig4/lostPredict.pdf")


### CORRELATION MATRIX
#Correlation matrix for supplemental figure 4
cormat <- round(cor(delta[,2:(ncol(delta)-1)]),2)
#Only get upper triangle: 
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
cormat <- get_lower_tri(cormat)
cormat
cormat <- melt(cormat, na.rm = TRUE)

#Plot 
ggplot(data = cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup4/corrmatrix.pdf")

#legend for corrmatrix: 
#set color palette
pal <- brewer.pal(8, "YlGnBu")
cols <- c("#8ed368", pal[4:8], pal[c(6,8)], "#71797E", "#D3D3D3")

key <- (strsplit(as.character(cormat$Var2), "_"))
key <- strsplit(as.character(unique(cormat$Var2)),)
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

gainedFeatures$x <- unique(key$x)
key <- merge(key, gainedFeatures)

##Generate key plot
ggplot(key, aes(x = x, y = variable, fill = value)) + 
  geom_tile(color = "white", alpha = 0.85) + coord_fixed() + geom_text(aes(label=value)) + 
  scale_fill_manual(values = cols) + 
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), 
                          axis.text.y=element_blank(), axis.ticks.y = element_blank(), legend.position = "") + 
  xlab("") + ylab("")
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup4/cormatKey.pdf")

