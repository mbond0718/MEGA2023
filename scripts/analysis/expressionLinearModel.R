## PREP DATA FOR LINEAR MODEL IN FIGURE 5
## This script performs linear modeling to predict changes in gene expression from changes in acetylation and looping 

## Load libraries -----------------------------------------
library(ggplot2)
library(mariner)
library(GenomicRanges)
library(InteractionSet)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## Read in data ------------------------------------------
## Loop classes: 
#loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds") #all loops
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopLFC.rds") #just filtered ones from deseq, use this 
loops$name <- rownames(loops)
loops <- as_ginteractions(loops)
loops <- swapAnchors(loops)
loops[loops$class=="gained" | loops$class=="lost"]$class <- "differential"
static <- loops[loops$class == "static",]
diff <- loops[loops$class == "differential",]
loopCounts <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds") |> as.data.frame()
rownames(loopCounts) <- paste0("loop", rownames(loopCounts))
loopCounts$counts0 <- rowSums(loopCounts[,21:30])
loopCounts$counts6 <- rowSums(loopCounts[,31:40])
loopCounts$counts72 <- rowSums(loopCounts[,41:50])
loopCounts$sumCounts <- rowSums(loopCounts[,21:50])
## Enhancers
enhCounts <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/H3K27Counts.rds")
enhCounts <- enhCounts[!enhCounts$`0`==enhCounts$`3`,] #Remove the peaks that have counts but are 100% STATIC 
enh <- GRanges(seqnames = Rle(enhCounts$chr), 
               ranges = IRanges(start = enhCounts$start, end = enhCounts$stop), name = rownames(enhCounts))

#normalize counts for each feature, rescale: 
enhCounts <- 
  cbind(enhCounts, 
        apply(enhCounts[,4:6], 2, scales::rescale, to = c(1, 100)))
loopCounts <- 
  cbind(loopCounts, 
        apply(loopCounts[,51:53], 2, scales::rescale, to = c(1, 100)))

## Genes
rnaCounts <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaCounts.rds")
rnaLFC <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/rnaLFC.rds")
rnaLFC <- rnaLFC[!is.na(rnaLFC$`4320`),]
#deg <- rnaLFC[rnaLFC$differential==TRUE,] |> rownames()
deg <- rnaLFC[rnaLFC$differential==TRUE & abs(rnaLFC$`4320`)>1,] |> rownames() #all deg
static <- rnaLFC[rnaLFC$differential==FALSE,] |> rownames()
genes <- loadDb("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/ref/hg38txdb.sqlite") |> genes() |> data.frame() #just 72h deg
genes$seqnames <- paste0("chr", genes$seqnames)
rownames(genes) <- genes$gene_id
rna <- merge(rnaCounts, genes, by=0)
rna <- rna[,c(17,12:13,1:9)]
rna$differential <- ""
rna[rna$gene_id %in% deg,]$differential <- "differential"
rna <- GRanges(seqnames = Rle(rna$seqnames), 
               ranges = IRanges(start = rna$start, end = rna$end), 
               name = rna$gene_id, differential = rna$differential)
promoters <- promoters(rna)

## Identify genes that are (1) differential (2) have acetylated promoters and (3) are looped to a distal acetylation peak
# 1. Differential Genes 
diffPromoters <- promoters[promoters$differential=="differential",]
# 2. Differential genes with acetylated promoters
kProm <- subsetByOverlaps(diffPromoters, enh)
# 3. Differential genes with acetylated promoters and looped distal acetylation peaks
loopedProm <- kProm[linkOverlaps(loops, kProm, enh)$subject1]$name |> unique()
# 4. Static gene controls
staticPromoters <- promoters[promoters$name %in% static]
staticK <- subsetByOverlaps(staticPromoters, enh)
staticControl <- staticK[linkOverlaps(loops, staticK, enh)$subject1]$name |> unique()
# 5. Ensure that static genes are matched for median expression 
staticConMatched <- 
  matchRanges(focal = rnaCounts[rownames(rnaCounts) %in% loopedProm,], 
            pool = rnaCounts[rownames(rnaCounts) %in% staticControl,], 
            covar = ~median, method = "stratified")
## Combine loopedProm with control 
#saveRDS(staticConMatched, file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/model/matchedStaticGenes_new.rds")
staticConMatched <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/model/matchedStaticGenes_new.rds")
loopedProm <- c(loopedProm, rownames(staticConMatched))

diffdata <- list()
for (i in 1:length(loopedProm)){
  print(i)
  gene <- loopedProm[i]
  #gene counts at 0h and 72h: 
  # gfc <- (rnaCounts[rownames(rnaCounts) == gene,])[c(1,8)]
  # gfc[gfc == 0] <- 1
  # gfc <- log2(gfc[2]/gfc[1]) #fold change 
  gfc <- rnaLFC[rownames(rnaLFC) == gene,]$`4320`
  
  #promoter k27 counts at 0h and 72h: 
  prom <- subsetByOverlaps(enh, promoters[promoters$name == gene])$name
  pfc <- (enhCounts[rownames(enhCounts) %in% prom,])[,c(4,6)]
  pfc[pfc == 0] <- 1
  if (nrow(pfc)>1){
    pfc <- log2(sum(pfc[,2])/sum(pfc[,1])) #fold change 
  } else if (nrow(pfc)==1){
    pfc <- log2(pfc[2]/pfc[1])
  } else if (nrow(pfc)==0){
    pfc <- 0
  }
  
  #distal k27 counts at 0h and 72h: 
  dist <- enh[linkOverlaps(loops, promoters[promoters$name == gene], enh)$subject2]$name
  dfc <- (enhCounts[rownames(enhCounts) %in% dist,])[,c(4,6)]
  dfc[dfc == 0] <- 1
  if (nrow(dfc)>1){
    dfc <- log2(sum(dfc[,2])/sum(dfc[,1]))
  } else if (nrow(dfc)==1){
    dfc <- log2(dfc[2]/dfc[1])
  } else if (nrow(dfc)==0){
    dfc <- 0
  }
  
  #loop counts at 0h and 72h (fold change): 
  l <- loops[linkOverlaps(loops, promoters[promoters$name == gene], enh)$query]$name
  lfc <- (loopCounts[rownames(loopCounts) %in% l,])[,c(51,53)]
  if (nrow(lfc)>1){
    lfc <- log2(sum(lfc[,2])/sum(lfc[,1]))
  } else if (nrow(lfc)==1){
    lfc <- log2(lfc[2]/lfc[1])
  } else if (nrow(lfc)==0){
    lfc <- 0
  }
  
  #differential loop; 
  diff <- loops[linkOverlaps(loops, promoters[promoters$name == loopedProm[i]], enh)$query]$class
  if (length(diff)>1 & "differential" %in% diff){
    diff <- "differential"
  } else if (length(diff)>1 & !"differential" %in% diff){
    diff <- "static"
  } else if (length(diff)==0){
    diff <- ""
  }
  
  #nearest unlooped enhancer: 
  exclude <- subsetByOverlaps(enh, loops)$name #remove the looped nearest enhancers
  near <- enh[nearest(promoters[promoters$name == gene], enh[!enh$name %in% exclude])]$name
  near <- enhCounts[rownames(enhCounts) %in% near,]
  near[near == 0] <- 1
  nfc <- log(near$`3`/near$`0`)
  if (length(nfc)==0){
    nfc <- 0
  }

  #abc:
  link <- linkOverlaps(loops, promoters[promoters$name == gene], enh)
  colnames(link) <- c("loop", "gene", "enh")
  if(nrow(link)>1 | nrow(link)==1){
    tmp <- list()
    for (x in 1:nrow(link)){
      l <- link[x,1]
      e <- link[x,3]
      
      l <- loops[l]
      e <- enh[e]
      
      l <- (loopCounts[rownames(loopCounts) %in% l$name,])[,c(55,57)]
      e <- (enhCounts[rownames(enhCounts) %in% e$name,])[,c(9,11)]
      
      l[l == 0] <- 1
      e[e == 0] <- 1
      
      tmp[[x]] <- l*e
    }
    tmp <- do.call(rbind, tmp)
    
    if (nrow(tmp)>1){
      abc <- log2(sum(tmp[,2])/sum(tmp[,1]))
    } else if (nrow(tmp)==1){
      abc <- log2(tmp[2]/tmp[1])
    } else if (nrow(tmp)==0){
      abc <- 0
    }
  } else if (nrow(link)==0){
    abc <- 0
  }
  
  diffdata[[i]] <- data.frame(gfc, pfc, nfc, dfc, lfc, abc, diff)
  
}

diffdata <- lapply(diffdata, setNames, c("geneFC", "promoterFC", "nearestFC", "distalFC", "loopFC", "abc", "class"))
diffdata <- do.call(rbind, diffdata) |> as.data.frame()
rownames(diffdata) <- loopedProm

#data <- rbind(diffdata, control) #not sure what this is, there is no object called control? seems like i may have originally done the loop twice--one for control and one for diff loops? 
data <- diffdata
data <- na.omit(data)

## Scale
#data[,1:6] <- scale(data[,1:6])

#Split into test and training datasets
#60% into training with an equal proportion of static and differential loops
#40% into testing with an equal proportion of static and differential loops 

set.seed(3) #seed for all 
set.seed(420) #seed for all
set.seed(4848) #seed for diff 


# st <- data[data$class=="static",]
# df <- data[data$class=="differential",]

df <- data[rownames(data) %in% deg,]
st <- data[rownames(data) %in% static,]


dfD <- df[df$class=="differential",]
dfS <- df[df$class=="static",]

stD <- st[st$class=="differential",]
stS <- st[st$class=="static",]

train <-
  rbind(dfD[rownames(dfD) %in% sample(rownames(dfD), nrow(dfD)*0.6),],
        dfS[rownames(dfS) %in% sample(rownames(dfS), nrow(dfS)*0.6),],
        stD[rownames(stD) %in% sample(rownames(stD), nrow(stD)*0.6),],
        stS[rownames(stS) %in% sample(rownames(stS), nrow(stS)*0.6),])

test <-
  rbind(dfD[!rownames(dfD) %in% rownames(train),],
        dfS[!rownames(dfS) %in% rownames(train),], 
        stD[!rownames(stD) %in% rownames(train),],
        stS[!rownames(stS) %in% rownames(train),])

# train <-
#   rbind(st[rownames(st) %in% sample(rownames(st), nrow(st)*0.6),],
#         df[rownames(df) %in% sample(rownames(df), nrow(df)*0.6),])
# test <-
#   rbind(st[!rownames(st) %in% rownames(train),],
#         df[!rownames(df) %in% rownames(train),])

## Define modeling function; 
features <- c("promoterFC")
features <- c("promoterFC", "nearestFC")
features <- c("promoterFC", "distalFC")
features <- c("promoterFC", "abc")

features <- list(c("promoterFC"), 
                 c("promoterFC", "nearestFC"), 
                 c("promoterFC", "distalFC"),
                 c("promoterFC", "abc"))

trainDF <- train
testDF <- test

modelData <- function(features, trainDF, testDF, modelName){
  #Generate formula
  if(length(features)>1){
    f <- paste(features, collapse = "+")
    formula <- paste0("geneFC ~ ", f)
  } else if (length(features)==1){
    formula <- paste0("geneFC ~ ", features)
  }
  #Model
  mod <- lm(formula, data = trainDF)
  summary(mod)
  testDF$predicted <- predict(mod, newdata = testDF)
  r2all <- (cor(testDF$geneFC, testDF$predicted, method = "pearson"))^2 |> round(3)
  r2diff <- (cor(testDF[testDF$class=="differential",]$geneFC, testDF[testDF$class=="differential",]$predicted, method = "pearson"))^2 |> round(3)
  
  ggplot(testDF, aes(x = geneFC, y = predicted)) +
    geom_point(data = testDF[rownames(testDF) %in% rownames(st),], col = "#c0c0c0", alpha = 0.75) +
    geom_point(data = testDF[rownames(testDF) %in% rownames(df),], col = "#181819", alpha = 0.75) +
    geom_smooth(method = "lm", color = "#c0c0c0", se = FALSE) +
    annotate("text", x = -1, y = 10, label = modelName) +
    annotate("text", x = -1, y = 9, label = paste0("R2 = ", r2all), color = "#c0c0c0") +
    theme_classic() + xlim(-6,13) + ylim(-6,13) +
    xlab("Actual Gene FC") + ylab("Predicted Gene FC") +
    geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3)
  # ggsave(last_plot(), filename = paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig5_new/", modelName, "_withControl.pdf"))
  
  #Plot
  # ggplot(testDF, aes(x = geneFC, y = predicted)) +
  #   geom_point(data = testDF[rownames(testDF) %in% rownames(st) & testDF$class=="differential",], col = "#c0c0c0", alpha = 0.75) +
  #   geom_point(data = testDF[rownames(testDF) %in% rownames(df) & testDF$class=="differential",], col = "#181819", alpha = 0.75) +
  #   geom_smooth(method = "lm", color = "#c0c0c0", se = FALSE) +
  #   annotate("text", x = -1, y = 10, label = modelName) +
  #   annotate("text", x = -1, y = 8, label = paste0("R2 = ", r2diff), color = "#181819") +
  #   theme_classic() + xlim(-5,14) + ylim(-5,10) +
  #   xlab("Actual Gene FC") + ylab("Predicted Gene FC") +
  #   geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3)
  # ggsave(last_plot(), filename = paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup5_new/", modelName, "_withControlsup.pdf"))


  

}

#Generate models: 
modelData(features = features[[1]], trainDF = train, testDF = test, modelName = "model1")
modelData(features = features[[2]], trainDF = train, testDF = test, modelName = "model2")
modelData(features = features[[3]], trainDF = train, testDF = test, modelName = "model3")
modelData(features = features[[4]], trainDF = train, testDF = test, modelName = "model4")

seeds <- c(0.351, 0.35, 0.381, 0.415)
seeds <- c(0.401, 0.401, 0.421, 0.467)

# Random 1000 permutations to get coefficient values for to test for significance: 
coef <- list()
for (i in 1:1000){
  print(i)
  c <- list()
  #Subset randomly train and test 
  subTrain <- rbind(st[sample.int(289, 173, replace = F),], 
                    df[sample.int(43, 25, replace = F),])
  subTest <- data[!rownames(data) %in% rownames(subTrain),]
  
  #Define function
  modelData <- function(features, trainDF, testDF, modelNum){
    #Generate formula
    if(length(features)>1){
      f <- paste(features, collapse = "+")
      formula <- paste0("geneFC ~ ", f)
    } else if (length(features)==1){
      formula <- paste0("geneFC ~ ", features)
    }
    #Model
    mod <- lm(formula, data = trainDF)
    c[[modelNum]] <<- mod$coefficients[-1]
  }
  
  #Model: 
  modelData(features = c("promoterFC"), trainDF = subTrain, testDF = subTest, modelNum = 1)
  modelData(features = c("promoterFC", "nearestFC"), trainDF = subTrain, testDF = subTest, modelNum = 2)
  modelData(features = c("promoterFC", "distalFC"), trainDF = subTrain, testDF = subTest, modelNum = 3)
  modelData(features = c("promoterFC", "abc"), trainDF = subTrain, testDF = subTest, modelNum = 4)
  
  coef[[i]] <- unlist(c)
  
}
coef <- do.call(rbind, coef) |> as.data.frame()
colnames(coef) <- as.character(1:7)
coef <- melt(coef)

coefData <- list()
for (i in 1:length(unique(coef$variable))){
  sub <- coef[coef$variable==i,]
  coefData[[i]] <- c(mean(sub$value), 
                 sd(sub$value))
}
coefData <- do.call(rbind, coefData) |> 
  as.data.frame()
coefData$x <- as.character(1:7)
colnames(coefData) <- c("Coefficient", "sd", "Variable")
coefData$feature <- c("Promoter", "Promoter", "Nearest", "Promoter", "Distal", "Promoter", "DistalxLoop")
coefData$feature <- factor(coefData$feature, levels = c("Promoter", "Nearest", "Distal", "DistalxLoop"))
coefData$x <- c("1", "2", "2", "3", "3", "4", "4")

# ggplot(data = coefData, aes(x = x, y = Coefficient, fill = feature)) + 
#   geom_bar(stat = "identity") + 
#   #geom_errorbar(aes(ymin=Coefficient-sd, ymax=Coefficient+sd), width=.2) +
#   #scale_fill_manual(values = c("#3C756B", "#44A391", "#44A391", "#638ECA", "#638ECA", "#5B718C", "#5B718C")) + 
#   scale_fill_manual(values = c("#3C756B", "#44A391", "#638ECA", "#5B718C")) + 
#   theme_classic() + theme(legend.position = "none")

ggplot(data = coefData, aes(x = x, y = Coefficient, fill = feature)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  #geom_errorbar(aes(ymin=Coefficient-sd, ymax=Coefficient+sd), width=.2) +
  #scale_fill_manual(values = c("#3C756B", "#44A391", "#44A391", "#638ECA", "#638ECA", "#5B718C", "#5B718C")) + 
  scale_fill_manual(values = c("#3C756B", "#44A391", "#638ECA", "#5B718C")) + ylim(-0.1,1) + ylab("Estimate") +  
  theme_classic() + theme(legend.position = "none") + 
  geom_errorbar(aes(ymin = Coefficient - sd, ymax = Coefficient + sd), width = 0.2, position = position_dodge(0.9))
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/fig5_new/Estimate.pdf")

# Random 1000 permutations to get R2 values for to test for significance: 
r2 <- list()
for (i in 1:1000){
  print(i)
  r2all <- list()
  r2diff <- list()
  #Subset randomly train and test 
  subTrain <- rbind(st[sample.int(289, 173, replace = F),], 
                    df[sample.int(43, 25, replace = F),])
  subTest <- data[!rownames(data) %in% rownames(subTrain),]
  #Define function
  modelData <- function(features, trainDF, testDF, modelNum){
    #Generate formula
    if(length(features)>1){
      f <- paste(features, collapse = "+")
      formula <- paste0("geneFC ~ ", f)
    } else if (length(features)==1){
      formula <- paste0("geneFC ~ ", features)
    }
    #Model
    mod <- lm(formula, data = trainDF)
    testDF$predicted <- predict(mod, newdata = testDF)
    r2all[modelNum] <<- (cor(testDF$geneFC, testDF$predicted, method = "pearson"))^2 |> round(3)
    r2diff[modelNum] <<- (cor(testDF[testDF$class=="differential",]$geneFC, testDF[testDF$class=="differential",]$predicted, method = "pearson"))^2 |> round(3)
  }
  
  #Model: 
  modelData(features = c("promoterFC"), trainDF = subTrain, testDF = subTest, modelNum = 1)
  modelData(features = c("promoterFC", "nearestFC"), trainDF = subTrain, testDF = subTest, modelNum = 2)
  modelData(features = c("promoterFC", "distalFC"), trainDF = subTrain, testDF = subTest, modelNum = 3)
  modelData(features = c("promoterFC", "abc"), trainDF = subTrain, testDF = subTest, modelNum = 4)
  
  r2[[i]] <- c(unlist(r2all), unlist(r2diff))
  
}

r2 <- do.call(rbind, r2) |> as.data.frame()
colnames(r2) <- c("allModel1", "allModel2", "allModel3", "allModel4", "diffModel1", "diffModel2", "diffModel3", "diffModel4")
r2 <- melt(r2)
r2all <- r2[r2$variable %in% c("allModel1", "allModel2", "allModel3", "allModel4"),]
r2diff <- r2[r2$variable %in% c("diffModel1", "diffModel2", "diffModel3", "diffModel4"),]

#add red highlight points
plotSeed <- 
  data.frame(variable = c("allModel1", "allModel2", "allModel3", "allModel4"), 
           value = seeds)

ggplot(data = r2all, aes(x = variable, y = value, fill = variable)) + 
  geom_jitter(alpha = 0.25, width = 0.2, color = "gray") +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
  geom_point(data = plotSeed, color = "red") + 
  scale_fill_manual(values = c("#85d1d8", "#64afd4", "#7584c1", "#9454a1")) + 
  theme_classic()
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/fig5.pdf")

plotSeed <- 
  data.frame(variable = c("diffModel1", "diffModel2", "diffModel3", "diffModel4"), 
             value = c(0.674, 0.673, 0.678, 0.763)) #for diff

ggplot(data = r2diff, aes(x = variable, y = value, fill = variable)) + 
  geom_jitter(alpha = 0.25, width = 0.2, color = "gray") +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) + 
  geom_point(data = plotSeed, color = "red") + 
  scale_fill_manual(values = c("#85d1d8", "#64afd4", "#7584c1", "#9454a1")) + 
  theme_classic()
ggsave(last_plot(), filename = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/fig5S.pdf")


#Stat tests: 
wilcox.test(r2[r2$variable=="allModel1",]$value, r2[r2$variable=="allModel2",]$value)
wilcox.test(r2[r2$variable=="allModel1",]$value, r2[r2$variable=="allModel3",]$value)
wilcox.test(r2[r2$variable=="allModel1",]$value, r2[r2$variable=="allModel4",]$value)
wilcox.test(r2[r2$variable=="allModel3",]$value, r2[r2$variable=="allModel4",]$value)


wilcox.test(r2[r2$variable=="diffModel1",]$value, r2[r2$variable=="diffModel2",]$value)
wilcox.test(r2[r2$variable=="diffModel1",]$value, r2[r2$variable=="diffModel3",]$value)
wilcox.test(r2[r2$variable=="diffModel1",]$value, r2[r2$variable=="diffModel4",]$value)
wilcox.test(r2[r2$variable=="diffModel3",]$value, r2[r2$variable=="diffModel4",]$value)


