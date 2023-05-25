### CTCF LOOPS###

#The purpose of this script is to identify which loops are CTCF bound (at one or both anchors) and which loops are NOT CTCF bound. 
#This is to address Reviewer 1's concerns about the nature of our loops: 

library(mariner)

#Lopps
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds") #all loops
loops$name <- rownames(loops)
loops <- as_ginteractions(loops)

#CTCF
ctcf <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/cnr/CTCFLFC.rds")

#CTCF 
ctcfGR <- GRanges(seqnames = Rle(ctcf$chr), 
                  ranges = IRanges(start = ctcf$start, end = ctcf$stop), 
                  lfc = ctcf$`4320`,
                  name = rownames(ctcf), 
                  differential = ctcf$differential)

#Intersect: 
subsetByOverlaps(anchors(loops, "first"), ctcfGR)
subsetByOverlaps(anchors(loops, "second"), ctcfGR)

subsetByOverlaps(loops, ctcfGR)



both <- subsetByOverlaps(loops, ctcfGR, use.region="both")
first <- subsetByOverlaps(loops, ctcfGR, use.region="first")
second <- subsetByOverlaps(loops, ctcfGR, use.region="second")

#in first not second
f1 <- first[!first$name %in% second$name]
#in second not first
s1 <- second[!second$name %in% first$name]
#in both
both <- first[first$name %in% second$name]
#in either
either <- unique(c(f1, s1))
neither <- loops[!loops$name %in% c(both$name, either$name)]


#plot: 
data <- 
  data.frame(category = c("both", "one", "neither"), 
             num = c(length(both), length(either), length(neither)))
data$category <- factor(data$category, levels = c("both", "one", "neither"))
data$freq <- ((data$num/33914)*100) |> round(2)
data$x <- "x"

ggplot(data, aes(x=x, y=freq, fill=category)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=20) + 
  scale_fill_manual(values=c("navy", "blue", "gray")) + 
  geom_text(aes(label = paste0(freq, "%")), position = position_stack(vjust=0.5), color = "white") +
  theme_void()
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/ctcfPie.pdf")

#ctcf loops
ctcfLoops <- c(both, either)
(table(ctcfLoops$class)/length(ctcfLoops))*100
write(ctcfLoops$name, file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/ctcfLoops.txt")
#non ctcf loops
otherLoops <- neither
(table(otherLoops$class)/length(otherLoops))*100
write(otherLoops$name, file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/nonCTCFLoops.txt")

### Understanding which loops are detected at which timepoints
loopFiles <- 
  list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/loops/", full.names = TRUE)

# Make list
loops <- lapply(loopFiles, read.table, header=T)
names(loops) <- c("0", "6", "72", "omega")

# Assign row names
for (i in 1:length(loops)){
  loops[[i]]$name <- paste0("loop", 1:nrow(loops[[i]]))
}

# Convert to ginteractions: 
loops <- lapply(loops, as_ginteractions)

## Merge loops and extract counts
mp <-
  mergePairs(x = loops,
             radius = 2,
             column = "APScoreAvg")

#get loops identified in at least each timepoint/omgea
in0 <- subsetBySource(mp, include = "0")
in6 <- subsetBySource(mp, include = "6")
in72 <- subsetBySource(mp, include = "72")
inOmega <- subsetBySource(mp, include = "omega")

data <- 
  data.frame(category = rep(c("CTCF", "nonCTCF"),4), 
             sample = c("0", "0", "6", "6", "72", "72", "omega", "omega"), 
             value = c(in0[in0$name %in% ctcfLoops$name] |> length(), 
                       in0[in0$name %in% otherLoops$name] |> length(), 
                       in6[in6$name %in% ctcfLoops$name] |> length(),
                       in6[in6$name %in% otherLoops$name] |> length(),
                       in72[in72$name %in% ctcfLoops$name] |> length(),
                       in72[in72$name %in% otherLoops$name] |> length(),
                       inOmega[inOmega$name %in% ctcfLoops$name] |> length(),
                       inOmega[inOmega$name %in% otherLoops$name] |> length()))

#make frequencies: 
data$freq <- c(data[1:2,3]/sum(data[1:2,3])*100,
               data[3:4,3]/sum(data[3:4,3])*100,
               data[5:6,3]/sum(data[5:6,3])*100,
               data[7:8,3]/sum(data[7:8,3])*100)

#values: 
ggplot(data, aes(x = sample, y = value, fill = category)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("navy", "gray")) + 
  theme_classic()
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/ctcfPerTimepoint.pdf")
#frequencies
ggplot(data, aes(x = sample, y = freq, fill = category)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("navy", "gray")) + 
  theme_classic()
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/ctcfPerTimepointFreq.pdf")

#proportion at differential loops
loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
loops$loop <- rownames(loops)
gainedLoops <- loops[loops$class=="gained",]
lostLoops <- loops[loops$class=="lost",]
staticLoops <- loops[loops$class=="static",]
diffLoops <- rbind(gainedLoops, lostLoops)

data <- data.frame(category = rep(c("CTCF", "nonCTCF"),3),
                   sample = c("gained", "gained", "lost", "lost", "static", "static"), 
                   value = c(gainedLoops[rownames(gainedLoops) %in% ctcfLoops$name,] |> nrow(),
                             gainedLoops[rownames(gainedLoops) %in% otherLoops$name,] |> nrow(),
                             lostLoops[rownames(lostLoops) %in% ctcfLoops$name,] |> nrow(),
                             lostLoops[rownames(lostLoops) %in% otherLoops$name,] |> nrow(),
                             staticLoops[rownames(staticLoops) %in% ctcfLoops$name,] |> nrow(),
                             staticLoops[rownames(staticLoops) %in% otherLoops$name,] |> nrow()))

#make frequencies: 
data$freq <- c(data[1:2,3]/sum(data[1:2,3])*100,
               data[3:4,3]/sum(data[3:4,3])*100,
               data[5:6,3]/sum(data[5:6,3])*100)

#values: 
ggplot(data, aes(x = sample, y = value, fill = category)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("navy", "gray")) + 
  theme_classic()
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/ctcfLoopClass.pdf")
#frequencies
ggplot(data, aes(x = sample, y = freq, fill = category)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("navy", "gray")) + 
  theme_classic()
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/ctcfLoopClassFreq.pdf")

## version with just differential: 
data <- data.frame(category = rep(c("CTCF", "nonCTCF"),2),
                   sample = c("differential", "differential", "static", "static"), 
                   value = c(diffLoops[rownames(diffLoops) %in% ctcfLoops$name,] |> nrow(),
                             diffLoops[rownames(diffLoops) %in% otherLoops$name,] |> nrow(),
                             staticLoops[rownames(staticLoops) %in% ctcfLoops$name,] |> nrow(),
                             staticLoops[rownames(staticLoops) %in% otherLoops$name,] |> nrow()))

#make frequencies: 
data$freq <- c(data[1:2,3]/sum(data[1:2,3])*100,
               data[3:4,3]/sum(data[3:4,3])*100)

#values: 
ggplot(data, aes(x = sample, y = value, fill = category)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("navy", "gray")) + 
  theme_classic()
#frequencies
ggplot(data, aes(x = sample, y = freq, fill = category)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("navy", "gray")) + 
  theme_classic()
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/ctcfDiffLoopFreq.pdf")

#fisher test: 
fisher.test(matrix(c(1394,109,23495,8916),nrow=2,ncol=2))$p.value #all diff
fisher.test(matrix(c(940,60,23495,8916),nrow=2,ncol=2))$p.value #just gained
fisher.test(matrix(c(454,49,23495,8916),nrow=2,ncol=2))$p.value #just lost 


#are we getting better predictions for ctcf loops or non ctcf loops? 
#this is using data from the all features diffLoopModel.R script: 

#ctcf loops
exclude <- diff6$V1[!diff6$V1 %in% diff72$V1]
cL <- read.delim("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/ctcfLoops.txt", 
                 header=F)
cL <- cL$V1
cL <- data[rownames(data) %in% cL,]
r2 <- (cor(cL$predict, cL$actual, method = "pearson"))^2 |> round(3)

ggplot(cL[!rownames(cL) %in% exclude,], aes(x = predict, y = actual, color = class)) +
  geom_point(alpha = 0.75, size = 1.5) + 
  scale_color_manual(values = c("gray50", "#227B7F", "#931C1E")) + 
  geom_smooth(inherit.aes = FALSE, aes(x = predict, y = actual), method = "lm", color = "gray50", se = FALSE) + 
  annotate("text", x = -1.5, y = 1.5, label = paste0("R2 = ", r2)) +
  theme_classic() +  theme(legend.position = "none") + xlab("Predicted Loop LFC") + ylab("Actual Loop LFC") + 
  xlim(-2,2.5) + ylim(-2,2.5) + 
  geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3) + 
  geom_abline(slope = 1, alpha = 0.35, lty = 2)
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/predictCTCFLoopExclude6h.pdf")

#other loops
oL <- read.delim("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/nonCTCFLoops.txt", 
                 header=F)
oL <- oL$V1
oL <- data[rownames(data) %in% oL,]
r2 <- (cor(oL$predict, oL$actual, method = "pearson"))^2 |> round(3)
ggplot(oL[!rownames(oL) %in% exclude,], aes(x = predict, y = actual, color = class)) +
  geom_point(alpha = 0.75, size = 1.5) + 
  scale_color_manual(values = c("gray50", "#227B7F", "#931C1E")) + 
  geom_smooth(inherit.aes = FALSE, aes(x = predict, y = actual), method = "lm", color = "gray50", se = FALSE) + 
  annotate("text", x = -1.5, y = 1.5, label = paste0("R2 = ", r2)) +
  theme_classic() +  theme(legend.position = "none") + xlab("Predicted Loop LFC") + ylab("Actual Loop LFC") + 
  xlim(-2,2.5) + ylim(-2,2.5) + 
  geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3) + 
  geom_abline(slope = 1, alpha = 0.35, lty = 2)
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/predictNonCTCFLoopExclude.pdf")






