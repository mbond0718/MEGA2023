## GENERATE MATCHES
## Generate matched static loops for gained and lost loops on distance and contact

## Load libraries -----------------------------------------
library(nullranges)

## Obtain loop matrix from mergeLoops script
loopMatrix <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds") |> as.data.frame()
loopMatrix <- loopMatrix[,c(1:3,6:8,21:50)]
#Abbreviate names and replace column names
sampleNames <- colnames(loopMatrix[7:(ncol(loopMatrix))])
sampleNames <- unlist(strsplit(sampleNames, "_inter.hic"))
colnames(loopMatrix)[7:(ncol(loopMatrix))] <- sampleNames
#Set loop ID names for each loop to later access the loop coordinates 
rownames(loopMatrix) <- paste0("loop", 1:nrow(loopMatrix))
#For use with DESeq, get just the counts for each loop 
countMatrix <- loopMatrix[,7:36]

#Merge technical replicates
mergedMatrix <- list()
for (i in 1:3){
  x=c(1,11,21)[i]
  sub <- matrix(apply((countMatrix[,x:(x+9)]), 1, sum))
  mergedMatrix[[i]] <- sub
}
mergedMatrix <- data.frame(do.call(cbind, mergedMatrix))
mergedMatrix <- (mergedMatrix)
#Set columnnames and rownames 
colnames(mergedMatrix) <- c("counts0", "counts360", "counts4320")
rownames(mergedMatrix) <- paste0("loop", 1:nrow(mergedMatrix))
mergedMatrix$meanCounts <- apply(mergedMatrix, 1, mean)
mergedMatrix <- cbind(loopMatrix[,1:6], mergedMatrix)
mergedMatrix$distance <- mergedMatrix$start2 - mergedMatrix$start1
mergedMatrix <- mergedMatrix[mergedMatrix$distance>0,]

## Read in loop lists: 
#loops <- readRDS("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diffLoops.rds")
loops <- read.delim("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/hicSummary.txt")
gained <- loops[loops$cluster == "gained",]$Row.names
lost <- loops[loops$cluster == "lost",]$Row.names
static <- loops[loops$cluster=="static" & !is.na(loops$log2FoldChange_4320),]$Row.names

#Split static randomly in two so gained matches are unique from lost matches: 
staticG <- sample(static, 16206)
staticL <- static[!static %in% staticG]

## Add count/distance information to loops; 
gained <- mergedMatrix[rownames(mergedMatrix) %in% gained,]
gained <- rbind(gained, gained)
lost <- mergedMatrix[rownames(mergedMatrix) %in% lost,]
lost <- rbind(lost, lost)
staticG <- mergedMatrix[rownames(mergedMatrix) %in% staticG,]
staticL <- mergedMatrix[rownames(mergedMatrix) %in% staticL,]

## Match static loops to gained: 
gainedMatch <- matchRanges(focal = gained, 
                           pool = staticG, 
                           covar = ~meanCounts + distance, 
                           replace = F)

lostMatch <- matchRanges(focal = lost, 
                         pool = staticL, 
                         covar = ~meanCounts + distance, 
                         replace = F)

## Compare distributions: 
plot(density(gained$distance), main="Match for distance")
lines(density(gainedMatch$distance), lty=3)
legend("topright", c("gained", "gainedMatched"),
       lty =c(1, 3))
plot(density(gained$meanCounts), main="Match for contact")
lines(density(gainedMatch$meanCounts), lty=3)
legend("topright", c("gained", "gainedMatched"),
       lty =c(1, 3))

plot(density(lost$distance), main="Match for distance")
lines(density(lostMatch$distance), lty=3)
legend("topright", c("lost", "lostMatch"),
       lty =c(1, 3))
plot(density(lost$meanCounts), main="Match for contact")
lines(density(lostMatch$meanCounts), lty=3)
legend("topright", c("lost", "lostMatch"),
       lty =c(1, 3))


## Save matched controls
saveRDS(gainedMatch, "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/gainedMatch.rds")
saveRDS(lostMatch, "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/lostMatch.rds")
