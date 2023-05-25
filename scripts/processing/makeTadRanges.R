## TAD PRE-PROCESSING:
## The arrowhead program provides TAD ranges
## The purpose of this script is to merge TAD ranges and identify cell-type specific TADs 

## Load libraries -----------------------------------------
library(mariner)

## Process ------------------------------------------------

# Read in data
tad0 <- 
  read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/tads/SCALE_0.bedpe")[,c(1:6,12)]
tad6 <- 
  read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/tads/SCALE_360.bedpe")[,c(1:6,12)]
tad72 <- 
  read.table("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/tads/SCALE_4320.bedpe")[,c(1:6,12)]

# Make list
tads <- list(tad0, tad6, tad72)

#Assign column names
tads <- lapply(tads, setNames, nm = c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "score"))

# Assign row names
for (i in 1:length(tads)){
  tads[[i]]$name <- paste0("TAD", 1:nrow(tads[[i]]))
}

# Convert to ginteractions: 
tads <- lapply(tads, as_ginteractions)

# Merge TAD ranges: 
merged <- mergePairs(x = tads, 
                     binSize = 10E3, 
                     radius = 2, 
                     column = "score")

# Identify "denovo" or cell type specific TADs
dn0 <- deNovo(merged)$`1`
dn6 <- deNovo(merged)$`2`
dn72 <- deNovo(merged)$`3`

dn <- list(dn0, dn6, dn72)

# Generate g_interactions where there is a column for unique TADs
for (i in 1:length(tads)){
  t <- tads[[i]] |> 
    as.data.frame()
  t$unique <- ""
  t[t$name %in% dn[[i]]$name,]$unique <- TRUE
  t[t$unique=="",]$unique <- FALSE
  t <- t[,c(1:3,6:8,12:13)]
  tads[[i]] <- as_ginteractions(t)
}

# Save as rds
saveRDS(object = tads, 
        file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/tads/tadRanges.rds")
