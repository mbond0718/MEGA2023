## LOOP PRE-PROCESSING: 
## SIP provides loop calls
## The purpose of this script is to merge loop calls and extract counts from all hic technical replicates

## Load libraries -----------------------------------------
remotes::install_github("EricSDavis/mariner" ,ref="main" ,auth_token = "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")
library(mariner)
library(magrittr)
library(hictoolsr)

## Process ------------------------------------------------
#List file paths:
hicFiles <- 
  list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/hic/techreps/", full.names = TRUE)
loopFiles <- 
  list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/loops/", full.names = TRUE)

# Make list
loops <- lapply(loopFiles, read.table, header=T)

# Assign row names
for (i in 1:length(loops)){
  loops[[i]]$name <- paste0("loop", 1:nrow(loops[[i]]))
}

# Convert to ginteractions: 
loops <- lapply(loops, as_ginteractions)


## Merge and extract counts: 
## Merge loops and extract counts
mp <-
  mergePairs(x = loops,
             #binSize = 10e03,
             radius = 2,
             column = "APScoreAvg")

counts <- 
  extractCounts(bedpe = mp, 
                hic = hicFiles,
                chroms = c(1:22, "X"),
                res = 10e3,
                norm = "NONE",
                matrix = "observed")

counts <- sort(counts)

# Identify "denovo" or cell type specific loops
dn0 <- deNovo(mp)$`1`
dn6 <- deNovo(mp)$`2`
dn72 <- deNovo(mp)$`3`
dnOm <- deNovo(mp)$`4`

dn <- list(dn0, dn6, dn72, dnOm)

# Generate g_interactions where there is a column for unique loop
for (i in 1:length(loops)){
  l <- loops[[i]] |> 
    as.data.frame()
  l$unique <- ""
  l[l$name %in% dn[[i]]$name,]$unique <- TRUE
  l[l$unique=="",]$unique <- FALSE
  l <- l[,c(1:3,6:8,20:21)]
  loops[[i]] <- as_ginteractions(l)
}

# Save as rds
# merged with unique column
saveRDS(object = loops, 
        file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/uniqueLoops.rds")
# counts
saveRDS(object = counts, 
        file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/loopCounts.rds")

