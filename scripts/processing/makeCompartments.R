## COMPARTMENT PRE-PROCESSING:
## The eigenvector code provides the compartment calls for all samples for each chromsome individaully. 
## The purpose of this script is to organize these files into a single dataframe: 

## Process ------------------------------------------------

# Read in data
comp0 <- 
  list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/compartments/", "MEGA_0", full.names = TRUE)[1:22]
comp6 <- 
  list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/compartments/", "MEGA_360", full.names = TRUE)[1:22]
comp72 <- 
  list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/raw/hic/compartments/", "MEGA_4320", full.names = TRUE)[1:22]

comps <- list(comp0, comp6, comp72)
chr <- c(1,10:19,2,20:22,3:9)
chrsizes <- read.delim("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/ref/hg38_chromSizes_filt.txt", header=F)

# Compile into one dataframe and assign bins
tplist <- list()
indlist <- list()
for (z in 1:3){
  c <- comps[[z]]
  for (i in 1:length(c)){
    x <- read.table(c[i])
    length <- (chrsizes[chrsizes$V1==chr[i],])$V2
    x$chr <- chr[i]
    x$start <- seq(1, length, by=10000)
    x$stop <- c(seq(10000, length, by=10000), length)
    x <- x[,c(2, 3, 4, 1)]
    colnames(x) <- c("chr", "start", "stop", "score")
    
    indlist[[i]] <- x
  }
  tplist[[z]] <- do.call(rbind, indlist)
  
}
data <- do.call(cbind, tplist)
data <- data[,c(1:3, 4, 8, 12)]
colnames(data) <- c("chr", "start", "end",  "score0", "score6", "score72")
data <- data[order(as.numeric(data$chr)),]

saveRDS(object = data, 
        file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/compartments/compartments.rds")