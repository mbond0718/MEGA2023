## Formatting results of makeModelDataframe.R to use in Figure4 

## Load Data
load("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA_shared/lm/megaLMdata.RData")
load("~/Desktop/MEGA/MEGAlmMax.RData")
maxBW <- data
data <- megaLM[,c(7:ncol(megaLM))]

#Adjust column names for bigwig to show that the max values came from bigwigs 
colnames(maxBW) <- str_replace_all(colnames(maxBW), "Max", "MaxBW")

#Merge original dataframe with the maxBW dataframe
data <- cbind(data, maxBW[,7:ncol(maxBW)])

#normalize int sums depending on loop dists : 
data[,grepl("Int_Sum",colnames(data))] <- ((data[,grepl("Int_Sum",colnames(data))]) / (data$L_Length/5000))

## Anchors - taking the sum or multiplying for each loop
anchors <- data[,grepl("Anc",colnames(data))]
anchors[anchors == 0] <- 1
anchors <- anchors[,order(colnames(anchors))]

#remove original anchors
remove <- colnames(data)[grep("Anc", colnames(data))]
data <- data[,-which(colnames(data) %in% remove)]

elements <- list()
generic <- list()
for (i in seq(1, 72, by=12)){
  #max0
  x <- i 
  generic[[1]] <- apply(anchors[,c(x,x+6)], 1, sum)
  generic[[2]] <-apply(anchors[,c(x,x+6)], 1, prod)
  #max72
  x <- i+1
  generic[[3]] <- apply(anchors[,c(x,x+6)], 1, sum)
  generic[[4]] <- apply(anchors[,c(x,x+6)], 1, prod)
  #maxBW0
  x <- i+2
  generic[[5]] <- apply(anchors[,c(x,x+6)], 1, sum)
  generic[[6]] <- apply(anchors[,c(x,x+6)], 1, prod)
  #maxBW72
  x <- i+3
  generic[[7]] <- apply(anchors[,c(x,x+6)], 1, sum)
  generic[[8]] <- apply(anchors[,c(x,x+6)], 1, prod)
  #sum0
  x <- i+4
  generic[[9]] <- apply(anchors[,c(x,x+6)], 1, sum)
  generic[[10]] <- apply(anchors[,c(x,x+6)], 1, prod)
  #sum72
  x <- i+5
  generic[[11]] <- apply(anchors[,c(x,x+6)], 1, sum)
  generic[[12]] <- apply(anchors[,c(x,x+6)], 1, prod)
  
  elements[[i]] <- generic
}

anchorTrans <- data.frame(bind_cols(elements[!unlist(lapply(elements, is.null))]))

generic <- c("Anchor_Max_sum_0", "Anchor_Max_prod_0", "Anchor_Max_sum_72", "Anchor_Max_prod_72", "Anchor_BWMax_sum_0", "Anchor_BWMax_prod_0", "Anchor_BWMax_sum_72", "Anchor_BWMax_prod_72", "Anchor_Sum_sum_0", "Anchor_Sum_prod_0", "Anchor_Sum_sum_72", "Anchor_Sum_prod_72")
elements <- c("A_", "C_", "J_", "K_", "R_", "T_")

names <- list()
for (i in 1:length(elements)){
  names[[i]] <- paste0(elements[i], generic)
}
names <- unlist(names)
colnames(anchorTrans) <- names

#add transformed anchors onto original
data <- cbind(data, anchorTrans)
cols <- colnames(data)
select <- c("L_counts_0", "L_counts_72", cols[grep("Int_Sum", cols)], cols[grep("Int_MaxBW", cols)], cols[grep("Anchor_BWMax_prod", cols)], cols[grep("Anchor_Sum_prod", cols)])

data <- data[colnames(data) %in% select]
colnames(data) <- str_replace_all(colnames(data), "MaxBW", "Max")
colnames(data) <- str_replace_all(colnames(data), "BWMax", "Max")
colnames(data) <- str_replace_all(colnames(data), "Sum_prod", "Sum")
colnames(data) <- str_replace_all(colnames(data), "Max_prod", "Max")

saveRDS(data, "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/model/lmData.rds")




