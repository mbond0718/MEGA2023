## SUPPLEMENT FIGURE 2: 
## GO and KEGG enrichment barplots

## Load libraries -----------------------------------------
library(ggplot2)

## Read in data: 
clusters <- (list.files("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/GO/"))[4:5]

## GO 
go <- list()
for (i in 1:length(clusters)){
  go[[i]] <- read.delim(paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/GO/", clusters[i], "/biological_process.txt"))
}
## KEGG
kegg <- list()
for (i in 1:length(clusters)){
  kegg[[i]] <- read.delim(paste0("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/rna/GO/", clusters[i], "/kegg.txt"))
}

## Convert p values for GO and KEGG 
for (i in 1:length(go)){
  g <- go[[i]]
  g <- g[,c(2,4)]
  g$log10pval <- -log10(exp(g$logP))
  go[[i]] <- g[,c(1,3)]
}

go <- do.call(rbind, go)
go <- aggregate(go$log10pval, by = list(go$Term), max)
colnames(go) <- c("Term", "log10pval")
go <- go[order(go$log10pval, decreasing = TRUE),]

ggplot(go[1:50,], aes(x = reorder(Term, +log10pval), y = log10pval)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + theme_classic() + xlab("")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup2/go.pdf", last_plot())

for (i in 1:length(kegg)){
  k <- kegg[[i]]
  k <- k[,c(2,4)]
  k$log10pval <- -log10(exp(k$logP))
  k <- unique(k)
  kegg[[i]] <- k[,c(1,3)]
}

kegg <- do.call(rbind, kegg)
kegg <- aggregate(kegg$log10pval, by = list(kegg$Term), max)
colnames(kegg) <- c("Term", "log10pval")
kegg <- kegg[order(kegg$log10pval, decreasing = TRUE),]

ggplot(kegg[1:50,], aes(x = reorder(Term, +log10pval), y = log10pval)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + theme_classic() + xlab("")
ggsave("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/sup2/kegg.pdf", last_plot())




















