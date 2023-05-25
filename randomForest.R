## Random Forest predicted vs actual loop LFC: 

rfData <- read.csv("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/hic/rownames.pred.L_counts.csv")


#Exclude loops: 
diff72 <- read.csv("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/loops/diff72Loops.txt", 
                   header=F)
diff6 <- read.csv("~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/data/processed/hic/diff6Loops.txt", 
                  header=F)
exclude <- diff6$V1[!diff6$V1 %in% diff72$V1]


#Plot
r2 <- (cor(rfData$pred.L_counts, rfData$L_counts, method = "pearson"))^2 |> round(3)
ggplot(rfData[!rfData$names %in% exclude,], aes(x = pred.L_counts, y = L_counts, color = class)) +
  geom_point(alpha = 0.75, size = 1.5) + 
  scale_color_manual(values = c("#227B7F", "#931C1E", "gray50")) + 
  geom_smooth(inherit.aes = FALSE, aes(x = pred.L_counts, y = L_counts), method = "lm", color = "gray50", se = FALSE) + 
  annotate("text", x = -1.5, y = 1.5, label = paste0("R2 = ", r2)) +
  theme_classic() +  theme(legend.position = "none") + xlab("Predicted Loop LFC") + ylab("Actual Loop LFC") + 
  xlim(-2,2.5) + ylim(-2,2.5) + 
  geom_hline(yintercept = 0, lty = 3) + geom_vline(xintercept = 0, lty = 3) + 
  geom_abline(slope = 1, alpha = 0.35, lty = 2)
ggsave(last_plot(), file = "~/Phanstiel Lab Dropbox/Marielle Bond/MEGA/plots/reviewFigs/randomForestPredict.pdf")


