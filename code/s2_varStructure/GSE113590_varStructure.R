# analysis of CD8 and PBMC samples 
load("./data/GSE113590_dat.RData")
source("./code/src/Functions_all.R")
library(LSD)
cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
unit <- c("count", "tpm", "fpkm", "cpm")
meta <- GSE113590_dat$meta
path = "./output/varStructure/"
# draw the correlation heatmap
clustering_method <- "complete"
p_list <- list()
for(u in 1:length(unit)){
  dat <- GSE113590_dat[[unit[u]]]
  cor_dat <- cor(dat, method = "spearman")
  p_list[[u]] <- heatmap_forDist(cor_dat, meta = meta, feature = c("phenotype", "tissue"), main = unit[u], clustering_method = clustering_method)[[4]]
}
p_concat <- cowplot::plot_grid(plotlist = p_list, ncol = 2)
file_name <- paste0(path, "GSE113590_distHeatmap_", clustering_method, ".pdf")
ggsave(file_name, p_concat, height = 6, width = 10)

# plot scatter plots for sample pairs with highest and lowest spearman correlations 
for(u in 1:length(unit)){
  dat <- GSE113590_dat[[unit[u]]]
  cor_dat <- cor(dat, method = "spearman")
  idx <- lower.tri(cor_dat, diag = TRUE)
  cor_dat[idx] <- NA
  low_ids <- index.arr.table(cor_dat,n = 10,lowVal = TRUE)
  high_ids <- index.arr.table(cor_dat,n = 10,lowVal = FALSE)
  png(paste0(path, "GSE113590_sampleScatter_spearman_", unit[u], ".png"), width = 15, height = 15, units = "in", res = 300, type = "cairo")
  layout.matrix <- matrix(c(21, 1:5,21, 6:10, rep(23, 6), 
                            22, 11:15, 22, 16:20), nrow = 5, byrow = TRUE)
  layout(mat = layout.matrix, widths = c(1.5,2,2,2,2,2), heights =  c(2,2, 0.5, 2,2))
  layout.show(23)
  for(j in 1:10){
    id_x <- rownames(cor_dat)[low_ids[j,'row']]
    id_y <- rownames(cor_dat)[low_ids[j,'col']]
    draw_scatterFromID(id_x = id_x, id_y = id_y, dat = GSE113590_dat[[unit[u]]], logTrans = TRUE) # plot data in the log scale 
  }
  for(j in 1:10){
    id_x <- rownames(cor_dat)[high_ids[j,'row']]
    id_y <- rownames(cor_dat)[high_ids[j,'col']]
    draw_scatterFromID(id_x = id_x, id_y = id_y, dat = GSE113590_dat[[unit[u]]], logTrans = TRUE) # plot data in the log scale 
  }
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, "lowCor", cex = 2, font = 2)
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, "highCor", cex = 2, font = 2)
  dev.off()
}

# calculate mean & variance & cv 
dat_mean <- list()
dat_sd <- list()
dat_cv <- list()
for(u in 1:length(unit)){
  dat <- GSE113590_dat[[unit[u]]]
  dat_mean[[u]] <- rowMeans(dat)
  dat_sd[[u]] <- apply(dat, MARGIN = 1, sd)
  dat_cv[[u]] <- apply(dat, MARGIN = 1, cv.calc)
  
}
names(dat_mean) <- names(dat_sd) <- names(dat_cv) <- unit

# Plot the mean and variance 
png(paste0(path, "GSE113590_meanVar.png"), width = 3.5, height = 4, units = "in", res = 300, type = "cairo")
x <- dat_mean[['count']]
y <- dat_sd[['count']]^2
heatscatter(x, y, xlab="log2(Mean)", ylab="log2(Variance)", 
            colpal = cols, pch =20, cex.lab = 1.3, main = "count", log = "xy")
abline(a = 0, b = 1, col = "red")
dev.off()
# plot the cv density 
p_list <- list()
for(u in 1:length(unit)){
  # generate data for ggplot2
  gg_cv <- melt(dat_cv[[u]])
  p <- ggplot(gg_cv, aes(x = value)) + geom_density(alpha = 0.7, fill = "red") + ggtitle(unit[u]) + xlab("cv")
  # draw distribution curves from different cell types into different facets 
  p_list[[u]] <- p + theme_bw()
}
gg_concat <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
ggsave(gg_concat, filename = paste0(path, "GSE113590_CV_density.pdf"), width = 7, height = 2)
# save the environment 
save.image(file = paste0(path,"GSE113590_varStructure.RData"))