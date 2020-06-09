# Perform analysis to illustrate the variance structure existed in the real data
# We picked GSE60424 as an example, with 134 samples with 6 blood samples derived from 7 disease status
# - for each cell type around 20 samples are available for vairnace structure analysis
# - seven disease status the variance extent

# 1. Correlation heatmap
# 2. Scatter plots
# 3. Mean-variance relationships and coefficient variation density
#####------------------------ source functions & load dataset ------------------------#####
load("./data/GSE60424_dat.RData")
source("./code/src/Functions_all.R")
library(LSD)
library(pheatmap)
library(reshape)
library(ggplot2)
library(gridExtra)
# set the path to store all results
path = "./output/varStructure/"
#####------------------------ dist heatmap ------------------------#####
unit <- names(GSE60424_dat)[-1]
meta <- GSE60424_dat$meta
# dist heatmap for all samples
# draw the correlation heatmap
clustering_method <- "complete"
p_list <- list()
for(u in 1:length(unit)){
  dat <- GSE60424_dat[[unit[u]]]
  cor_dat <- cor(dat, method = "spearman")
  p_list[[u]] <- heatmap_forDist(cor_dat, meta = meta, feature = c("cell_type", "diseasestatus"), main = unit[u], clustering_method = clustering_method)[[4]]
}
p_concat <- cowplot::plot_grid(plotlist = p_list, ncol = 2)
file_name <- paste0(path, "GSE60424_distHeatmap_", clustering_method, ".pdf")
ggsave(file_name, p_concat, height = 8, width = 12)

#####------------------------ scatter plots of highly correlated and lowly correlated sample pairs ------------------------#####
# Categorize samples based on the cell type
celltype_list <- unique(meta[,'cell_type'])
  for(k in 1:length(celltype_list)){
    # calculate correlation
    dat <- GSE60424_dat[['count']]
    cor_dat <- cor(dat, method = "spearman")
    # pick samples
    sample_selector <- which(meta[,'cell_type'] == celltype_list[k])
    sample_id <- paste0(meta[sample_selector,'Run'], "_", "count")
    # cormat subsetting
    cormat_sub <- cor_dat[sample_id,sample_id]
    # eliminate redundant cor values
    idx <- lower.tri(cormat_sub, diag = TRUE)
    cormat_sub[idx] <- NA

    # retrieve sample ids
    low_ids <- index.arr.table(cormat_sub,n=10,lowVal = TRUE)
    high_ids <- index.arr.table(cormat_sub,n=10,lowVal = FALSE)

    png(paste0(path, "GSE60424_sampleScatter_spearman_", celltype_list[k], "_count.png"), width = 15, height = 15, units = "in", res = 300, type = "cairo")
    layout.matrix <- matrix(c(21, 1:5,21, 6:10, rep(23, 6),
                              22, 11:15, 22, 16:20), nrow = 5, byrow = TRUE)
    layout(mat = layout.matrix, widths = c(1.5,2,2,2,2,2), heights =  c(2,2, 0.5, 2,2))
    layout.show(23)
    for(j in 1:10){
      id_x <- rownames(cormat_sub)[low_ids[j,'row']]
      id_y <- rownames(cormat_sub)[low_ids[j,'col']]
      draw_scatterFromID(id_x = id_x, id_y = id_y, dat = dat, logTrans = TRUE) # plot data in the log scale
    }
    for(j in 1:10){
      id_x <- rownames(cormat_sub)[high_ids[j,'row']]
      id_y <- rownames(cormat_sub)[high_ids[j,'col']]
      draw_scatterFromID(id_x = id_x, id_y = id_y, dat = GSE60424_dat[['count']], logTrans = TRUE) # plot data in the log scale
    }
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, "lowCor", cex = 2, font = 2)
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, "highCor", cex = 2, font = 2)
    dev.off()
  }
#####------------------------  mean-variance relationships and coefficient variation density ------------------------#####
# calculate mean & variance & cv
dat_mean <- list()
dat_sd <- list()
dat_cv <- list()
# dat_sd <- list()
for(i in 1:length(unit)){
  dat <- GSE60424_dat[[unit[i]]]
  dat_mean[[i]] <- matrix(NA, nrow=nrow(dat), ncol=length(celltype_list))
  dat_sd[[i]] <- matrix(NA, nrow=nrow(dat), ncol=length(celltype_list))
  dat_cv[[i]] <- matrix(NA, nrow=nrow(dat), ncol=length(celltype_list))
  for(k in 1:length(celltype_list)){
    # pick samples
    sample_selector <- which(meta[, 'cell_type'] == celltype_list[k])
    sample_id <- paste0(meta[sample_selector,'Run'], "_",unit[i])
    # profile subsetting
    dat_sub <- dat[, sample_id]
    # calculate mean&sd
    dat_mean[[i]][,k] <- rowMeans(dat_sub)
    dat_sd[[i]][,k] <- apply(dat_sub, MARGIN = 1, sd)
    dat_cv[[i]][,k] <-  apply(dat_sub, MARGIN = 1, cv.calc)

  }
  colnames(dat_mean[[i]]) <- celltype_list
  colnames(dat_sd[[i]]) <- celltype_list
  colnames(dat_cv[[i]]) <- celltype_list
}
names(dat_mean) <- names(dat_sd) <- names(dat_cv) <- unit
# Plot the mean and variance
png(paste0(path, "GSE60424_meanVar_count.png"), width = 10, height = 6.5, units = "in", res = 300, type = "cairo")
par(mfrow=c(2,4))
for(k in 1:length(celltype_list)){
  library(LSD)
  cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  x <- dat_mean[['count']][, k]
  y <- dat_sd[['count']][, k]^2
  heatscatter(x, y, xlab="log2(Mean)", ylab="log2(Variance)",
              colpal = cols, pch =20, cex.lab = 1.3, main = celltype_list[k], log = "xy")
  abline(a=0,b=1,col="red")
  }
dev.off()

# plot the cv density
p_concat <- list()
for(u in 1:length(unit)){
  dat <- dat_cv[[unit[u]]]
  p_list <- list()
  for(k in 1:length(celltype_list)){
  # generate data for ggplot2
  gg_cv <- melt(dat[,k])
  p_list[[k]] <- ggplot(gg_cv, aes(x = value)) + geom_density(alpha = 0.7, fill = "red") + xlab("cv") + theme_bw() + ggtitle(celltype_list[k])
  }
  title <- cowplot::ggdraw() + cowplot::draw_label(unit[u], fontface = "bold")
  p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  p_concat[[u]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.5,7))
}
gg_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
ggsave(gg_concat_concat, filename = paste0(path, "GSE60424_CV_density.pdf"), width = 15, height = 9)
# save the environment
save.image(file = paste0(path,"GSE60424_varStructure.RData"))