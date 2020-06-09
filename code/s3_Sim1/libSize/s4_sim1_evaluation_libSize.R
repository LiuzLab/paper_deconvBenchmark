load("./output/Sim1/libSize/sim1_libSize_estW.RData")
load('./output/Sim1/libSize/sim1_mix_libSize.RData')
source("./code/src/Generics_all.R")
source("./code/src/Functions_all.R")
source("./code/src/Methods_evaluation.R")
library(pheatmap)
library(gridExtra)
# Define parameters 
unit <- sim1_params_ob@unit
path <- "./output/Sim1/libSize/"
nGrid <- length(sim1_params_ob@grid)
nDataset <- nMarker <- length(sim1_params_ob@dataset_name)
metric_col <- c("all", "T", "B", "Mono")
line_color_manual <- c("#999999", "#FF9933", "#56B4E9", "#009E73", "#CC00FF", "#0072B2", "#993300", "#99CC33", "#000000", "#FF0000")
methods <- c("DSA", "MMAD", "CAMmarker", "EPIC", "DeconRNASeq", "CIBERSORT", "TIMER", "MuSiC", "CAMfree", "LinSeed") 
for(u in 1:length(unit)){
    for(m in 1:length(methods)){
      estW_name <- paste("sim1", "estW", methods[m], "libSize", unit[u], sep = "_")
      
      cellcor_name <- paste("sim1", "cellcorList", methods[m], "libSize", unit[u], sep = "_")
      cellmad_name <- paste("sim1", "cellmadList", methods[m], "libSize", unit[u], sep = "_")
      
      assign(cellcor_name, deconv_evaluation("cellcor", sim1_params_ob, get(estW_name), sim1_W))
      assign(cellmad_name, deconv_evaluation("cellmad", sim1_params_ob, get(estW_name), sim1_W))
    }
}

# summarize eval_metrics
for(u in 1:length(unit)){
    sum_name_cellcor <- paste("sim1", "sum", "cellcor", "libSize", unit[u], sep = "_")
    sum_name_cellmad <- paste("sim1", "sum", "cellmad", "libSize", unit[u], sep = "_")
    
    assign(sum_name_cellcor, wrap_eval_summary(sim_param = sim1_params_ob, eval_type = "cellcor", all_methods = methods, sim_model_val = "libSize", unit_val = unit[u]))
    assign(sum_name_cellmad, wrap_eval_summary(sim_param = sim1_params_ob, eval_type = "cellmad", all_methods = methods, sim_model_val = "libSize", unit_val = unit[u]))
}

# Draw Figures 
# draw spectrum of celltype_corList 
p_list <- list() 
for(u in 1:length(unit)){
    spectrum_mat <- matrix(0, nrow = length(methods), ncol = nGrid*nMarker*nDataset)
    for(i in 1:length(methods)){
      list_name <- paste("sim1", "cellcorList", methods[i], "libSize", unit[u], sep = "_")
      tmp_list <- list() 
      for(k in 1:length(celltype)){
        tmp_list[[k]] <- evaluation_mat_celltype(get(list_name), sim_param = sim1_params_ob, celltype_val = celltype[k])
      }
      spectrum_mat[i,] <- rowMeans(do.call(cbind, tmp_list), na.rm = TRUE)
    }
    rownames(spectrum_mat) <- methods
    p_list[[u]] <- pheatmap(spectrum_mat, cluster_cols = FALSE, cluster_rows = FALSE, breaks = seq(0, 1, length = 100), main = paste(unit[u]))[[4]]
}
p_concat <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
file_name = paste0(path, "sim1_libSize_spectrum_cellcor.pdf")
ggsave(p_concat, filename = file_name ,width = 18,height=5)

# draw spectrum of celltype_madList 
p_list <- list() 
for(u in 1:length(unit)){
  spectrum_mat <- matrix(0, nrow = length(methods), ncol = nGrid*nMarker*nDataset)
  for(i in 1:length(methods)){
    list_name <- paste("sim1", "cellmadList", methods[i], "libSize", unit[u], sep = "_")
    tmp_list <- list() 
    for(k in 1:length(celltype)){
      tmp_list[[k]] <- evaluation_mat_celltype(get(list_name), sim_param = sim1_params_ob, celltype_val = celltype[k])
    }
    spectrum_mat[i,] <- rowMeans(do.call(cbind, tmp_list), na.rm = TRUE)
  }
  rownames(spectrum_mat) <- methods
  p_list[[u]] <- pheatmap(spectrum_mat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100), cluster_cols = FALSE, cluster_rows = FALSE, main = paste(unit[u]))[[4]]
}
p_concat <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
file_name = paste0(path, "sim1_libSize_spectrum_cellmad.pdf")
ggsave(p_concat, filename = file_name ,width = 18,height=5)



# draw modified main figure 
# cor heatmap 
library(pheatmap)
tmp <- matrix(0, nrow = length(methods), ncol = length(unit))
for(u in 1:length(unit)){
    sum_name_cellcor <- paste("sim1", "sum", "cellcor", "libSize", unit[u], sep = "_")
    tmp[,u] <- get(sum_name_cellcor)[['all']][,'all']
  }
  rownames(tmp) <- methods
  colnames(tmp) <- unit
  p <- pheatmap(tmp, cluster_cols = FALSE, cluster_rows = FALSE, 
                display_numbers = T, 
                breaks = seq(0,1, length = 100),
                angle_col = 45, number_color = "black")[[4]]
file_name = paste0(path,"sim1_libSize_heatmap_cellcor_all_mean.pdf")
ggsave(p, filename = file_name ,width = 4,height=4)

# mad heatmap 
library(pheatmap)
tmp <- matrix(0, nrow = length(methods), ncol = length(unit))
for(u in 1:length(unit)){
  sum_name_cellmad <- paste("sim1", "sum", "cellmad", "libSize", unit[u], sep = "_")
  tmp[,u] <- get(sum_name_cellmad)[['all']][,'all']
}
rownames(tmp) <- methods
colnames(tmp) <- unit
p <- pheatmap(tmp, cluster_cols = FALSE, cluster_rows = FALSE, 
              display_numbers = T, 
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100),
              angle_col = 45, number_color = "black")[[4]]
file_name = paste0(path,"sim1_libSize_heatmap_cellmad_all_mean.pdf")
ggsave(p, filename = file_name ,width = 4,height=4)

# cor rank
tmp <- matrix(0, nrow = length(methods), ncol = length(unit))
for(u in 1:length(unit)){
  sum_name_cellcor <- paste("sim1", "sum", "cellcor", "libSize", unit[u], sep = "_")
  tmp[,u] <- rank(-get(sum_name_cellcor)[['all']][,'all'])
}
rownames(tmp) <- methods
colnames(tmp) <- unit

rank_mat_gg <- melt(tmp)
colnames(rank_mat_gg) <- c("methods", "unit", "rank")
rank_mat_gg$unit <- factor(rank_mat_gg$unit, levels = unit)
rank_mat_gg$rank <- factor(rank_mat_gg$rank)
rank_mat_gg$methods <- factor(rank_mat_gg$methods, levels = methods)
nt <- theme(legend.position='none')
p <- ggplot(data = rank_mat_gg, aes(x = unit, y = rank, group = methods)) + 
  geom_line(aes(color = methods), size = 1) + geom_point(aes(color = methods), size = 2) + 
  theme_bw() + scale_color_manual(values = line_color_manual) + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) , axis.title = element_text(size = 15), text = element_text(size = 15))
file_name = paste0(path,"sim1_libSize_rank_cellcor_all_mean.pdf")
ggsave(p, filename = file_name ,width = 4,height=4)

# mad rank
tmp <- matrix(0, nrow = length(methods), ncol = length(unit))
for(u in 1:length(unit)){
  sum_name_cellmad <- paste("sim1", "sum", "cellmad", "libSize", unit[u], sep = "_")
  tmp[,u] <- rank(get(sum_name_cellmad)[['all']][,'all'])
}
rownames(tmp) <- methods
colnames(tmp) <- unit

rank_mat_gg <- melt(tmp)
colnames(rank_mat_gg) <- c("methods", "unit", "rank")
rank_mat_gg$unit <- factor(rank_mat_gg$unit, levels = unit)
rank_mat_gg$rank <- factor(rank_mat_gg$rank)
rank_mat_gg$methods <- factor(rank_mat_gg$methods, levels = methods)
nt <- theme(legend.position='none')
p <- ggplot(data = rank_mat_gg, aes(x = unit, y = rank, group = methods)) + 
  geom_line(aes(color = methods), size = 1) + geom_point(aes(color = methods), size = 2) + 
  theme_bw() + scale_color_manual(values = line_color_manual) + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) , axis.title = element_text(size = 15), text = element_text(size = 15))
file_name = paste0(path,"sim1_libSize_rank_cellmad_all_mean.pdf")
ggsave(p, filename = file_name ,width = 4,height=4)

save.image("./output/sim1_evaluation_libSize_submission.RData")