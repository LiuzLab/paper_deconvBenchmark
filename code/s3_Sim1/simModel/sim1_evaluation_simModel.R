load("./output/Sim1/simModel/sim1_simModel_estW.RData") # Please set the path to the targeted path. 
load('./output/Sim1/simModel/sim1_mix_simModel.RData')
source("./code/src/Generics_all.R")
source("./code/src/Functions_all.R")
source("./code/src/Methods_evaluation.R")
library(pheatmap)
library(gridExtra)
# Define parameters 
unit <- sim1_params_ob@unit
sim_model <- sim1_params_ob@sim_model
path <- "./output/Sim1/simModel/" # Please set the path to the targeted path. 
methods <- c("DSA", "MMAD", "CAMmarker", "EPIC", "DeconRNASeq", "CIBERSORT", "CIBERSORTx", "TIMER", "MuSiC", "CAMfree", "LinSeed") 
nGrid <- length(sim1_params_ob@grid)
nDataset <- nMarker <- length(sim1_params_ob@dataset_name)
celltype <- sim1_params_ob@celltype
line_color_manual <- c("#999999", "#FF9933", "#56B4E9", "#009E73", "#CC00FF", "#0072B2", "#0064ff", "#993300", "#99CC33", "#000000", "#FF0000")
metric_col <- c("all", "T", "B", "Mono")
# Calculate evaluation metrics
for(u in 1:length(unit)){
  for(s in 1:length(sim_model)){
    for(m in 1:length(methods)){
      estW_name <- paste("sim1", "estW", methods[m], sim_model[s], unit[u], sep = "_")
 #     cor_name <- paste("sim1", "corList", methods[m], sim_model[s], unit[u], sep = "_")
 #      mad_name <- paste("sim1", "madList", methods[m], sim_model[s], unit[u], sep = "_")
      cellcor_name <- paste("sim1", "cellcorList", methods[m], sim_model[s], unit[u], sep = "_")
      cellmad_name <- paste("sim1", "cellmadList", methods[m], sim_model[s], unit[u], sep = "_")
      
  #    assign(cor_name, deconv_evaluation("cor",sim1_params_ob, get(estW_name), sim1_W))
  #    assign(mad_name, deconv_evaluation("mad",sim1_params_ob, get(estW_name), sim1_W))
      assign(cellcor_name, deconv_evaluation("cellcor", sim1_params_ob, get(estW_name), sim1_W))
      assign(cellmad_name, deconv_evaluation("cellmad", sim1_params_ob, get(estW_name), sim1_W))
    }
  }
}

# summarize eval_metrics
for(u in 1:length(unit)){
  for(s in 1:length(sim_model)){
    sum_name_cellcor <- paste("sim1", "sum", "cellcor", sim_model[s], unit[u], sep = "_")
    sum_name_cellmad <- paste("sim1", "sum", "cellmad", sim_model[s], unit[u], sep = "_")
 
    assign(sum_name_cellcor, wrap_eval_summary(sim_param = sim1_params_ob, eval_type = "cellcor", all_methods = methods, sim_model_val = sim_model[s], unit_val = unit[u]))
    assign(sum_name_cellmad, wrap_eval_summary(sim_param = sim1_params_ob, eval_type = "cellmad", all_methods = methods, sim_model_val = sim_model[s], unit_val = unit[u]))
  }
}

# Draw Figures 
# Draw spectrum of celltype_corList 
p_concat <- list() 
for(u in 1:length(unit)){
  p_list <- list() 
  for(s in 1:length(sim_model)){
    spectrum_mat <- matrix(0, nrow = length(methods), ncol = nGrid*nMarker*nDataset)
    for(i in 1:length(methods)){
      list_name <- paste("sim1", "cellcorList", methods[i], sim_model[s], unit[u], sep = "_")
      tmp_list <- list() 
      for(k in 1:length(celltype)){
        tmp_list[[k]] <- evaluation_mat_celltype(get(list_name), sim_param = sim1_params_ob, celltype_val = celltype[k])
      }
      tmp_all <- do.call(cbind, tmp_list)
      spectrum_mat[i, ] <- rowMeans(tmp_all)
    }
    rownames(spectrum_mat) <- methods
    p_list[[s]] <- pheatmap(spectrum_mat, breaks = seq(0,1, length = 100), cluster_cols = FALSE, cluster_rows = FALSE, main = paste(sim_model[s]))[[4]]
  }
  title <- cowplot::ggdraw() + cowplot::draw_label(unit[u], fontface = "bold")
  p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  p_concat[[u]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
file_name = paste0(path, "sim1_simModel_spectrum_cellcor.pdf")
ggsave(p_concat_concat, filename = file_name ,width = 15,height=20)



# draw spectrum of celltype_madList 
p_concat <- list() 
for(u in 1:length(unit)){
    p_list <- list() 
    for(s in 1:length(sim_model)){
      spectrum_mat <- matrix(0, nrow = length(methods), ncol = nGrid*nMarker*nDataset)
      for(i in 1:length(methods)){
        list_name <- paste("sim1", "cellmadList", methods[i], sim_model[s], unit[u], sep = "_")
        tmp_list <- list() 
        for(k in 1:length(celltype)){
          tmp_list[[k]] <- evaluation_mat_celltype(get(list_name), sim_param = sim1_params_ob, celltype_val = celltype[k])
        }
        tmp_all <- do.call(cbind, tmp_list)
        spectrum_mat[i, ] <- rowMeans(tmp_all)
      }
      rownames(spectrum_mat) <- methods
      p_list[[s]] <- pheatmap(spectrum_mat,color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100), cluster_cols = FALSE, cluster_rows = FALSE, main = paste(sim_model[s]))[[4]]
    }
    title <- cowplot::ggdraw() + cowplot::draw_label(unit[u], fontface = "bold")
    p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
    p_concat[[u]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
file_name = paste0(path, "sim1_simModel_spectrum_cellmad.pdf")
ggsave(p_concat_concat, filename = file_name ,width = 15,height=20)

# cor heatmap 
library(pheatmap)
p_list <- list() 
for(u in 1:length(unit)){
  tmp <- matrix(0, nrow = length(methods), ncol = length(sim_model))
  for(s in 1:length(sim_model)){
    sum_name_cellcor <- paste("sim1", "sum", "cellcor", sim_model[s], unit[u], sep = "_")
    tmp[,s] <- get(sum_name_cellcor)[['all']][,'all']
    
  }
  rownames(tmp) <- methods
  colnames(tmp) <- sim_model 
  p_list[[u]] <- pheatmap(tmp, cluster_cols = FALSE, cluster_rows = FALSE, 
                          display_numbers = T, main = unit[u], breaks = seq(0,1, length = 100),
                          angle_col = 45, number_color = "black")[[4]]
}
p_concat <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
file_name = paste0(path,"sim1_simModel_heatmap_cellcor_all_mean.pdf")
ggsave(p_concat, filename = file_name ,width = 12,height=4)

# mad heatmap 
library(pheatmap)
p_list <- list() 
for(u in 1:length(unit)){
  tmp <- matrix(0, nrow = length(methods), ncol = length(sim_model))
  for(s in 1:length(sim_model)){
    sum_name_cellmad <- paste("sim1", "sum", "cellmad", sim_model[s], unit[u], sep = "_")
    tmp[,s] <- get(sum_name_cellmad)[['all']][,'all']
    
  }
  rownames(tmp) <- methods
  colnames(tmp) <- sim_model 
  p_list[[u]] <- pheatmap(tmp, color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100),
                          cluster_cols = FALSE, cluster_rows = FALSE, 
                          display_numbers = T, main = unit[u],
                          angle_col = 45, number_color = "black")[[4]]
}
p_concat <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
file_name = paste0(path,"sim1_simModel_heatmap_cellmad_all_mean.pdf")
ggsave(p_concat, filename = file_name ,width = 12,height=4)

# cor rank
library(pheatmap)
p_list <- list() 
for(u in 1:length(unit)){
  tmp <- matrix(0, nrow = length(methods), ncol = length(sim_model))
  for(s in 1:length(sim_model)){
    
    sum_name_cellcor <- paste("sim1", "sum", "cellcor", sim_model[s], unit[u], sep = "_")
    tmp[,s] <- rank(-get(sum_name_cellcor)[['all']][,'all'])
    
  }
  rownames(tmp) <- methods
  colnames(tmp) <- sim_model 
  
  rank_mat_gg <- melt(tmp)
  colnames(rank_mat_gg) <- c("methods", "simModel", "rank")
  rank_mat_gg$simModel <- factor(rank_mat_gg$simModel, levels = c("nb", "lognormal", "normal"))
  rank_mat_gg$rank <- factor(rank_mat_gg$rank)
  rank_mat_gg$methods <- factor(rank_mat_gg$methods, levels = methods)
  nt <- theme(legend.position='none')
  p_list[[u]] <- ggplot(data = rank_mat_gg, aes(x = simModel, y = rank, group = methods)) + 
    geom_line(aes(color = methods), size = 1) + geom_point(aes(color = methods), size = 2) + 
    theme_bw() + scale_color_manual(values = line_color_manual) + 
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 15), text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(unit[u])

}
p_concat <- grid_arrange_common_legend(p_list = p_list, common_title = "", col_ratio = c(0,6,1), legend_position = "right", row_number = 1)
file_name = paste0(path,"sim1_simModel_rank_cellcor_all_mean.pdf")
ggsave(p_concat, filename = file_name ,width = 12,height=4)

# mad rank 
library(pheatmap)
p_list <- list() 
for(u in 1:length(unit)){
  tmp <- matrix(0, nrow = length(methods), ncol = length(sim_model))
  for(s in 1:length(sim_model)){
    
    sum_name_cellmad <- paste("sim1", "sum", "cellmad", sim_model[s], unit[u], sep = "_")
    tmp[,s] <- rank(get(sum_name_cellmad)[['all']][,'all'])
    
  }
  rownames(tmp) <- methods
  colnames(tmp) <- sim_model 
  
  rank_mat_gg <- melt(tmp)
  colnames(rank_mat_gg) <- c("methods", "simModel", "rank")
  rank_mat_gg$simModel <- factor(rank_mat_gg$simModel, levels = c("nb", "lognormal", "normal"))
  rank_mat_gg$rank <- factor(rank_mat_gg$rank)
  rank_mat_gg$methods <- factor(rank_mat_gg$methods, levels = methods)
  nt <- theme(legend.position='none')
  p_list[[u]] <- ggplot(data = rank_mat_gg, aes(x = simModel, y = rank, group = methods)) + 
    geom_line(aes(color = methods), size = 1) + geom_point(aes(color = methods), size = 2) + 
    theme_bw() + scale_color_manual(values = line_color_manual) + 
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 15), text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(unit[u])
  
}
p_concat <- grid_arrange_common_legend(p_list = p_list, common_title = "", col_ratio = c(0,6,1), legend_position = "right", row_number = 1)
file_name = paste0(path,"sim1_simModel_rank_cellmad_all_mean.pdf")
ggsave(p_concat, filename = file_name ,width = 12,height=4)

# Averaged correlation evaluation heatmap ordered by noise level 

p_concat <- list() 
for(u in 1:length(unit)){
  p_list <- list() 
  for(s in 1:length(sim_model)){
    spectrum_mat <- matrix(0, nrow = length(methods), ncol = nGrid)
    for(i in 1:length(methods)){
      list_name <- paste("sim1", "cellcorList", methods[i], sim_model[s], unit[u], sep = "_")
      tmp_list <- list() 
      for(k in 1:length(celltype)){
        tmp <- evaluation_mat_celltype(get(list_name), sim_param = sim1_params_ob, celltype_val = celltype[k])
        p_seq <- seq(1, 90, 9)
        tmp_list[[k]] <- sapply(p_seq, function(i) {mean(tmp[i:(i+8)])})
      }
      tmp_all <- do.call(cbind, tmp_list)
      spectrum_mat[i, ] <- rowMeans(tmp_all)
    }
    rownames(spectrum_mat) <- methods
    colnames(spectrum_mat) <- paste0("P",1:10)
    p_list[[s]] <- pheatmap(spectrum_mat, display_numbers = TRUE, breaks = seq(0,1, length = 100), cluster_cols = FALSE, cluster_rows = FALSE, main = paste(sim_model[s]))[[4]]
  }
  title <- cowplot::ggdraw() + cowplot::draw_label(unit[u], fontface = "bold")
  p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  p_concat[[u]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
file_name = paste0(path, "sim1_simModel_heatmap_cellcor_noise.pdf")
ggsave(p_concat_concat, filename = file_name ,width = 15,height=20)



# Averaged mAD evaluation heatmap ordered by noise level 
p_concat <- list() 
for(u in 1:length(unit)){
  p_list <- list() 
  for(s in 1:length(sim_model)){
    spectrum_mat <- matrix(0, nrow = length(methods), ncol = nGrid)
    for(i in 1:length(methods)){
      list_name <- paste("sim1", "cellmadList", methods[i], sim_model[s], unit[u], sep = "_")
      tmp_list <- list() 
      for(k in 1:length(celltype)){
        tmp <- evaluation_mat_celltype(get(list_name), sim_param = sim1_params_ob, celltype_val = celltype[k])
        p_seq <- seq(1, 90, 9)
        tmp_list[[k]] <- sapply(p_seq, function(i) {mean(tmp[i:(i+8)])})
      }
      tmp_all <- do.call(cbind, tmp_list)
      spectrum_mat[i, ] <- rowMeans(tmp_all)
    }
    rownames(spectrum_mat) <- methods
    colnames(spectrum_mat) <- paste0("P",1:10)
    p_list[[s]] <- pheatmap(spectrum_mat,display_numbers = TRUE, color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100), cluster_cols = FALSE, cluster_rows = FALSE, main = paste(sim_model[s]))[[4]]
  }
  title <- cowplot::ggdraw() + cowplot::draw_label(unit[u], fontface = "bold")
  p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  p_concat[[u]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
file_name = paste0(path, "sim1_simModel_heatmap_cellmad_noise.pdf")
ggsave(p_concat_concat, filename = file_name ,width = 15,height=20)

save.image("./output/Sim1/sim1_evaluation_simModel_submission.RData") # Please set the path to the targeted path. 
