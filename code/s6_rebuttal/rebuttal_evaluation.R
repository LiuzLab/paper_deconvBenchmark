load("./output/rebuttal/rebuttal_estW.RData")
load("./output/rebuttal/rebuttal_mix.RData")
source("./code/src/Generics_all.R")
source("./code/src/Functions_all.R")
source("./code/src/Methods_evaluation.R")
library(pheatmap)

unit <- sim2_params_ob@unit
mix_type <- c("a", "b")
methods <- c("DSA", "MMAD", "CAMmarker", "EPIC", "DeconRNASeq", "CIBERSORT", "CIBERSORTx", "TIMER", "MuSiC", "LinSeed") 

line_color_manual <- c("#999999", "#FF9933", "#56B4E9", "#009E73", "#CC00FF", "#0072B2", "#993300", "#99CC33", "#FF0000")
path <- "./output/rebuttal/"
# calculate evaluation metrics 
for(u in 1:length(unit)){
  for(s in 1:length(mix_type)){
    for(m in 1:length(methods)){
      estW_name <- paste("sim2", "estW", methods[m], unit[u], mix_type[s], sep = "_")
      
      #   cor_name <- paste("sim2", "corList", methods[m], mix_type[s], unit[u], sep = "_")
      #    mad_name <- paste("sim2", "madList", methods[m], mix_type[s], unit[u], sep = "_")
      cellcor_name <- paste("sim2", "cellcorList", methods[m], mix_type[s], unit[u], sep = "_")
      cellmad_name <- paste("sim2", "cellmadList", methods[m], mix_type[s], unit[u], sep = "_")
      cellad_name <- paste("sim2", "celladList", methods[m], mix_type[s], unit[u], sep = "_")
      
      truth_name <- paste("weight", mix_type[s], sep ="_")
      
      # assign(cor_name, deconv_evaluation("cor",sim2_params_ob, get(estW_name), get(truth_name)))
      #  assign(mad_name, deconv_evaluation("mad",sim2_params_ob, get(estW_name), get(truth_name)))
      
      assign(cellcor_name, deconv_evaluation("cellcor", sim2_params_ob, get(estW_name), get(truth_name)))
      assign(cellmad_name, deconv_evaluation("cellmad", sim2_params_ob, get(estW_name), get(truth_name)))
      assign(cellad_name, deconv_evaluation("cellad", sim2_params_ob, get(estW_name), get(truth_name)))
    }
  }
}

# summarize eval_metrics by the method name 
for(u in 1:length(unit)){
  for(s in 1:length(mix_type)){
    sum_name_cellcor <- paste("sim2", "sum", "cellcor", mix_type[s], unit[u], sep = "_")
    sum_name_cellmad <- paste("sim2", "sum", "cellmad", mix_type[s], unit[u], sep = "_")
    
    assign(sum_name_cellcor, wrap_eval_summary(sim_param = sim2_params_ob, eval_type = "cellcor", all_methods = methods, unit_val = unit[u], mix_type_val = mix_type[s]))
    assign(sum_name_cellmad, wrap_eval_summary(sim_param = sim2_params_ob, eval_type = "cellmad", all_methods = methods, unit_val = unit[u], mix_type_val = mix_type[s]))
  }
}

n_comp <- sim2_params_ob@n_comp 

####### cor related figures #######

# heatmap of cellcor for all cell types (Supp. Figures)
library(pheatmap)
for(u in 1:length(unit)){
  p_concat <- list() 
  for(s in 1:length(mix_type)){
    sum_name <- paste("sim2", "sum", "cellcor", mix_type[s], unit[u], sep = "_")
    sum_list <- get(sum_name)
    p_list <- list() 
    for(i in 1:length(n_comp)){
      sum_mat <- sum_list[[i]]
      p_list[[i]] <- pheatmap(sum_mat, cluster_cols = FALSE, cluster_rows = FALSE,
                              display_numbers = T, main = paste("Comp", n_comp[i]), 
                              breaks = seq(0, 1, length = 100), angle_col = 45, number_color = "black")[[4]]
    }
    title <- cowplot::ggdraw() + cowplot::draw_label(mix_type[s], fontface = "bold")
    p_tmp <- cowplot::plot_grid(plotlist = p_list, ncol = 1)
    p_concat[[s]] <- cowplot::plot_grid(title, p_tmp, ncol = 1, rel_heights = c(0.2, 6))
  }
  p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 2)
  file_name = paste0(path,"sim2_heatmap_cellcor_", unit[u], ".pdf")
  ggsave(p_concat_concat, filename = file_name ,width = 10,height=30)
}

# heatmap of cellcor (main figure)
p_concat <- list() 
for(u in 1:length(unit)){
  p_list <- list() 
  for(s in 1:length(mix_type)){
    sum_name <- paste("sim2", "sum", "cellcor", mix_type[s], unit[u], sep = "_")
    sum_list <- get(sum_name)
    sum_mat <- matrix(0, nrow = length(methods), ncol = length(n_comp))
    for(i in 1:length(n_comp)){
      sum_mat[,i] <- sum_list[[i]][,'all']
    }
    rownames(sum_mat) <- methods
    colnames(sum_mat) <- paste("Comp", n_comp)
    
    p_list[[s]] <- pheatmap(sum_mat, cluster_cols = FALSE, cluster_rows = FALSE,
                            display_numbers = T, main = mix_type[s],
                            breaks = seq(0, 1, length = 100), angle_col = 45,
                            number_color = "black")[[4]]
  }
  title <- cowplot::ggdraw() + cowplot::draw_label(unit[u], fontface = "bold")
  p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  p_concat[[u]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
file_name <- paste0(path,"sim2_heatmap_cellcor_all.pdf")
ggsave(p_concat_concat, filename = file_name ,width = 12,height = 18)



####### mad related figures #######
# heatmap of cellmad for all cell types (Supp. Figures)
library(pheatmap)
for(u in 1:length(unit)){
  p_concat <- list() 
  for(s in 1:length(mix_type)){
    sum_name <- paste("sim2", "sum", "cellmad", mix_type[s], unit[u], sep = "_")
    sum_list <- get(sum_name)
    p_list <- list() 
    for(i in 1:length(n_comp)){
      sum_mat <- sum_list[[i]]
      p_list[[i]] <- pheatmap(sum_mat, cluster_cols = FALSE, cluster_rows = FALSE,
                              color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100),
                              display_numbers = T, main = paste("Comp", n_comp[i]), 
                              angle_col = 45, number_color = "black")[[4]]
    }
    title <- cowplot::ggdraw() + cowplot::draw_label(mix_type[s], fontface = "bold")
    p_tmp <- cowplot::plot_grid(plotlist = p_list, ncol = 1)
    p_concat[[s]] <- cowplot::plot_grid(title, p_tmp, ncol = 1, rel_heights = c(0.2, 6))
  }
  p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 2)
  file_name = paste0(path,"sim2_heatmap_cellmad_", unit[u], ".pdf")
  ggsave(p_concat_concat, filename = file_name ,width = 10,height=30)
}

# heatmap of cellmad (main figure)
p_concat <- list() 
for(u in 1:length(unit)){
  p_list <- list() 
  for(s in 1:length(mix_type)){
    sum_name <- paste("sim2", "sum", "cellmad", mix_type[s], unit[u], sep = "_")
    sum_list <- get(sum_name)
    sum_mat <- matrix(0, nrow = length(methods), ncol = length(n_comp))
    for(i in 1:length(n_comp)){
      sum_mat[,i] <- sum_list[[i]][,'all']
    }
    rownames(sum_mat) <- methods
    colnames(sum_mat) <- paste("Comp", n_comp)
    
    p_list[[s]] <- pheatmap(sum_mat, cluster_cols = FALSE, cluster_rows = FALSE,
                            color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100),
                            display_numbers = T, main = mix_type[s],
                            angle_col = 45, number_color = "black")[[4]]
  }
  title <- cowplot::ggdraw() + cowplot::draw_label(unit[u], fontface = "bold")
  p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  p_concat[[u]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
file_name <- paste0(path,"sim2_heatmap_cellmad_all.pdf")
ggsave(p_concat_concat, filename = file_name ,width = 12,height = 18)

# sample scatter 
# scatter plots
for(u in 1:length(unit)){
  for(s in 1:length(mix_type)){
    p_cocat <- list() 
    for(m in 1:length(methods)){
      estW_name <- paste("sim2", "estW", methods[m], unit[u], mix_type[s], sep = "_")
      truth_name <- paste("weight", mix_type[s], sep ="_")
      p_list <- list()
      p_list <- deconv_plot_sampleScatter(sim_param = sim2_params_ob, 
                                          est_W = get(estW_name), 
                                          truth = get(truth_name))
      p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
      title <- cowplot::ggdraw() + cowplot::draw_label(methods[m], fontface = "bold")
      p_concat[[m]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
    }
    p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
    file_name <- paste0(path, "sim2_sampleScatter_", mix_type[s], "_", unit[u], ".pdf")
    ggsave(filename = file_name, p_concat_concat, width = 30, height = 30)
  }
}
# scatter plots (main Figure)
point_color_manual <- c("#999999", "#FF9933", "#56B4E9", "#009E73", "#CC00FF", "#0072B2", "#993300", "#99CC33", "#000000", "#FF0000")
p_concat <- list()
for(s in 1:length(mix_type)){
  p_list <- list()
  for(m in 1:length(methods)){
    estW_name <- paste("sim2", "estW", methods[m], "tpm", mix_type[s], sep = "_")
    estW <- get(estW_name)
    truth_name <- paste("weight", mix_type[s], sep ="_")
    truth <- get(truth_name)
    
    tmp <- estW[[6]]
    celltype <- rownames(truth[[6]])
    if(is.logical(tmp)){
      est_val <- NA 
      est_val_gg <- data.frame(NA)
      p_list[[s]] <- ggplot(data = est_val_gg) + geom_abline(slope = 1, colour = "red") + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 15), text = element_text(size = 15)) +
        ggtitle(paste0("Comp ", n_comp[i]))
    }else{
      if(nrow(tmp) == 10){
        est_val <- as.matrix(t(tmp))
      }
      else if(ncol(tmp) == 10){
        est_val <- as.matrix(tmp)
      }
      colnames(est_val) <- celltype 
      est_val_gg <- cbind(melt(est_val), melt(t(truth[[6]]))$value)
      colnames(est_val_gg) <- c("sample_id", "celltype", "est", "truth")
      p_list[[m]] <- ggplot(data = est_val_gg, aes(x = truth, y = est, group = celltype)) + 
        geom_point(aes(color = celltype), alpha = 0.7) + geom_abline(slope = 1, colour = "red") + 
        scale_x_continuous(limits=c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
        scale_y_continuous(limits=c(min(est_val_gg$est), max(est_val_gg$est)), breaks = seq(from = round(min(est_val_gg$est), 1), to = round(max(est_val_gg$est), 1), length.out = 6)) +
        scale_color_manual(values = point_color_manual) + theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 10), text = element_text(size = 10, face = "bold")) +
        ggtitle(methods[m])
    }
  }
  p_no_legend <- lapply(p_list, function(x) x + theme(legend.position = "none"))
  legend <- cowplot::get_legend(p_list[[1]] + theme(legend.position = "right"))
  title <- cowplot::ggdraw() + cowplot::draw_label(mix_type[s], fontface = "bold")
  p_grid <- cowplot::plot_grid(plotlist = p_no_legend, nrow = 3)
  p_concat[[s]] <- cowplot::plot_grid(title, p_grid, legend, nrow = 1, rel_widths = c(1,7,2))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
file_name <- paste0(path, "sim2_sampleScatter_", mix_type[s], "_", "tpm", "_", "Comp10", ".pdf")
ggsave(filename = file_name, p_concat_concat, width = 10, height = 13)


for(s in 1:length(mix_type)){
  for(i in 1:length(n_comp)){
    p_list <- list()
    for(m in 1:length(methods)){
      eval_list_name <- paste("sim2", "celladList", methods[m], mix_type[s], unit[u], sep = "_")
      eval_list <- get(eval_list_name)
      
      boxplot(t(eval_list[[i]]))
      
      eval_gg <- data.frame()
      eval_gg <- melt(eval_list[[i]])
      colnames(eval_gg) <- c("celltype", "mixture_idx", "ad")
      p_list[[m]] <- ggplot(eval_gg, aes(x = celltype, y = ad, fill = celltype)) + 
        geom_boxplot(alpha = 0.6) +
        stat_summary(fun.y = mean, geom = "point", shape=20, size=2, color="red", fill="red") +
        theme_bw() + scale_fill_brewer(palette="Set3") + xlab("") + ggtitle(methods[m]) + 
        theme(plot.title = element_text(hjust = 0.5), 
              axis.title = element_text(size = 15), 
              axis.text.x = element_text(angle = 45, vjust = 0.6),
              text = element_text(size = 15), legend.position = "none")
      
    }
    p_grid <- cowplot::plot_grid(plotlist = p_list, nrow = 5)
    file_name <- paste0(path, "sim2_ad_boxPlot_", mix_type[s], "_", "tpm", "_", "Comp", n_comp[i], ".pdf")
    ggsave(filename = file_name, p_grid, width = 14, height = 20)
  }
}

save.image("./output/rebuttal/rebuttal_evaluation_submission.RData")