load("./output/Sim3/sim3_estW.RData")
load("./output/Sim3/sim3_mix.RData")
source("./code/src/Generics_all.R")
source("./code/src/Functions_all.R")
source("./code/src/Methods_evaluation.R")
library(pheatmap)

unit <- sim3_params_ob@unit
methods <- c("DSA", "MMAD", "CAMmarker", "EPIC", "EPICabsolute", "DeconRNASeq",  "CIBERSORT", "TIMER","TIMERtumor", "MuSiC", "LinSeed") 
n_comp <- sim3_params_ob@n_comp
mix_type <- sim3_params_ob@mix_type
path <- "./output/Sim3/"
line_color_manual <- c("#999999", "#FF9933", "#56B4E9", "#009E73", "#00ebab", "#CC00FF", "#0072B2", "#993300", "#e64d00", "#99CC33", "#FF0000")
# manipulate epic weight estimation for further evaluation 
# calculate evaluation metrics 
weight_type <- c("absolute", "relative")
for(u in 1:length(unit)){
  for(s in 1:length(mix_type)){
    for(j in 1:length(tumor_content)){
      for(m in 1:length(methods)){
        for(w in 1:length(weight_type)){
          estW_name <- paste("sim3", "estW", methods[m], unit[u], mix_type[s], tumor_content[j], sep = "_")
          
          cellcor_name <- paste(weight_type[w], "sim3", "cellcorList", methods[m], mix_type[s], tumor_content[j], unit[u], sep = "_")
          cellmad_name <- paste(weight_type[w], "sim3", "cellmadList", methods[m], mix_type[s], tumor_content[j], unit[u], sep = "_")
          
          estW <- get(estW_name)
          if(methods[m] == "EPICabsolute"){
            estW <- lapply(estW, function(x) x[[1]])
          }
          
          if(weight_type[w] == "absolute"){
            truth_name <- paste("sim3", "W", mix_type[s], tumor_content[j], sep ="_")
          }else if(weight_type[w] == "relative"){
            truth_name <- paste("sim3", "W", mix_type[s], sep = "_")
          }
          assign(cellcor_name, deconv_evaluation("cellcor", sim3_params_ob, estW, get(truth_name), weight_type = weight_type[w]))
          assign(cellmad_name, deconv_evaluation("cellmad", sim3_params_ob, estW, get(truth_name), weight_type = weight_type[w]))
          
        }
      }
    }
  }
}

# summarize eval_metrics by the method name    
for(u in 1:length(unit)){
  for(s in 1:length(mix_type)){
    for(j in 1:length(tumor_content)){
      for(w in 1:length(weight_type)){
        sum_name_cellcor <- paste(weight_type[w], "sim3", "sum", "cellcor", mix_type[s], tumor_content[j], unit[u], sep = "_")
        sum_name_cellmad <- paste(weight_type[w], "sim3", "sum", "cellmad", mix_type[s], tumor_content[j], unit[u], sep = "_")
        
        assign(sum_name_cellcor, wrap_eval_summary(sim_param = sim3_params_ob, eval_type = "cellcor", all_methods = methods, unit_val = unit[u], mix_type_val = mix_type[s], tumor_content_val= tumor_content[j], weight_type_val = weight_type[w]))
        assign(sum_name_cellmad, wrap_eval_summary(sim_param = sim3_params_ob, eval_type = "cellmad", all_methods = methods, unit_val = unit[u], mix_type_val = mix_type[s], tumor_content_val = tumor_content[j], weight_type_val = weight_type[w]))
      }
    }
  }
}

# Draw figures
#--- summarized by tumor content, heatmap, cor ---#
p_concat <- list() 
for(w in 1:length(weight_type)){
  p_list <- list()
  for(s in 1:length(mix_type)){
    sum_mat <- list()
    for(c in 1:length(n_comp)){
      sum_mat[[c]] <- matrix(0, nrow = length(methods), ncol = length(tumor_content))
      for(j in 1:length(tumor_content)){
        sum_name <- paste(weight_type[w], "sim3", "sum", "cellcor", mix_type[s], tumor_content[j], "tpm", sep = "_")
        sum_list <- get(sum_name)
        sum_mat[[c]][,j] <- sum_list[[c]][,'all']
      }
      rownames(sum_mat[[c]]) <- methods
      colnames(sum_mat[[c]]) <- tumor_content
    }
    p_list[[s]] <- pheatmap(Reduce("+", sum_mat) / length(n_comp), cluster_cols = FALSE, cluster_rows = FALSE,
                            display_numbers = T, main = mix_type[s],
                            breaks = seq(0, 1, length = 100), angle_col = 45,
                            number_color = "black")[[4]]
  }
  title <- cowplot::ggdraw() + cowplot::draw_label(weight_type[w], fontface = "bold")
  p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  p_concat[[w]] <-  cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
file_name = paste0(path,"sim3_heatmap_cellcor_sum_all_tpm", ".pdf")
ggsave(p_concat_concat, filename = file_name ,width = 8,height = 8)

#--- summarized by tumor content, heatmap, mad ---#
p_concat <- list() 
for(w in 1:length(weight_type)){
  p_list <- list()
  for(s in 1:length(mix_type)){
    sum_mat <- list()
    for(c in 1:length(n_comp)){
      sum_mat[[c]] <- matrix(0, nrow = length(methods), ncol = length(tumor_content))
      for(j in 1:length(tumor_content)){
        sum_name <- paste(weight_type[w], "sim3", "sum", "cellmad", mix_type[s], tumor_content[j], "tpm", sep = "_")
        sum_list <- get(sum_name)
        sum_mat[[c]][,j] <- sum_list[[c]][,'all']
      }
      rownames(sum_mat[[c]]) <- methods
      colnames(sum_mat[[c]]) <- tumor_content
    }
    p_list[[s]] <- pheatmap(Reduce("+", sum_mat) / length(n_comp), cluster_cols = FALSE, cluster_rows = FALSE,
                            display_numbers = T, main = mix_type[s],
                            angle_col = 45,
                            color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100),
                            number_color = "black")[[4]]
  }
  title <- cowplot::ggdraw() + cowplot::draw_label(weight_type[w], fontface = "bold")
  p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  p_concat[[w]] <-  cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
file_name = paste0(path,"sim3_heatmap_cellmad_sum_all_tpm", ".pdf")
ggsave(p_concat_concat, filename = file_name ,width = 8,height = 8)

#--- scatter plots, comp 5 ---#
point_color_manual <- c("#FF9933", "#56B4E9", "#009E73", "#CC00FF", "#0072B2")
# mosaic.small.large.absolute and mosaic.small.large.relative 
for( w in 1:length(weight_type)){
  p_concat <- list() 
  for(j in 1:length(tumor_content)){
    p_list <- list()
    for(m in 1:length(methods)){
      estW_name <- paste("sim3", "estW", methods[m], "tpm", "real", tumor_content[j], sep = "_")
      
      estW <- get(estW_name)
      if(methods[m] == "EPICabsolute"){
        estW <- lapply(estW, function(x) x[[1]])
      }
      
      if(weight_type[w] == "absolute"){
        truth_name <- paste("sim3", "W", "real", tumor_content[j], sep ="_")
      }else if(weight_type[w] == "relative"){
        truth_name <- paste("sim3", "W", "real", sep = "_")
      }
      truth <- get(truth_name)[[1]][1:5,]
      tmp <- estW[[1]]
      celltype <- rownames(truth)
      if(is.logical(tmp)){
        est_val <- NA 
        est_val_gg <- data.frame(NA)
        p_list[[s]] <- ggplot(data = est_val_gg) + geom_abline(slope = 1, colour = "red") + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 15), text = element_text(size = 15)) +
          ggtitle(methods[m])
      }else{
        if(nrow(tmp) == 5){
          est_val <- as.matrix(t(tmp))
        }else if(ncol(tmp) == 5){
          est_val <- as.matrix(tmp)
        }
        colnames(est_val) <- c("T", "B", "Mono", "Neutro", "NK")
        est_val_gg <- cbind(melt(t(est_val)), melt(truth)$value)
        colnames(est_val_gg) <- c("celltype", "sample_id", "est", "truth")
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
    title <- cowplot::ggdraw() + cowplot::draw_label(tumor_content[j], fontface = "bold")
    p_tmp <- cowplot::plot_grid(plotlist = p_no_legend, nrow = 3)
    p_concat[[j]] <- cowplot::plot_grid(title, p_tmp, legend, nrow = 1, rel_widths = c(0.8,6))
  }
  p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
  file_name <- paste0(path, "sim3_sampleScatter_", "real", "_", weight_type[w], "_", "tpm", ".pdf")
  ggsave(filename = file_name, p_concat_concat, width = 10, height = 15)
}

#--- summarized by comp number, heatmap, cor ---#
for(u in 1:length(unit)){
  for(w in 1:length(weight_type)){
    p_concat <- list() 
    for(s in 1:length(mix_type)){
      p_list <- list() 
      for(j in 1:length(tumor_content)){
        sum_name <- paste(weight_type[w], "sim3", "sum", "cellcor", mix_type[s], tumor_content[j], unit[u], sep = "_")
        sum_list <- get(sum_name)
        sum_mat <- matrix(0, nrow = length(methods), ncol = length(n_comp))
        for(i in 1:length(n_comp)){
          sum_mat[,i] <- sum_list[[i]][,'all']
        }
        rownames(sum_mat) <- methods
        colnames(sum_mat) <- paste("Comp", n_comp)
        
        
        p_list[[j]] <- pheatmap(sum_mat, cluster_cols = FALSE, cluster_rows = FALSE,
                                display_numbers = T, main = tumor_content[j],
                                breaks = seq(0, 1, length = 100), angle_col = 45,
                                number_color = "black")[[4]]
        
      }
      title <- cowplot::ggdraw() + cowplot::draw_label(mix_type[s], fontface = "bold")
      p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
      p_concat[[s]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
    }
    title <- cowplot::ggdraw() + cowplot::draw_label(weight_type[w], fontface = "bold")
    p_tmp <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
    p_concat_concat <- cowplot::plot_grid(title, p_tmp, ncol = 1, rel_heights = c(0.5, 6))
    file_name = paste0(path,"sim3_heatmap_cellcor_all_", weight_type[w], "_", unit[u], ".pdf")
    ggsave(p_concat_concat, filename = file_name ,width = 15,height = 10)
  }
}

#--- summarized by comp number, heatmap, mad ---#
for(u in 1:length(unit)){
  for(w in 1:length(weight_type)){
    p_concat <- list() 
    for(s in 1:length(mix_type)){
      p_list <- list() 
      for(j in 1:length(tumor_content)){
        sum_name <- paste(weight_type[w], "sim3", "sum", "cellmad", mix_type[s], tumor_content[j], unit[u], sep = "_")
        sum_list <- get(sum_name)
        sum_mat <- matrix(0, nrow = length(methods), ncol = length(n_comp))
        for(i in 1:length(n_comp)){
          sum_mat[,i] <- sum_list[[i]][,'all']
        }
        rownames(sum_mat) <- methods
        colnames(sum_mat) <- paste("Comp", n_comp)
        
        
        p_list[[j]] <- pheatmap(sum_mat, cluster_cols = FALSE, cluster_rows = FALSE,
                                colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(100),
                                display_numbers = T, main = tumor_content[j],
                                angle_col = 45,
                                number_color = "black")[[4]]
        
      }
      title <- cowplot::ggdraw() + cowplot::draw_label(mix_type[s], fontface = "bold")
      p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
      p_concat[[s]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8,6))
    }
    title <- cowplot::ggdraw() + cowplot::draw_label(weight_type[w], fontface = "bold")
    p_tmp <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
    p_concat_concat <- cowplot::plot_grid(title, p_tmp, ncol = 1, rel_heights = c(0.5, 6))
    file_name = paste0(path,"sim3_heatmap_cellmad_all_", weight_type[w], "_", unit[u], ".pdf")
    ggsave(p_concat_concat, filename = file_name ,width = 15,height = 10)
  }
}

save.image("./output/Sim3/sim3_evaluation_submission.RData")