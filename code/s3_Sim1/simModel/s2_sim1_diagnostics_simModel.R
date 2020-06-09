load("./output/Sim1/simModel/sim1_mix_simModel.RData")
source("./code/src/Functions_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Methods_evaluation.R")
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(gridExtra)
path <- paste("./output/Sim1/simModel/")
unit <- sim1_params_ob@unit
sim_model <- sim1_params_ob@sim_model
gg_theme <- theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                               axis.title = element_text(size = 15), 
                               text = element_text(size = 15))
# Structure of sim1_M_nb_count 
# 2 layers of lists 
# first layer: derived datasets (length = 3)
# second layer: noise level (length = 10) 

# calculate eval metrics
for(i in 1:length(sim_model)){
  for( u in 1:length(unit)){
    mix_name <- paste("sim1", "M", sim_model[i], unit[u], sep = "_")
    cross_cor_name <- paste(mix_name, "cross_cor", sep = "_")
    within_cor_name <- paste(mix_name, "within_cor", sep = "_")
    mean_name <- paste(mix_name, "mean", sep = "_")
    var_name <- paste(mix_name, "var", sep = "_")
    cv_name <- paste(mix_name, "cv", sep = "_")
    
    assign(cross_cor_name, sim_eval_metrics(sim_param = sim1_params_ob,
                                      dat_list = get(mix_name),
                                      eval_type = "cross_perturb_cor"))
    assign(within_cor_name, sim_eval_metrics(sim_param = sim1_params_ob,
                                            dat_list = get(mix_name),
                                            eval_type = "within_perturb_cor"))
    assign(mean_name, sim_eval_metrics(sim_param = sim1_params_ob,
                                       dat_list = get(mix_name),
                                       eval_type = "mean"))
    assign(var_name, sim_eval_metrics(sim_param = sim1_params_ob,
                                      dat_list = get(mix_name),
                                      eval_type = "var"))
    assign(cv_name, sim_eval_metrics(sim_param = sim1_params_ob,
                                     dat_list = get(mix_name),
                                     eval_type = "cv"))
  }
}

# draw sampleScatter_within
for(u in 1:length(unit)){
  png(paste0(path, "sim1_simModel_sampleScatter_within_", unit[u], ".png"), width = 13, height = 18, units = "in", res = 300, type = "cairo")
  layout.matrix <- matrix(c(31, 1:5,31, 6:10, rep(34, 6), 
                            32, 11:15, 32, 16:20, rep(35,6),  
                            33, 21:25, 33, 26:30), nrow = 8, byrow = TRUE)
  layout(mat = layout.matrix, widths = c(1.7,2,2,2,2,2), heights =  c(2,2, 0.5, 2,2,0.5, 2,2))
  layout.show(35)
  for(s in 1:length(sim_model)){
    mix_name <- paste("sim1", "M", sim_model[s], unit[u], sep = "_")
    sample_scatter(sim_param = sim1_params_ob, mix_list = get(mix_name), scatter_type = "within_perturb")
    
  }
  for(s in 1:length(sim_model)){
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, sim_model[s], cex = 2, font = 2)
  }
  dev.off()
}

# figures to illustrate mean-variance relationships of simulated mixtures 
u = 1
for(i in 1:length(sim_model)){
    mix_name <- paste("sim1", "M", sim_model[i], unit[u], sep = "_")
#    pdf(paste0(path, sim_model[i], "_", unit[u], "_meanVar.pdf"), width = 13, height = 6.5)
    png(paste0(path, sim_model[i], "_", unit[u], "_meanVar.png"), width = 13, height = 6.5, units = "in", res = 300, type = "cairo")
    par(mfrow = c(2,5))
    mean_name <- paste(mix_name, "mean", sep = "_")
    var_name <- paste(mix_name, "var", sep = "_")
    # draw the mean-var plot for each unit  
    meanVar_scatter(sim_param = sim1_params_ob, mean_mat = get(mean_name)[[1]], var_mat = get(var_name)[[1]]) 
    dev.off()  
}

# heatmap_within_cor
p_concat <- list( )
for(s in 1:length(sim_model)){
  p_list <- list() 
  for(u in 1:length(unit)){
    mix_name <- paste("sim1", "M", sim_model[s], unit[u], sep = "_")
    within_cor_name <- paste(mix_name, "within_cor", sep = "_")
    p_list[[u]] <- pheatmap(get(within_cor_name), display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE,  main = unit[u], breaks = seq(0,1,length = 100), clustering_method = "single")[[4]] 
  }
  title <- cowplot::ggdraw() + cowplot::draw_label(sim_model[s], fontface = "bold")
  p_tmp <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
  p_concat[[s]] <- cowplot::plot_grid(title, p_tmp, nrow = 1, rel_widths = c(0.8, 6))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, ncol = 1)
ggsave(p_concat_concat, filename = paste0(path, "sim1_simModel_mix_cor_within.pdf"), width = 15, height = 8)

# cv density plot 
u = 1
p_concat <- list() 
for(s in 1:length(sim_model)){
  p_list <- list() 
  mix_name <- paste("sim1", "M", sim_model[s], unit[u], sep = "_")
  cv_name <- paste(mix_name, "cv", sep = "_")
  cv_mat <- get(cv_name)[[1]]
  for(i in 1:ncol(cv_mat)){
    gg_cv <- melt(cv_mat[,i])
    p_list[[i]] <- ggplot(gg_cv, aes(x = value)) + geom_density(alpha = 0.7, fill = "red") + xlab("cv") + ggtitle(paste0("P", i)) + theme_bw()
  }
  title <- cowplot::ggdraw() + cowplot::draw_label(sim_model[s], fontface = "bold")
  p_tmp <- cowplot::plot_grid(plotlist = p_list, ncol = 1)
  p_concat[[s]] <- cowplot::plot_grid(title, p_tmp, ncol = 1, rel_heights = c(0.5,10))
}
p_concat_concat <- cowplot::plot_grid(plotlist = p_concat, nrow = 1)
file_name = paste0(path,"sim1_simModel_cv_D1_", unit[u], ".pdf")
ggsave(p_concat_concat, filename = file_name ,width = 5,height=20)

# Draw marker gene and signature gene plots (remember to load ref list)
p_marker <- list()
p_sigGene <- list() 
for(i in 1:length(sim1_params_ob@dataset_name)){
  ref <- sim1_ref_countNorm[[i]]
  marker <- unlist(sim1_marker_countNorm[[i]][[2]])
  sigGene <- sim1_sigGene_count[[i]]
  
  p_marker[[i]] <- expr_heatmap(sim1_params_ob, ref, marker, main = paste0("D",i))
  p_sigGene[[i]] <- expr_heatmap(sim1_params_ob, ref, sigGene, main = paste0("D",i))
}
p_marker_concat <- grid.arrange(arrangeGrob(grobs = p_marker, nrow = 1))
p_filename <- paste0(path, "sim1_marker_heatmap.pdf")
ggsave(p_marker_concat, filename = p_filename, width = 15, height = 5)

p_sigGene_concat <- grid.arrange(arrangeGrob(grobs = p_sigGene, nrow = 1))
p_filename <- paste0(path, "sim1_sigGene_heatmap.pdf")
ggsave(p_sigGene_concat, filename = p_filename, width = 15, height = 5)

save.image(file = "./output/Sim1/simModel/sim1_diagnostics_simModel.RData")