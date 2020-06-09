# 0: General functions and methods 
scatter_refSim <- function(x,y,logTrans=TRUE, xlab = "Profile 1", ylab = "Profile 2" ){
  library(LSD)
  library(scales)
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  cor_xy <- round(cor(x,y,method = "spearman"),digits=2)
  euclidean_xy <- scientific_format(digits = 1)(round(dist(rbind(x,y)),digits=0))
  if(logTrans){
    heatscatter(x, y, xlab = xlab,ylab= ylab,colpal = cols, pch=20,cex.lab=1.3,main=paste0("r = ",cor_xy," , ","d = ",euclidean_xy), log = "xy")
    abline(a=0,b=1,col="red")
  }
  else{
    heatscatter(x,y, xlab=x_index,ylab=y_index,colpal = cols, pch=20,cex.lab=1.3,main=paste0("r = ",cor_xy," , ","d = ",euclidean_xy))
    abline(a=0,b=1,col="red")
    
  }
}
# 
cv_density <- function(cv_mat){
  library(reshape2)
  gg_cv <- melt(cv_mat)
  colnames(gg_cv) <- c("gene","sample","cv")
  p <- ggplot(gg_cv, aes(x = cv)) + geom_density(alpha = 0.7, fill = "red")
  p <- p + facet_grid(. ~ sample) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title = element_text(size = 15), 
          text = element_text(size = 15))
  return(p)
}

##################################################################################################################
#-------------------------------- sim1 related methods ------------------------------------------------------#
##################################################################################################################
#@ evaluation metrics 
#@ sim1 
setMethod(f = "sim_eval_metrics",
          signature = c("sim1_params", "list", "character"),
          definition = function(sim_param, dat_list, eval_type = c("cross_perturb_cor","within_perturb_cor", "cv", "mean", "var")){
            #@ sim_param: sim1_params ob 
            #@ dat_list: list of simulated expressions 
            #@ eval_type: cross_perturb_cor - spearman's correlation of cross perturb simulations 
            #@            within_perturb_cor - spearman's correlation of inter perturb simulations 
            #@            cv - coefficient of variation, mean - mean, var - variance 
            #@ output: eval_result with selected evaluation metrics 
            nsim <- sim_param@nsim
            gene_id <- sim_param@gene_id
            nDataset <- length(sim_param@dataset_name)
            grid <- sim_param@grid
            nGrid <- length(grid)
            
            if(eval_type == "cross_perturb_cor"){
              eval_result <- list() 
              for(i in 1:nDataset){
                tmp <- list() 
                for(s in 1:nsim){
                  tmp[[s]] <- matrix(NA, nrow = nGrid, ncol = nGrid)
                  for(j in 1:nGrid){
                    for(k in 1:nGrid){
                      tmp[[s]][j, k] <- cor(dat_list[[i]][[j]][, s], dat_list[[i]][[k]][, s], method = "spearman")
                    }
                  }
                }
                # average out the cross_perturb_cor across simulation profiles
                tmp_arr <- array(unlist(tmp), dim = c(nGrid,nGrid,nsim))
                eval_result[[i]] <- apply(tmp_arr, 1:2, mean, na.rm = TRUE)
                rownames(eval_result[[i]]) <- paste0("D", i, "P", 1:nGrid)
                colnames(eval_result[[i]]) <- paste0("D", i, "P", 1:nGrid)
              }
            }
            else if(eval_type == "within_perturb_cor"){
              eval_result <- matrix(0, nrow = nDataset, ncol = nGrid)
              for(i in 1:nDataset){
                for(j in 1:nGrid){
                  pair_id <- combn(c(1:nsim), 2)
                  cor_val <- vector()
                  for(p in 1:ncol(pair_id)){
                    pair <- pair_id[,p]
                    cor_val[p] <- cor(dat_list[[i]][[j]][, pair[1]], dat_list[[i]][[j]][, pair[2]], method = "spearman")
                  }
                  eval_result[i,j] <- sum(cor_val)/length(cor_val)
                  rownames(eval_result) <- paste0("D", 1:nDataset)
                  colnames(eval_result) <- paste0("P", 1:nGrid)
                }
              }
            }
            else if(eval_type == "cv" | eval_type == "mean" | eval_type == "var"){
              func_option <- function(x){
                switch(x,
                       cv = "cv.calc",
                       mean = "mean",
                       var = "var",
                       stop("Invalid eval_type value"))
              }
              func_name <- func_option(eval_type)
              eval_result <- list() 
              for(i in 1:nDataset){
                eval_result[[i]] <- matrix(NA, nrow = length(gene_id), ncol = nGrid)
                for(j in 1:nGrid){
                  eval_result[[i]][, j] <- apply(dat_list[[i]][[j]], 1, func_name)
                }
                colnames(eval_result[[i]]) <- paste0("D", i, "P", 1:nGrid)
                rownames(eval_result[[i]]) <- gene_id
              }
            
            }
            return(eval_result)
          }
)

#@ draw the scatter plots of one pair of expression profiles 
#@ sim1
setMethod(f = "sample_scatter",
          signature = c("sim1_params", "list", "character"),
          definition = function(sim_param, mix_list, scatter_type = c("cross_perturb", "within_perturb")){
            #@ sim_param: sim1_params object 
            #@ mix_list: simulated mixtures 
            #@ scatter_type: scatter plot type, all samples will be derived from D1 simulation 
            #@               cross_perturb: sample pairs are simulation 1 in P1 vs. simulation 1 in P1-10 
            #@               within_perturb: sample pairs are simulation 1 vs. simulation 2 in the same perturbation  level
            #@ output: scatter plot of selected scatter_type 
            library(LSD)
            cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
            nDataset <- length(sim_param@dataset_name)
            nGrid <- length(sim_param@grid)
            # for each row 
            if(scatter_type == "cross_perturb"){
              for(j in 1:nGrid){
                x <- mix_list[[1]][[1]][,1]
                y <- mix_list[[1]][[j]][,1]
                scatter_refSim(x, y, logTrans = TRUE, xlab = "P1", ylab = paste0("P",j))
              }
            }
            else if(scatter_type == "within_perturb"){
              for(j in 1:nGrid){
                x <- mix_list[[1]][[j]][,1]
                y <- mix_list[[1]][[j]][,2]
                scatter_refSim(x, y, logTrans = TRUE, xlab = paste0("P",j, "S1"), ylab = paste0("P",j, "S2"))
              }
            }
            
          }
)

#@ draw mean-variance scatterplots 
#@ sim1
setMethod(f = "meanVar_scatter",
          signature = c("sim1_params", "matrix", "matrix"),
          definition = function(sim_param, mean_mat, var_mat){
            #@ sim_param: sim_params object 
            #@ mean_mat: (gene) mean for each simulation set 
            #@ var_mat: (gene) variance for each simulation set 
            #@ output: mean-variance scatter plot 
            library(LSD)
            cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
            nGrid <- length(sim_param@grid)
            # for each row 
            for( i in 1:nGrid){
              x <- mean_mat[, i]
              y <- var_mat[, i]
              heatscatter(log2(x), log2(y), xlab = "log2(mean)", ylab = "log2(var)", colpal = cols, pch = 20, cex.lab = 1.3, main = paste0("P",i))
              abline(a=0,b=1,col="red")
            }
          }
)

#@ expr heatmap 
#@ sim1 
setMethod(f = "expr_heatmap",
          signature = c("sim1_params", "ANY", "vector"),
          definition = function(sim_param, ref, targeted_gene, main){
            #@ sim_param: sim_params object 
            #@ ref: reference expression profile for the heatmap 
            #@ targeted_gene: marker gene or signature gene list 
            #@ main: main title for the figure 
            #@ output: heatmap for the marker gene expression 
            library(LSD)
            library(limma)
            
            expr <- log2(ref[targeted_gene, ] + 1)
            celltype_anno <- strsplit2(colnames(ref), "_")[,2]
            celltype <- sim_param@celltype
            n_celltype <- length(celltype)
            
            dataset_anno <- data.frame(celltype = celltype_anno,
                                       row.names = colnames(ref))
            
            col2 <- c(brewer.pal(n_celltype, "Set3"))
            names(col2) <- celltype
            anno_col <- list(celltype = col2)
            
            p <- pheatmap(expr, 
                          annotation_col = dataset_anno,
                          annotation_colors = anno_col,
                          cluster_col = FALSE, 
                          cluster_row = FALSE, 
                          show_rownames = FALSE, 
                          show_colnames = FALSE,
                          main = main)[[4]]
            return(p)
          }
)



setMethod(f = "deconv_plot_sampleScatter",
          signature = c("sim1_params", "list", "matrix"),
          definition = function(sim_param, est_W, truth, noise_index){
            
            library(plyr)
            library(cowplot)
            
            nMarker <- nDataset <- length(sim_param@dataset_name)
            nGrid <- length(sim_param@grid)
            nCelltype <- length(sim_param@celltype)
            celltype <- sim_param@celltype
            
            p_list <- list()
            for(i in 1:length(noise_index)){
              tmp <- est_W[[1]][[1]][[noise_index[i]]]
              if(is.logical(tmp)){
                est_val <- NA
                est_val_gg <- data.frame(NA)
                p_list[[i]] <- ggplot(data = est_val_gg) + geom_abline(slope = 1, colour = "red") + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 15), text = element_text(size = 15)) +
                ggtitle(paste0("Noise level ", noise_index[i]))
              } 
              else{ 
                if(nrow(tmp) == nCelltype){
                  est_val <- as.matrix(t(tmp))
                }
                else if(ncol(tmp) == nCelltype){
                  est_val <- as.matrix(tmp)
                }
                colnames(est_val) <- celltype 
                est_val_gg <- cbind(melt(est_val), melt(truth)$value)
                colnames(est_val_gg) <- c("sample_id", "celltype", "est", "truth")
                p_list[[i]] <- ggplot(data = est_val_gg, aes(x = truth, y = est, group = celltype)) + 
                  geom_point(aes(color = celltype)) + geom_abline(slope = 1, colour = "red") + 
                  scale_x_continuous(limits=c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
                  scale_y_continuous(limits=c(min(est_val_gg$est), max(est_val_gg$est)), breaks = seq(from = round(min(est_val_gg$est), 1), to = round(max(est_val_gg$est), 1), length.out = 6)) +
                  scale_color_brewer(palette="Set1") + theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 15), text = element_text(size = 15)) +
                  ggtitle(paste0("Noise level ", noise_index[i]))
              }
              }
            return(p_list)
          })


##################################################################################################################
#-------------------------------- comp_sim (sim2 and sim3) related methods ------------------------------------------------------#
##################################################################################################################
#@ evaluation metrics 
#@ sim2 and sim3 
setMethod(f = "sim_eval_metrics",
          signature = c("comp_sim_params", "list", "character"),
          definition = function(sim_param, dat_list, eval_type = c("cor", "cv", "mean", "var")){
            #@ sim_param: sim1_params ob 
            #@ dat_list: list of simulated expressions 
            #@ eval_type: cross_perturb_cor - spearman's correlation of cross perturb simulations 
            #@            within_perturb_cor - spearman's correlation of inter perturb simulations 
            #@            cv - coefficient of variation, mean - mean, var - variance 
            #@ output: eval_result with selected evaluation metrics 
            nsim <- sim_param@nsim
            n_comp <- sim_param@n_comp
            gene_id <- sim_param@gene_id
            
            if(eval_type == "cross_perturb_cor"){
              cor_mat <- list()
              for(s in 1:nsim){
                cor_mat[[s]] <- matrix(NA, nrow = length(n_comp), ncol = length(n_comp))
                for(j in 1:length(n_comp)){
                  for(k in 1:length(n_comp)){
                    cor_mat[[s]][j, k] <- cor(dat_list[[j]][, s], dat_list[[k]][, s], method = "spearman")
                    }
                }
                rownames(cor_mat[[s]]) <- paste("Comp", n_comp, sep = "_")
                colnames(cor_mat[[s]]) <- paste("Comp", n_comp, sep = "_")
              }
              eval_result <- rowMeans(do.call(cbind, cor_mat), na.rm = TRUE)
            }
            else if(eval_type == "within_perturb_cor"){
              eval_result <- vector()
              for(i in 1:length(n_comp)){
                pair_id <- combn(c(1:nsim),2)
                cor_val <- vector()
                for(p in 1:ncol(pair_id)){
                  pair <- pair_id[, p]
                  cor_val[p] <- cor(dat_list[[i]][,pair[1]], dat_list[[i]][,pair[2]], method = "spearman")
                }
                eval_result[i] <- sum(cor_val)/length(cor_val)
              }
              names(eval_result) <- paste("Comp", n_comp)
              }
            else if(eval_type == "cv" | eval_type == "mean" | eval_type == "var"){
              func_option <- function(x){
                switch(x,
                       cv = "cv.calc",
                       mean = "mean",
                       var = "var",
                       stop("Invalid eval_type value"))
              }
              func_name <- func_option(eval_type)
              eval_result <- matrix(NA, nrow = length(gene_id), ncol = length(dat_list))
              for(i in 1:length(dat_list)){
                eval_result[,i] <- apply(dat_list[[i]], 1, func_name)
              }
              rownames(eval_result) <- gene_id
              colnames(eval_result) <- paste(n_comp, "Comp Mix")
            }
            return(eval_result)
          }
)

#@ draw the scatter plots of one pair of expression profiles 
#@ sim2 and sim3 
setMethod(f = "sample_scatter",
          signature = c("comp_sim_params", "list", "character"),
          definition = function(sim_param, mix_list, scatter_type = c("cross_perturb", "within_perturb")){
            #@ sim_param: sim1_params object 
            #@ mix_list: simulated mixtures 
            #@ scatter_type: scatter plot type, all samples will be derived from D1 simulation 
            #@               cross_perturb: sample pairs are simulation 1 in P1 vs. simulation 1 in P1-10 
            #@               within_perturb: sample pairs are simulation 1 vs. simulation 2 in the same perturbation  level
            #@ output: scatter plot of selected scatter_type 
            library(LSD)
            cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
            n_comp <- sim_param@n_comp
        
            if(scatter_type == "cross_perturb"){
              for(i in 1:length(n_comp)){
                x <- mix_list[[1]][,1]
                y <- mix_list[[i]][,1]
                scatter_refSim(x, y, logTrans = TRUE, xlab = "Comp 1", ylab = paste0("Comp",i))
              }
            }
            else if(scatter_type == "within_perturb"){
              for(i in 1:length(n_comp)){
                x <- mix_list[[i]][,1]
                y <- mix_list[[i]][,2]
                scatter_refSim(x, y, logTrans = TRUE, xlab = paste0("Comp",i, "S1"), ylab = paste0("Comp",i, "S2"))
              }
            }
            
          }
)

#@ draw mean-variance scatterplots 
#@ sim2 and sim3
setMethod(f = "meanVar_scatter",
          signature = c("comp_sim_params", "matrix", "matrix"),
          definition = function(sim_param, mean_mat, var_mat){
            #@ sim_param: sim_params object 
            #@ mean_mat: (gene) mean for each simulation set 
            #@ var_mat: (gene) variance for each simulation set 
            #@ output: mean-variance scatter plot 
            library(LSD)
            cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
            n_comp <- sim_param@n_comp
            # for each row 
            for( i in 1:length(n_comp)){
              x <- mean_mat[, i]
              y <- var_mat[, i]
              heatscatter(log2(x), log2(y), xlab = "log2(mean)", ylab = "log2(var)", colpal = cols, pch = 20, cex.lab = 1.3, main = paste(n_comp[i], "Comp Mix"))
              abline(a=0,b=1,col="red")
            }
          }
)

#@ draw expr heatmap 
#@ sim2 and sim3 
setMethod(f = "expr_heatmap",
          signature = c("comp_sim_params", "ANY", "vector"),
          definition = function(sim_param, ref, targeted_gene, main){
            #@ sim_param: sim_params object 
            #@ ref: reference expression profile for the heatmap 
            #@ targeted_gene: marker gene or signature gene list 
            #@ main: main title for the figure 
            #@ output: heatmap for the marker gene expression 
            library(LSD)

            expr <- log2(ref[targeted_gene, ] + 1)
            celltype_anno <- strsplit2(colnames(ref), "_")[,1]
            celltype <- names(sim_param@extract_pattern_immune)
            n_celltype <- length(celltype)
            
            dataset_anno <- data.frame(celltype = celltype_anno,
                                       row.names = colnames(ref))
            
            col2 <- c(brewer.pal(n_celltype, "Set3"))
            names(col2) <- celltype
            anno_col <- list(celltype = col2)
            
            p <- pheatmap(expr, 
                          annotation_col = dataset_anno,
                          annotation_colors = anno_col,
                          cluster_col = FALSE, 
                          cluster_row = FALSE, 
                          show_rownames = FALSE, 
                          show_colnames = FALSE,
                          main = main, fontsize = 7)[[4]]
            return(p)
          }
)
#@ evaluation of deconvolution results 
#@ sim1 
setMethod(f = "deconv_evaluation",
          signature = c("character", "sim1_params", "list", "matrix"),
          definition = function(eval_type = c("cor", "mad", "cellcor", "cellmad"), sim_param, est_W, truth){
            #@ eval_type: evalution type 
            #@ sim_param: sim_params_ob 
            #@ est_W: estimated weights 
            #@ truth: ground truth
            #@ output: evaluation metric  
            library(plyr)
            
            nMarker <- nDataset <- length(sim_param@dataset_name)
            nGrid <- length(sim_param@grid)
            nCelltype <- length(sim_param@celltype)
            celltype <- sim_param@celltype
            
            if(eval_type == "cellcor"){
              w_eval <- list() 
              for(m in 1:nMarker){
                tmp <- list()
                for(j in 1:nDataset){
                  tmp[[j]] <- laply(est_W[[m]][[j]], calc.cor_celltype, truth, 3)
                  colnames(tmp[[j]]) <- celltype
                }
                names(tmp) <- paste0("D", 1:3)
               w_eval[[m]] <- tmp
              }
              names(w_eval) <- paste0("M", 1:3)
            }
            else if(eval_type == "cellmad"){
              w_eval <- list() 
              for(m in 1:nMarker){
                tmp <- list()
                for(j in 1:nDataset){
                  tmp[[j]] <- laply(est_W[[m]][[j]], calc.mad_celltype, truth, 3)
                  colnames(tmp[[j]]) <- celltype
                }
                names(tmp) <- paste0("D", 1:3)
                w_eval[[m]] <- tmp
              }
              names(w_eval) <- paste0("M", 1:3)
            }
            else if(eval_type == "cor"){
              w_eval <- list() 
              for(m in 1:nMarker){
                tmp <- matrix(NA, nrow = nDataset, ncol = nGrid)
                for(j in 1:nDataset){
                  tmp[j, ] <- laply(est_W[[m]][[j]], calc.cor, truth)
                }
                rownames(tmp) <- paste0("D",1:3)
                w_eval[[m]] <- tmp
              }
              names(w_eval) <- paste0("M", 1:3)
            }
            else if(eval_type == "mad"){
              w_eval <- list() 
              for(m in 1:nMarker){
                tmp <- matrix(NA, nrow = nDataset, ncol = nGrid)
                for(j in 1:nDataset){
                  tmp[j, ] <- laply(est_W[[m]][[j]], calc.mad, truth)
                }
                rownames(tmp) <- paste0("D",1:3)
                w_eval[[m]] <- tmp
              }
              names(w_eval) <- paste0("M", 1:3)
            }
            
            return(w_eval)
          }
)

#@ summary of evaluation derive from deconv_evaluation
#@ sim1 
setMethod(f = "deconv_eval_summary",
          signature = c("character", "character", "sim1_params", "list"),
          definition = function(eval_type = c("cor", "mad", "cellcor", "cellmad"), robustness = c("all", "high_noise", "ref_noise"), sim_param, eval_metric){
            #@ eval_type: evaluation type of input 
            #@ robustness: summary type 
            #@             all - all datasets 
            #@             high_noise - noise level 6 - 10 
            #@             ref_noise - datasets with different combination of ref and source 
            #@ sim_param: sim1_params_ob 
            #@ eval_metric: evaluation metric (output of deconv_evaluation)
            #@ output: summary of deconv evaluation metric 
            library(plyr)
            
            nMarker <- nDataset <- length(sim_param@dataset_name)
            nGrid <- length(sim_param@grid)
            nCelltype <- length(sim_param@celltype)
            celltype <- sim_param@celltype
            
            if(eval_type %in% c("cellcor", "cellmad")){
              if(robustness == "all"){
                tmp_all <- do.call(cbind, lapply(eval_metric, function(x) {do.call(cbind, x)}))
                sum_eval <- vector()
                
                for(i in 1:nCelltype){
                  tmp_celltype <- as.numeric(tmp_all[,grep(celltype[i], colnames(tmp_all))])
                  sum_eval[i] <- mean(tmp_celltype)
                }
                
                tmp_all <- as.numeric(tmp_all)
                
                sum_eval[nCelltype + 1] <- mean(tmp_all)
                names(sum_eval) <- c(celltype, "all")
              }
              else if(robustness == "high_noise"){
                tmp_all <- do.call(cbind, lapply(eval_metric, function(x) {do.call(cbind, x)}))
                sum_eval <- vector()
                for(i in 1:nCelltype){
                  tmp_celltype <- as.numeric(tmp_all[6:10,grep(celltype[i], colnames(tmp_all))])
                  sum_eval[i] <- sum(tmp_celltype, na.rm = TRUE)/length(tmp_celltype)
                }
                
                tmp_all <- as.numeric(tmp_all[6:10,])
                sum_eval[4] <- sum(tmp_all, na.rm = TRUE)/length(tmp_all)
                names(sum_eval) <- c(celltype, "all")
              }
              else if(robustness == "ref_noise"){
                eval_metric$M1$D1 <- eval_metric$M2$D2 <- eval_metric$M3$D3 <- NULL 
                tmp_all <- do.call(cbind, lapply(eval_metric, function(x) {do.call(cbind, x)}))
                sum_eval <- vector()
                for(i in 1:nCelltype){
                  tmp_celltype <- as.numeric(tmp_all[,grep(celltype[i], colnames(tmp_all))])
                  sum_eval[i] <- sum(tmp_celltype, na.rm = TRUE)/length(tmp_celltype)
                }
                
                tmp_all <- as.numeric(tmp_all)
                sum_eval[4] <- sum(tmp_all, na.rm = TRUE)/length(tmp_all)
                names(sum_eval) <- c(celltype, "all")
              }
              return(sum_eval)
            }
            else if(eval_type %in% c("cor", "mad")){
              if(robustness == "all"){
                tmp <- unlist(eval_metric)
            #    tmp[is.na(tmp)] <- 0 
                sum_eval <- sum(tmp, na.rm = TRUE)/length(tmp)
              }
              else if(robustness == "high_noise"){
                tmp <- lapply(eval_metric, function(x) return(x[,5:10]))
                tmp <- unlist(tmp)
              #  tmp[is.na(tmp)] <- 0 
                sum_eval <- sum(tmp, na.rm = TRUE)/length(tmp)
              }
              else if(robustness == "ref_noise"){
                tmp <- list()
                for(i in 1:length(eval_metric)){
                  tmp[[i]] <- eval_metric[[i]][-i,]
                }
                tmp <- unlist(tmp)
              #  tmp[is.na(tmp)] <- 0 
                sum_eval <- sum(tmp, na.rm = TRUE)/length(tmp)
              }
            return(sum_eval)
            }
          }
)

#@ wrap function for deconv_eval_summary w.r.t methods 
#@ sim1
setMethod(f = "wrap_eval_summary",
          signature = c("sim1_params", "character", "character", "character"),
          definition = function(sim_param, eval_type = c("cor", "mad", "cellcor", "cellmad"), all_methods, unit_val, sim_model_val){
            #@ sim_param: sim1_params_ob 
            #@ eval_type: evaluation type (same with eval_type in deconv_eval_summary)
            #@ all_methods: all of the deconvolution methods tested 
            #@ sim_model_val: sim_model value 
            #@ unit_val: unit value 
            #@ output: aggregated results deconv_eval_summary w.r.t all_methods
            library(plyr)
            celltype <- sim_param@celltype
            all_robustness <- c("all", "high_noise", "ref_noise")
            if(eval_type %in% c("cellcor", "cellmad")){
              sum_list <- list()
              for(i in 1:length(all_robustness)){
                sum_list[[i]] <- matrix(data = NA, nrow = length(all_methods), ncol = 4)
                for(m in 1:length(all_methods)){
                  eval_metric_name <- paste("sim1", paste0(eval_type, "List"), all_methods[m], sim_model_val, unit_val, sep = "_")
                  eval_metric <- get(eval_metric_name)
                  sum_list[[i]][m, ] <- deconv_eval_summary(eval_type = eval_type, robustness = all_robustness[i], sim1_params_ob, eval_metric)
                }
                colnames(sum_list[[i]]) <- c(celltype, "all")
                rownames(sum_list[[i]]) <- all_methods
              }
              names(sum_list) <- all_robustness
              
              return(sum_list)
            }
            else if(eval_type %in% c('cor', "mad")){
              
              sum_mat <- matrix(data = NA, nrow = length(all_methods), ncol = 3)
              for(m in 1:length(all_methods)){
                eval_metric_name <- paste("sim1", paste0(eval_type, "List"), all_methods[m], sim_model_val, unit_val, sep = "_")
                eval_metric <- get(eval_metric_name)
                for(i in 1:length(all_robustness)){
                  sum_mat[m, i] <- deconv_eval_summary(eval_type = eval_type, robustness = all_robustness[i], sim1_params_ob, eval_metric)
                }
              }
              colnames(sum_mat) <- all_robustness
              rownames(sum_mat) <- all_methods
              return(sum_mat)
            }
          }
)

# sim2 related functions 
setMethod(f = "deconv_evaluation",
          signature = c("character", "sim2_params", "list", "list"),
          definition = function(eval_type = c("cor", "mad", "cellcor", "cellmad"), sim_param, est_W, truth){
            #@ eval_type: evalution type 
            #@ sim_param: sim_params_ob 
            #@ est_W: estimated weights 
            #@ truth: ground truth
            #@ output: evaluation metric  
            library(plyr)
            n_comp <- sim_param@n_comp 
            
            if(eval_type == "cellcor"){
              w_eval <- list() 
              for(i in 1:length(est_W)){
                w_eval[[i]] <- calc.cor_celltype(est_W[[i]], truth[[i]], n_comp[i])
              }
            }
            else if(eval_type == "cellmad"){
              w_eval <- list() 
              for(i in 1:length(est_W)){
                w_eval[[i]] <- calc.mad_celltype(est_W[[i]], truth[[i]], n_comp[i])
              }
            }
            else if(eval_type == "cor"){
              w_eval <- vector()
              for(i in 1:length(est_W)){
                w_eval[i] <- calc.cor(est_W[[i]], truth[[i]])
              }
            }
            else if(eval_type == "mad"){
              w_eval <- vector()
              for(i in 1:length(est_W)){
                w_eval[i] <- calc.mad(est_W[[i]], truth[[i]])
              }
            }
            return(w_eval)
          }
)

#@ wrap function for deconv_eval_summary w.r.t methods 
#@ sim2
setMethod(f = "wrap_eval_summary",
          signature = c("sim2_params", "character", "character", "character"),
          definition = function(sim_param, eval_type = c("cor", "mad", "cellcor", "cellmad"), all_methods, unit_val, mix_type_val){
            #@ sim_param: sim1_params_ob 
            #@ eval_type: evaluation type (same with eval_type in deconv_eval_summary)
            #@ all_methods: all of the deconvolution methods tested 
            #@ sim_model_val: sim_model value 
            #@ unit_val: unit value 
            #@ output: aggregated results deconv_eval_summary w.r.t all_methods
            library(plyr)
            n_comp <- sim_param@n_comp 
            sum_list <- list() 
            for(i in 1:length(n_comp)){
              sum_list[[i]] <- matrix(data = NA, nrow = length(all_methods), ncol = n_comp[i] + 1)
              for(m in 1:length(all_methods)){
                eval_metric_name <- paste("sim2", paste0(eval_type, "List"), all_methods[m], mix_type_val, unit_val, sep = "_")
                eval_metric <- get(eval_metric_name)[[i]]
                sum_list[[i]][m,1:n_comp[i]] = eval_metric
               # eval_metric[is.na(eval_metric)] <- 0
                sum_list[[i]][m, n_comp[i] + 1] = sum(eval_metric)/length(eval_metric)
              }
              colnames(sum_list[[i]]) <- c(names(eval_metric), "all")
              rownames(sum_list[[i]]) <- all_methods
            }
            return(sum_list)
}
)

# draw sampled scatter 
# sim2 
setMethod(f = "deconv_plot_sampleScatter",
          signature = c("sim2_params", "list", "list"),
          definition = function(sim_param, est_W, truth){
            
            library(plyr)
            library(cowplot)
            
            n_comp <- sim_param@n_comp
            line_color_manual <- c("#999999", "#FF9933", "#56B4E9", "#009E73", "#CC00FF", "#0072B2", "#993300", "#99CC33", "#000000", "#FF0000")
            
            
            p_list <- list()
            for(i in 1:length(n_comp)){
              tmp <- est_W[[i]]
              celltype <- rownames(truth[[i]])
              if(is.logical(tmp)){
                est_val <- NA 
                est_val_gg <- data.frame(NA)
                p_list[[i]] <- ggplot(data = est_val_gg) + geom_abline(slope = 1, colour = "red") + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 15), text = element_text(size = 15)) +
                  ggtitle(paste0("Comp ", n_comp[i]))
              }
              else{
                if(nrow(tmp) == n_comp[i]){
                  est_val <- as.matrix(t(tmp))
                }
                else if(ncol(tmp) == n_comp[i]){
                  est_val <- as.matrix(tmp)
                }
                colnames(est_val) <- celltype 
                est_val_gg <- cbind(melt(est_val), melt(t(truth[[i]]))$value)
                colnames(est_val_gg) <- c("sample_id", "celltype", "est", "truth")
                p_list[[i]] <- ggplot(data = est_val_gg, aes(x = truth, y = est, group = celltype)) + 
                  geom_point(aes(color = celltype)) + geom_abline(slope = 1, colour = "red") + 
                  scale_x_continuous(limits=c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
                  scale_y_continuous(limits=c(min(est_val_gg$est), max(est_val_gg$est)), breaks = seq(from = round(min(est_val_gg$est), 1), to = round(max(est_val_gg$est), 1), length.out = 6)) +
                  scale_color_manual(values = line_color_manual) + theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 15), text = element_text(size = 15)) +
                  ggtitle(paste0("Comp ", n_comp[i]))
                }
            }
            return(p_list)
          }
          )


setMethod(f = "deconv_evaluation",
          signature = c("character", "sim3_params", "list", "list"),
          definition = function(eval_type = c("cor", "mad", "cellcor"), sim_param, est_W, truth, weight_type){
            n_comp <- sim_param@n_comp 
            if(weight_type == "absolute"){
              tmp <- list()
              for(i in 1:length(est_W)){
                tmp[[i]] <- truth[[i]][-(n_comp[i]+1),]
              }
              truth <- tmp
            }

            if(eval_type == "cellcor"){
              wCor <- list() 
              for(i in 1:length(est_W)){
                wCor[[i]] <- calc.cor_celltype(est_W[[i]], truth[[i]], n_comp[i]) 
              }
            }
            else if(eval_type == "cellmad"){
              wCor <- list() 
              for(i in 1:length(est_W)){
                wCor[[i]] <- calc.mad_celltype(est_W[[i]], truth[[i]], n_comp[i]) 
              }
            }
            else if(eval_type == "cor"){
              wCor <- vector() 
              for(i in 1:length(est_W)){
                wCor[i] <- calc.cor(est_W[[i]], truth[[i]])
              }
            }
            else if(eval_type == "mad"){
              wCor <- vector() 
              for(i in 1:length(est_W)){
                wCor[i] <- calc.mad(est_W[[i]], truth[[i]])
              }
            }
            return(wCor)
          }
)


#@ wrap function for deconv_eval_summary w.r.t methods 
#@ sim3
setMethod(f = "wrap_eval_summary",
          signature = c("sim3_params", "character", "character", "character"),
          definition = function(sim_param, eval_type = c("cor", "mad", "cellcor", "cellmad"), all_methods, unit_val, mix_type_val, tumor_content_val, weight_type_val){
            #@ sim_param: sim1_params_ob 
            #@ eval_type: evaluation type (same with eval_type in deconv_eval_summary)
            #@ all_methods: all of the deconvolution methods tested 
            #@ sim_model_val: sim_model value 
            #@ unit_val: unit value 
            #@ output: aggregated results deconv_eval_summary w.r.t all_methods
            library(plyr)
            n_comp <- sim_param@n_comp 
            sum_list <- list() 
            for(i in 1:length(n_comp)){
              sum_list[[i]] <- matrix(data = NA, nrow = length(all_methods), ncol = n_comp[i] + 1)
              for(m in 1:length(all_methods)){
                eval_metric_name <- paste(weight_type_val, "sim3", paste0(eval_type, "List"), all_methods[m], mix_type_val, tumor_content_val, unit_val, sep = "_")
                eval_metric <- get(eval_metric_name)[[i]]
                sum_list[[i]][m,1:n_comp[i]] = eval_metric
                # indicate NA 
               # eval_metric[is.na(eval_metric)] <- 0
                sum_list[[i]][m, n_comp[i] + 1] = sum(eval_metric)/length(eval_metric)
              }
              colnames(sum_list[[i]]) <- c(names(eval_metric), "all")
              rownames(sum_list[[i]]) <- all_methods
            }
            return(sum_list)
          }
)

# draw sampled scatter 
# sim3
setMethod(f = "deconv_plot_sampleScatter",
          signature = c("sim3_params", "list", "list"),
          definition = function(sim_param, est_W, truth, weight_type){
            
            library(plyr)
            library(cowplot)
            
            n_comp <- sim_param@n_comp
            line_color_manual <- c("#999999", "#FF9933", "#56B4E9", "#009E73", "#CC00FF", "#0072B2", "#993300", "#99CC33", "#000000", "#FF0000")
            
            
            if(weight_type == "absolute"){
              tmp <- list()
              for(i in 1:length(est_W)){
                tmp[[i]] <- truth[[i]][-(n_comp[i]+1),]
              }
              truth <- tmp
            }
            
            p_list <- list()
            for(i in 1:length(n_comp)){
              tmp <- est_W[[i]]
              celltype <- rownames(truth[[i]])
              if(is.logical(tmp)){
                est_val <- NA 
                est_val_gg <- data.frame(NA)
                p_list[[i]] <- ggplot(data = est_val_gg) + geom_abline(slope = 1, colour = "red") + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 15), text = element_text(size = 15)) +
                  ggtitle(paste0("Comp ", n_comp[i]))
              }else{
                if(nrow(tmp) == n_comp[i]){
                  est_val <- as.matrix(t(tmp))
                }
                else if(ncol(tmp) == n_comp[i]){
                  est_val <- as.matrix(tmp)
                }
                colnames(est_val) <- celltype 
                est_val_gg <- cbind(melt(est_val), melt(t(truth[[i]]))$value)
                colnames(est_val_gg) <- c("sample_id", "celltype", "est", "truth")
                p_list[[i]] <- ggplot(data = est_val_gg, aes(x = truth, y = est, group = celltype)) + 
                  geom_point(aes(color = celltype)) + geom_abline(slope = 1, colour = "red") + 
                  scale_x_continuous(limits=c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
                  scale_y_continuous(limits=c(min(est_val_gg$est), max(est_val_gg$est)), breaks = seq(from = round(min(est_val_gg$est), 1), to = round(max(est_val_gg$est), 1), length.out = 6)) +
                  scale_color_manual(values = line_color_manual) + theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 15), text = element_text(size = 15)) +
                  ggtitle(paste0("Comp ", n_comp[i]))
              }
            }
            return(p_list)
          }
)


setMethod( f = "evaluation_mat",
            signature = c("list", "sim1_params"),
            definition = function(eval_list, sim_param){
              nGrid = length(sim_param@grid)
              nMarker = length(sim_param@dataset_name)
              eval_val <- c() 
              for(g in 1:nGrid){
                for(m in 1:nMarker){
                    eval_val <- c(eval_val, as.numeric(eval_list[[m]][,g]))
                }
              }
              return(eval_val)
            }
)

setMethod( f = "evaluation_mat_celltype",
           signature = c("list", "sim1_params", "character"),
           definition = function(eval_list, sim_param, celltype_val){
             nGrid = length(sim_param@grid)
             nMarker = length(sim_param@dataset_name)
             eval_val <- c() 
             for(g in 1:nGrid){
               for(m in 1:nMarker){
                 for(d in 1:nDataset){
                   eval_val <- c(eval_val, as.numeric(eval_list[[m]][[d]][g, celltype_val]))
                 }
               }
             }
             return(eval_val)
           }
)



