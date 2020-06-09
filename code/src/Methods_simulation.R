# 0: General functions and methods 
#@ marker gene filter for each gene 
markerGene.filter <- function(ref_rowExpr, s_sample_id, max_val, min_val, ns_max, ns_min){
  #@ ref_rowExpr: gene expression across all samples 
  #@ s_sample_id: sample ids for the targeted cell type 
  #@ max_val: threshold for expression (targeted cell type)
  #@ min_val: threshold for expression (other cell types) 
  #@ ns_max: number of samples pass the test (targeted cell type)
  #@ ns_min: number of samples pass the test (other cell types)
  #@ output: TRUE if it is a marker gene (otherwise FALSE)
  if(sum(ref_rowExpr[s_sample_id] >= max_val) >= ns_max && sum(ref_rowExpr[-s_sample_id] <= min_val) >= ns_min){ # use && for length comparison (first element)
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

#@ extract marker genes from reference 
#@ sim1, sim2 and sim3 
setMethod(f = "marker_gene_from_ref",
          signature = c("sim_params", "marker_params", "ANY", "character"),
          definition = function(sim_param, marker_param, ref, celltype){
            #@ sim_param: sim_params ob 
            #@ marker_param: marker_params ob 
            #@ ref: reference profile from create_ref() 
            #@ output: marker gene list 
            #@         marker_set[[1]]: marker gene in the table format 
            #@         marker_set[[2]]: marker gene in the list format 
            #@         marker_set[[3]]: sample number threshold 
            max_val <- marker_param@max_val 
            min_val <- marker_param@min_val
            p_initial <- marker_param@p_initial
            step <- marker_param@step 
            gene_id <- sim_param@gene_id
            
            marker <- data.frame()
            p_val <- vector()
            # for each cellType k, extract marker genes 
            for(k in 1:length(celltype)){
              pattern = paste0(celltype[k],"_")
              s_sample_id <- grep(pattern = pattern, colnames(ref)) # extract sample ids for the targeted cell type 
              marker_k <- vector()
              p <- p_initial
              while(length(marker_k) < 2){ # iterate over p until we have more than 2 marker genes
                ns_max <- round(p*length(s_sample_id), 0)   # therehold for number of samples pass the test (target cell type)
                ns_min <- round(p*(ncol(ref) - length(s_sample_id)), 0) # threshold for number of samples pass the test (other cell types)
                marker_k_filter <- apply(ref,1,markerGene.filter,s_sample_id, max_val, min_val, ns_max, ns_min) 
                # marker gene for targeted cell 
                marker_k <- gene_id[marker_k_filter]
                # you need to lower the threshold if less than 2 marker genes are identified 
                p <- p - step
              }
              p_val[k] <- p + step
              # combine all marker genes from all cell types, C1 - marker gene, C2 - cell type identity 
              marker <- rbind(marker,cbind(markers = marker_k, celltype = rep(celltype[k],length(marker_k))))
            }
            marker <- mapply(marker,FUN = as.character) # transform every element to the character class 
            # eliminate duplicated genes, matrix version of marker gene can be used as the input of DSA
            dup_genes <- marker[duplicated(marker[,1]),1]
            marker <- marker[which(!marker[,1] %in% dup_genes),]
            
            # convert it to the list (input of EstimateWeight function)
            marker_list <- list()
            for( i in 1 : length(celltype)){
              marker_list[[i]] <- marker[marker[,2] == celltype[i],1]
            }
            names(marker_list) <- celltype
            # wrap marker gene in both matrix version and list version
            marker_set <- list()
            marker_set[[1]] <- marker
            marker_set[[2]] <- marker_list
            marker_set[[3]] <- p_val
            return(marker_set)
          })


#@ extract signature gene from matrix (union of DE genes)
#@ sim1, sim2 and sim3 
setMethod(f = "sig_gene_from_ref",
          signature = c("sim_params", "ANY", "numeric", "numeric", "data.frame"),
          definition = function(sim_param, ref, alpha, nFold, meta){
            #@ sim_param: sim_params object
            #@ ref: reference for signature gene extraction 
            #@ alpha: adj p-val threshold for sigGene filter 
            #@ nFold: log fold change threshold for sigGene filter 
            #@ output: sigGene from reference 
            library(DESeq2)
            
            # round ref to the discrete data 
            ref <- round(ref, digits = 0)
            # Use DESeq2 extract overexpressed genes with adjusted p value < 0.01 
            sample_name <- colnames(ref)
            
            if(all(rownames(meta) == sample_name)){
              dds <- DESeqDataSetFromMatrix(countData = ref,
                                            colData = meta,
                                            design = ~ celltype)
            }
            
            # DE analysis 
            dds <- DESeq(dds)
            res <- results(dds)
            
            group <- combn(unique(meta$celltype),2)
            de_genes <- vector() 
            for( g in 1:ncol(group)){
              # filter out genes that have adj.p value larger than 0.01
              res_tmp <- data.frame(results(dds, contrast = c("celltype", as.character(group[,g]))))
              res_tmp <- na.omit(res_tmp)
              filter <- res_tmp$padj <= alpha & abs(res_tmp$log2FoldChange) >= nFold
              if(is.na(sum(filter))){
                print(paste0(" Constract group ", group[,g], " doesn't contain genes that passed the DE criteria!"))
              }
              else{
                de_genes <- c(de_genes ,rownames(res_tmp[filter,]))
              }
            }
            return(unique(de_genes))
          })

#@ set parameters for nb simulation 
setMethod(f = "set_nb_params",
          signature = "nb_params",
          definition = function(nb_param, nsim_val, bv_val, libsize_val, sigma_val){
            #@ nb_param: nb_params object 
            #@ nsim_val: number of simulations (number of samples in one simulation dataset)
            #@ bv_val: bv parameter for nb simulation 
            #@ libsize_val: library size for simulation
            #@ sigma_val: sigma parameter for nb simulation 
            #@ nb_param: nb_params for nb simulations 
            nb_param@nsim = nsim_val
            nb_param@bv = bv_val
            nb_param@libsize = libsize_val
            nb_param@sigma = sigma_val
            return(nb_param)
          }
)
#@ generate one set of nb simulation 
setMethod(f = "sim_gammaPoisson",
          signature = c("nb_params", "matrix"),
          definition = function(nb_param, reference){
            #@ nb_param: nb_params object from set_nb_params() 
            #@ reference: reference matrix to extract proportions of genomic features 
            #@ output: one set of nb simulation (one matrix, row - gene, col - sim_sample) 
            library(edgeR)
            # Define the paramters for the simulation 
            # 1. prop - Calculate proportions of genomics features in the library 
            prop <- est.prop(reference)
            ngenes <- length(prop)
            # 2. mu0 - expected count value (same for all reference files)
            mu0 <- matrix(prop, ngenes,1) %*% matrix(nb_param@libsize, 1, nb_param@nsim)
            # 3. bcv0 - set the biological variation 
            bcv0 <- nb_param@bv +  1/sqrt(mu0) # source of the perturbation grid in the simulation 
            # 4. bcv - set the variance to the variance 
            # use log - normal to add additional perturbation on the bcv0 -> this simulation is derived from Lu & Voom paper 
            # * Alternatively you can use chi-square distribution to determine bcv0 -> this simulation is derived from Voom & Splatter paper: code --- df.BCV = 40 , bcv <- bcv0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV)) --- 
            bcv <- bcv0*exp(rnorm(ngenes, mean = 0, sd = nb_param@sigma)/2) 
            # 5. shape & scale parameters for gamma sampling 
            # perform gamma sampling to give variance for the mu0 
            shape <- 1/bcv^2
            scale <- mu0/shape
            mu <- matrix(rgamma(ngenes*nb_param@nsim, shape = shape, scale = scale), ngenes, nb_param@nsim)
            mu[is.na(mu)] <- 0 
            # 6. Poisson sampling 
            sim_count <- matrix(rpois(ngenes*nb_param@nsim, lambda = mu), ngenes, nb_param@nsim)
            rownames(sim_count) <- names(prop)
            colnames(sim_count) <- paste0(nb_param@prefix, 1:nb_param@nsim)
            return(sim_count)
          }
)

# 1. Methods for sim1 
#@ create reference
#@ sim1 
setMethod(f = "create_ref",
          signature = c("sim1_params", "character"),
          definition = function(sim_param, unit_val){
            #@ sim_param: sim_params object 
            #@ unit_val: targeted unit value 
            #@ output: reference profile with all targeted cell type profiles, 1 - ref, 2 - ref_mdian 
            celltype <- sim_param@celltype
            dataset_name <- sim_param@dataset_name
            
            ref <- list()
            ref_median <- list() 
            
            for(i in 1:length(dataset_name)){
              sample_name <- vector()
              tmp_ref <- list() 
              tmp_ref_median <- list() 
              for(k in 1:length(celltype)){
                profile_name_input <- paste(celltype[k], unit_val, sep = "_")
                profile <- get(profile_name_input)[[i]]
                tmp_ref[[k]] <- profile 
                sample_name <- c(sample_name, paste0("D", i, "_", celltype[k],"_", 1:ncol(profile)))
                tmp_ref_median[[k]] <- apply(profile, 1, median)
              }
              ref[[i]] <- do.call(cbind, tmp_ref)
              ref_median[[i]] <- do.call(cbind, tmp_ref_median)
              colnames(ref[[i]]) <- sample_name
              colnames(ref_median[[i]]) <- celltype
            }
            return(list(ref, ref_median))
          }
)

#@ create cibersort specific ref anno
#@ sim1
setMethod(f = "ref_anno_cibersort",
          signature = c("sim1_params", "list"),
          definition = function(sim_param, ref_list){
            #@ sim_param: sim_params ob 
            #@ ref_list: reference profile created by create_ref() 
            #@ output: cibersort specific ref anno
            dataset_name <- sim_param@dataset_name
            celltype <- sim_param@celltype
            cibersort_anno <- list()
            for(i in 1:length(dataset_name)){
              dat_tmp <- ref_list[[i]]
              
              cibersort_anno[[i]] <- matrix(2, nrow = length(celltype), ncol = ncol(dat_tmp))
              for(k in 1:length(celltype)){
                select_id <- grep(pattern = paste0("_", celltype[k],"_"), colnames(dat_tmp))
                cibersort_anno[[i]][k, select_id] <- 1
              }
              rownames(cibersort_anno[[i]]) <- celltype
              colnames(cibersort_anno[[i]]) <- colnames(dat_tmp)
            }
            return(cibersort_anno)  
          }
)

#@ create ref anno 
#@ sim1 
setMethod(f = "ref_anno_others",
          signature = c("sim1_params", "list"),
          definition = function(sim_param, ref_list){
            #@ sim_param: sim_params_ob 
            #@ ref_list: reference profile created by create_ref() 
            #@ output: ref anno for non-cibersort methods (i.e. timer)
            others_anno <- list()
            dataset_name <- sim_param@dataset_name
            for(i in 1:length(dataset_name)){
              ref <- ref_list[[i]]
              others_anno[[i]] <- cbind(colnames(ref), strsplit2(colnames(ref), "_")[,2])
              colnames(others_anno[[i]]) <- c("sample", "celltype")
              rownames(others_anno[[i]]) <- colnames(ref)
              others_anno[[i]] <- data.frame(others_anno[[i]])
            }
            return(others_anno)
          }
)

#@ create reference for epic algorithm 
#@ sim1
setMethod(f = "create_epic_ref",
          signature = c("sim1_params", "list","list","list"),
          definition = function(sim_param, sig_gene, ref, median_ref){
            #@ sim_param: sim_params ob 
            #@ sig_gene: signature gene from sig_gene_from_ref()
            #@ ref: reference expression - for IQR extraction 
            #@ median_ref: median reference expression (of each cell type) 
            #@ output: reference for epic 
            #@         refProfiles - median ref
            #@         refProfiles.var = interquartile range of gene (across profiles from the same cell type)
            #@         sigGenes = signature gene 
            gene_id <- sim_param@gene_id
            celltype <- sim_param@celltype
            dataset_name <- sim_param@dataset_name
            ref_IQR <- list()
            epic_ref <- list()
            for(i in 1:length(dataset_name)){
              ref_IQR[[i]] <- matrix(NA, nrow = length(gene_id), ncol = length(celltype))
              for(k in 1:length(celltype)){
                sample_selector <- grep(pattern = celltype[k], colnames(ref[[i]]))
                ref_IQR[[i]][,k] <- apply(ref[[i]][, sample_selector], 1, IQR)
              }
              rownames(ref_IQR[[i]]) <- rownames(ref[[i]]) 
              colnames(ref_IQR[[i]]) <- celltype
              epic_ref[[i]] <- list(refProfiles = median_ref[[i]], refProfiles.var = ref_IQR[[i]], sigGenes = sig_gene[[i]])
            }
            return(epic_ref)
          }
)

#@ generate nb simulations 
#@ sim1
setMethod(f = "wrap_gammaPoisson",
          signature = c("sim1_params", "nb_params", "list", "numeric"),
          definition = function(sim_param, nb_param, expr_list, interval){
            #@ sim_param: sim_params object 
            #@ nb_param: nb_params object 
            #@ expr_list: expression matrix from 3 datasets 
            #@ interval: scale the noise level 
            #@ output: generate simulations of 3 datasets - for each dataset, 10 simulation matrix for 10 gradients of noise 
            sim_expr <- list() 
            for(i in 1:length(expr_list)){
              sim_expr_tmp <- list()
              grid <- sim_param@grid
              nb_param@prefix = paste0("D",i,"_",celltype,"_") # set the smaple name of simulations 
              
              for(j in 1:length(grid)){
                nb_param@bv = 0.1*j*interval
                sim_expr_tmp[[j]] <- sim_gammaPoisson(nb_param, round(expr_list[[i]],0))
              }
              sim_expr[[i]] <- sim_expr_tmp
            }
            return(sim_expr)
          }
)

#@ generate log normal noise 
#@ sim1  
setMethod(f = "sim_lognormal_noise",
          signature = c("list", "sim1_params", "numeric"),
          definition = function(clean_mix_list, sim_param, sigma){
          #@ clean_mix_list: weighted sum of purified profiles 
          #@ sim_param: sim_params object 
          #@ sigma: preset base variance 
          #@ output: list of simulations derived from 3 datasets, each with 10 levels of log-normal noise 
          grid <- sim_param@grid
          
          noise_mix_list <- list()
          for(i in 1:length(clean_mix_list)){
            mix_tmp <- list() 
            for(j in 1:length(grid)){
              mix_tmp[[j]] <- apply(clean_mix_list[[i]], 2, function(x) x+2^rnorm(length(x),0,grid[j]*sigma))
            }
            noise_mix_list[[i]] <- mix_tmp
          }
          return(noise_mix_list)
          }
          )

#@ generate normal nosie 
#@ sim1 
setMethod(f = "sim_normal_noise",
          signature = c("list", "sim1_params", "numeric"),
          definition = function(clean_mix_list, sim_param, sigma){
            #@ clean_mix_list: weighted sum of purified profiles 
            #@ sim_param: sim_params object 
            #@ sigma: preset base variance 
            #@ output: list of simulations derived from 3 datasets, each with 10 levels of normal noise 
            grid <- sim_param@grid
            
            noise_mix_list <- list()
            for(i in 1:length(clean_mix_list)){
              mix_tmp <- list() 
              for(j in 1:length(grid)){
                mix_tmp[[j]] <- apply(clean_mix_list[[i]], 2, function(x) return(2^(log2(x+1)+rnorm(length(x),0,grid[j]*sigma))))
              }
              noise_mix_list[[i]] <- mix_tmp
            }
            return(noise_mix_list)
          }
          )

#@ generate simulated mixture (3 simulation models: nb, log-normal, normal)
#@ sim1 
setMethod(f = "sim_mix",
          signature = c("sim1_params", "ANY"),
          definition = function(sim_param, W, ref_list, model_name){
            #@ sim_param: sim_params object 
            #@ W: pre-defined weights for mixture (ground-truth)
            #@ ref_list: reference used to generate simulations 
            #@ model_name: targeted simulation model 
            #@ output: simulated mixtures generated from targeted simulation model, 3 datasets, 10 perturbation, each perturbation dataset with 20 profiles 
          dataset_name <- sim_param@dataset_name
          grid <- sim_param@grid
          nsim <- sim_param@nsim
          gene_id <- sim_param@gene_id
          
          mix_list <- list()
          if(model_name == "nb"){
            for(i in 1:length(dataset_name)){
              mix_list_tmp <- list()
              for(j in 1:length(grid)){
                mix_list_tmp[[j]] <- matrix(0, nrow = length(gene_id), ncol = nsim)
                  for( s in 1:nsim){
                    S_tmp <- cbind( T = get(paste("T", "count", "nb", sep = "_"))[[i]][[j]][,s],
                                    B = get(paste("B", "count", "nb", sep = "_"))[[i]][[j]][,s],
                                    Mono = get(paste("Mono", "count", "nb", sep = "_"))[[i]][[j]][,s])
                    mix_list_tmp[[j]][,s] <- as.matrix(S_tmp) %*% (W[s,])
                  } # s
                rownames(mix_list_tmp[[j]]) <- gene_id
              } # j
              mix_list[[i]] <- mix_list_tmp
            } # i 
            return(mix_list)
          }
          else if(model_name == "lognormal" | model_name == "normal"){
            # create clean mixtures 
            for(i in 1:length(dataset_name)){
              mix_list[[i]] <- as.matrix(ref_list[[i]]) %*% as.matrix(t(W))
              rownames(mix_list[[i]]) <- gene_id
            }
            # add noise term to the mixture 
            if(model_name == "lognormal"){
              noise_mix_list <- sim_lognormal_noise(clean_mix_list = mix_list,
                                                     sim_param = sim1_params_ob,
                                                     sigma = 10)
            }else if(model_name == "normal"){
              noise_mix_list <- sim_normal_noise(clean_mix_list = mix_list,
                                                  sim_param = sim1_params_ob,
                                                  sigma = 10)
            }
            return(noise_mix_list)
          }
          
          }
)

#@ generate nb simulated mixtures with disprepant library sizes 
#@ sim1 
setMethod(f = "sim_libsize_mix",
          # this returns the mixtures in the count unit 
          signature = c("sim1_params", "character","character", "matrix"),
          definition = function(sim_param, suffix_lib1, suffix_lib2, W){
            #@ sim_param: sim_params object 
            #@ suffix_lib1: suffix(except cell type name) of cell-type specific reference object of lib1 
            #@ suffix_lib2: suffix(except cell type name) of cell-type specific reference object of lib2
            #@ W: pre-defined weights for mixture (ground-truth)
            #@ output: simulated mixtures with discrepant library sizes
            gene_id <- sim_param@gene_id
            nsim <- sim_param@nsim 
            nDataset <- length(sim_param@dataset_name)
            nGrid <- length(sim_param@grid)
            
            mix_count <- list()
            for(i in 1:nDataset){
              count_tmp <- list() 
              for(j in 1:nGrid){
                count_tmp[[j]] <- matrix(0, nrow = length(gene_id), ncol = nsim)
                for(s in 1:nsim){
                  
                  if(s <= 10){
                    S_tmp <- cbind( T = get(paste("T", suffix_lib1, sep = "_"))[[i]][[j]][,s],
                                    B = get(paste("B", suffix_lib1, sep = "_"))[[i]][[j]][,s],
                                    Mono = get(paste("Mono", suffix_lib1, sep = "_"))[[i]][[j]][,s])
                    count_tmp[[j]][,s] <- as.matrix(S_tmp) %*% as.matrix((W[s,]))
                  }else{
                    S_tmp <- cbind( T = get(paste("T", suffix_lib2, sep = "_"))[[i]][[j]][,s-10],
                                    B = get(paste("B", suffix_lib2, sep = "_"))[[i]][[j]][,s-10],
                                    Mono = get(paste("Mono", suffix_lib2, sep = "_"))[[i]][[j]][,s-10])
                    count_tmp[[j]][,s] <- as.matrix(S_tmp) %*% as.matrix((W[s,]))
                  } # if-else
                  
                } # s
                rownames(count_tmp[[j]]) <- gene_id
              } # j 
              mix_count[[i]] <- count_tmp
            } # i 
            return(mix_count)
          }
)

#@ unit transformation from count unit 
#@ sim1 
setMethod(f = "sim_mix_unit_transform",
          signature = c("sim1_params", "list", "character", "list"),
          definition = function(sim_param, mix_list, unit_val, eff_length){
            #@ sim_param: sim_params obejct 
            #@ mix_list: simulated mixture  lists 
            #@ unit_val: targeted output unit
            #@ eff_length: effective length of gene for tpm transformation 
            #@ output: mixtures with transformed targeted unit
            library(edgeR)
            library(scater)
            nDataset <- length(sim_param@dataset_name)
            mix_transform <- list()
            for(i in 1:nDataset){
              if(unit_val == "countNorm"){
                mix_transform[[i]] <- lapply(mix_list[[i]], function(x) apply(x, MARGIN = 2, median_libSize_norm, median(colSums(x))))
              }else if(unit_val == "cpm"){
                mix_transform[[i]] <- lapply(mix_list[[i]], function(x) edgeR::cpm(x))
              }
              else if(unit_val == "tpm"){
                mix_transform[[i]] <- lapply(mix_list[[i]], function(x) calculateTPM(x, effective_length = eff_length[[i]][, 'effective_length']))
              }
            } #i 
            return(mix_transform)
          }
)

#2. Methods for sim2 
#@ create reference
#@ sim2 and sim3 
setMethod(f = "create_ref",
          signature = c("comp_sim_params", "character"),
          definition = function(sim_param, unit_val, celltype_list, sim_name){
            #@ sim_param: sim_params object 
            #@ unit_val: targeted unit value 
            #@ celltype_list: a list of cell types to be added to the reference list (5 comp - 10 comp)
            #@ output: reference profile with all targeted cell type profiles, 1 - ref, 2 - ref_mdian 
            ref <- list() 
            ref_median <- list() 
            
            for( i in 1:length(celltype_list)){
              ref[[i]] <- data.frame(NA)
              ref_median[[i]] <- data.frame(NA)
              for(j in 1:length(celltype_list[[i]])){
                celltype_expr <- get(paste(sim_name,celltype_list[[i]][j], unit_val, sep = "_"))
                colnames(celltype_expr) <- paste(celltype_list[[i]][j], 1:ncol(celltype_expr),sep = "_") 
                median_expr <- apply(celltype_expr,1,median)
                
                ref[[i]] <- cbind(ref[[i]], celltype_expr)
                ref_median[[i]] <- cbind(ref_median[[i]], median_expr)
              }
              # eliminate the first column (NA values due to initialization)
              ref[[i]] <- ref[[i]][,-1]
              ref_median[[i]] <- ref_median[[i]][,-1]
              # set colnames
              colnames(ref_median[[i]]) <- celltype_list[[i]]
            }
            return(list(ref, ref_median))
          })

#@ create cibersort specific ref anno
#@ sim2 and sim3 
setMethod( f = "ref_anno_cibersort",
           # for sim2 
           signature = c("comp_sim_params", "list"),
           definition = function(sim_param, ref_list){
             #@ sim_param: sim_params ob 
             #@ ref_list: reference profile created by create_ref() 
             #@ output: cibersort specific ref anno
             library(limma)
             cibersort_anno <- list()
             for(i in 1:length(ref_list)){
               ref <- ref_list[[i]]
               celltype <- unique(strsplit2(colnames(ref), "_")[,1])
               cibersort_anno[[i]] <- matrix(2, nrow = length(celltype), ncol = ncol(ref))
               for(k in 1:length(celltype)){
                 pattern <- celltype[k]
                 select_id <- grep(pattern = paste0(pattern,"_"), colnames(ref))
                 cibersort_anno[[i]][k, select_id] <- 1
               }
               rownames(cibersort_anno[[i]]) <- celltype
               colnames(cibersort_anno[[i]]) <- colnames(ref)
             }
             return(cibersort_anno)
           }
)

#@ create ref anno 
#@ sim2 and sim3 
setMethod(f = "ref_anno_others",
          signature = c("comp_sim_params", "list"),
          definition = function(sim_param, ref_list){
            #@ sim_param: sim_params_ob 
            #@ ref_list: reference profile created by create_ref() 
            #@ output: ref anno for non-cibersort methods (i.e. timer)
            library(limma)
            others_anno <- list()
            for(i in 1:length(ref_list)){
              ref <- ref_list[[i]]
              others_anno[[i]] <- cbind(colnames(ref), strsplit2(colnames(ref), "_")[,1])
              colnames(others_anno[[i]]) <- c("sample", "celltype")
              rownames(others_anno[[i]]) <- colnames(ref)
              others_anno[[i]] <- data.frame(others_anno[[i]])
            }
            return(others_anno)
          }
)




#@ create reference for epic algorithm 
#@ sim2 and sim3
setMethod(f = "create_epic_ref",
          signature = c("comp_sim_params", "list","list","list"),
          definition = function(sim_param, sig_gene, ref, median_ref){
            #@ sim_param: sim_params ob 
            #@ sig_gene: signature gene from sig_gene_from_ref()
            #@ ref: reference expression - for IQR extraction 
            #@ median_ref: median reference expression (of each cell type) 
            #@ output: reference for epic 
            #@         refProfiles - median ref
            #@         refProfiles.var = interquartile range of gene (across profiles from the same cell type)
            #@         sigGenes = signature gene 
            gene_id <- sim_param@gene_id
            
            ref_IQR <- list()
            epic_ref <- list()
            for(i in 1:length(ref)){
              celltype <- colnames(median_ref[[i]])
              ref_IQR[[i]] <- matrix(NA, nrow = length(gene_id), ncol = length(celltype))
              
              for(k in 1:length(celltype)){
                sample_selector <- grep(pattern = celltype[k], colnames(ref[[i]]))
                ref_IQR[[i]][,k] <- apply(ref[[i]][, sample_selector], 1, IQR)
              }
              rownames(ref_IQR[[i]]) <- rownames(ref[[i]]) 
              colnames(ref_IQR[[i]]) <- celltype
              epic_ref[[i]] <- list(refProfiles = median_ref[[i]], refProfiles.var = ref_IQR[[i]], sigGenes = sig_gene[[i]])
            }
            return(epic_ref)
          }
)

#@ generate simulated mixture (6 sets of mixtures with a diverse component number)
#@ sim2 and sim3 
setMethod(f = "sim_mix",
          signature = c("comp_sim_params", "ANY"),
          definition = function(sim_param, W, sim_celltype, sim_name){
            #@ sim_param: sim_params object 
            #@ W: pre-defined weights for mixture (ground-truth)
            #@ ref_list: reference used to generate simulations 
            #@ model_name: targeted simulation model 
            #@ output: simulated mixtures with component number from 5 to 10. 
            nsim <- sim_param@nsim
            gene_id <- sim_param@gene_id
            
            mix_list <- list()
            for(i in 1:length(sim_celltype)){
              mix_list[[i]] <- matrix(0, nrow = length(gene_id), ncol = nsim)
              for(s in 1:nsim){
                S_tmp <- matrix(0, nrow = length(gene_id), ncol = length(sim_celltype[[i]]))
                for(k in 1:length(sim_celltype[[i]])){
                  S_tmp[, k] <- get(paste(sim_name, sim_celltype[[i]][k], "nb", "count", sep = "_"))[,s]
                }
                mix_list[[i]][ ,s] <- as.matrix(S_tmp %*% W[[i]][, s])
              }
              rownames(mix_list[[i]]) <- gene_id
            }
            
            
            return(mix_list)
          }
)




#@ unit transformation from count unit 
#@ sim2 and sim3
setMethod(f = "sim_mix_unit_transform",
          signature = c("comp_sim_params", "list", "character", "data.frame"),
          definition = function(sim_param, mix_list, unit_val, eff_length){
            #@ sim_param: sim_params obejct 
            #@ mix_list: simulated mixture  lists 
            #@ unit_val: targeted output unit
            #@ eff_length: effective length of gene for tpm transformation 
            #@ output: mixtures with transformed targeted unit
            library(edgeR)
            library(scater)
            mix_transform <- list()
            if(unit_val == "countNorm"){
              mix_transform <- lapply(mix_list, function(x) apply(x, MARGIN = 2, median_libSize_norm, median(colSums(x))))
            }else if(unit_val == "cpm"){
              mix_transform <- lapply(mix_list, function(x) edgeR::cpm(x))
            }else if(unit_val == "tpm"){
              mix_transform <- lapply(mix_list, function(x) calculateTPM(x, effective_length = eff_length[,'effective_length']))
            }
            return(mix_transform)
          }
          )









setMethod( f = "ref_anno_cibersort_oneMat",
           # for sim2 
           signature = c("ANY", "character"),
           definition = function(ref, celltype_identity){
             celltype <- unique(celltype_identity)
             cibersort_anno <- matrix(2, nrow = length(celltype), ncol = ncol(ref))
            
             for(k in 1:length(celltype)){
               pattern <- celltype[k]
               select_id <- grep(pattern = paste0(pattern, "_"), colnames(ref))
               cibersort_anno[k, select_id] <- 1
             }
             rownames(cibersort_anno) <- celltype
             colnames(cibersort_anno) <- colnames(ref)
             
             return(cibersort_anno)
           }
)







# oneMat 
setMethod(f = "create_epic_ref_oneMat",
          signature = c("character","ANY","ANY"),
          definition = function(sig_gene, ref, median_ref){
            
            celltype <- colnames(median_ref)
            ref_IQR <- matrix(NA, nrow = nrow(median_ref), ncol = length(celltype))
            
            for(k in 1:length(celltype)){
              sample_selector <- grep(pattern = celltype[k], colnames(ref))
              ref_IQR[,k] <- apply(ref[, sample_selector], 1, IQR)
              }
            rownames(ref_IQR) <- rownames(ref) 
            colnames(ref_IQR) <- celltype
            epic_ref <- list(refProfiles = median_ref, refProfiles.var = ref_IQR, sigGenes = sig_gene)
        
            return(epic_ref)
          }
)








setMethod(f = "marker_gene_from_ref_oneMat", 
          signature = c("marker_params", "ANY", "character"),
          definition = function(marker_param, ref, celltype){
            max_val <- marker_param@max_val
            min_val <- marker_param@min_val
            
            p_initial <- marker_param@p_initial
            step <- marker_param@step
            
            gene_id <- rownames(ref)
            marker <- data.frame()
            p_val <- vector() 
            # for each celltype k, extract marker genes 
            for(k in 1:length(celltype)){
              pattern = paste0(celltype[k], "_")
              s_sample_id <- grep(pattern = pattern, colnames(ref))
              marker_k <- vector()
              p <- p_initial
              while(length(marker_k) < 2){
                ns_max <- round(p*length(s_sample_id), 0)
                ns_min <- round(p*(ncol(ref) - length(s_sample_id)), 0)
                
                marker_k_filter <- apply(ref, 1, markerGene.filter, s_sample_id, max_val, min_val, ns_max, ns_min)
                
                marker_k <- gene_id[marker_k_filter]
                p <- p-step
              }
              p_val[k] <- p
              # combine all marker genes from all cell types 
              marker <- rbind(marker, cbind(markers = marker_k, celltype = rep(celltype[k], length(marker_k))))
            }
            
            marker <- mapply(marker, FUN = as.character)
            # eliminate duplication 
            marker <- marker[!duplicated(marker[,1]), ]
            # convert the marker from the matrix form to the list form  
            marker_list <- list() 
            for( k in 1:length(celltype)){
              marker_list[[k]] <- marker[marker[,2] == celltype[k], 1]
            }
            names(marker_list) <- celltype 
            
            # wrap marker gene in both matrix version and list version 
            marker_set <- list() 
            marker_set[[1]] <- marker 
            marker_set[[2]] <- marker_list 
            marker_set[[3]] <- p_val + step 
            names(marker_set) <- c("matrix_ver", "list_ver", "sample_p")
            return(marker_set)
          }
)  

setMethod(f = "marker_gene_from_refMean_oneMat",
          signature = c("marker_params", "ANY", "character"),
          definition = function(marker_param, ref_mean, celltype){
            
            gene_id <- rownames(ref_mean)
            max_val <- marker_param@max_val
            min_val <- marker_param@min_val
            
            marker <- data.frame()
            for(k in 1:length(celltype)){
              s_sample_id <- k 
              marker_k <- vector()
              marker_k_filter <- apply(ref_mean, 1, markerGene_mean.filter, s_sample_id, max_val, min_val)
              marker_k <- gene_id[marker_k_filter]
              marker <- rbind(marker, cbind(markers = marker_k, celltype = rep(celltype[k], length(marker_k))))
            }
            
            marker <- mapply(marker, FUN = as.character)
            # eliminate duplication 
            marker <- marker[!duplicated(marker[,1]), ]
            # convert the marker from the matrix form to the list form  
            marker_list <- list() 
            for( k in 1:length(celltype)){
              marker_list[[k]] <- marker[marker[,2] == celltype[k], 1]
            }
            names(marker_list) <- celltype 
            
            # wrap marker gene in both matrix version and list version 
            marker_set <- list() 
            marker_set[[1]] <- marker 
            marker_set[[2]] <- marker_list 
            names(marker_set) <- c("matrix_ver", "list_ver")
            return(marker_set)
          }
)  
          






setMethod(f = "sig_gene_from_ref_oneMat",
          signature = c("ANY", "character", "character", "numeric", "numeric"),
          definition = function(ref, celltype, identity, alpha, nFold){
            library(DESeq2)
            
            # Use DESeq2 extract overexpressed genes with adjusted p value < 0.01 
            sample_name <- colnames(ref)
            meta <- data.frame(sample = sample_name, 
                               celltype = identity)
            rownames(meta) <- meta[,1] # set column 1 as the rownames 
            
            if(all(rownames(meta) == sample_name)){
              dds <- DESeqDataSetFromMatrix(countData = ref,
                                            colData = meta,
                                            design = ~ celltype)
            }
            
            # DE analysis 
            dds <- DESeq(dds)
            res <- results(dds)
            
            group <- combn(celltype,2)
            de_genes <- vector() 
            for( g in 1:ncol(group)){
              # filter out genes that have adj.p value larger than 0.01
              res_tmp <- data.frame(results(dds, contrast = c("celltype", group[,g])))
              res_tmp <- na.omit(res_tmp)
              filter <- res_tmp$padj <= alpha & abs(res_tmp$log2FoldChange) >= nFold
              if(is.na(sum(filter))){
                print(paste0(" Constract group ", group[,g], " doesn't contain genes that passed the DE criteria!"))
              }
              else{
                de_genes <- c(de_genes ,rownames(res_tmp[filter,]))
              }
            }
            return(unique(de_genes))
          })


setMethod(f = "tumor_weight",
          signature = c("list", "character"),
          definition = function(known_weight, tumorContent = c("small", "large", "mosaic")){
            # generate weight of the tumor 
            if(tumorContent == "small"){
              t_prop <- runif(20, min = 0.2, max = 0.3)
            }
            else if(tumorContent == "large"){
              t_prop <- runif(20, min = 0.7, max = 0.8)
            }
            else if(tumorContent == "mosaic"){
              t_prop <- runif(20, min = 0.05, max = 0.95)
            }
            weight <- list()
            for(i in 1:length(known_weight)){
              # add weight of the tumor to the predefined immune weight matrix 
              tmp <- rbind(known_weight[[i]], HCT = t_prop)
              idx <- nrow(tmp) # get the index of HCT 
              # row normalization so that the mixing ratio of all ccomponents sums to 1 for each sample 
              weight[[i]] <- apply(tmp, 2, function(x) c(x[-idx]*(1-x[idx]), x[idx]))
            }
            return(weight)
          })

setMethod(f = "orthog_weight",
          signature = c("comp_sim_params", "numeric", "list"),
          definition = function(sim_param, N, sim_celltype)
          { 
            n_comp <- sim_param@n_comp
            nsim <- sim_param@nsim
            
            orthog_w_list <- list()
            for(i in 1:length(n_comp)){
              orthog_w_list[[i]] <- w_min_condition(nsim, n_comp[i], N = N)
              rownames(orthog_w_list[[i]]) <- sim_celltype[[i]]
            }
            return(orthog_w_list)
          }
)

w_min_condition <- function(nsim, n_comp_val, N){
  library(plyr)
  w_list <- list()
  kappa_value <- vector() 
  
  for(n in 1:N){
    # create initial probability based on uniform distribution 
    w <- laply(as.list(rep(20, n_comp_val)), runif)
    # normalize probability so that it satisfies sum-to-1 rule 
    w_norm <- apply(w, MARGIN = 2, function(x) x/sum(x))
    # save the weight 
    w_list[[n]] <- w_norm
    # save the condition number 
    kappa_value[n] <- kappa(w_norm)
  }
  orthog_w <- w_list[[which.min(kappa_value)]]
  
  return(orthog_w)
}

setMethod(f = "real_weight",
          signature = c("comp_sim_params", "numeric", "list"),
          definition = function(sim_param, N, sim_prop_mat)
          { 
            n_comp <- sim_param@n_comp
            nsim <- sim_param@nsim
            
            real_w_list <- list() 
            for(i in 1:length(n_comp)){
              real_w_list[[i]] <- w_real_weight(nsim, n_comp[i], N = N, sim_prop_mat[[i]])
              while(w_real_test(real_w_list[[i]], sim_prop_mat[[i]])){
                print("Warning! Simulated proportions doesn't satisfy predefined proportion range!")
                print("Re-simulate weight set", i)
                real_w_list[[i]] <- w_real_weight(nsim, n_comp[i], N = N, sim_prop_mat[[i]])
              }
              rownames(real_w_list[[i]]) <- rownames(sim_prop_mat[[i]])
            }
            return(real_w_list)
          }
)

w_real_weight <- function(nsim, n_comp_val, N, sim_prop_mat_val){
  w_sample <- matrix(0, nrow = n_comp_val, ncol = N)
  
  for(k in 1:n_comp_val){
    w_sample[k, ] <- runif(N, min = sim_prop_mat_val[k, 1], max = sim_prop_mat_val[k, 2])
  }
  # rank the weights according to the colSums, find the one that is closest to 100
  rank_p <- order(abs(colSums(w_sample) - 100))
  # pick the first 20 sampled weights and normalize weights to make it a valid probability 
  w_real <- apply(w_sample[, rank_p[1:20]], 2, function(x) x/sum(x))
  
  return(w_real)
}

w_real_test <- function(w_real, sim_prop_mat_val){
  for(k in 1:nrow(w_real)){
    if(range(w_real)[1] < sim_prop_mat_val[k, 1]/100 | range(w_real)[2] > sim_prop_mat_val[k, 2]/100){
      return(FALSE)
    }else{
      return(TRUE)
    }
  }
}

