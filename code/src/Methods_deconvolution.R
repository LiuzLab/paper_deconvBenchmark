##### I: General functions ########

library(limma)
# deconv_run is write for different groups of deconvolution algorithms that 
# require different numbers and types of inputs 
setMethod( f = "deconv_run_marker", 
           # deconv_run for DSA, CAMmarker, deconv, ssKL, ssFrobenius
           signature = c("matrix", "character", "list"),
           definition = function(mix, method_name, marker){
             if(method_name == "DSA"){
               library(DSA)
               est_w = EstimateWeight(mix, marker)$weight
             }
             else if(method_name == "CAMmarker"){
               library(CAMTHC)
               # note CAMmarker requires one cell type has at least two markers 
               est_w <- tryCatch({AfromMarkers(mix, marker, scaleRecover = TRUE)},
                                 error = function(e) NA)
             }
#             else if(method_name %in% c("deconf", "ssKL", "ssFrobenius")){
#               library(CellMix)
#               # transform mix and marker 
#               mix <- ExpressionMix(mix + 0.1)
#               marker <- MarkerList(marker)
               
#               cellmix_wrap <- function(expr_mat, marker_mat, method_name){
                 # retrieve the weights and process failed result 
#                 tryCatch({
#                   return(coef(ged(expr_mat, marker_mat, method_name, log = FALSE, rng = 1234, nrun = 10)))
#                 }, error = function(e) NA)
#               }
#               
#               est_w <- cellmix_wrap(mix, marker, method_name)
#             }

             return(est_w)
           }
)

setMethod( f = "deconv_run_refSig", 
           # deconv_run for EPIC and DeconRNASeq
           signature = c("matrix", "character", "ANY"),
           definition = function(mix, method_name, ref, ref_anno = "NA", sig = "NA"){
             if(method_name == "EPIC"){
               library(EPIC) # the ref for the EPIC is in the form of a list: contains sigGene, reference and variance of the .. 
               est_w <- tryCatch({EPIC(mix, ref, withOtherCells = FALSE)$mRNAProportion},
                                 error = function(e) NA)
             }
             else if(method_name == "EPICabsolute"){
               library(EPIC)
               tmp_w <- tryCatch({EPIC(mix, ref)$mRNAProportion},
                                 error = function(e) NA)
               n <- ncol(ref[[1]])
               est_w <- list()
               est_w[[1]] <- tmp_w[,1:n]
               est_w[[2]] <- tmp_w[,n+1]
             }
             else if(method_name == "DeconRNASeq"){
               library(DeconRNASeq)
               mix <- data.frame(mix)
               ref <- data.frame(ref)
               est_w <- DeconRNASeq(mix, ref, fig = FALSE)$out.all
             }
             else if(method_name %in% c("TIMER", "TIMERtumor")){
               source("./code/src/TIMER_codes.R")
               est_w <- TIMER_deconv(mix, ref, ref_anno, sig)
             }
             else if(method_name == "MuSiC"){
               library(MuSiC)
               library(xbioc)
               # ref_eset create 
               sample_name <- colnames(ref)
               meta <- data.frame(sample_name = sample_name, celltype = ref_anno[,2], sample_type = "ref", sampleID = rep(1, ncol(ref)))
               rownames(meta) <- sample_name
               assaydata <- as.matrix(ref)
               phenodata <- new("AnnotatedDataFrame",data = meta)
               ref_eset <- ExpressionSet(assayData = assaydata, phenoData = phenodata)
               # mix_eset create 
               colnames(mix) <- paste0('mix',1:ncol(mix))
               sample_name <- colnames(mix)
               meta <- data.frame(sample_name = sample_name, sample_type = "mix", sampleID = 1:ncol(mix))
               rownames(meta) <- sample_name
               assaydata <- as.matrix(mix)
               phenodata <- new("AnnotatedDataFrame", data = meta)
               mix_eset <- ExpressionSet(assayData = assaydata, phenoData = phenodata)
               # run MuSiC
               music_result <- music_prop(mix_eset, ref_eset, clusters = "celltype", samples = "sample_name", select.ct = unique(ref_anno[,2]), verbose = F)
               est_w <- music_result$Est.prop.weighted
             }
             
             return(est_w)
             }
)







setMethod( f = "deconv_run_free",
           # deconv_run for CAMfree
           signature = c("matrix", "character", "matrix", "numeric"),
           definition = function(mix, method_name, ref, K){
             if(method_name == "CAMfree"){
               library(CAMTHC)
               # run the camfree 
               camfree_result = CAMfree_run(mix, K = K)
               if(class(camfree_result) == "logical"){
                 est_w = NA
               }else{
                 # assign the cell type
                 marker <- MGsforA(camfree_result, K = K) # if the result turn too bed, try to add one more step for further filterign of marker genes 
                 est_s <- Smat(camfree_result,K) 
                 celltype_assign <- spcor_assign(est_marker = marker, est_s = est_s, ref_median = ref)
                 # reorder estimated W
                 est_w <- freeW_reorder(est_w = t(Amat(camfree_result, K)),
                                          assign_result = celltype_assign,
                                          celltype_order = colnames(ref))
                 }
             }
             else if(method_name == "LinSeed"){
                 library(linseed)
                 # run the LinSeed
                 LinSeed_result <- LinSeed_run(mix, K = K)

                 # calculate the fc for cell type assignment 
                 fc <- fc_calc(ref)
                 celltype_assign <- fgsea_assign(LinSeed_result$markers, ref, fc)
                 if(class(celltype_assign) == "logical"){
                   celltype_assign <- zscore_assign(LinSeed_result$markers, ref)
                 }
                 # reorder estimated W 
                 est_w <- freeW_reorder(est_w = LinSeed_result$proportions,
                                        assign_result = celltype_assign,
                                        celltype_order = colnames(ref))
             }
             return(est_w)
           }
)


#@ run deconvolution methods from other platform (MATLAB and Java)
#@ sim1, sim2 and sim3
setMethod(f = "deconv_run_crossPlatform",
          # deconv_run for MMAD and CIBERSORT
          signature = c("character", "sim_params", "character", "character", "character", "character"),
          definition = function(method_name, sim_param, mix_name_prefix, ref_name_prefix, dat_path, sim_prefix){
            #@ method: deconvolution method 
            #@ sim_param: sim_params object 
            #@ mix_name_prefix: mixture object name 
            #@ ref_name_prefix: reference object name 
            #@ dat_path: path of directory that have mixture and reference files 
            #@ sim_prefix: should be one of sim1,sim2 and sim3 (this corresponds to the wrap-up deconvolution script name) 
            #@ output: result of deconvolution will be directly stored in .txt files in the output_path directory
            # create output path if it does no exists 
            output_path <- paste0(dat_path, "/out")
            if (!dir.exists(output_path)){
              dir.create(file.path(output_path))
            }
            
            logfile_name = paste(tail(strsplit2(dat_path, split = "/")[1, ], n = 1), method_name, sep = "_")
            
            if(method_name == "MMAD"){
              # run the sim1_wrapMMAD
              # function sim1_wrapMMAD(mix_name_prefix,ref_name_prefix,nMarker, nDataset, nGrid, dat_path, output_path)
              MMAD_command <- paste(paste0('cd /mnt/data/haijing/simDeconv/paper_deconvBenchmark/code/otherPlatform/MMAD/; nohup matlab -nosplash -nodisplay -nodesktop -r "try; ', sim_prefix, '_wrapMMAD('),
                                    mix_name_prefix, ',', ref_name_prefix, ',', dat_path, ',', output_path, paste0(');catch;end;quit" > ', logfile_name, '.out &'),sep = '\'') 
              
              system(MMAD_command)
            }
            else if(method_name == "CIBERSORT"){
              CIBERSORT_command <- paste0("cd /mnt/data/haijing/simDeconv/paper_deconvBenchmark/code/otherPlatform/CIBERSORT/; nohup bash ", sim_prefix, "_wrapCIBERSORT.sh ", dat_path," ", output_path, " ", mix_name_prefix, " ", ref_name_prefix, " > ", logfile_name, ".out &")
              system(CIBERSORT_command)
            }
            
            else if(method_name == "CIBERSORTx"){
              CIBERSORTx_command <- paste0("cd /mnt/data/haijing/simDeconv/paper_deconvBenchmark/code/otherPlatform/CIBERSORTx/; nohup bash ", sim_prefix, "_wrapCIBERSORTx.sh ", dat_path," ", output_path, " ", mix_name_prefix, " ", ref_name_prefix, " > ", logfile_name, ".out &")
              system(CIBERSORTx_command)
            }

            
          })

setMethod(f = "deconv_run_crossPlatform_oneMat",
          # deconv_run for MMAD and CIBERSORT
          signature = c("character", "character", "character", "character"),
          definition = function(method_name, mix_name, ref_name, dat_path){
            # create output path if it does no exists 
            output_path <- paste0(dat_path, "/out")
            if (!dir.exists(output_path)){
              dir.create(file.path(output_path))
            }
            
            logfile_name = paste(tail(strsplit2(dat_path, split = "/")[1, ], n = 1), method_name, sep = "_")
            
            if(method_name == "MMAD"){
              # run the sim1_wrapMMAD
              # function sim1_wrapMMAD(mix_name_prefix,ref_name_prefix,nMarker, nDataset, nGrid, dat_path, output_path)
              MMAD_command <- paste(paste0('cd /mnt/data/haijing/simDeconv/otherPlatform/MMAD/; nohup matlab -nosplash -nodisplay -nodesktop -r "try; mmad_oneMat('),
                                    mix_name, ',', ref_name, ',', dat_path, ',', output_path, paste0(');catch;end;quit" > ', logfile_name, '.out &'),sep = '\'') 
              system(MMAD_command)
            }
            else if(method_name == "CIBERSORT"){
              ref_anno_name <- paste0(ref_name, "_anno_cibersort")
              CIBERSORT_command <- paste0("cd /mnt/data/haijing/simDeconv/otherPlatform/CIBERSORT/; nohup bash cibersort_oneMat.sh ", dat_path," ", output_path, " ", mix_name, " ", ref_name, " ", ref_anno_name, " > ", logfile_name, ".out &")
              system(CIBERSORT_command)
            }
          })



setMethod(f = "deconv_write",
          # function to write data to the write_path for cross-platform deconvolution algorithms
          signature = c("matrix", "character", "character", "character", "character"),
          definition = function(dat, method_name, file_type, file_name, write_path){
            # first check if the directory exists 
            if (!dir.exists(write_path)){
              dir.create(file.path(write_path))
            }
            # 
            if(method_name == "MMAD"){
              if(file_type == "mix"){
                write.table(log2(dat + 1), file = paste0(write_path, "/", file_name), sep = "\t")
              }
              else if(file_type == "ref"){
                write.table(dat, file = paste0(write_path, "/", file_name), sep = "\t")
              }
              
            }
            else if(method_name %in% c("CIBERSORT", "CIBERSORTx")){
              
              if(file_type == "mix"){
                first_row <- c("!Sample_title", paste0("Mixture", 1:ncol(dat)))
                dat <- CIBERSORT.transform(dat, first_row)
                write.table(dat, file = paste0(write_path, "/", file_name), row.names = FALSE, quote = FALSE, sep = "\t")
              }
              else if(file_type == "ref"){
                first_row <- c("!Sample_title", colnames(dat))
                dat <- CIBERSORT.transform(dat, first_row)
                write.table(dat, file = paste0(write_path, "/", file_name), row.names = FALSE, quote = FALSE, sep = "\t")
              }
              else if(file_type == "ref_anno"){
                write.table(dat, file = paste0(write_path, "/", file_name), row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
              }
              
            }
            
            
          }
)



setMethod(f = "deconv_read",
          signature = c("character", "character", "character", "numeric", "numeric"),
          definition = function(file_name, method_name, read_path, nsim, nComp){
            
            full_file_name = paste0(read_path, file_name)
            if (method_name == "MMAD") {
              est_W <- read.table(full_file_name, header = FALSE, sep = "\t")[, 1:nsim]
            } else if (method_name == "CIBERSORT") {
              est_W <- tryCatch({read.table(full_file_name, header = FALSE, sep = "\t", skip = 10)[, 1:nComp + 1]
              }, error = function(e) NA)
            } else if (method_name == "CIBERSORTx"){
              est_W <- tryCatch({read.table(full_file_name, header = TRUE, sep = "\t")[, 1:nComp + 1]
              }, error = function(e) NA)
            }
            return(est_W)
          }
)


##### II: sim1 specific functions  ########
setMethod( f = "deconv_analysis", 
           signature = c("character", "deconv_params", "sim1_params", "list"),
           definition = function(method_name, deconv_param, sim_param, mix_list){
           
           nMarker <- length(sim_param@dataset_name)
           nDataset <- length(sim_param@dataset_name)
           gene_id <- sim_param@gene_id
           
           est_W <- list() 
           for(m in 1:nMarker){
             tmp <- list() 
             for( j in 1:nDataset){
               M <- mix_list[[j]]
            #   Marker <- marker_list[[m]][[2]]
               
               if(method_name %in% c("DSA", "CAMmarker", "deconf", "ssKL", "ssFrobenius")){
                 marker <- deconv_param@marker_list[[m]][[2]]
                 tmp[[j]] <- lapply(M, deconv_run_marker, method_name = method_name, marker = marker)
               }
               else if(method_name %in% c("EPIC", "EPICabsolute")){
                 ref <- deconv_param@epic_ref_list[[m]]
                 tmp[[j]] <- lapply(M, deconv_run_refSig, method_name = method_name, ref = ref)
               }
               else if(method_name == "DeconRNASeq"){
                 ref <- deconv_param@ref_median_list[[m]][deconv_param@sig_list[[m]], ]
                 tmp[[j]] <- lapply(M, deconv_run_refSig, method_name = method_name, ref = ref)
               }
               else if(method_name == "TIMER"){
                 ref <- as.matrix(deconv_param@ref_list[[m]])
                 anno <- deconv_param@ref_anno_others[[m]][,2]
                 sig <- deconv_param@sig_list[[m]]
                 ref_anno <- colnames(ref)
                 names(ref_anno) <- anno
                 
                 tmp[[j]] <- lapply(M, deconv_run_refSig, method_name = method_name, ref_anno = ref_anno, sig = sig, ref = ref)
               }
               else if(method_name == "MuSiC"){
                 ref <- deconv_param@ref_list[[m]]
                 ref_anno <- deconv_param@ref_anno_others[[m]]
                 tmp[[j]] <- lapply(M, deconv_run_refSig, method_name = method_name, ref = ref, ref_anno = ref_anno)
               }
               else if(method_name %in% c("CAMfree", "LinSeed")){
                 # retreive result 
                 ref <- deconv_param@ref_median_list[[m]]
                 tmp[[j]] <- lapply(M, deconv_run_free, method_name, ref = ref, K = ncol(ref))
               }
             }
             est_W[[m]] <- tmp
           }
           return(est_W)
           }
)

# 
setMethod(f = "wrap_mix_write",
          signature = c("list", "character", "sim1_params", "character", "character", "character"),
          definition = function(mix_list, prefix, sim_param, file_type, method_name, write_path){
            
            nDataset<- length(sim_param@dataset_name)
            nGrid <- length(sim_param@grid)
            
            # if the destined directory does not exit, create one
            if (!dir.exists(write_path)){
              dir.create(file.path(write_path))
            }
            for(i in 1:nDataset){
              for(j in 1:nGrid){
                file_name = paste0(prefix, "_D", i, "P", j, ".txt")
                deconv_write(mix_list[[i]][[j]], method_name = method_name, file_type = file_type,file_name = file_name, write_path = write_path)
              }
            }
          }
)

setMethod(f = "wrap_mix_write",
          signature = c("list", "character", "comp_sim_params", "character", "character", "character"),
          definition = function(mix_list, prefix, sim_param, file_type, method_name, write_path){
            # if the destined directory does not exit, create one
            if (!dir.exists(write_path)){
              dir.create(file.path(write_path))
            }
            
            for(i in 1:length(mix_list)){
              file_name = paste0(prefix, "_C", i, ".txt")
              deconv_write(mix_list[[i]], method_name = method_name, file_type = file_type, file_name = file_name, write_path = write_path)
            }

          }
)



setMethod(f = "wrap_ref_write",
          signature = c("list", "character", "sim1_params", "character", "character"),
          definition = function(ref_list, prefix, sim_param, file_type, method_name, write_path){
            
            if(method_name == "MMAD"){
              gene_id <- sim_param@gene_id
              for(i in 1:length(ref_list)){
                ref <- MMAD_markerTransform(ref_list[[i]][[1]], gene_id)
                file_name = paste0(prefix, "_D", i, ".txt")
                deconv_write(ref, method_name = method_name, file_type = file_type, file_name = file_name, write_path = write_path)
              }
            } else if (method_name %in% c("CIBERSORT", "CIBERSORTx")){
              gene_id <- sim_param@gene_id
              for(i in 1:length(ref_list)){
                ref <- as.matrix(ref_list[[i]])
                file_name = paste0(prefix, "_D", i, ".txt")
                deconv_write(ref, method_name = method_name, file_type = file_type, file_name = file_name, write_path = write_path)
              }
            }
          }
)

setMethod(f = "wrap_ref_write",
          signature = c("list", "character", "comp_sim_params", "character", "character"),
          definition = function(ref_list, prefix, sim_param, file_type, method_name, write_path){
            
            if(method_name == "MMAD"){
              gene_id <- sim_param@gene_id
              for(i in 1:length(ref_list)){
                ref <- MMAD_markerTransform(ref_list[[i]][[1]], gene_id)
                file_name = paste0(prefix, "_C", i, ".txt")
                deconv_write(ref, method_name = method_name, file_type = file_type, file_name = file_name, write_path = write_path)
              }
            }else if(method_name %in% c("CIBERSORT", "CIBERSORTx")){
              gene_id <- sim_param@gene_id
              for(i in 1:length(ref_list)){
                ref <- as.matrix(ref_list[[i]])
                file_name = paste0(prefix, "_C", i, ".txt")
                deconv_write(ref, method_name = method_name, file_type = file_type, file_name = file_name, write_path = write_path)
              }
            }
          }
)

#@ read deconvolution results 
#@ sim1 
setMethod(f = "wrap_result_read",
          signature = c("character", "sim1_params", "character", "character"),
          definition = function(method_name, sim_param, prefix, read_path){
            #@ method_name: deconvolution method name
            #@ sim_param: sim_params object 
            #@ prefix: result prefix (corresponds to the prefix of results .txt files name)
            #@ read_path: path to read results 
            #@ output: est_W[[m]][[d]], lists of deconvolution results, [[m]] for Marker, [[d]] for Datasets
            nDataset<- length(sim_param@dataset_name)
            nMarker <- length(sim_param@dataset_name)
            nsim <- sim_param@nsim
            # sim_param is deconv_read needed   
            est_W <- list()
            for(m in 1:nMarker){
              tmp <- list()
              for(d in 1:nDataset){
                file_names = paste0("Results_", prefix, "_D", d, "P", 1:10, "M", m, ".txt")
                tmp[[d]] <- lapply(file_names, deconv_read, method_name, read_path, nsim, nComp = 3)
              }
              est_W[[m]] <- tmp
            }
            return(est_W)
          }
)

#@ read deconvolution results 
#@ sim2 and sim3 
setMethod(f = "wrap_result_read",
          signature = c("character", "comp_sim_params", "character", "character"),
          definition = function(method_name, sim_param, prefix, read_path){
            
            n_comp <- sim_param@n_comp
            nsim <- sim_param@nsim
            # sim_param is deconv_read needed   
            est_W <- list()
            for(i in 1:length(n_comp)){
              file_name = paste0("Results_", prefix, "_C", i, ".txt")
              est_W[[i]] <- deconv_read(file_name = file_name, method_name = method_name, read_path = read_path, nsim = nsim, nComp = n_comp[i])
            }
            return(est_W)
          }
)


#@ wrap-up deconvolution analysis function
#@ sim2 and sim3 
setMethod( f = "deconv_analysis", 
           signature = c("character", "deconv_params", "comp_sim_params", "list"),
           definition = function(method_name, deconv_param, sim_param, mix_list){
             #@ method_name: targed deconvolution method 
             #@ deconv_param: deconv_params object 
             #@ sim_param: sim_params object 
             #@ mix_list: lists of mixtures 
             gene_id <- sim_param@gene_id
             n_comp <- sim_param@n_comp
             est_W <- list()
             for(i in 1:length(mix_list)){
               if(method_name %in% c("DSA", "CAMmarker")){
                 marker <- deconv_param@marker_list[[i]][[2]]
                 est_W[[i]] <- deconv_run_marker(mix = mix_list[[i]], method_name = method_name, marker = marker)
               }
               else if(method_name %in% c("EPIC", "EPICabsolute")){
                 ref <- deconv_param@epic_ref_list[[i]]
                 est_W[[i]] <- deconv_run_refSig(mix = mix_list[[i]], method_name = method_name, ref = ref)
               }
               else if(method_name == "DeconRNASeq"){
                 ref <- deconv_param@ref_median_list[[i]][deconv_param@sig_list[[i]], ]
                 est_W[[i]] <- deconv_run_refSig(mix = mix_list[[i]], method_name = method_name, ref = ref)
               }
               else if(method_name == "TIMER"){
                 ref <- as.matrix(deconv_param@ref_list[[i]])
                 anno <- deconv_param@ref_anno_others[[i]][,2]
                 sig <- deconv_param@sig_list[[i]]
                 ref_anno <- colnames(ref)
                 names(ref_anno) <- anno
                 
                 est_W[[i]] <- deconv_run_refSig(mix = mix_list[[i]], method_name = method_name, ref_anno = ref_anno, sig = sig, ref = ref)
               }
               else if(method_name == "TIMERtumor"){
                 ref <- as.matrix(deconv_param@ref_list[[i]])
                 anno <- deconv_param@ref_anno_others[[i]][,2]
                 #  TIMER_tumor is specific for sim3 where we add additional filter to only include genes that are negatively correlates to the tumor purity 
                 sig <- intersect(deconv_param@sig_list[[i]], timer_purity_genes)
                 ref_anno <- colnames(ref)
                 names(ref_anno) <- anno
                 
                 est_W[[i]] <- deconv_run_refSig(mix = mix_list[[i]], method_name = method_name, ref_anno = ref_anno, sig = sig, ref = ref)
               }
               else if(method_name == "MuSiC"){
                 ref <- deconv_param@ref_list[[i]]
                 ref_anno <- deconv_param@ref_anno_others[[i]]
                 est_W[[i]] <- deconv_run_refSig(mix = mix_list[[i]], method_name = method_name, ref = ref, ref_anno = ref_anno)
               }
               else if(method_name %in% c("CAMfree", "LinSeed")){
                 ref <- as.matrix(deconv_param@ref_median_list[[i]])
                 est_W[[i]] <- deconv_run_free(mix = mix_list[[i]], method_name = method_name, ref = ref, K = n_comp[i])
               }
             }
             return(est_W)
           }
)

# oneMat specific functions 
setMethod( f = "deconv_analysis_oneMat", 
           signature = c("character", "deconv_params", "ANY"),
           definition = function(method_name, deconv_param_oneMat, mix){
             gene_id <- rownames(mix)
             mix <- as.matrix(mix)
             
             est_W <- list()
             if(method_name %in% c("DSA", "CAMmarker")){
               marker <- deconv_param_oneMat@marker_list[[2]]
               est_W <- deconv_run_marker(mix = mix, method_name = method_name, marker = marker)
             }
             else if(method_name == "EPIC"){
               ref <- deconv_param_oneMat@epic_ref_list
               est_W <- deconv_run_refSig(mix = mix, method_name = method_name, ref = ref)
             }
             else if(method_name == "DeconRNASeq"){
               ref <- deconv_param_oneMat@ref_median[deconv_param_oneMat@sigGene, ]
               est_W <- deconv_run_refSig(mix = mix, method_name = method_name, ref = ref)
             }
             else if(method_name == "TIMER"){
               ref <- as.matrix(deconv_param_oneMat@ref)
               anno <- deconv_param_oneMat@ref_anno_others[,2]
               sig <- deconv_param_oneMat@sigGene
               ref_anno <- colnames(ref)
               names(ref_anno) <- anno
               
               est_W <- deconv_run_refSig(mix = mix, method_name = method_name, ref_anno = ref_anno, sig = sig, ref = ref)
             }
             else if(method_name == "MuSiC"){
               ref <- deconv_param_oneMat@ref
               ref_anno <- deconv_param_oneMat@ref_anno_others
               est_W <- deconv_run_refSig(mix = mix, method_name = method_name, ref = ref, ref_anno = ref_anno)
             }
             else if(method_name %in% c("CAMfree", "LinSeed")){
               ref <- as.matrix(deconv_param_oneMat@ref_median)
               est_W <- deconv_run_free(mix = mix, method_name = method_name, ref = ref, K = ncol(ref))
             }
             return(est_W)
           }
)













