# for deconvolution of sim1 simModel - methods include all R platform based deconvolution algorithms 
# We recommend to run this script in background as it takes a while to finish. 
setwd("/abs_path/") # Please set the path to the absolute path of your working directory. 
load("./output/Sim3/sim3_mix.RData")
load("./output/Sim3/sim3_ref.RData")
source("./code/src/Functions_all.R")
source("./code/src/Classes_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Methods_deconvolution.R")
# initiate deconv param object
deconv_params_ob <- deconv_params()
deconv_params_ob@marker_list <- sim3_marker_countNorm
deconv_params_ob@sig_list <- sim3_sigGene_count
deconv_params_ob@ref_anno_cibersort <- sim3_ref_anno_cibersort
deconv_params_ob@ref_anno_others <- sim3_ref_anno_others
timer_purity_genes <- readRDS("./output/Sim3/timer_purity_genes.rds")
# start to run the deconv 
unit <- sim3_params_ob@unit
# run first set of methods, based on R platform and quick 
methods <- c("DSA", "CAMmarker", "EPIC", "EPICabsolute", "DeconRNASeq","TIMER", "TIMERtumor", "MuSiC", "LinSeed")
mix_type <- sim3_params_ob@mix_type
tumor_content <- sim3_params_ob@tumor_content
sim3_time <- list() 
for(u in 1:length(unit)){
  deconv_params_ob@epic_ref_list <- get(paste("sim3", "epic", "ref", unit[u], sep = "_"))
  deconv_params_ob@ref_list <- get(paste("sim3", "ref", unit[u], sep = "_"))
  deconv_params_ob@ref_median_list <- get(paste("sim3", "ref", "median", unit[u], sep = "_"))
  sim3_time[[u]] <- matrix(0, nrow = length(methods), ncol = length(mix_type)*length(tumor_content))
  tmp_name <- vector()
  for(s in 1:length(mix_type)){
    for(j in 1:length(tumor_content)){
      for(m in 1:length(methods)){
        start_time = Sys.time()
        result_name <- paste("sim3", "estW", methods[m], unit[u], mix_type[s], tumor_content[j], sep = "_")
        
        input_mix_name <- paste("sim3", "M", unit[u], mix_type[s], tumor_content[j], sep = "_")
        assign(result_name, deconv_analysis(method_name = methods[m],
                                            deconv_param = deconv_params_ob,
                                            sim_param = sim3_params_ob,
                                            mix_list = get(input_mix_name)))
        
        end_time = Sys.time()
        sim3_time[[u]][m, length(tumor_content)*(s-1) + j] <- end_time - start_time
        tmp_name[length(tumor_content)*(s-1) + j] <- paste(mix_type[s], tumor_content[j], sep = "_")
      }
    }
  }
  colnames(sim3_time[[u]]) <- tmp_name
  rownames(sim3_time[[u]]) = methods
}
save(list = ls(pattern = "estW"), file = "/abs_path/output/Sim3/sim3_estW.RData") # Please set the path to the targeted directory. 