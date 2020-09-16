# for deconvolution of sim1 simModel - methods include all R platform based deconvolution algorithms 
# We recommend to run this script in background as it takes a while to finish. 
setwd("/abs_path/") # Please set the path to the absolute path of your working directory. 
load("./output/Sim2/sim2_mix.RData")
load("./output/Sim2/sim2_ref.RData")
source("./code/src/Functions_all.R")
source("./code/src/Classes_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Methods_deconvolution.R")
# initiate deconv param object
deconv_params_ob <- deconv_params()
deconv_params_ob@marker_list <- sim2_marker_countNorm
deconv_params_ob@sig_list <- sim2_sigGene_count
deconv_params_ob@ref_anno_cibersort <- sim2_ref_anno_cibersort
deconv_params_ob@ref_anno_others <- sim2_ref_anno_others
# start to run the deconv 
unit <- sim2_params_ob@unit
# run first set of methods, based on R platform and quick 

methods <- c("DSA", "CAMmarker", "EPIC", "DeconRNASeq","TIMER", "MuSiC", "LinSeed")
mix_type <- sim2_params_ob@mix_type
sim2_time <- list() 
for(u in 1:length(unit)){
  deconv_params_ob@epic_ref_list <- get(paste("sim2", "epic", "ref", unit[u], sep = "_"))
  deconv_params_ob@ref_list <- get(paste("sim2", "ref", unit[u], sep = "_"))
  deconv_params_ob@ref_median_list <- get(paste("sim2", "ref", "median", unit[u], sep = "_"))
  sim2_time[[u]] <- matrix(0, nrow = length(methods), ncol = length(mix_type))
  for(s in 1:length(mix_type)){
    for(m in 1:length(methods)){
      start_time = Sys.time()
      result_name <- paste("sim2", "estW", methods[m], unit[u], mix_type[s], sep = "_")
      
      input_mix_name <- paste("sim2", "M", unit[u], mix_type[s], sep = "_")
      assign(result_name, deconv_analysis(method_name = methods[m], 
                                          deconv_param = deconv_params_ob,
                                          sim_param = sim2_params_ob,
                                          mix_list = get(input_mix_name)))
      
      end_time = Sys.time()
      sim2_time[[u]][m,s] <- end_time - start_time
    }
  }
  rownames(sim2_time[[u]]) = methods
  colnames(sim2_time[[u]]) = mix_type
}
save(list = ls(pattern = "estW"), file = "/abs_path/output/Sim2/sim2_estW.RData") # Please set the path to the targeted directory. 