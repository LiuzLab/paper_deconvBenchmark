# for deconvolution of sim1 simModel - methods include all R platform based deconvolution algorithms 
# This scripts takes a while to finish. Therefore, running the scripts in background is recommended. 
# Since all path will be absolute, please make sure to reset the path accordingly. 
setwd("/abs_path/") # Please set the path to the absolute path of your working directory. 
load("./output/Sim1/simModel/sim1_mix_simModel.RData")
source("./code/src/Generics_all.R")
source("./code/src/Classes_all.R")
source("./code/src/Methods_deconvolution.R")
source("./code/src/Functions_all.R")
# initiate deconv param object
deconv_params_ob <- deconv_params()
deconv_params_ob@marker_list <- sim1_marker_countNorm
deconv_params_ob@sig_list <- sim1_sigGene_count
deconv_params_ob@ref_anno_cibersort <- sim1_ref_anno_cibersort
deconv_params_ob@ref_anno_others <- sim1_ref_anno_others
# start to run the deconv 
unit <- sim1_params_ob@unit
# run first set of methods, based on R platform and quick 

methods <- c("DSA", "CAMmarker", "EPIC", "DeconRNASeq","TIMER", "MuSiC", "CAMfree", "LinSeed") 
sim_model <- sim1_params_ob@sim_model

sim1_time <- list() 
for(u in 1:length(unit)){
  deconv_params_ob@epic_ref_list <- get(paste("sim1", "epic", "ref", unit[u], sep = "_"))
  deconv_params_ob@ref_list <- get(paste("sim1", "ref", unit[u], sep = "_"))
  deconv_params_ob@ref_median_list <- get(paste("sim1", "ref", "median", unit[u], sep = "_"))
  sim1_time[[u]] <- matrix(0, nrow = length(methods), ncol = length(sim_model))
  for(s in 1:length(sim_model)){
    for(m in 1:length(methods)){
      start_time = Sys.time()
      result_name <- paste("sim1", "estW", methods[m], sim_model[s], unit[u], sep = "_")
      
      input_mix_name <- paste("sim1", "M", sim_model[s], unit[u], sep = "_")
      assign(result_name, deconv_analysis(method_name = methods[m], 
                                          deconv_param = deconv_params_ob,
                                          sim_param = sim1_params_ob,
                                          mix_list = get(input_mix_name)))

      end_time = Sys.time()
      sim1_time[[u]][m,s] <- end_time - start_time
    }
  }
}

save(list = ls(pattern = "estW"), file = "/abs_path/output/Sim1/simModel/sim1_simModel_estW.RData") # Please set the path to the targeted directory.