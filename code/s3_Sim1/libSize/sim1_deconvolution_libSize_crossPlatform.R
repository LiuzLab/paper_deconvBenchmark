load("./output/Sim1/libSize/sim1_mix_libSize.RData")
load("./output/Sim1/libSize/sim1_libSize_estW.RData")
load('./output/sim1_env.RData')
source("./code/src/Functions_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Methods_deconvolution.R")

unit <- sim1_params_ob@unit
sim_model <- sim1_params_ob@sim_model

#library(Rserve)
#Rserve(args="--no-save")

# have finished running all the crossPlatform methods, please skip this step and go directly for the result reading 
methods <- "CIBERSORTx" # Please run the method one by one (CIBERSORTx, CIBERSORT and MMAD) to avoid exaustion of computing resource. 
for(u in 1:length(unit)){
    # set the name of inputs for analysis 
    mix_name <- paste("sim1", "M","libSize", unit[u], sep = "_")
    ref_name <- paste("sim1", "ref", unit[u], sep = "_")
    ref_anno_name <- paste("sim1", "ref", "anno", "cibersort", sep = "_")
    marker_name <- paste("sim1", "marker", "countNorm", sep = "_")
    for(m in 1:length(methods)){
      write_path <- paste0("/abs_path/output/Sim1/libSize/", methods[m]) # Please set the path to the targeted path. 
      # Step 1: write mix
      wrap_mix_write(mix_list = get(mix_name), prefix = mix_name, sim_param = sim1_params_ob, file_type = "mix", method_name = methods[m], write_path = write_path)
      # Step 2: write reference
      if(methods[m] == "MMAD"){
        wrap_ref_write(ref_list = get(marker_name), prefix = marker_name, sim_param = sim1_params_ob, file_type = "ref", method_name = methods[m], write_path = write_path)
        ref_name_prefix = marker_name
      }
      else if(methods[m] %in% c("CIBERSORT", "CIBERSORTx")){
        wrap_ref_write(ref_list = get(ref_name), prefix = ref_name, sim_param = sim1_params_ob, file_type = "ref", method_name = methods[m], write_path = write_path)
        wrap_ref_write(ref_list = get(ref_anno_name), prefix = ref_anno_name, sim_param = sim1_params_ob, file_type = "ref_anno", method_name = methods[m], write_path = write_path)
        ref_name_prefix = ref_name
      }
      
      # Step 3: run the deconvolution 
      deconv_run_crossPlatform(method_name = methods[m], 
                               sim_param = sim1_params_ob, 
                               mix_name_prefix = mix_name, 
                               ref_name_prefix = ref_name_prefix, 
                               dat_path = write_path,
                               sim_prefix = "sim1")
      
    }
}


methods <- c("CIBERSORT", "CIBERSORTx", "MMAD")
# result read 
for(u in 1:length(unit)){
  for(s in 1:length(sim_model)){
    mix_name <- paste("sim1", "M", "libSize", unit[u], sep = "_")
    
    for( m in 1:length(methods)){
      read_path <- paste0("/abs_path/output/Sim1/libSize/", methods[m],"/out/") # Please set the path to the targeted path. 
      
      result_name <- paste("sim1", "estW", methods[m], "libSize", unit[u], sep = "_")
      assign(result_name, wrap_result_read(method_name = methods[m],
                                           prefix = mix_name, 
                                           sim_param = sim1_params_ob,
                                           read_path = read_path))
    }
  }
}
# save results 
save(list = ls(pattern = "estW"), file = "./output/Sim1/libSize/sim1_libSize_estW.RData") # Please set the path to the targeted path. 
