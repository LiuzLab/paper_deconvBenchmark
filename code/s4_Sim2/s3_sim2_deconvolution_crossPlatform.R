load("./output/Sim2/sim2_mix.RData")
load("./output/Sim2/sim2_ref.RData")
load("./output/Sim2/sim2_estW.RData")
source("./code/src/Generics_all.R")
source("./code/src/Methods_deconvolution.R")
source("./code/src/Functions_all.R")
unit <- sim2_params_ob@unit
mix_type <- sim2_params_ob@mix_type
# library(Rserve)
# Rserve(args="--no-save")
methods <- "MMAD"
for(u in 1:length(unit)){
  for(s in 1:length(mix_type)){
    # set the name of inputs for analysis 
    mix_name <- paste("sim2", "M", unit[u], mix_type[s], sep = "_")
    ref_name <- paste("sim2", "ref", unit[u], sep = "_")
    ref_anno_name <- paste("sim2", "ref", "anno", "cibersort", sep = "_")
    marker_name <- paste("sim2", "marker", "countNorm", sep = "_")
    for(m in 1:length(methods)){
      write_path <- paste0("/mnt/data/haijing/simDeconv/paper_deconvBenchmark/output/Sim2/", methods[m])
      # Step 1: write mix
      wrap_mix_write(mix_list = get(mix_name), prefix = mix_name, sim_param = sim2_params_ob, file_type = "mix", method_name = methods[m], write_path = write_path)
      # Step 2: write reference
      if(methods[m] == "MMAD"){
        wrap_ref_write(ref_list = get(marker_name), prefix = marker_name, sim_param = sim2_params_ob, file_type = "ref", method_name = methods[m], write_path = write_path)
        ref_name_prefix = marker_name
      }
      else if(methods[m] == "CIBERSORT"){
        wrap_ref_write(ref_list = get(ref_name), prefix = ref_name, sim_param = sim2_params_ob, file_type = "ref", method_name = methods[m], write_path = write_path)
        wrap_ref_write(ref_list = get(ref_anno_name), prefix = ref_anno_name, sim_param = sim2_params_ob, file_type = "ref_anno", method_name = methods[m], write_path = write_path)
        ref_name_prefix = ref_name
      }
      
      # Step 3: run the deconvolution 
      deconv_run_crossPlatform(method_name = methods[m], 
                               sim_param = sim2_params_ob, 
                               mix_name_prefix = mix_name, 
                               ref_name_prefix = ref_name_prefix, 
                               dat_path = write_path,
                               sim_prefix = "sim2")
      
    }
  }
}
methods <- c("CIBERSORT", "MMAD")
# result read 
for(u in 1:length(unit)){
  for(s in 1:length(mix_type)){
    mix_name <- paste("sim2", "M", unit[u], mix_type[s], sep = "_")
    
    for( m in 1:length(methods)){
      read_path <- paste0("/mnt/data/haijing/simDeconv/paper_deconvBenchmark/output/Sim2/", methods[m],"/out/")
      
      result_name <- paste("sim2", "estW", methods[m], unit[u], mix_type[s], sep = "_")
      assign(result_name, wrap_result_read(method_name = methods[m],
                                           prefix = mix_name, 
                                           sim_param = sim2_params_ob,
                                           read_path = read_path))
    }
  }
}
# save results 
save.image(file = "./output/Sim2/sim2_deconvolution_crossPlatform.RData")
save(list = ls(pattern = "estW"), file = "./output/Sim2/sim2_estW.RData")