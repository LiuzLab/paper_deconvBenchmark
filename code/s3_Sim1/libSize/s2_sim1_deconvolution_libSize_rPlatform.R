# for deconvolution of sim1 libSize - methods include all R platform based deconvolution algorithms 
setwd("/mnt/data/haijing/simDeconv/paper_deconvBenchmark")
load("./output/Sim1/libSize/sim1_mix_libSize.RData")
source("./code/src/Generics_all.R")
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


for(u in 1:length(unit)){
  deconv_params_ob@epic_ref_list <- get(paste("sim1", "epic", "ref", unit[u], sep = "_"))
  deconv_params_ob@ref_list <- get(paste("sim1", "ref", unit[u], sep = "_"))
  deconv_params_ob@ref_median_list <- get(paste("sim1", "ref", "median", unit[u], sep = "_"))
  for(m in 1:length(methods)){
      result_name <- paste("sim1", "estW", methods[m], "libSize", unit[u], sep = "_")
      
      input_mix_name <- paste("sim1", "M", "libSize", unit[u], sep = "_")
      assign(result_name, deconv_analysis(method_name = methods[m], 
                                          deconv_param = deconv_params_ob,
                                          sim_param = sim1_params_ob,
                                          mix_list = get(input_mix_name)))
    }
}
save.image(file = "/mnt/data/haijing/simDeconv/paper_deconvBenchmark/output/Sim1/sim1_libSize_deconvolution_rPlatform.RData")
save(list = ls(pattern = "estW"), file = "/mnt/data/haijing/simDeconv/paper_deconvBenchmark/output/Sim1/libSize/sim1_libSize_estW.RData")