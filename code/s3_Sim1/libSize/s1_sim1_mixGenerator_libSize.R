# generate mixture for Sim1 libSize
load("./output/Sim1/sim1_ref.RData")
source("./code/src/Generics_all.R")
source("./code/src/Methods_simulation.R")
source("./code/src/Functions_all.R")
unit <- sim1_params_ob@unit
#------------ read predefined weight 
sim1_W <- as.matrix(read.table("./data/sim1_W.txt", sep = "\t"))
#------------ simulation with negative binomial distribution 
# set the negative binomial parameters
nb_params_ob <- nb_params()
nb_params_ob_lib1 <- set_nb_params(nb_params_ob, nsim_val = 10, libsize_val = 12*10^6)
nb_params_ob_lib2 <- set_nb_params(nb_params_ob, nsim_val = 10, libsize_val = 24*10^6 )
# add noise to the original cell type specific profiles - create perturbed references 
for(k in sim1_params_ob@celltype){
  # set the output name of signature profiles 
  # note that negative binomial simulation take the countNorm units and output count units 
  profile_name_input <- paste(k, "countNorm", sep = "_")
  profile_name_output1 <- paste(k, "count", "nb", "lib1", sep = "_")
  profile_name_output2 <- paste(k, "count", "nb", "lib2", sep = "_")
  
  assign(profile_name_output1, wrap_gammaPoisson(sim_param = sim1_params_ob, 
                                                      nb_param = nb_params_ob_lib1, 
                                                      expr_list = get(profile_name_input),
                                                      interval = 1.8))
  assign(profile_name_output2, wrap_gammaPoisson(sim_param = sim1_params_ob, 
                                                      nb_param = nb_params_ob_lib2, 
                                                      expr_list = get(profile_name_input),
                                                      interval = 1.8))
}
# in silico mixing with perturbed references 
mix_name <- paste("sim1","M", "libSize","count", sep = "_")
assign(mix_name, sim_libsize_mix(sim_param = sim1_params_ob, 
                                 suffix_lib1 = "count_nb_lib1",
                                 suffix_lib2 = "count_nb_lib2",
                                 W = sim1_W))
# unit transform from count to countNorm, cpm and tpm 
for(u in 2:length(unit)){
  mix_name_out <- paste("sim1", "M", "libSize", unit[u], sep = "_")
  assign(mix_name_out, sim_mix_unit_transform(sim_param = sim1_params_ob,
                                              mix_list = sim1_M_libSize_count,
                                              unit_val = unit[u],
                                              eff_length = sim1_eff_length))
  
}
save.image(file = "./output/Sim1/libSize/sim1_mix_libSize.RData")