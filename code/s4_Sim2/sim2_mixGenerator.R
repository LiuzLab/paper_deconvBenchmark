load("./output/Sim2/sim2_env.RData")
#####------------------------ I: source functions & load dataset ------------------------#####
# load processed data for simulation 
source("./code/src/Classes_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Functions_all.R")
source("./code/src/Methods_simulation.R")
# ------------------------ simulate weight 
# generate orthogonal weight matrix 
sim2_W_orthog <- orthog_weight(sim2_params_ob, N = 1000, sim_celltype = sim2_cellType)
sim2_W_orthog_kappa <- laply(sim2_W_orthog, kappa)
saveRDS(sim2_W_orthog, "./output/Sim2/sim2_W_orthog.rds")
# generate real weight matrix 
sim2_W_real <- real_weight(sim2_params_ob, N = 1000, sim_prop_mat = sim2_prop_mat)
saveRDS(sim2_W_real, "./output/Sim2/sim2_W_real.rds")
# ------------------------ create perturbed nb simulations 
# initiate nb simulation paramters 
nb_params_ob <- nb_params()
nb_params_ob <- set_nb_params(nb_params_ob, nsim_val = 20, libsize_val = 12*10^6)
celltype <- names(sim2_params_ob@extract_pattern_immune)
for(k in 1:length(celltype)){
  nb_params_ob@prefix = celltype[k]
  input_name <- paste("sim2", celltype[k], "countNorm", sep = "_")
  output_name <- paste("sim2", celltype[k], "nb", "count", sep = "_")
  
  assign(output_name, sim_gammaPoisson(nb_param = nb_params_ob,
                                       reference = as.matrix(get(input_name))))
}
#------------------------ create nb mixtures - noTumor 
# first create nb mixture in count unit 
sim2_M_count_orthog <- sim_mix(sim_param = sim2_params_ob, W = sim2_W_orthog, sim_celltype = sim2_cellType, sim_name = "sim2")
sim2_M_count_real <- sim_mix(sim_param = sim2_params_ob,  W = sim2_W_real, sim_celltype = sim2_cellType, sim_name = "sim2")
# unit transform: from count to countNorm, cpm and tpm
eff_length <- sim2_eff_length_mean
unit <- sim2_params_ob@unit
for(u in 2:length(unit)){
  input_name1 <- paste("sim2", "M","count", "orthog", sep = "_")
  input_name2 <- paste("sim2", "M", "count", "real", sep = "_")
  output_name1 <- paste("sim2", "M", unit[u], "orthog", sep = "_")
  output_name2 <- paste("sim2", "M", unit[u], "real", sep = "_")
  assign(output_name1, sim_mix_unit_transform(sim_param = sim2_params_ob, 
                                              mix_list = get(input_name1), 
                                              unit_val = unit[u], 
                                              eff_length = eff_length))
  assign(output_name2, sim_mix_unit_transform(sim_param = sim2_params_ob, 
                                              mix_list = get(input_name2), 
                                              unit_val = unit[u],
                                              eff_length = eff_length))
}
save.image("./output/Sim2/sim2_mix.RData")