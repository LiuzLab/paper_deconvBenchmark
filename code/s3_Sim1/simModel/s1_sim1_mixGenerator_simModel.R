# generate mixture for Sim1 simModel 
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
nb_params_ob <- set_nb_params(nb_params_ob, nsim_val = 20, libsize_val = 12*10^6)
# add noise to the original cell type specific profiles - create perturbed references 
for(k in sim1_params_ob@celltype){
  # set the output name of signature profiles 
  # note that negative binomial simulation take the countNorm units and output count units 
  input_signature_name <- paste(k, "countNorm", sep = "_")
  output_simulation_name <- paste(k, "count", "nb", sep = "_")
  
  assign(output_simulation_name, wrap_gammaPoisson(sim_param = sim1_params_ob, 
                                                   nb_param = nb_params_ob,
                                                   expr_list = get(input_signature_name),
                                                   interval = 1.8))
}
# in silico mixing with perturbed references 
mix_name <- paste("sim1","M", "nb","count", sep = "_")
assign(mix_name, sim_mix(sim_param = sim1_params_ob,
                          W = sim1_W,
                          ref_list = list(),
                          model_name = "nb"))
# unit transform from count to countNorm, cpm and tpm 
for(u in 2:length(unit)){
  mix_name_out <- paste("sim1", "M", "nb", unit[u], sep = "_")
  assign(mix_name_out, sim_mix_unit_transform(sim_param = sim1_params_ob,
                                           mix_list = sim1_M_nb_count,
                                           unit_val = unit[u],
                                           eff_length = sim1_eff_length))
  
}
#------------ simulation with log-normal distribution
# note that log-normal simulation first perform mixing process than add the noise 
for(u in 1:length(unit)){
  ref_list <- paste("sim1", "ref", "median", unit[u], sep = "_")
  mix_name <- paste("sim1", "M", "lognormal", unit[u], sep = "_")
  assign(mix_name, sim_mix(sim_param = sim1_params_ob,
                      ref_list = get(ref_list),
                      W = sim1_W,
                      model_name = "lognormal"))

}
#------------ simulation with normal distribution
for(u in 1:length(unit)){
  ref_list <- paste("sim1", "ref", "median", unit[u], sep = "_")
  mix_name <- paste("sim1", "M", "normal", unit[u], sep = "_")
  assign(mix_name, sim_mix(sim_param = sim1_params_ob,
                      ref_list = get(ref_list),
                      W = sim1_W,
                      model_name = "normal"))
  
}  
save.image(file = "./output/Sim1/simModel/sim1_mix_simModel.RData")