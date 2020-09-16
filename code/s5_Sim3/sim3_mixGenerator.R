# Parameter set for sim3
# Load libraries and scripts 
load("./output/Sim3/sim3_env.RData")
source("./code/src/Functions_all.R")
source("./code/src/Classes_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Methods_params.R")
source("./code/src/Methods_simulation.R")
# ------------------------ simulate weight 
# generate orthogonal weight matrix 
library(plyr)
# use the same set of relative proportions of immune cell type of sim2 
sim3_W_orthog <- readRDS("./output/Sim2/sim2_W_orthog.rds")
sim3_W_real <- readRDS("./output/Sim2/sim2_W_real.rds")
# generate tumor weight 
mix_type <- sim3_params_ob@mix_type
tumor_content <- sim3_params_ob@tumor_content
for(s in 1:length(mix_type)){
  known_weight <- get(paste("sim3", "W", mix_type[s], sep = "_"))
  for(j in 1:length(tumor_content)){
    weight_name <- paste("sim3", "W", mix_type[s], tumor_content[j], sep = "_")
    
    assign(weight_name, tumor_weight(known_weight = known_weight,
                                     tumorContent = tumor_content[j]))
  }
}
# ------------------------ create perturbed nb simulations 
# initiate nb simulation paramters 
nb_params_ob <- nb_params()
nb_params_ob <- set_nb_params(nb_params_ob, nsim_val = 20, libsize_val = 12*10^6)
celltype <- names(sim3_params_ob@extract_pattern)
for(k in 1:length(celltype)){
  nb_params_ob@prefix = celltype[k]
  input_name <- paste("sim3", celltype[k], "countNorm", sep = "_")
  output_name <- paste("sim3", celltype[k], "nb", "count", sep = "_")
  
  assign(output_name, sim_gammaPoisson(nb_param = nb_params_ob,
                                       reference = as.matrix(get(input_name))))
}
# 
#------------------------ create nb mixtures - withTumor 
# first create nb mixture in count unit 
for( s in 1:length(mix_type)){
  for(j in 1:length(tumor_content)){
    mix_name <- paste("sim3", "M", "count", mix_type[s], tumor_content[j], sep = "_")
    W_name <- paste("sim3", "W", mix_type[s], tumor_content[j], sep = "_")
    assign(mix_name, sim_mix(sim_param = sim3_params_ob,
                             W = get(W_name),
                             sim_celltype = sim3_cellType,
                             sim_name = "sim3"))
    
  }
}
# unit transform: from count to countNorm, cpm and tpm
eff_length <- sim3_eff_length_mean
unit <- sim3_params_ob@unit
for(u in 2:length(unit)){
  for(s in 1:length(mix_type)){
    for(j in 1:length(tumor_content)){
      input_name <- paste("sim3", "M", "count", mix_type[s], tumor_content[j], sep = "_")
      output_name <- paste("sim3", "M", unit[u], mix_type[s], tumor_content[j], sep = "_")
      assign(output_name, sim_mix_unit_transform(sim_param = sim3_params_ob,
                                                 mix_list = get(input_name),
                                                 unit_val = unit[u],
                                                 eff_length = eff_length))
      
    }
  }
}
save.image("./output/Sim3/sim3_mix.RData")