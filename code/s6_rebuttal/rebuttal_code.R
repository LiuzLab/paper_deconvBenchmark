# To my knowledge, the statement on lines 179-182 is actually due to the fact that smaller proportions are more difficult to predict accurately.
# Two new scenarios must be evaluated:
# a) cases with 5 to 10 components where n-1 cell types have all small proportions and only one (always the same cell type) is big; 
# b) cases with 5 to 10 components with equal proportions across all cell types.

load("./data/sim2_env.RData")
#####------------------------ I: source functions & load dataset ------------------------#####
# load processed data for simulation 
source("./code/src/Classes_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Functions_all.R")
source("./code/src/Methods_simulation.R")
# ------------------------ simulate weight 
n_comp <- sim2_params_ob@n_comp
nsim <- sim2_params_ob@nsim
# weight a 
weight_a <- list() 
for(i in 1:length(n_comp)){
  weight_a[[i]] <- matrix(data = 0, nrow = n_comp[i], ncol = 20)
  for(j in 1:nsim){
    major_comp <- runif(1, min = 0.9, max = 0.99)
    m_n <- n_comp[i] - 1
    minor_comp <- runif(n = m_n, min = ((1-major_comp)/m_n), max = ((1-major_comp)/m_n+0.01))
    tmp <- c(major_comp, minor_comp)
    # eliminate non-negative number 
    v <- tmp/sum(tmp)
    weight_a[[i]][,j] <- v
  }
  rownames(weight_a[[i]]) <- sim2_cellType[[i]]
}
# weight b: 
weight_b <- list() 
for(i in 1:length(n_comp)){
  weight_b[[i]] <- matrix(data = 0, nrow = n_comp[i], ncol = 20)
  for(j in 1:nsim){
    tmp <- runif(n = n_comp[i], min = (1/n_comp[i]), max = (1/n_comp[i]+0.04))
    v <- tmp/sum(tmp)
    weight_b[[i]][,j] <- v
  }
  rownames(weight_b[[i]]) <- sim2_cellType[[i]]
}

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
sim2_M_count_a <- sim_mix(sim_param = sim2_params_ob, W = weight_a, sim_celltype = sim2_cellType, sim_name = "sim2")
sim2_M_count_b <- sim_mix(sim_param = sim2_params_ob,  W = weight_b, sim_celltype = sim2_cellType, sim_name = "sim2")
# unit transform: from count to countNorm, cpm and tpm
eff_length <- sim2_eff_length_mean
unit <- sim2_params_ob@unit
for(u in 2:length(unit)){
  input_name1 <- paste("sim2", "M","count", "a", sep = "_")
  input_name2 <- paste("sim2", "M", "count", "b", sep = "_")
  output_name1 <- paste("sim2", "M", unit[u], "a", sep = "_")
  output_name2 <- paste("sim2", "M", unit[u], "b", sep = "_")
  assign(output_name1, sim_mix_unit_transform(sim_param = sim2_params_ob, 
                                              mix_list = get(input_name1), 
                                              unit_val = unit[u], 
                                              eff_length = eff_length))
  assign(output_name2, sim_mix_unit_transform(sim_param = sim2_params_ob, 
                                              mix_list = get(input_name2), 
                                              unit_val = unit[u],
                                              eff_length = eff_length))
}
save.image("./output/rebuttal/rebuttal_mix.RData")