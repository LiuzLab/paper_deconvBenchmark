# Generate reference, reference annotation and marker genes for Sim1
load("./output/Sim1/sim1_params.RData")
source("./code/src/Functions_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Methods_simulation.R")
#------------------create reference profiles with different units for sim1 
# unit includes count, countNorm, cpm and tpm 
# sim1_ref_<unit> will be created
# sim1_ref_median_<unit> will be created
for(u in unit){
  ref_tmp <- create_ref(sim1_params_ob, u)
  profile_name_output1 <- paste("sim1", "ref", u, sep = "_")
  profile_name_output2 <- paste("sim1", "ref", "median", u, sep = "_")
  
  assign(profile_name_output1, ref_tmp[[1]])
  assign(profile_name_output2, ref_tmp[[2]])
}
#------------------create reference annotation 
sim1_ref_anno_cibersort <- ref_anno_cibersort(sim1_params_ob, sim1_ref_count)
sim1_ref_anno_others <- ref_anno_others(sim1_params_ob, sim1_ref_count)
#------------------create marker gene list - use ref_countNorm
marker_params_ob <- marker_params()
sim1_marker_countNorm <- list()
celltype <- sim1_params_ob@celltype
for(i in 1:length(sim1_params_ob@dataset_name)){
  marker_params_ob@max_val <- quantile(as.matrix(sim1_ref_countNorm[[i]]), 0.8)
  marker_params_ob@min_val <- quantile(as.matrix(sim1_ref_countNorm[[i]]), 0.5)
  
  sim1_marker_countNorm[[i]] <- marker_gene_from_ref(sim_param = sim1_params_ob,
                                                    marker_param = marker_params_ob,
                                                    ref = sim1_ref_countNorm[[i]],
                                                    celltype = celltype)
}
#------------------create signature gene list - use ref_count
sim1_sigGene_count <- list()
for( i in 1:length(sim1_params_ob@dataset_name)){
  sim1_sigGene_count[[i]] <- sig_gene_from_ref(sim_param = sim1_params_ob,
                                               ref = sim1_ref_count[[i]],
                                               alpha = 0.01,
                                               nFold = 10,
                                               meta = sim1_ref_anno_others[[i]])
}
#------------------create ref list for epic 
for(u in 1:length(unit)){
  input_name1 = paste("sim1", "ref",unit[u], sep = "_")
  input_name2 = paste("sim1", "ref", "median", unit[u], sep = "_")
  output_name = paste("sim1", "epic", "ref", unit[u], sep = "_")
  
  assign(output_name, create_epic_ref(sim_param = sim1_params_ob,
                                      sig_gene = sim1_sigGene_count,
                                      ref = get(input_name1),
                                      median_ref = get(input_name2)))

  }
save.image("./output/Sim1/sim1_ref.RData")