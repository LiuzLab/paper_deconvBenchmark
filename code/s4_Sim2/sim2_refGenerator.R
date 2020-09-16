# Generate reference, reference annotation and marker genes for Sim2
# We recommend to run this script in background as it takes a while to finish. 
setwd("/abs_path/") # Please set the path to the absolute path of your working directory. 
load("./output/Sim2/sim2_env.RData")
source("./code/src/Classes_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Methods_params.R")
source("./code/src/Methods_simulation.R")
#------------------ create reference profiles with different units for sim2 
# unit includes count, countNorm, cpm and tpm 
# Step 1: create sim2_ref_raw_<unit> profiles first, note this list includes all samples from 4 datasets, and samples are sorted according to the dataset
# Step 2: extract cell type specific profiles from sim2_ref_raw_<unit> so that samples are sorted according to the cell type identity
# Step 3: sim2_ref_<unit> and sim2_ref_median_<unit> will be created based on the cell type specific profiles created in the Step2
unit <- sim2_params_ob@unit
gene_id <- sim2_params_ob@gene_id
for(u in unit){
  all_ref_name <- paste("sim2", "ref_raw", u, sep = "_")
  if(u == "countNorm"){
    tmp <- cbind(GSE51984_dat_raw[['count']][gene_id,],
                 GSE64655_dat_raw[['count']][gene_id,],
                 GSE60424_dat_raw[['count']][gene_id,],
                 GSE115736_dat_raw[['count']][gene_id,])
    tmp <- apply(tmp,2,median_libSize_norm, median(colSums(tmp))) 
    colnames(tmp) <- gsub("count", "countNorm", colnames(tmp))
  }
  else{
    tmp <- cbind(GSE51984_dat_raw[[u]][gene_id,],
                 GSE64655_dat_raw[[u]][gene_id,],
                 GSE60424_dat_raw[[u]][gene_id,],
                 GSE115736_dat_raw[[u]][gene_id,])
  }
  assign(all_ref_name, tmp)
  }
# cell type specific profiles will be used to create reference for deconvolution where samples are sorted by the cell type identity
extract_pattern_immune <- sim2_params_ob@extract_pattern_immune
for(u in unit){
  all_ref_name <- paste("sim2", "ref_raw", u, sep = "_")
  for(i in 1:length(extract_pattern_immune)){
    celltype_ref_name <- paste("sim2", "ref", names(extract_pattern_immune)[i], u, sep = "_")
    
    assign(celltype_ref_name, extract_celltype_profiles(sim_param = sim2_params_ob,
                                                        unit_val = u,
                                                        pattern_val = extract_pattern_immune[i],
                                                        sample_anno = sim2_sampleAnno,
                                                        dat = get(all_ref_name)))
  }
}
# create ref and ref_median profiles for deconvolution 
for(u in unit){
  ref_tmp <- create_ref(sim2_params_ob, unit_val = u, celltype_list = sim2_cellType, sim_name = "sim2_ref")
  profile_name_output1 <- paste("sim2", "ref", u, sep = "_")
  profile_name_output2 <- paste("sim2", "ref", "median", u, sep = "_")
  
  assign(profile_name_output1, ref_tmp[[1]])
  assign(profile_name_output2, ref_tmp[[2]])
}
#------------------ create reference annotation 
sim2_ref_anno_cibersort <- ref_anno_cibersort(sim2_params_ob, sim2_ref_count)
sim2_ref_anno_others <- ref_anno_others(sim2_params_ob, sim2_ref_count)

#------------------ marker gene extraction using countNorm data 
marker_params_ob <- marker_params()
sim2_marker_countNorm <- list()
for(i in 1:length(sim2_ref_countNorm)){
  marker_params_ob@max_val <- quantile(as.matrix(sim2_ref_countNorm[[i]]), 0.8)
  marker_params_ob@min_val <- quantile(as.matrix(sim2_ref_countNorm[[i]]), 0.5)
  # modify the marker_gene_from_ref 
  sim2_marker_countNorm[[i]] <- marker_gene_from_ref(sim_param = sim2_params_ob,
                                                     marker_param = marker_params_ob,
                                                     ref = sim2_ref_countNorm[[i]],
                                                     celltype = sim2_cellType[[i]])
}
#------------------  signature gene extraction using count data
sim2_sigGene_count <- list()
for( i in 1:length(sim2_ref_countNorm)){
  sim2_sigGene_count[[i]] <- sig_gene_from_ref(sim_param = sim2_params_ob,
                                               ref = sim2_ref_count[[i]],
                                               alpha = 0.01,
                                               nFold = 10,
                                               meta = sim2_ref_anno_others[[i]])
}
#------------------ signature gene for epic data 
for(u in 1:length(unit)){
  input_name1 = paste("sim2", "ref",unit[u], sep = "_")
  input_name2 = paste("sim2", "ref", "median", unit[u], sep = "_")
  output_name = paste("sim2", "epic", "ref", unit[u], sep = "_")
  
  
  assign(output_name, create_epic_ref(sim_param = sim2_params_ob,
                                      sig_gene = sim2_sigGene_count,
                                      ref = get(input_name1),
                                      median_ref = get(input_name2)))
  
}
save.image("./output/Sim2/sim2_ref.RData")