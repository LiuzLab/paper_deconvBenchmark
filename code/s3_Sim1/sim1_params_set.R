# Parameter set for Sim1
# Load libraries and scripts 
source("./code/src/Functions_all.R")
source("./code/src/Classes_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Methods_params.R")
load("./data/GSE51984_dat.RData")
load("./data/GSE60424_dat.RData")
load("./data/GSE64655_dat.RData")
#--------------------------------set a sim1 class object 
# initialize a sim1 class object
# most of the parameters have been set prior in the class definition
sim1_params_ob <- sim1_params()
# set gene_id parameter
# consensus gene sets for 3 datasets 
sim1_params_ob@gene_id <- set_gene_id(sim1_params_ob, "_dat")
length(sim1_params_ob@gene_id) # 15123 consensus genes been identified 
# set the extract_pattern 
# extract_pattern parameter will be used to extract samples derived from corresponding cell types from 3 datasets
sim1_params_ob@extract_pattern_immune <- c(T = "T-cells|CD4|CD8|T cells",B = "B-cells|B cells",  Mono = "Monocytes|monocytes")
#--------------------------------set additional feature variable for simulation 1 
# set the eff_length (it will take a while to finish)
# returns a 2-column matrix, the first is the length and the second column is the effective length
sim1_eff_length <- set_eff_length(sim1_params_ob, path = "/mnt/data/haijing/simDeconv/Raw/")
# set the celltype_anno 
sim1_celltype_anno <- get_celltype_anno(sim1_params_ob,"_dat")
#--------------------------------extract all cell type specific profiles for simulation 
extract_pattern_immune <- sim1_params_ob@extract_pattern_immune
celltype <- sim1_params_ob@celltype
# check if the squence of celltype and extract pattern are the same 
names(extract_pattern_immune) == celltype
unit <- sim1_params_ob@unit
# creates lists of cell type specific profiles 
# T_count, B_count and Mono_count (applies same for other units countNorm, cpm and tpm)
# Each list consists of 3 expression matrices, each matrix is derived from one dataset 
# The order of datasets is the same with the order in the sim1_params_ob@dataset_name
for(k in celltype){
  # define the cell type 
  for(u in unit){
    profile_name = paste(k, u, sep = "_")
    if(grepl(pattern = "countNorm", profile_name)){
      assign(profile_name, extract_celltype_profiles(sim1_params_ob,
                                                     pattern_val = extract_pattern_immune[k],
                                                     unit_val = "count",
                                                     medianNorm = TRUE))
    }else{
      assign(profile_name, extract_celltype_profiles(sim1_params_ob,
                                                     pattern_val = extract_pattern_immune[k],
                                                     unit_val = u,
                                                     medianNorm = FALSE))
    }
  }
}
# 6. save the environment for the simulation 1. 
save.image(file = "./output/Sim1/sim1_env.RData")