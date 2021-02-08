# Parameter set for Sim2
# Load libraries and scripts 
source("./code/src/Functions_all.R")
source("./code/src/Classes_all.R")
source("./code/src/Generics_all.R")
source("./code/src/Methods_params.R")
source("./code/src/Methods_simulation.R")
load("./data/GSE51984_dat_raw.RData")
load("./data/GSE60424_dat_raw.RData")
load("./data/GSE64655_dat_raw.RData")
load("./data/GSE115736_dat_raw.RData")
#--------------------------------set a sim2 class object 
# extract count data 
sim2_params_ob <- sim2_params()
# extract consensus ids across 5 datasets
sim2_params_ob@all_id <- set_gene_id(sim2_params_ob, "_dat_raw")
all_id <- sim2_params_ob@all_id
dataset_name <- sim2_params_ob@dataset_name
# reorder the rows of 5 datasets
sim2_all_expr_count <- dataset_concat(dataset_name, suffix = "_dat_raw", unit = "count", id = all_id)
# 342 samples with 57189 (combination of samples derive from 4 datasets)
# perform filtering 
sim2_all_expr_count <- dat.filter(sim2_all_expr_count,1,10) # 22508
sim2_params_ob@gene_id <- rownames(sim2_all_expr_count)
#--------------------------------set additional feature variable for simulation 2 
# set the eff_length (it will take a while to finish)
# returns a 2-column matrix, the first is the length and the second column is the effective length
sim2_eff_length <- set_eff_length(sim2_params_ob, path = "./Raw")
sim2_eff_length_mean <- Reduce("+", sim2_eff_length)/length(sim2_eff_length)
# third extract sample annotations from each dataset  
sim2_anno_inf <- get_celltype_anno(sim2_params_ob, suffix = "_dat_raw")
# you'll find each dataset uses slightly different names for cell types 
# we need to replace these names to a consensus version 
# ----------------------------- construct a consensus annotation for all datasets 
# cell type annotation for GSE51984
anno <- GSE51984_dat_raw[['meta']][,'cell_type']
anno <- gsub("T-cells","T cells", anno)
anno <- gsub("B-cells","B cells", anno) 
anno <- gsub("Total white blood cells","PBMC", anno) 
names(anno) <- GSE51984_dat_raw[['meta']][,'Run']
GSE51984_dat_raw[['anno']] <- anno
# Cell type annotation for GSE64655
anno <- GSE64655_dat_raw[['meta']][,'cell_type']
anno <- gsub("primary human ","", anno)
anno <- gsub("monocytes","Monocytes", anno)
anno <- gsub("neutrophils","Neutrophils", anno)
anno <- gsub("myeloid DC","Myeloid DC", anno)
names(anno) <- GSE64655_dat_raw[['meta']][,'Run']
GSE64655_dat_raw[['anno']] <- anno
# Cell type annotation for GSE60424
anno <- GSE60424_dat_raw[['meta']][,'cell_type']
anno <- gsub("B-cells","B cells", anno)
anno <- gsub("CD4","CD4 T cells", anno)
anno <- gsub("CD8","CD8 T cells", anno)
anno <- gsub("NK","NK cells", anno)
anno <- gsub("Whole Blood","PBMC", anno)
names(anno) <- GSE60424_dat_raw[['meta']][,'Run']
GSE60424_dat_raw[['anno']] <- anno
# Cell type annotation for GSE115736
anno <- GSE115736_dat_raw[['meta']][,'cell_type']
anno <- gsub("Cell","cells", anno)
anno <- gsub("Monocyte","Monocytes", anno)
anno <- gsub("monocyte","Monocytes", anno)
anno <- gsub("Neutrophil","Neutrophils", anno)
anno <- gsub("Eosinophil","Eosinophils", anno)
anno <- gsub("Myeloid DC CD123+","Myeloid DC", anno,fixed = TRUE)
anno <- gsub(":+","", anno, fixed = TRUE)
names(anno) <- GSE115736_dat_raw[['meta']][,'Run']
GSE115736_dat_raw[['anno']] <- anno
#----------------------------- store all sample-wise annotation to a list
sim2_sampleAnno <- list()
for(i in 1:length(dataset_name)){
  anno_dat <- get(paste0(dataset_name[i], "_dat_raw"))[['anno']]
  sim2_sampleAnno[[i]] <- anno_dat
}

#----------------- create celltype annotation 
sim2_cellType <- list()
sim2_cellType[[1]] <- c("T", "B", "Mono","Neutro", "NK")
sim2_cellType[[2]] <- c("T", "B", "Mono","Neutro","NK", "Eosino")
sim2_cellType[[3]] <- c("T", "B", "Mono","Neutro","NK", "Eosino","MyeloidDC")
sim2_cellType[[4]] <-  c("T", "B", "Mono","Neutro","NK", "Eosino","MyeloidDC","HSC")
sim2_cellType[[5]] <- c("CD4T","CD8T", "B", "Mono","Neutro","NK", "Eosino","MyeloidDC","HSC")
sim2_cellType[[6]] <- c("CD4T","CD8T", "NaiveB", "MemoryB","Mono","Neutro","NK", "Eosino","MyeloidDC","HSC")
#----------------- create list for extract pattern 
sim2_extractPattern <- list() 
sim2_extractPattern[[1]] <- c("T cells", "B cells", "Monocytes", "Neutrophils", "NK cells")
sim2_extractPattern[[2]] <- c("T cells", "B cells", "Monocytes", "Neutrophils", "NK cells", "Eosinophils")
sim2_extractPattern[[3]] <- c("T cells", "B cells", "Monocytes", "Neutrophils", "NK cells", "Eosinophils", "Myeloid DC")
sim2_extractPattern[[4]] <- c("T cells", "B cells", "Monocytes", "Neutrophils", "NK cells", "Eosinophils", "Myeloid DC", "CD34+")
sim2_extractPattern[[5]] <- c("CD4 T cells", "CD8 T cells", "B cells", "Monocytes", "Neutrophils", "NK cells", "Eosinophils", "Myeloid DC", "CD34+")
sim2_extractPattern[[6]] <- c("CD4 T cells", "CD8 T cells", "Naive B", "Memory B", "Monocytes", "Neutrophils", "NK cells", "Eosinophils", "Myeloid DC", "CD34+")
#---------------- extract pattern immune ----- match the cell type to the etract pattern 
all_patterns <- unique(unlist(sim2_extractPattern))
names(all_patterns) <- unique(unlist(sim2_cellType))
sim2_params_ob@extract_pattern_immune <- all_patterns 
# component range 
# apply simulation based on hierarchical probability model 
# Lymphocytes: 14 - 47
#+ T: 7 - 24
#++ CD4T: 4 - 20 *We use the same proportion as T to investigate impact of reference collineraty without chaning the range (CD4T 7 - 24 in actual simulation)
#++ CD8T: 2 - 11
#+ B: 1 - 7
#++ Naive B: 0.7 - 4.9
#++ Memory B: 0.2 - 1.7
#+ NK: 1 - 6
# Myloid cells: 53 - 86
#+ Granulocytes: 35 - 80
#++ Neutrophils: 30 - 80
#++ Eosinophils: 0 - 7 
#+ Monocytes: 2 - 12
#+ DC: 0.3 - 0.9 *We used myeloidDC and use proportion of DC to compensate other missing cell types in the blood (myeloidDC 0.3 - 0.9 in actual simulation)
# HSC: 0.03 - 0.06
sim2_prop_mat <- list()
sim2_prop_mat[[1]] <- rbind(T = c(7, 24), B = c(1, 7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6))
sim2_prop_mat[[2]] <- rbind(T = c(7, 24), B = c(1, 7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6), Eosino = c(0,7))
sim2_prop_mat[[3]] <- rbind(T = c(7, 24), B = c(1, 7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6), Eosino = c(0,7), MyeloidDC = c(0.3, 0.9))
sim2_prop_mat[[4]] <- rbind(T = c(7, 24), B = c(1, 7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6), Eosino = c(0,7), MyeloidDC = c(0.3, 0.9), HSC = c(0.03, 0.06))
sim2_prop_mat[[5]] <- rbind(CD4T = c(7, 24), CD8T = c(2, 11), B = c(1, 7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6), Eosino = c(0,7), MyeloidDC = c(0.3, 0.9), HSC = c(0.03, 0.06))
sim2_prop_mat[[6]] <- rbind(CD4T = c(7, 24), CD8T = c(2, 11), NaiveB = c(0.7, 4.9), MemoryB = c(0.2, 1.7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6), Eosino = c(0,7), MyeloidDC = c(0.3, 0.9), HSC = c(0.03, 0.06))
#----------------------- create celltype specific countNorm expr for nb simulation 
# normalize count to countNorm (by median library size)
sim2_all_expr_countNorm <- apply(sim2_all_expr_count, 2, median_libSize_norm, median(colSums(sim2_all_expr_count)))
colnames(sim2_all_expr_countNorm) <- gsub("count", "countNorm", colnames(sim2_all_expr_countNorm))
# create cell type specific expr profiles 
extract_pattern_immune <- sim2_params_ob@extract_pattern_immune
for(i in 1:length(extract_pattern_immune)){
  celltype_ref_name <- paste("sim2", names(extract_pattern_immune)[i], "countNorm", sep = "_")
  
  assign(celltype_ref_name, extract_celltype_profiles(sim_param = sim2_params_ob,
                                                      unit_val = "countNorm",
                                                      pattern_val = extract_pattern_immune[i],
                                                      sample_anno = sim2_sampleAnno,
                                                      dat = sim2_all_expr_countNorm))
}
#----------------------- generate dataset and celltype annotations 
# note that this annotation has nothing to do with the concatenated expression matrix 
celltype_anno <- unlist(sim2_sampleAnno)
dataset_anno <- vector()
for(i in 1:length(sim2_sampleAnno)){
  anno <- rep(dataset_name[i], length(sim2_sampleAnno[[i]]))
  names(anno) <- names(sim2_sampleAnno[[i]])
  dataset_anno <- c(dataset_anno, anno)
}
all.equal(names(dataset_anno), names(celltype_anno))
save.image(file = "./output/Sim2/sim2_env.RData")
