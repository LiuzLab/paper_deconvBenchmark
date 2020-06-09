# Parameter set for sim3
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
load("./data/GSE118490_dat_raw.RData")
#--------------------------------set a sim3 class object 
# extract count data 
sim3_params_ob <- sim3_params()
# extract consensus ids across 5 datasets
sim3_params_ob@all_id <- set_gene_id(sim3_params_ob, "_dat_raw")
all_id <- sim3_params_ob@all_id
dataset_name <- sim3_params_ob@dataset_name
# reorder the rows of 5 datasets
sim3_all_expr_count <- dataset_concat(dataset_name, suffix = "_dat_raw", unit = "count", id = all_id)
# 345 samples with 57189 genes (combination of samples derive from 5 datasets)
# perform filtering 
sim3_all_expr_count <- dat.filter(sim3_all_expr_count,1,10) # 23324
sim3_params_ob@gene_id <- rownames(sim3_all_expr_count)
#--------------------------------set additional feature variable for simulation 2 
# set the eff_length (it will take a while to finish)
# returns a 2-column matrix, the first is the length and the second column is the effective length
sim3_eff_length <- set_eff_length(sim3_params_ob, path = "/mnt/data/haijing/simDeconv/Raw/")
sim3_eff_length_mean <- Reduce("+", sim3_eff_length)/length(sim3_eff_length)
# third extract sample annotations from each dataset  
sim3_anno_inf <- get_celltype_anno(sim3_params_ob, suffix = "_dat_raw")
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
# Cell type annotation for GSE118490
anno <- as.character(GSE118490_dat_raw[['meta']][,'cell_type'])
names(anno) <- GSE118490_dat_raw[['meta']][,'Run']
GSE118490_dat_raw[['anno']] <- anno
#----------------------------- store all sample-wise annotation to a list
# this annotation is set for the purpose of extracting samples from targeted cell types 
sim3_sampleAnno <- list() #! fix this using dataset_name 
for(i in 1:length(dataset_name)){
  anno_dat <- get(paste0(dataset_name[i], "_dat_raw"))[['anno']]
  sim3_sampleAnno[[i]] <- anno_dat
}
#----------------- create celltype annotation 
# this annotation is set for simplifying the naming for  data manipulation and compacted figure legends (a shorter version of _sampleAnno)
sim3_immune_cellType <- list()
sim3_immune_cellType[[1]] <- c("T", "B", "Mono","Neutro", "NK")
sim3_immune_cellType[[2]] <- c("T", "B", "Mono","Neutro","NK", "Eosino")
sim3_immune_cellType[[3]] <- c("T", "B", "Mono","Neutro","NK", "Eosino","MyeloidDC")
sim3_immune_cellType[[4]] <-  c("T", "B", "Mono","Neutro","NK", "Eosino","MyeloidDC","HSC")
sim3_immune_cellType[[5]] <- c("CD4T","CD8T", "B", "Mono","Neutro","NK", "Eosino","MyeloidDC","HSC")
sim3_immune_cellType[[6]] <- c("CD4T","CD8T", "NaiveB", "MemoryB","Mono","Neutro","NK", "Eosino","MyeloidDC","HSC")

sim3_cellType <- list()
sim3_cellType[[1]] <- c("T", "B", "Mono","Neutro", "NK", "HCT")
sim3_cellType[[2]] <- c("T", "B", "Mono","Neutro","NK", "Eosino", "HCT")
sim3_cellType[[3]] <- c("T", "B", "Mono","Neutro","NK", "Eosino","MyeloidDC", "HCT")
sim3_cellType[[4]] <-  c("T", "B", "Mono","Neutro","NK", "Eosino","MyeloidDC","HSC", "HCT")
sim3_cellType[[5]] <- c("CD4T","CD8T", "B", "Mono","Neutro","NK", "Eosino","MyeloidDC","HSC", "HCT")
sim3_cellType[[6]] <- c("CD4T","CD8T", "NaiveB", "MemoryB","Mono","Neutro","NK", "Eosino","MyeloidDC","HSC", "HCT")
#----------------- create list for extract pattern 
sim3_extractPattern <- list() 
sim3_extractPattern[[1]] <- c("T cells", "B cells", "Monocytes", "Neutrophils", "NK cells", "HCT")
sim3_extractPattern[[2]] <- c("T cells", "B cells", "Monocytes", "Neutrophils", "NK cells", "Eosinophils", "HCT")
sim3_extractPattern[[3]] <- c("T cells", "B cells", "Monocytes", "Neutrophils", "NK cells", "Eosinophils", "Myeloid DC", "HCT")
sim3_extractPattern[[4]] <- c("T cells", "B cells", "Monocytes", "Neutrophils", "NK cells", "Eosinophils", "Myeloid DC", "CD34+", "HCT")
sim3_extractPattern[[5]] <- c("CD4 T cells", "CD8 T cells", "B cells", "Monocytes", "Neutrophils", "NK cells", "Eosinophils", "Myeloid DC", "CD34+", "HCT")
sim3_extractPattern[[6]] <- c("CD4 T cells", "CD8 T cells", "Naive B", "Memory B", "Monocytes", "Neutrophils", "NK cells", "Eosinophils", "Myeloid DC", "CD34+", "HCT")
all_patterns <- unique(unlist(sim3_extractPattern))
names(all_patterns) <- unique(unlist(sim3_cellType))
sim3_params_ob@extract_pattern <- all_patterns 
sim3_params_ob@extract_pattern_immune <- all_patterns[-which(all_patterns == 'HCT')]
# component range 
# apply simulation based on hierarchical probability model 
# Lymphocytes: 14 - 47
#+ T: 7 - 24
#++ CD4T: 4 - 20
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
#+ DC: 0.3 - 0.9
# HSC: 0.03 - 0.06
sim3_prop_mat <- list()
sim3_prop_mat[[1]] <- rbind(T = c(7, 24), B = c(1, 7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6))
sim3_prop_mat[[2]] <- rbind(T = c(7, 24), B = c(1, 7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6), Eosino = c(0,7))
sim3_prop_mat[[3]] <- rbind(T = c(7, 24), B = c(1, 7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6), Eosino = c(0,7), MyeloidDC = c(0.3, 0.9))
sim3_prop_mat[[4]] <- rbind(T = c(7, 24), B = c(1, 7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6), Eosino = c(0,7), MyeloidDC = c(0.3, 0.9), HSC = c(0.03, 0.06))
sim3_prop_mat[[5]] <- rbind(CD4T = c(7, 24), CD8T = c(2, 11), B = c(1, 7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6), Eosino = c(0,7), MyeloidDC = c(0.3, 0.9), HSC = c(0.03, 0.06))
sim3_prop_mat[[6]] <- rbind(CD4T = c(7, 24), CD8T = c(2, 11), NaiveB = c(0.7, 4.9), MemoryB = c(0.2, 1.7), Mono = c(2, 12), Neutro = c(30, 80), NK = c(1, 6), Eosino = c(0,7), MyeloidDC = c(0.3, 0.9), HSC = c(0.03, 0.06))
#----------------------- create celltype specific countNorm expr for nb simulation 
# normalize count to countNorm (by median library size)
sim3_all_expr_countNorm <- apply(sim3_all_expr_count, 2, median_libSize_norm, median(colSums(sim3_all_expr_count)))
colnames(sim3_all_expr_countNorm) <- gsub("count", "countNorm", colnames(sim3_all_expr_countNorm))
# create cell type specific expr profiles 
extract_pattern <- sim3_params_ob@extract_pattern
for(i in 1:length(extract_pattern)){
  celltype_ref_name <- paste("sim3", names(extract_pattern)[i], "countNorm", sep = "_")
  
  assign(celltype_ref_name, extract_celltype_profiles(sim_param = sim3_params_ob,
                                                      unit_val = "countNorm",
                                                      pattern_val = extract_pattern[i],
                                                      sample_anno = sim3_sampleAnno,
                                                      dat = sim3_all_expr_countNorm))
}
# generate dataset annotation 
#----------------------- create celltype specific countNorm expr for nb simulation 
celltype_anno <- unlist(sim3_sampleAnno)
dataset_anno <- vector()
for(i in 1:length(sim3_sampleAnno)){
  anno <- rep(dataset_name[i], length(sim3_sampleAnno[[i]]))
  names(anno) <- names(sim3_sampleAnno[[i]])
  dataset_anno <- c(dataset_anno, anno)
}
all.equal(names(dataset_anno), names(celltype_anno))
save.image(file = "./output/Sim3/sim3_env.RData")