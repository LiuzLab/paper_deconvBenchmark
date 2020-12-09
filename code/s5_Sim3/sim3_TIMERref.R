# TIMER ref gene extraction 
load("./output/Sim3/sim3_env.RData")
source("./code/src/TIMER_codes.R")
# extract HCT and PBMC profiles
extract_pattern_TIMER <- c("PBMC", "HCT")
names(extract_pattern_TIMER) <- c("PBMC", "HCT")
for(i in 1:length(extract_pattern_TIMER)){
  celltype_ref_name <- paste("sim3", names(extract_pattern_TIMER)[i], "countNorm", sep = "_")
  
  assign(celltype_ref_name, extract_celltype_profiles(sim_param = sim3_params_ob,
                                                      unit_val = "countNorm",
                                                      pattern_val = extract_pattern_TIMER[i],
                                                      sample_anno = sim3_sampleAnno,
                                                      dat = sim3_all_expr_countNorm))
}
# build HCT mixtures with a grid HCT proportions from 0% to 100% 
# each level with 3 replicates 

hct_prop <- seq(0, 1, length.out = 11)
n_PBMC <- length(sim3_PBMC_countNorm)
n_HCT <- length(sim3_HCT_countNorm)

tmp_outer <- list()
for(i in 1:length(hct_prop)){
  # sample 3 profiles each time 
  id_PBMC <- sample(1:n_PBMC, 3)
  id_HCT <- sample(1:n_HCT, 3)
  tmp_inner <- list() 
  for(j in 1:3){
    tmp_inner[[j]] <- hct_prop[i]*sim3_HCT_countNorm[,id_HCT[j]] + (1-hct_prop[i])*sim3_PBMC_countNorm[,id_PBMC[j]]
  }
  tmp_outer[[i]] <- do.call(cbind, tmp_inner)
  colnames(tmp_outer[[i]]) <- paste("HCT", hct_prop[i], 1:3, sep = "_")
}
hct_timer_mix <- do.call(cbind, tmp_outer)
rownames(hct_timer_mix) <- rownames(sim3_HCT_countNorm)
# get purity gene 
timer_purity_genes <- getPurityGenes(dd = hct_timer_mix, purity = rep(hct_prop, each = 3))

# 
saveRDS(timer_purity_genes, file = "./output/Sim3/timer_purity_genes.rds")

