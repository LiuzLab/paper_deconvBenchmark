#@ set the consensus gene id across all units and datasets for the simulation 
setMethod(f = "set_gene_id",
          signature = c("sim_params", "character"),
          definition = function(sim_param, suffix)
          #@ sim_param: sim_params objects
          #@ suffix: suffix of data object after GEO query id (i.e. _dat or _dat_raw)
          #@ output: consensus gene ids across datasets and units 
          {
            unit = sim_param@dat_unit 
            dataset_name = sim_param@dataset_name
            
            id_list = list()
            for(u in 1:length(unit)){
              for(i in 1:length(dataset_name)){
                tmp_dat = get(paste0(dataset_name[i],suffix))[[unit[u]]] # get the gene expression matrix of dataset[[i]] and unit[u]
                id_list[[ length(dataset_name)*(u-1) + i]] = rownames(tmp_dat)
              }
            }
            consensus_id = Reduce(intersect, id_list)
            return(consensus_id)
          }
)

#@ set the effective length of across datasets 
setMethod(f = "set_eff_length",
          signature = "sim_params",
          definition = function(sim_param, path)
            #@ sim_param: sim_params object 
            #@ path: path of raw files that store the length data (RSEM results)
            #@ output: matrix, first column - length, second column - effective length, rownames - gene symbol 
          {
            eff_length <- list() 
            
            dataset_name = sim_param@dataset_name
            gene_id = sim_param@gene_id
            id_name_table <- id.name.db(version = 95)
            
            rownames(id_name_table) <- id_name_table[,1]
            for(i in 1:length(dataset_name)){
              file_path <- c(paste0(path, dataset_name[i], "/RSEM_quant/"))
              file <- list.files(path = file_path, pattern = "genes.results")[1]
              tmp_length <- read.table(paste0(file_path, file), sep="\t", header=TRUE)[,c(1,3,4)] # extract column of ensembl id, length and effective length
              tmp_length <- eliminate.duplicate(tmp_length)
              # convert ensembl ID to gene symbol 
              eff_length[[i]] <- id.to.name(tmp_length, id_name_table = id_name_table)[gene_id,]
            }
            return(eff_length)
          }
)

#@ extract levels of cell type annotation 
setMethod(f = "get_celltype_anno",
          signature = "sim_params",
          definition = function(sim_param, suffix)
            #@ sim_params obejct 
            #@ output: levels of cell type annotation across all dataset 
          { 
            dataset_name <- sim_param@dataset_name 
            
            celltype_anno_val <- list() 
            for(i in 1:length(dataset_name)){
              celltype_anno_val[[i]] <- unique(get(paste0(dataset_name[i],suffix))[['meta']][,'cell_type'])
            }
            
            
            return(celltype_anno_val)
            }
)

#@ extract cell type specific profiles and store them in a list 
#@ this Method is specific for sim1 
setMethod( f = "extract_celltype_profiles",
           signature = c("sim1_params", "character", "character"),
           definition = function(sim_param, pattern_val, unit_val, medianNorm){
             #@ sim_param: sim_params obejct 
             #@ pattern_val: pattern val corresponds to the targeted cell type 
             #@ unit_val: unit of expression matrix 
             #@ medianNorm: TRUE if the unit is countNorm 
             #@ output: cell type specific profile for targeted cell type and unit of 3 datasets 
             dataset_name = sim_param@dataset_name
             gene_id = sim_param@gene_id
             
             celltype_expr <- list() 
             for(i in 1:length(dataset_name)){
               # pick the right sample 
               dataset_list <- get(paste0(dataset_name[i],"_dat"))
               dat_tmp <- dataset_list[[unit_val]][gene_id,]
               # normalize data by the median library size if medianNorm = TRUE
               if(medianNorm == TRUE){
                 median_val <- median(colSums(dat_tmp))
                 dat_expr <- apply(dat_tmp, 2, median_libSize_norm, median_val)
               }else{
                 dat_expr <- dat_tmp
               }
               
               # select sample that belongs to the targeted cell type 
               sample_selector <- grep(pattern_val, dataset_list[['meta']][,'cell_type'])
               # select the right sample 
               if(length(sample_selector) == 0){
                 celltype_expr[[i]] <- NA
               }
               else{
                 sample_id <- paste0(dataset_list[['meta']][sample_selector,'Run'],"_",unit_val)
                 celltype_expr[[i]] <- dat_expr[,sample_id]  
               }
             }
             names(celltype_expr) <- dataset_name
             return(celltype_expr)
           }
)

#@ concatenate expr profiles derived from different datasets 
dataset_concat <- function(dataset_name, suffix, unit, id){
  #@ dataset_name: geo id
  #@ suffix: object name suffix
  #@ unit: unit of datastet: count 
  #@ id: consensus gene ids of profiles to be extracted  
  tmp <- data.frame(NA)
  for(i in 1:length(dataset_name)){
    dat_expr <- get(paste0(dataset_name[i], suffix))[[unit]][id, ]
    tmp <- cbind(tmp, dat_expr)
  }
  dat_concat <- tmp[,-1]
  return(dat_concat)
}

#@ extract cell type specific profiles and store them in a list 
#@ this Method is specific for sim2
setMethod( f = "extract_celltype_profiles",
           signature = c("comp_sim_params", "character", "character"),
           definition = function(sim_param, pattern_val, unit_val, sample_anno, dat){
             #@ sim_param: sim_parmas object 
             #@ pattern_val: patterns match to the meta data cell type annotation 
             #@ unit_val: value of the unit 
             #@ sample_anno: annotation of the sample 
             #@ dat: orignal concatenated expression matrix 
             
             dataset_name = sim_param@dataset_name
             gene_id = sim_param@gene_id
             
             expr <- data.frame(NA)
             
             for(i in 1:length(dataset_name)){
               # pick the right sample 
               # for each dataset, pick the right sample 
               anno <- sample_anno[[i]]
               run_id <- names(anno[grep(pattern_val, anno)])
               
               if(length(run_id) == 0){
                 expr == expr
               }
               else{
                 expr <- cbind(expr, dat[gene_id, paste(run_id, unit_val, sep = "_")])
               }
               
             }
             
             return(expr[,-1])
           }
           
)

