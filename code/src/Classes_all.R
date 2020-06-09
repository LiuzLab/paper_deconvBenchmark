# marked with * need to be defined later 
sim1_params <- setClass(
  # set the name for the class
  "sim1_params",
  
  # define the slots 
  slots = c(
    dataset_name = "character",
    dat_unit = "character",
    unit = "character",
    nsim = "numeric",
    grid = "numeric",
    celltype = "character",
    gene_id = "character", #* 
    extract_pattern_immune = "character", #*
    sim_model = "character",
    analysis_type = "character"
  ),
  
  prototype = list(
    dataset_name = c("GSE51984", "GSE64655", "GSE60424"),
    dat_unit = c("count", "cpm", "tpm"),
    unit = c("count", "countNorm", "cpm", "tpm"),
    nsim = 20,
    grid = seq(0, 0.9, by = 0.1),
    celltype = c("T", "B", "Mono"),
    sim_model = c("nb", "lognormal", "normal"),
    analysis_type = c("simModel", "libSize")
  )
)

# marked with * need to be defined later 
nb_params <- setClass(
  # set the name for the class
  "nb_params",
  
  # define the slots 
  slots = c(
    prefix = "character", 
    nsim = "numeric", #*
    bv = "numeric", #* 
    libsize = "numeric", #* 
    sigma = "numeric"
  ),
  
  prototype = list(
    prefix = "sim",
    nsim = NA_integer_,
    bv = NA_real_,
    libsize = NA_integer_,
    sigma = NA_real_
  )
)

# class for marker gene params 
marker_params <- setClass(
  # set the name for the class
  "marker_params",
  
  # define the slots 
  slots = c(
    max_val = "numeric", 
    min_val = "numeric",
    p_initial = "numeric",
    step = "numeric"
  ),
  
  prototype = list(
    max_val = 1000, 
    min_val = 10,
    p_initial = 0.95,
    step = 0.03 
  ) 
)


sim2_params <- setClass(
  # set the name for the class
  "sim2_params",
  
  # define the slots 
  slots = c(
    dataset_name = "character",
    dat_unit = "character",
    unit = "character",
    nsim = "numeric",
    all_id = "character",
    gene_id = "character", #* 
    extract_pattern_immune = "character", #*
    n_comp = "numeric",
    mix_type = "character",
    tumor_content = "character"
  ),
  
  prototype = list(
    dataset_name = c("GSE51984", "GSE64655", "GSE60424", "GSE115736"),
    dat_unit = c("count", "cpm", "tpm"),
    nsim = 20,
    n_comp = c(5:10),
    unit = c("count", "countNorm", "cpm", "tpm"),
    mix_type = c("orthog", "real")
  )
)



# marked with * need to be defined later 
sim3_params <- setClass(
  # set the name for the class
  "sim3_params",
  
  # define the slots 
  slots = c(
    dataset_name = "character",
    dat_unit = "character",
    unit = "character",
    nsim = "numeric",
    all_id = "character", #* 
    gene_id = "character",
    extract_pattern = "character", #*
    extract_pattern_immune = "character",
    n_comp = "numeric",
    mix_type = "character",
    tumor_content = "character"
  ),
  
  prototype = list(
    dataset_name = c("GSE51984", "GSE64655", "GSE60424", "GSE115736","GSE118490"),
    dat_unit = c("count", "cpm", "tpm"),
    nsim = 20,
    n_comp = c(5:10),
    unit = c("count", "countNorm", "cpm", "tpm"),
    mix_type = c("orthog", "real"),
    tumor_content = c("small", "large", "mosaic")
  )
)

deconv_params <- setClass(
  "deconv_params",
  
  slots = c(
    marker_list = "list",
    sig_list = "list",
    ref_anno_cibersort = "list",
    ref_anno_others = "list",
    ref_list = "list",
    ref_median_list = "list",
    epic_ref_list = "list"
  )
)

deconv_params_oneMat <- setClass(
  "deconv_params_oneMat",
  
  slots = c(
    marker_list = "list",
    sigGene = "character",
    ref_anno_cibersort = "matrix",
    ref_anno_others = "matrix",
    ref = "data.frame",
    ref_median = "data.frame",
    epic_ref_list = "list"
  )
)



# create a method to assign the value of velocity
setClassUnion("sim_params", c("sim1_params", "sim2_params", "sim3_params"))
setClassUnion("comp_sim_params", c("sim2_params", "sim3_params"))