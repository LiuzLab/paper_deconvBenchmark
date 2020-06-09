# sim1_params_generics 
setGeneric(name = "set_gene_id",
           def = function(sim_param, suffix)
           {
             standardGeneric("set_gene_id")
           }
)

setGeneric(name = "set_eff_length",
           def = function(sim_param, path)
           {
             standardGeneric("set_eff_length")
           }
)

setGeneric(name = "get_celltype_anno",
           def = function(sim_param, suffix)
           {
             standardGeneric("get_celltype_anno")
           }
)


# sim1 ref extract generics 
setGeneric(name = "extract_celltype_profiles",
           def = function(sim_param, pattern_val, unit_val, ...)
           {
             standardGeneric("extract_celltype_profiles")
           }
)




# sim1 nb simulation generics 
setGeneric(name = "set_nb_params",
           def = function(nb_param, nsim_val = 20, bv_val = 0.2, libsize_val = 12*10^6, sigma_val = 0.25)
           {
             standardGeneric("set_nb_params")
           }
)

setGeneric(name = "sim_gammaPoisson",
           def = function(nb_param, reference)
           {
             standardGeneric("sim_gammaPoisson")
           }
)

setGeneric(name = "wrap_gammaPoisson",
           def = function(sim_param, nb_param, expr_list, interval)
           {
             standardGeneric("wrap_gammaPoisson")
           }
)

setGeneric(name = "sim_libsize_mix",
           def = function(sim_param, suffix_lib1, suffix_lib2, W)
           {
             standardGeneric("sim_libsize_mix")
           }
)

setGeneric(name = "sim_mix_unit_transform",
          def = function(sim_param, mix_list, unit_val, eff_length)
            {
            standardGeneric("sim_mix_unit_transform")
            }
)

setGeneric(name = "create_ref",
           def = function(sim_param, unit_val, ...)
           {
             standardGeneric("create_ref")
           }
)

setGeneric(name = "ref_anno_cibersort",
           def = function(sim_param, ref_list)
           {
             standardGeneric("ref_anno_cibersort")
           }
)

setGeneric(name = "ref_anno_cibersort_oneMat",
           def = function(ref, celltype_identity)
             {
             standardGeneric("ref_anno_cibersort_oneMat")
           })

setGeneric(name = "ref_anno_others",
           def = function(sim_param, ref_list)
             {
             standardGeneric("ref_anno_others")
           }
           )

setGeneric(name = "marker_gene_from_ref",
           def = function(sim_param, marker_param, ref, celltype)
             {
             standardGeneric("marker_gene_from_ref")
           }
)

setGeneric(name = "marker_gene_from_ref_oneMat",
           def = function(marker_param, ref, celltype)
           {
             standardGeneric("marker_gene_from_ref_oneMat")
           }
)

setGeneric(name = "marker_gene_from_refMean_oneMat",
           def = function(marker_param, ref_mean, celltype)
             {
             standardGeneric("marker_gene_from_refMean_oneMat")
           })


setGeneric(name = "sig_gene_from_ref",
           def = function(sim_param, ref, alpha, nFold,meta)
             {
             standardGeneric("sig_gene_from_ref")
           })

setGeneric(name = "sig_gene_from_ref_oneMat", 
           def = function(ref, celltype, identity, alpha, nFold)
           {
             standardGeneric("sig_gene_from_ref_oneMat")
           })

setGeneric(name = "create_epic_ref",
           def = function(sim_param, sig_gene, ref, median_ref)
           {
             standardGeneric("create_epic_ref")
           }
)

setGeneric(name = "create_epic_ref_oneMat",
           def = function(sig_gene, ref, median_ref)
             {
             standardGeneric("create_epic_ref_oneMat")
           })

setGeneric(name = "sim_mix",
           def = function(sim_param, W, ...)
          {
             standardGeneric("sim_mix")
           }
)

setGeneric(name = "sim_lognormal_noise", 
           def = function(clean_mix_list, sim_param, sigma){
             standardGeneric("sim_lognormal_noise")
           }
)

setGeneric(name = "sim_normal_noise", 
           def = function(clean_mix_list, sim_param, sigma){
             standardGeneric("sim_normal_noise")
           }
)

# Generics for deconvolution analysis 
setGeneric(name = "deconv_run_marker",
           # function for dsa, CAMmarker, deconv, ssKL, ssFrobenius
           def = function(mix, method_name, marker)
           {
             standardGeneric("deconv_run_marker")
           }
)

setGeneric(name = "deconv_run_refSig",
           # deconv_run for epic and DeconRNASeq
           def = function(mix, method_name, ref, ...){
             standardGeneric("deconv_run_refSig")
             }
           )

setGeneric(name = "deconv_run_free",
           def = function(mix, method_name, ref, K){
             standardGeneric("deconv_run_free")
           }
           )

             
setGeneric(name = "deconv_run_timer",
           # deconv_run for timer 
           def = function(mix, method_name, ref, anno, sig)
           {
             standardGeneric("deconv_run_timer")
           }
)

setGeneric(name = "deconv_run_CAMfree",
           # generic function for CAMfree 
           def = function(mix, method_name, K, ref)
           {
             standardGeneric("deconv_run_CAMfree")
           }
)

setGeneric(name = "deconv_run_crossPlatform",
           # generic function for mmad and cibersort
           def = function(method_name, sim_param, mix_name_prefix, ref_name_prefix, dat_path, sim_prefix)
           {
             standardGeneric("deconv_run_crossPlatform")
           }
)

setGeneric(name = "deconv_run_crossPlatform_oneMat",
           # generic function for mmad and cibersort 
           def = function(method_name, mix_name, ref_name, dat_path)
           {
            standardGeneric("deconv_run_crossPlatform_oneMat") 
           }
)



setGeneric(name = "deconv_analysis",
           def = function(method_name, deconv_param, sim_param, mix_list)
           {
             standardGeneric("deconv_analysis")
           }
)

setGeneric(name = "deconv_analysis_oneMat",
           def = function(method_name, deconv_param_oneMat, mix)
           {
             standardGeneric("deconv_analysis_oneMat")
           }
)

setGeneric(name = "deconv_write",
           def = function(dat, method_name, file_type, file_name, write_path)
           {
             standardGeneric("deconv_write")
           }
)

setGeneric(name = "deconv_read",
           def = function(file_name, method_name, read_path, nsim, nComp)
           {
             standardGeneric("deconv_read")
           }
)

setGeneric(name = "wrap_mix_write",
           def = function(mix_list, prefix, sim_param, file_type, method_name, write_path)
           {
             standardGeneric("wrap_mix_write")
           }
)

setGeneric(name = "wrap_ref_write",
           def = function(ref_list, prefix, sim_param, file_type, method_name, write_path)
           {
             standardGeneric("wrap_ref_write")
           }
)

setGeneric(name = "wrap_result_read",
           def = function(method_name, sim_param, prefix, read_path)
           {
             standardGeneric("wrap_result_read")
           }
)


setGeneric(name = "deconv_evaluation",
           def = function(eval_type = c("cor", "mad", "cellcor", "cellmad"), sim_param, est_W, truth, ...)
           {
             standardGeneric("deconv_evaluation")
           }
)

setGeneric(name = "deconv_eval_summary", 
          def = function(eval_type = c("cor", "mad", "cellcor", "cellmad"), robustness = c("all", "high_noise", "ref_noise"), sim_param, eval_metric)
          {
            standardGeneric("deconv_eval_summary")
          }
)

setGeneric(name = "wrap_eval_summary",
           def = function(sim_param, eval_type = c("cor", "mad", "cellcor", "cellmad"), all_methods, unit_val, ...)
           {
             standardGeneric("wrap_eval_summary")
           }
)


setGeneric(name = "deconv_plot_sampleScatter",
           def = function(sim_param, est_W, truth, ...)
           {
             standardGeneric("deconv_plot_sampleScatter")
           }
)
            
setGeneric(name = "evaluation_mat",
           def = function(eval_list, sim_param)
           {
             standardGeneric("evaluation_mat")
           }
           )

setGeneric(name = "evaluation_mat_celltype",
           def = function(eval_list, sim_param, celltype_val)
             {
             standardGeneric("evaluation_mat_celltype")
           })
            
setGeneric("tumor_weight", 
           def = function(known_weight, tumorContent = c("small", "large", "mosaic"))
             {
             starndardGeneric("tumor_weight")
           }
           )
setGeneric("orthog_weight",
           def=function(sim_param, N, sim_celltype)
           {
             standardGeneric("orthog_weight")
           }
)


setGeneric("real_weight",
           def=function(sim_param, N, sim_prop_mat)
           {
             standardGeneric("real_weight")
           }
)

setGeneric(name = "sim_eval_metrics",
           def = function(sim_param, dat_list, eval_type = c("cor", "cv", "mean", "var"))
           {
             standardGeneric("sim_eval_metrics")
           }
)

setGeneric("sample_scatter",
           def = function(sim_param, mix_list, scatter_type)
           {
             standardGeneric("sample_scatter")
           }
)

setGeneric("meanVar_scatter",
           def = function(sim_param, mean_mat, var_mat)
           {
             standardGeneric("meanVar_scatter")
           }
)


setGeneric("expr_heatmap",
           def = function(sim_param, ref, targeted_gene, ...)
           {
             standardGeneric("expr_heatmap")
           }
)
