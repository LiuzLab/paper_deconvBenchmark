Please execute scripts in the following sequence: 

1. sim1_params_set.R 
2. sim1_refGenerator.R

simModel 
1. sim1_mixGenerator_simModel.R 
2. sim1_diagnostics_simModel.R
3. sim1_deconvolution_simModel_rPlatform.R 
4. sim1_deconvolution_simModel_crossPlatform.R 
5. sim1_evaluation_simModel.R

libSize 
1. sim1_mixGenerator_libSize.R 
2. sim1_deconvolution_libSize_rPlatform.R 
3. sim1_deconvolution_libSize_crossPlatform.R 
4. sim1_evaluation_libSize.R