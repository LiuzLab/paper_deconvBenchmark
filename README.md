# paper_deconvBenchmark

We have included all the analysis of the deconvolution benchmark paper in the manuscript Jin et al., "A benchmark for RNA-seq deconvolution analysis under dynamic testing environments". 

## Data
All the data described in the manuscript will be available in the following figshare link, https://figshare.com/projects/Data_for_the_manuscript_A_benchmark_for_RNA-seq_deconvolution_analysis_under_dynamic_testing_environments_/96908

Raw quantification data are stored as Raw.tar.gz

To free yourself from path reassignment, please store all raw data under ./Raw folder and intermediate .Rdata under ./data folder.

Detailed data description are in the file Data_description.txt

## Library

All required libraries are listed in the Libraries.R file. Please source this file before you start any analysis.  

## Execution order 
We organized the folder based on the analysis. 

The recommended execution order is indicated in the prefix(s1, s2, ..., s5) of the folder name. The execution order of the scripts under each subfolder (i.e. s1) is in the *_README.md file. 

You can start the analysis at any step after loading relevant intermediate .RData to the ./data folder. 



## Output
We suggest downloading the ./output folder so that you don't need to suffer from non-existent-folder error. 


#
Please contact us at haijing.jin@bcm.edu, if you have any questions. 
