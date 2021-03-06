# preprocess GSE113590 
source("./code/src/Functions_all.R")
GSE113590_path = "./Raw/GSE113590/RSEM_quant/"
# read expr from .txt file (output of RSEM)
# unit: count, tpm and fpkm 
gene_file <- list.files(path=GSE113590_path,pattern="genes.results")
dat_count <- file.reader(GSE113590_path, gene_file, feature_index = 1, column_index=5, header=TRUE, unit="count")
dat_tpm <- file.reader(GSE113590_path, gene_file, feature_index = 1, column_index=6, header=TRUE, unit="tpm")
dat_fpkm <- file.reader(GSE113590_path, gene_file, feature_index = 1, column_index=7, header=TRUE, unit="fpkm")
# transfrom count to cpm
library(edgeR)
dat_cpm <- data.frame(gene_id = dat_count[,'gene_id'],edgeR::cpm(dat_count[,-1])) # beware of the first column is gene name thus should be elminated when using cpm() transformation 
colnames(dat_cpm) <- c('gene_id',paste0(substr(colnames(dat_cpm)[-1],1,10),"_cpm"))
# eliminate rows with duplication ids 
dat_count_p1 <- eliminate.duplicate(dat_count)
dat_tpm_p1 <- eliminate.duplicate(dat_tpm)
dat_fpkm_p1 <- eliminate.duplicate(dat_fpkm)
dat_cpm_p1 <- eliminate.duplicate(dat_cpm)
# convert ensembl id to gene symbol
# id_name_table is a two-column dataframe, column1 - gene id, column2 - gene name 
id_name_table <- id.name.db(version = 95) # define the version of ensembl annotation 
rownames(id_name_table) <- id_name_table[,1]
# id transform 
dat_count_p2 <- round(id.to.name(dat_count_p1,id_name_table = id_name_table))
dat_tpm_p2 <- id.to.name(dat_tpm_p1,id_name_table = id_name_table)
dat_fpkm_p2 <- id.to.name(dat_fpkm_p1,id_name_table = id_name_table)
dat_cpm_p2 <- id.to.name(dat_cpm_p1,id_name_table = id_name_table)
# filter out genes that are lowly expressed - morderate filtering 
dat_count_p3 <- dat.filter(dat_count_p2,4,10) # 13256 genes remained
dat_tpm_p3 <- dat.filter(dat_tpm_p2,4,1) # 14684 genes remained 
dat_fpkm_p3 <- dat.filter(dat_fpkm_p2,4,1) # 14441 genes remained 
dat_cpm_p3 <- dat.filter(dat_cpm_p2,4,1) #14134 genes remained 
# read meta data 
meta <- read.table("./Raw/GSE113590/GSE113590_meta.txt",header=TRUE,sep=",")
# store all processed data into a list 
GSE113590_dat <- list(meta = meta ,count = dat_count_p3,tpm = dat_tpm_p3,fpkm = dat_fpkm_p3,cpm = dat_cpm_p3)
GSE113590_dat_raw <- list(meta = meta ,count = dat_count_p2,tpm = dat_tpm_p2,fpkm = dat_fpkm_p2,cpm = dat_cpm_p2)
# save all preprocessed data and meta data for downstream analysis 
save(GSE113590_dat,file="./Data/GSE113590_dat.RData")
save(GSE113590_dat_raw, file="./Data/GSE113590_dat_raw.RData")