
#@ Read expression data from quantification results
file.reader <- function(path,filelist,feature_index,column_index,header,unit="expr"){
  #@ path: folder that has quantification results stored
  #@ filelist: quantification files you want to read from the path
  #@ feature_index: gene ids or transcript ids
  #@ column_index: the exact column you want to extract from the file
  #@ output: a gene expression matrix with feature ids in the first row
  for(i in 1:length(filelist)){
    if(i == 1){
      dat <- read.table(paste0(path,filelist[i]),header=header,sep="\t")[,c(feature_index,column_index)]
    }
    else if(i > 1){
      tmp_dat <- read.table(paste0(path,filelist[i]),header=header,sep="\t")[,column_index]
      dat <- cbind(dat,tmp_dat)
      rm(tmp_dat)
    }
  }
  # assign column names according to the filelist
  # extract the first 10 character, concatenate it with the unit of the quantification results (if no unit is predetermined, it will just add expr)

  sample_name <- paste(substr(filelist,1,10),unit,sep="_")
  first_col <- colnames(dat)[1]
  colnames(dat) <- c(first_col,sample_name)
  return(dat)
}

#@ Eliminate row with duplicated ids (ENSEMBL Gene IDs)
eliminate.duplicate <- function(dat){
  #@ dat: a gene expression matrix with ensembl ID in the first row (row for feature id, column for sample name)
  #@ output: return a matrix with duplicated feature ids removed & use filtered feature ids in the first row as the row name
  gene_names = dat[,1]
  if(sum(duplicated(gene_names))==0){
    dat_tr <- dat[,-1]
    rownames(dat_tr) <- dat[,1]
  }
  else{
    duplicated_genes <- gene_names[duplicated(gene_names)]
    del_id <- vector()
    for(i in 1:length(duplicated_genes)){
      tmp_dupRow <- dat[dat[,1] == duplicated_genes[i],-1]
      # keep the rows which has maximal mean count value across samples
      keep <- names(which.max(rowMeans(tmp_dupRow)))
      # store the index of rows to delete
      del_id <- c(del_id,setdiff(rownames(tmp_dupRow),keep))
    }
    # store the transformed expression matrix
    dat_tr <- dat[-as.numeric((del_id)),-1]
    # set the row names
    rownames(dat_tr) <- dat[-as.numeric(del_id),1]

  }

  return(dat_tr)
}

#@ Perform filtering based on the threshold set according to the gene expression value
dat.filter <- function(dat,n,threshold){
  # This function picks genes that are expressed higher than a presetted threshold across all samples
  #@ dat: preprocessed gene expression matrix (row: gene symbols / column: sample names)
  #@ n: number of samples passed threshold to have the feature picked
  #@ threshold: gene expression threshold
  # output: filtered gene expression matrix (row:gene symbols column:sample names)
  library(genefilter)
  f1 <- kOverA(n, threshold)
  ffun <- filterfun(f1)
  wh1 <- genefilter(dat, ffun)
  p <- round(sum(wh1)*100/nrow(dat),digits = 2)
  print(paste0("Input gene number: ",nrow(dat)))
  print(paste0("Output gene number: ",sum(wh1)))
  print(paste0(p,"% of genes remained"))
  filtered_dat <- dat[wh1,]
  return(filtered_dat)
}

#@ Extract Ensembl ID and Gene Symbol Mapping Table
id.name.db <- function(version=NA, dataset="hsapiens_gene_ensembl"){
  #@ version: ensembl release version
  #           In default, the newest release of ensembl annotation will be extracted
  library(biomaRt)
  if(!is.na(version)){
    ensembl <- useEnsembl(biomart = "ensembl", version = version)
  }
  else{
    ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL") # just use most recent release
  }
  ensembl=useDataset(dataset,mart=ensembl)
  id_name_anno <- getBM(attributes = c("ensembl_gene_id","external_gene_name"), mart = ensembl )
  return(id_name_anno)
}

#@ Convert Ensembl Id to Gene Symbol
id.to.name <- function(dat,id_name_table=id_name_table){
  # match ensembl ids in dat with corresponding gene names
  #@ dat: gene expression matrix, rownames are ensembl ids
  #@ id_name_table: id_name_table is the output of id.name.db
  #@ output: return gene expression matrix with gene symbol as feature id
  ids <- rownames(dat)
  gene_names <- vector()
  gene_names <- sapply(ids,function(x) id_name_table[x,2])

  # remove gene names with NA values
  gene_names <- na.omit(gene_names)
  # remove rows of dat with non-identified gene names
  dat_tmp <- dat[names(gene_names),]

  duplicated_gene_names <- unique(gene_names[duplicated(gene_names)]) # extract duplicated
  del_id <- vector()
  for(i in 1:length(duplicated_gene_names)){
    sub_id <- names(gene_names[gene_names == duplicated_gene_names[i]])
    keep <- names(which.max(rowMeans(dat[sub_id,])))
    del_id <- c(del_id,setdiff(sub_id,keep))
  }
  # remove duplicated rows
  dat_tmp <- dat_tmp[!rownames(dat_tmp) %in% del_id,]
  # extract gene names from id_name_table
  dat_gene_names <- id_name_table[rownames(dat_tmp),2]
  # assign extracted gene names to the row names of dat_count_prep
  rownames(dat_tmp) <- dat_gene_names
  dat = dat_tmp
  return(dat)
}

#@ Draw sample-sample heatmap with specified sample annotations
heatmap_forDist <- function(cor_dat, meta, feature, main = "", breaks, clustering_method){
  #@ cor_dat: sample-sample correlation value - spearman's correlation coefficient is recommended
  #@ meta: meta data of samples used to pheatmap annotations
  #@ feature: sample feature added to the pheatmap annotation (should be one of column names of meta)
  #@ output: pheatmap plot
  library(pheatmap)
  library(RColorBrewer)
  id = substr(colnames(cor_dat), 1, 10)
  # build up the annotation according to the meta and feature
  rownames(meta) <- meta[, 'Run']
  annotation <- as.data.frame(matrix(data = NA, nrow = ncol(cor_dat), ncol = length(feature)))
  annotation_colors <- list()
  for(i in 1:length(feature)){
    annotation[, i] <- meta[id, feature[i]]
    tmp <- brewer.pal(n = length(unique(annotation[, i])), name = paste0("Set",i))
    names(tmp) <- unique(annotation[, i])
    annotation_colors[[i]] <- tmp
  }
  colnames(annotation) <- feature
  names(annotation_colors) <- feature
  rownames(annotation) <- rownames(cor_dat)
  # draw the heatmap
  if(missing(breaks)) {
    p <- pheatmap(cor_dat,
                  labels_row = rep("",ncol(cor_dat)),
                  labels_col=rep("",ncol(cor_dat)),
                  annotation_col=annotation,
                  annotation_colors = annotation_colors,
                  clustering_method = clustering_method,
                  treeheight_row = 0,
                  treeheight_col = 20,
                  main = main)
  } else {
    p <- pheatmap(cor_dat,
                  labels_row = rep("",ncol(cor_dat)),
                  labels_col = rep("",ncol(cor_dat)),
                  annotation_col = annotation,
                  annotation_colors = annotation_colors,
                  clustering_method = clustering_method,
                  treeheight_row = 0,
                  treeheight_col = 20,
                  main = main,
                  breaks=breaks)

  }

  return(p)
}


#@ draw scatter plot based on ID name
#@ correlation and euclidean distance between samples will be printed in the output figure
draw_scatterFromID <- function(id_x,id_y,dat,logTrans=TRUE){
  #@ id_x: sample name on the x axis
  #@ id_y: sample name on the y axis
  #@ logTrans: perform log transformation or not? (recommended)
  #@ output: sample-sample scatter plot, correlation (spearman) and euclidian distance will be printed out in the figure
  library(LSD)
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  x <- dat[,id_x]
  y <- dat[,id_y]
  cor_xy <- round(cor(x,y, method = "spearman"), digits = 2)
  euclidean_xy <- scientific_format(digits = 1)(round(dist(rbind(x,y)), digits = 0))
  if(logTrans){
    heatscatter(x, y,
                xlab=paste0("log2(",id_x,")"),
                ylab = paste0("log2(",id_y,")"),
                colpal = cols,
                pch = 20,
                cex.lab = 1.3,
                main = paste0("r = ", cor_xy, " , ","d = ", euclidean_xy),
                log = "xy")
    abline(a=0,b=1,col="red")
  }
  else{
    heatscatter(x, y,
                xlab = id_x,
                ylab = id_y,
                colpal = cols,
                pch = 20,
                cex.lab = 1.3,
                main=paste0("r = ", cor_xy, " , ", "d = ", euclidean_xy))
    abline(a = 0, b = 1, col = "red")
  }
}

#@ calculate coefficient of variance
cv.calc <- function(x){
  #@ x: expression of a gene across samples
  #@ output: coefficient of variation
  if(!is.na(x) && mean(x) != 0){
    cv_val <- sd(x)/mean(x)
  }
  else{
    cv_val = NA
  }
  return(cv_val)
}


#@ extract samples pair ids with lowest/highest correlations
index.arr.table <- function(dat,n,lowVal=TRUE){
  #@ dat: sample correlation matrix
  #@ n: number of sample pairs you want to extract
  #@ lowVal: TRUE if you want to extract sample pairs with lowest correlation
  if(lowVal == TRUE){
    low_values <- unique(sort(dat)[1:n])
    low_ids <- data.frame()
    for(i in 1:n){
      low_ids <- rbind(low_ids,which(dat == low_values[i],arr.ind=TRUE))
    }
    return(low_ids)
  }
  else{
    high_values <- unique(sort(dat,decreasing=TRUE)[1:n])
    high_ids <- data.frame()
    for(i in 1:length(high_values)){
      high_ids <- rbind(high_ids,which(dat == high_values[i], arr.ind=TRUE))
    }
    return(high_ids)
  }
}


#@ normalize each dataset by the median library size
median_libSize_norm <- function(x,M){
  #@ x: one expresion profile
  #@ M: median library size
  #@ output: normalized profile (by median library size)
  return(round(x*M/sum(x), digits = 0))
}

grid_arrange_shared_legend <- function(...) {
  library(ggplot2)
  library(gridExtra)
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

est.prop <- function(count){
  library(edgeR)
  p = goodTuringProportions(count)
  pooled_p = rowMeans(p, na.rm = TRUE)
  # remove na values
  pooled_p = na.omit(pooled_p)
  nGene = nrow(count) # number of genes in the count matrix, also the number of quantiles that is going to be generated based on the smoothed empirical CDF
  sorted_count <- sort(rowMeans(count)) # sort gene expression
  sorted_gene <- names(sorted_count) # list name according to the sorting result
  # generate smoothed emprifical CDF distribution and quantiles
  sim_quantile = quantile(ecdf(pooled_p),1:nGene/nGene)
  # attribute each sim_quantile value with sorted gene name
  names(sim_quantile) <- sorted_gene
  # return the proportion with original order
  est_prop <- sim_quantile[rownames(count)]
  return(est_prop)
}

# median.calc <- function(x){
#   return(apply(x,1,median))
# }

# 9/10/2019 deconvolution related functions
CAMfree_run <- function(mix, K = K){
  library(CAMTHC)
  tryCatch({
    return(CAM(mix, K))
  }, error = function(e) NA)
}

LinSeed_run <- function(mix, K = K){
  library(linseed)
  lo <- LinseedObject$new(mix, topGenes = 10000)
  lo$calculatePairwiseLinearity()
  lo$calculateSpearmanCorrelation()
  lo$calculateSignificanceLevel(100)
  lo$setCellTypeNumber(K)
  lo$filterDatasetByPval(0.01)
  lo$project("filtered")
  lo$smartSearchCorners(dataset = "filtered", error = "norm")
  lo$deconvolveByEndpoints()
  lo$proportions
  lo$selectGenes(100)
  return(lo)
}

zscore_assign <- function(est_marker, ref_median){
  K <- ncol(ref_median)
  celltype <- colnames(ref_median)
  assign_result <- rep("celltype", length.out = K)
  ref_z <- scale(ref_median, center = TRUE, scale = TRUE)
  for(k in 1:K){
    marker <- est_marker[[k]]
    ref_z_sub <- ref_z[marker, ]
    sum_z <- colSums(ref_z_sub)
    assign_result[k] <- celltype[which.max(sum_z)]
  }
  return(assign_result)
}

fc_calc <- function(ref_median){
  K <- ncol(ref_median)
  celltype <- colnames(ref_median)
  fc <- matrix(0, nrow = nrow(ref_median), ncol = ncol(ref_median))
  for(k in 1:K){
    # calculate fold change across for all other cell types vs. targeted cell type
    fc_all <- matrix(rep(ref_median[,k], K-1), ncol = K - 1) - ref_median[, -k]
    # pick the valu
    fc[,k] <- apply(fc_all, MARGIN = 1, function(x) return(x[which.min(x)]))
  }
  rownames(fc) <- rownames(ref_median)
  return(fc)
}

spcor_assign <- function(est_marker, est_s, ref_median){
  # calculate the spearman's correlation
  K <- ncol(ref_median)
  celltype <- colnames(ref_median)
  assign_result <- rep("celltype", length.out = K)
  for(k in 1:K){
    ref_sub <- ref_median[est_marker[[k]],]
    est_sub <- est_s[est_marker[[k]],k]
    spcor <- apply(ref_sub, 2, function(x) return(cor(est_sub, x, method = "spearman")))
    assign_result[k] <- celltype[which.max(spcor)]
  }
  return(assign_result)
}

fgsea_assign <- function(est_marker, ref_median, fc){
  library(fgsea)
  K <- ncol(ref_median)
  celltype <- colnames(ref_median)
  assign_result <- rep("celltype", length.out = K)
  for(k in 1:K){
    # run the test
    fgsea_result <- fgsea(pathways = est_marker,
                          stats = rank(fc[,k], ties.method = "random"),
                          minSize = 2,
                          maxSize = 500,
                          nperm = 10000)
    id_corner <- which.min(fgsea_result$padj)
    min_padj <- fgsea_result$padj[id_corner]
    if(min_padj == 1){
      #
      print("fgsea assignment failed! Unknown components detected")
      return(NA)
    }
    else{
      assign_result[id_corner] <- celltype[k]
    }
  }
  # if two or more celltypes assign to the same corner
  if(sum(grep("celltype", assign_result))){
    print("Overlapping assignments!")
    return(NA)
  }
  return(assign_result)
}

freeW_reorder <- function(est_w, assign_result, celltype_order){
  if(nrow(est_w) != length(assign_result)){
    est_w <- t(est_w)
  }
  if(sum(duplicated(assign_result)) > 0){
    # first check if there is any duplicated elements
    reorder_w <- matrix(0, nrow = nrow(est_w), ncol = ncol(est_w))
    rownames(reorder_w) <- celltype_order
    duplicated <- assign_result[duplicated(assign_result)]
    nonduplicated <- intersect(setdiff(celltype_order, duplicated), unique(assign_result))
    for(i in duplicated){
      # get the index of duplicates
      idx <- which(assign_result == i)
      # sum over the duplicated rows
      reorder_w[i,] = colSums(est_w[idx,])
    }
    for(i in nonduplicated){
      idx <- which(assign_result == i)
      reorder_w[i,] = est_w[idx,]
    }
  }
  else if(sum(duplicated(assign_result)) == 0){
    # set the rownames of the weight matrix as the assign_result

    rownames(est_w) <- assign_result
    # re-order the row with pre-defined order of the ground-truth
    reorder_w <- est_w[celltype_order,]
  }
  return(reorder_w)
}


MMAD_markerTransform <- function(marker, gene_id){
  K = length(table(marker[,2]))
  marker_mat <- matrix(0, ncol = K, nrow = length(gene_id))
  colnames(marker_mat) <- unique(marker[,2])
  rownames(marker_mat) <- gene_id
  for( k in 1:K){
    row_idx <- as.character(marker[marker[,2] == colnames(marker_mat)[k],1])
    marker_mat[row_idx,k] = 1
  }
  return(marker_mat)
}


CIBERSORT.transform <- function(X,first_row){
  X <- cbind(rownames(X),X)
  colnames(X) <- first_row
  return(X)
}

calc.cor <- function(x,y){
  x <- as.matrix(x)
  y <- as.matrix(y)

  if(dim(x)[1] == dim(y)[1] && length(x) == length(y)){
    x_t <- as.numeric(x)
    x_t[is.na(x_t)] <- 0
    y_t <- as.numeric(y)
    y_t[is.na(y_t)] <- 0
    result <- cor(x_t, y_t)
  }
  else if(dim(x)[1] == dim(y)[2] && length(x) == length(y)){
    x_t <- as.numeric(x)
    x_t[is.na(x_t)] <- 0
    y_t <- as.numeric(t(y))
    y_t[is.na(y_t)] <- 0
    result <- cor(x_t, y_t)
  }
  else{
    result <- NA
  }

  return(result)
}


calc.mad <- function(x,y){
  x <- as.matrix(x)
  y <- as.matrix(y)
  if(dim(x)[1] == dim(y)[1] && length(x) == length(y)){
    x_t <- as.numeric(x)
    x_t[is.na(x_t)] <- 0
    y_t <- as.numeric(y)
    y_t[is.na(y_t)] <- 0
    result <- sum(abs(x_t - y_t))/length(x_t)
  }
  else if(dim(x)[1] == dim(y)[2] && length(x) == length(y)){
    x_t <- as.numeric(x)
    x_t[is.na(x_t)] <- 0
    y_t <- as.numeric(t(y))
    y_t[is.na(y_t)] <- 0
    result <- sum(abs(x_t - y_t))/length(x_t)
  }
  else{
    result <- NA
  }
  return(result)
}

calc.cor_celltype <- function(est, truth,n){
  est = as.matrix(est)
  truth = as.matrix(truth)
  if(nrow(truth) != n){
    truth <- t(truth)
  }
  cor_result <- vector()

  if(is.na(est) || !length(est) == length(truth)){
    cor_result <- rep(-1,n)
  }
  else if(dim(est)[1] == n && length(est) == length(truth)){
    for(i in 1:n){
      x <- as.numeric(est[i,])
      x[is.na(x)] <- -1
      y <- as.numeric(truth[i,])
      y[is.na(y)] <- -1
      cor_result[i] <- cor(x, y)
    }
    names(cor_result) <- rownames(truth)
  }
  else if(dim(est)[2] == n && length(est) == length(truth)){
    for(i in 1:n){
      x <- as.numeric(est[,i])
      x[is.na(x)] <- -1
      y <- as.numeric(truth[i,])
      y[is.na(y)] <- -1
      cor_result[i] <- cor(x, y)
    }
    names(cor_result) <- rownames(truth)

  }

  cor_result[is.na(cor_result)] <- -1
  return(cor_result)
}

calc.mad_celltype <- function(est, truth,n){
  est = as.matrix(est)
  truth = as.matrix(truth)
  if(nrow(truth) != n){
    truth <- t(truth)
  }
  mad_result <- vector()

  if(is.na(est) || !length(est) == length(truth)){
    mad_result <- rep(1,n)
  }
  else if(dim(est)[1] == n && length(est) == length(truth)){
    for(i in 1:n){
      x <- as.numeric(est[i,])
      x[is.na(x)] <- 1
      y <- as.numeric(truth[i,])
      y[is.na(y)] <- 1
      mad_result[i] <- sum(abs(x - y))/length(x)
    }
    names(mad_result) <- rownames(truth)
  }
  else if(dim(est)[2] == n && length(est) == length(truth)){
    for(i in 1:n){
      x <- as.numeric(est[,i])
      x[is.na(x)] <- 1
      y <- as.numeric(truth[i,])
      y[is.na(y)] <- 1
      mad_result[i] <- sum(abs(x - y))/length(x)
    }
    names(mad_result) <- rownames(truth)

  }
  mad_result[is.na(mad_result)] <- 1
  return(mad_result)
}

cpm.HJ <- function(x_count){
  # This function is slower than the cpm from edgeR package
  # Suggests for double-check with edgeR::cpm() results
x_cpm <- 10^6*x_count/sum(x_count)
return(x_cpm)
}

grid_arrange_common_legend <- function(p_list, common_title, col_ratio, legend_position, row_number){
  library(cowplot)
  p_no_legend <- lapply(p_list, function(x) x + theme(legend.position = "none"))
  legend <- cowplot::get_legend(p_list[[1]] + theme(legend.position = legend_position))
  title <- cowplot::ggdraw() + cowplot::draw_label(common_title, fontface = "bold")
  p_grid <- cowplot::plot_grid(plotlist = p_no_legend, nrow = row_number)
  p_arrange <- cowplot::plot_grid(title, p_grid, legend, nrow = row_number, rel_widths = col_ratio)
  return(p_arrange)
}

