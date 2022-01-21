

normalize_seurat <- function(dat, all = T, PCA = F){
  
  test_set <- dat
  n <- nrow(test_set)
  d <- ncol(test_set)
  
  # data to sparse matrix
  sdata <- as(as.matrix(t(test_set)), "dgCMatrix")
  
  # Initialize the Seurat object with the raw (non-normalized data).
  sdata <- Seurat::CreateSeuratObject(counts = sdata)
  sdata[["percent.mt"]] <- Seurat::PercentageFeatureSet(sdata, pattern = "^MT-")
  sdata <- Seurat::NormalizeData(object = sdata, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # 2000 variable genes were determined using the Seurat FindVariableFeatures function
  sdata <- Seurat::FindVariableFeatures(sdata, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(sdata)
  
  # Seurat ScaleData function was used to scale the data
  if(all){
    sdata <- Seurat::ScaleData(sdata, features = all.genes)
  }
  if(!all){
    sdata <- Seurat::ScaleData(sdata)
  }
  
  if (!PCA){
    return(t(sdata[["RNA"]]@scale.data))}
  
  if(PCA){
    sdata <- Seurat::RunPCA(object = sdata,  features = VariableFeatures(object = sdata))
    sdata <- Seurat::FindNeighbors(sdata, reduction = "pca", dims = 1:15)
    sdata <- Seurat::FindClusters(sdata, resolution = 0.8)
    sdata <- Seurat::RunUMAP(sdata, dims = 1:15)
    sdata <- Seurat::RunTSNE(object = sdata, dims.use = 1:15, do.fast = TRUE)
    
    tsne_data <- Seurat::Embeddings(object = sdata[["tsne"]])
    umap_data <- Seurat::Embeddings(object = sdata[["umap"]])
    
    return(list(t(sdata[["RNA"]]@scale.data), tsne_data, umap_data, Seurat::Embeddings(sdata, reduction = "pca")[,1:15]))
  }
}

para_est_new <- function(test_set){
  
  test_set <- test_set[, which(colSums(test_set) > 0)]
  n <- nrow(test_set)
  d <- ncol(test_set)
  N <- n+d+1
  
  mu <- log(sum(test_set) / (n*d))
  w <- log(rowMeans(test_set)) -mu
  r <- log(colMeans(test_set)) - mu
  
  para <- c(mu, w, r)
  
  return(para)
}


adj_CDF_logit <- function(dat, change = 1e-10){
  test_set <- dat
  test_set <- test_set[, which(colSums(test_set) > 0)]
  n <- nrow(test_set)
  d <- ncol(test_set)
  para <- para_est_new(test_set)
  mu <- para[1]
  w <- para[2:(n+1)]
  r <- para[(n+2):(n+d+1)]
  log_Pe = mu + w %*% matrix(rep(1, len = d), nrow = 1) + matrix(rep(1, len = n), ncol = 1) %*% r
  Pe <- exp(log_Pe)
  rownames(Pe) <- rownames(test_set)
  colnames(Pe) <- colnames(test_set)
  
  # adjust CDF
  F_adj <- function(a, b){
    Fadj <- (ppois(a,b) + ppois(a-1,b)) / 2 
    return(Fadj)
  }
  
  #compare with adjusted CDF
  l_adj_test <- list(a = as.list(as.vector(as.matrix(test_set))), b = as.list(Pe))
  cdf_adj_test <- purrr::pmap(l_adj_test, function(a,b) F_adj(a,b)) %>% 
    matrix(ncol = ncol(test_set))
  mcdf_adj_test <- unlist(cdf_adj_test) %>% matrix(ncol = ncol(test_set))
  rownames(mcdf_adj_test) <- rownames(test_set)
  colnames(mcdf_adj_test) <- colnames(test_set)
  mcdf_adj_test[mcdf_adj_test < change] <- change
  mcdf_adj_test[mcdf_adj_test > (1 - change)] <- 1 - change
  return(logit(mcdf_adj_test))
}
