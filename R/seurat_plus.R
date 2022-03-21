#' Louvain clustering using departure as data representation
#'
#' logit_Seurat_clustering returns a list with elements useful to check and compare cell clustering
#'
#' This is a function used to get cell clustering using Louvain clustering implemented in Seurat package
#'
#' @param test_set a UMI count data frame or matirx with cells as rows and genes as columns
#' @param pdat a matrix representated by model departure for each entry
#' @param PCA whether apply PCA before Louvain clustering, default is \code{TRUE}
#' @param N the number of pricipal components included for further clustering, default = 15
#' @param pres the resolution parameter in Louvain clustering, default = 0.8
#'
#' @return a list with the following elements:
#' \itemize{
#' \item{\code{sdata}}: {a Seurat object}
#' \item{\code{tsne_data}}: {a matrix containing t-SNE dimensionality reduction results,
#' with cells as rows, and first two t-SNE dimension as columns}
#' \item{\code{umap_data}}: {a matrix containing UMAP dimensionality reduction results,
#' with cells as rows, and first two UMAP dimension as columns}
#' \item{\code{sdata$seurat_clusters}}: {a vector with clustering index for each cell}
#' }
#'
#' @examples
#'
#' test_set <- matrix(rpois(500, 0.5), nrow = 10)
#' logit_Seurat_clustering(test_set)
#'
#' @references
#' \itemize{
#'     \item Stuart, T. et al. Comprehensive integration of single-cell data. Cell 177, 1888 (2019)
#' }
#'
#'
#' @export
LouvainDepart <- function(data, pdat = NULL, PCA = T,
                          N = 15, pres = 0.8,
                          tsne = F, umap = F, ...) UseMethod("LouvainDepart")
#' @export
#' @retuns scppp

LouvainDepart.scppp <- function(data, pdat = NULL, PCA = T,
                                     N = 15, pres = 0.8,
                                     tsne = F, umap = F){

  test_set <- data[["data"]]
  stopifnot('Require a matrix or data frame as input' = is.matrix(test_set))
  data$clust_results[["Lclust"]] <- LouvainDepart.matrix(test_set, pdat, PCA,
                       N, pres,
                       tsne, umap)
  return(data)
}
#' @export
#' @returns scppp_hclust_results

LouvainDepart.matrix <- function(data, pdat = NULL, PCA = T,
                          N = 15, pres = 0.8,
                          tsne = F, umap = F){

  test_set <- data
  stopifnot('Require a matrix or data frame as input' = is.matrix(test_set))
  sdat <- as(as.matrix(t(test_set)), "dgCMatrix")
  sdata <- Seurat::CreateSeuratObject(counts = sdat)
  if(is.null(pdat)){
    pdat <- adj_CDF_logit(test_set)
  }

  sdata[["scppp"]] <- Seurat::CreateAssayObject(counts = sdat)
  sdata[["scppp"]] <- SeuratObject::SetAssayData(sdata[["scppp"]], slot = "scale.data", new.data = t(pdat))
  .n <- ncol(t(pdat)) - 1L
  N <- min(c(N, .n - 1L))
  if(PCA){
    sdata <- Seurat::RunPCA(object = sdata, assay = "scppp",
                            reduction.name = "scppp_pca",
                            features = rownames(sdata),
                            npcs = min(c(50, .n)))
    sdata <- Seurat::FindNeighbors(sdata, reduction = "scppp_pca",
                                   dims = 1:N)
    sdata <- Seurat::FindClusters(sdata, graph.name = "scppp_snn", resolution = pres)
    if(umap){
      sdata <- Seurat::RunUMAP(sdata, reduction = "scppp_pca", dims = 1:N,
                               assay = "scppp",
                               reduction.name = "scppp_umap")
    }

    if(tsne){
      sdata <- Seurat::RunTSNE(object = sdata, reduction="scppp_pca",
                               assay = "scppp",
                               reduction.name = "scppp_tsne",
                               dims = 1:N, do.fast = TRUE, check_duplicates = FALSE)
    }

  } else {
    #do FindNeighbors on pdata instead of PCA cell embeddings
    neighbor.graphs <- Seurat::FindNeighbors(pdat, dims = 1:ncol(pdat))
    graph.name <- paste0("scppp_", names(neighbor.graphs))
    for (ii in seq_along(graph.name)) {
      if (inherits(x = neighbor.graphs[[ii]], what = "Graph")) {
        DefaultAssay(object = neighbor.graphs[[ii]]) <- "scppp"
      }
      sdata[[graph.name[[ii]]]] <- neighbor.graphs[[ii]]
    }
    sdata <- Seurat::FindClusters(sdata, graph.name = "scppp_snn", resolution = pres)
    if(umap){
      sdata[["scppp_umap"]] <- RunUMAP(pdat[, 1:N], assay = "scppp")
    }
    if(tsne){
      sdata[["scppp_tsne"]] <- RunTSNE(pdat[, 1:N], assay = "scppp")
    }

  }



  if(tsne){
    tsne_data <- Seurat::Embeddings(object = sdata[["scppp_tsne"]])
  }
  else tsne_data <- NULL
  if(umap){
    umap_data <- Seurat::Embeddings(object = sdata[["scppp_umap"]])
  }
  else umap_data <- NULL

  res_clust <- data.frame(names = names(sdata$seurat_clusters), cluster = sdata$seurat_clusters)

  return(list(sdata,
              tsne_data, umap_data,
              res_clust))
}


new_scppp_Seurat <- function(sdata, pdata) {

  obj <- Seurat::CreateSeuratObject(counts = t(sdata))
  assay_sc <- Seurat::CreateAssayObject(counts = t(sdata))
  assay_sc@scale.data <- t(pdata)

  obj@assays[["scppp"]] <- assay_sc
  return(obj)
}
