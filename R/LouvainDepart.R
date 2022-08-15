#' Louvain clustering using departure as data representation
#'
#' This function returns a list with elements useful to check and compare cell clustering.
#'
#' This is a function used to get cell clustering using Louvain clustering algorithm implemented in the Seurat package.
#'
#' @param data A UMI count matrix with genes as rows and cells as columns or an S3 object for class 'scppp'.
#' @param pdat A matrix used as input for cell clustering. If not specify, the departure matrix will be calculated within the function.
#' @param PCA A logic value specifying whether apply PCA before Louvain clustering, default is \code{TRUE}.
#' @param N A numeric value specifying the number of principal components included for further clustering (default 15).
#' @param pres A numeric value specifying the resolution parameter in Louvain clustering (default 0.8)
#' @param tsne A logic value specifying whether t-SNE dimension reduction should be applied for visualization.
#' @param umap A logic value specifying whether UMAP dimension reduction should be applied for visualization.
#' @param ... not used.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{sdata}}: {a Seurat object}
#' \item{\code{tsne_data}}: {a matrix containing t-SNE dimension reduction results,
#' with cells as rows, and first two t-SNE dimensions as columns; NULL if \code{tsne = FALSE}.}
#' \item{\code{umap_data}}: {a matrix containing UMAP dimension reduction results,
#' with cells as rows, and first two UMAP dimensions as columns; NULL if \code{tsne = FALSE}.}
#' \item{\code{res_clust}}: {a data frame contains two columns: names (cell names) and clusters (cluster label)}
#' }
#'
#' @examples
#'
#' set.seed(1234)
#' test_set <- matrix(rpois(500, 2), nrow = 20)
#' rownames(test_set) <- paste0("gene", 1:nrow(test_set))
#' colnames(test_set) <- paste0("cell", 1:ncol(test_set))
#' LouvainDepart(test_set)
#'
#' @references
#' \insertRef{Seurat}{scpoisson}
#'
#'
#' @export
LouvainDepart <- function(data, pdat = NULL, PCA = TRUE,
                          N = 15, pres = 0.8,
                          tsne = FALSE, umap = FALSE, ...) UseMethod("LouvainDepart")
#' @export
#' @return scppp
LouvainDepart.scppp <- function(data, pdat = NULL, PCA = TRUE,
                                     N = 15, pres = 0.8,
                                     tsne = FALSE, umap = FALSE, ...){

  test_set <- data[["data"]]
  stopifnot('Require a matrix or data frame as input' = is.matrix(test_set))
  data$clust_results[["Lclust"]] <- LouvainDepart.matrix(test_set, pdat, PCA,
                       N, pres,
                       tsne, umap)
  return(data)
}
#' @export
#' @return scppp_Lclust_results
LouvainDepart.matrix <- function(data, pdat = NULL, PCA = TRUE,
                          N = 15, pres = 0.8,
                          tsne = FALSE, umap = FALSE, ...){

  test_set <- t(data)
  stopifnot('Require a matrix' = is.matrix(test_set))
  sdat <- methods::as(as.matrix(t(test_set)), "dgCMatrix")
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
        Seurat::`DefaultAssay<-`(object = neighbor.graphs[[ii]], "scppp")
      }
      sdata[[graph.name[[ii]]]] <- neighbor.graphs[[ii]]
    }
    sdata <- Seurat::FindClusters(sdata, graph.name = "scppp_snn", resolution = pres)
    if(umap){
      sdata[["scppp_umap"]] <- Seurat::RunUMAP(pdat[, 1:N], assay = "scppp")
    }
    if(tsne){
      sdata[["scppp_tsne"]] <- Seurat::RunTSNE(pdat[, 1:N], assay = "scppp")
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
