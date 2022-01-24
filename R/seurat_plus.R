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
logit_Seurat_clustering <- function(test_set, pdat = NULL, PCA = T, N = 15, pres = 0.8){

  sdata <- as(as.matrix(t(test_set)), "dgCMatrix")
  sdata <- Seurat::CreateSeuratObject(counts = sdata)
  if(is.null(pdat)){
    pdat <- adj_CDF_logit(test_set)
  }

  sdata@assays$RNA@scale.data=t(pdat)
  sdata[["RNA"]]@scale.data = t(pdat)

  if(PCA){
    sdata <- Seurat::RunPCA(object = sdata,  features = rownames(sdata))
    sdata <- Seurat::FindNeighbors(sdata, reduction = "pca", dims = 1:N)
  }

  if(!PCA){
    sdata <- Seurat::RunPCA(object = sdata,  features = rownames(sdata))
    sdata@reductions$pca@cell.embeddings <- pdat
    sdata <- Seurat::FindNeighbors(sdata, reduction = "pca", dims = 1:ncol(pdat))
  }

  sdata <- Seurat::FindClusters(sdata, resolution = pres)
  sdata <- Seurat::RunUMAP(sdata, dims = 1:N)
  sdata <- Seurat::RunTSNE(object = sdata, dims.use = 1:N, do.fast = TRUE, check_duplicates = FALSE)

  tsne_data <- Seurat::Embeddings(object = sdata[["tsne"]])
  umap_data <- Seurat::Embeddings(object = sdata[["umap"]])

  return(list(sdata,
              tsne_data, umap_data,
              sdata$seurat_clusters))
}
