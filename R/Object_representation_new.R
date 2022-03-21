

obj <- structure(.Data = list(data = NULL, clust_results = NULL,
), class = "scppp")


new_result <- function(.Data = list(), ..., .subclass = NULL) {
  structure(.Data = .Data, ..., class = c(.subclass ,"scppp_results"))
}

new_departure_result <- function(data, ...) {
  new_result(.Data, ..., .subclass = "scppp_departure_results")
}

new_seurat_result <- function(data, ...) {
  new_result(.Data, ..., .subclass = "scppp_seurat_results")
}

