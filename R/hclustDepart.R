#' Significance for first split using sigclust2
#'
#' This function returns a list with elements mainly generated from sigclust2.
#'
#' This is a function used to calculate the significance level of the first split from hierarchical clustering
#' based on euclidean distance and Ward's linkage.
#'
#' @param test_dat A UMI count data matrix with samples to cluster as rows and features as columns.
#' @param minSize A numeric value specifying the minimal allowable cluster size (the number of cells for the smallest cluster, default 10).
#' @param sim A numeric value specifying the number of simulations during the Monte Carlo simulation procedure (default 100).
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{p}}: {p-value for the first split}
#' \item{\code{z}}: {z-score for the first split}
#' \item{\code{shc_result}}: {a \code{shc} S3-object as defined in sigclust2 package}
#' \item{\code{clust2}}: {a vector with group index for each cell}
#' \item{\code{clust_dat}}: {a matrix of data representation used as input for hierarchical clustering}
#' }
#'
#'
#' @references
#' \insertRef{kimes2017statistical}{scpoisson}
#' \insertRef{Rclusterpp}{scpoisson}
#'
#'
sigp <- function(test_dat, minSize = 10, sim = 100){

  test_dat <- test_dat
  stopifnot('Require a matrix' = is.matrix(test_dat))
  stopifnot('Remove columns with only zero values' = min(colSums(test_dat)) > 0)
  stopifnot('Remove rows with only zero values' = min(rowSums(test_dat)) > 0)
  if(nrow(test_dat) <= minSize) {
    return(list(NA,NA, NA, NA))
  }
  else {
    test_dat <- test_dat[, which(colSums(test_dat) > 0)]
  }
  clust_dat <- adj_CDF_logit(test_dat)

  shc_result <- shc_test(as.matrix(clust_dat), metric="euclidean", linkage="ward.D2", icovest = 2, n_sim = sim)
  p <- round(shc_result$p_norm[1], 7)
  m_idx <- colMeans(as.matrix(shc_result$ci_sim[1, , ]))
  s_idx <- apply(as.matrix(shc_result$ci_sim[1, , ]), 2, sd)
  z <- (shc_result$ci_dat[1, ] - m_idx) / s_idx

  clust2 <- cutree(shc_result$hc_dat, k = 2)
  clust2 <- clust2[match(rownames(test_dat), names(clust2))]
  return(list(p, z, shc_result,
              clust2, clust_dat))
}

#' Cluster cells in a recursive way
#'
#' This function returns a list with clustering results.
#'
#' This is a function used to get cell clustering results in a recursive way.
#' At each step, the two-way approximation is re-calculated again within each subcluster,
#' and the potential for further splitting is calculated using sigclust2.
#' A non significant result suggests cells are reasonably homogeneous
#' and may come from the same cell type. In addition, to avoid over splitting,
#' the maximum allowable number of splitting steps \code{maxSplit}
#' (default is 10, which leads to at most \eqn{2^{10} = 1024} total number of clusters) and
#' minimal allowable cluster size \code{minSize}
#' (the number of cells in a cluster allowed for further splitting, default is 10)
#' may be set beforehand.
#' Thus the process is stopped when any of the conditions
#' is satisfied: (1) the split is no longer statistically significant;
#' (2) the maximum allowable number of splitting steps is reached;
#' (3) any current cluster has less than 10 cells.
#'
#' @param data A UMI count matrix with genes as rows and cells as columns or an S3 object for class 'scppp'.
#' @param maxSplit A numeric value specifying the maximum allowable number of splitting steps (default 10).
#' @param minSize A numeric value specifying the minimal allowable cluster size (the number of cells for the smallest cluster, default 10).
#' @param sim A numeric value specifying the number of simulations during the Monte Carlo simulation procedure for statistical significance test, i.e. n_sim argument when apply sigclust2 (default = 100).
#' @param ... not used.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{res2}}: {a data frame contains two columns: names (cell names) and clusters (cluster label)}
#' \item{\code{sigclust_p}}: {a matrix with cells to cluster as rows, split index as columns,
#' the entry in row \code{i} and column \code{j} denoting the p-value
#' for the cell \code{i} at split step \code{j}}
#' \item{\code{sigclust_z}}: {a matrix with cells to cluster as rows, split index as columns,
#' the entry in row \code{i} and column \code{j} denoting the z-score
#' for the cell \code{i} at split step \code{j}}
#' }
#' If the input is an S3 object for class 'scppp', clustering result will be stored in object scppp under "clust_results".
#'
#' @examples
#'
#' test_set <- matrix(rpois(500, 0.5), nrow = 10)
#' HclustDepart(test_set)
#'
#' @export
HclustDepart <- function(data, maxSplit = 10, minSize = 10, sim = 100, ...) UseMethod("HclustDepart")

#' @export
#' @return scppp
HclustDepart.scppp <- function(data, maxSplit = 10, minSize = 10, sim = 100, ...) {

  test_dat <- data[["data"]]
  data$clust_results[["Hclust"]] <- HclustDepart.matrix(test_dat,
                                                        maxSplit,
                                                        minSize,
                                                        sim)
  return(data)
}

#' @export
#' @return scppp_hclust_results
HclustDepart.matrix <- function(data, maxSplit = 10, minSize = 10, sim = 100, ...){

  test_dat <- t(data)
  stopifnot('Remove columns with only zero values' = min(rowSums(test_dat)) > 0)
  stopifnot('Remove rows with only zero values' = min(colSums(test_dat)) > 0)
  test_dat <- test_dat[, which(colSums(test_dat) > 0)]

  S <- minSize
  J <- maxSplit
  Env <- new.env()

  clustonly <- function(test_dat, j = j){
    sigclust_obj <- sigp(test_dat, S, sim)
    if(j > J |
       cluster_size(test_dat) <= S |
       sigclust_obj[[1]] > 0.05) {
      warning("finish required split or cluster homogeneous enough")

      if(is.null(Env$clust1all[[j]]) & is.null(Env$clust2all[[j]])){
        Env$sigclust_p[, j] <- sigclust_obj[[1]]
        Env$sigclust_z[, j] <- sigclust_obj[[2]]
      }
      if(!is.null(Env$clust1all[[j]]) & all(rownames(test_dat) %in% Env$clust1all[[j]])){
        Env$sigclust_p[which(rownames(Env$res_split) %in% Env$clust1all[[j]]), j] <- sigclust_obj[[1]]
        Env$sigclust_z[which(rownames(Env$res_split) %in% Env$clust1all[[j]]), j] <- sigclust_obj[[2]]
      }

      if(!is.null(Env$clust2all[[j]]) & all(rownames(test_dat) %in% Env$clust2all[[j]])){
        Env$sigclust_p[which(rownames(Env$res_split) %in% Env$clust2all[[j]]), j] <- sigclust_obj[[1]]
        Env$sigclust_z[which(rownames(Env$res_split) %in% Env$clust2all[[j]]), j] <- sigclust_obj[[2]]
      }

      if(is.null(dim(test_dat))) {
        test_dat <- matrix(test_dat, nrow = 1)
        test_dat <- test_dat[, which(colSums(test_dat) > 0)]
        test_dat <- matrix(test_dat, nrow = 1)
      }
      test_dat <- test_dat[, which(colSums(test_dat) > 0)]
      #Env$res[[i]] <- test_dat
      Env$i <- Env$i + 1
    }
    else if (is.null(dim(test_dat))){
      Env$sigclust_p[, j] <- sigclust_obj[[1]]
      Env$sigclust_z[, j] <- sigclust_obj[[2]]
      test_dat <- matrix(test_dat, nrow = 1)
      test_dat <- test_dat[, which(colSums(test_dat) > 0)]
      test_dat <- matrix(test_dat, nrow = 1)
      #Env$res[[i]] <- test_dat
      Env$i <- Env$i + 1
    }
    else if (nrow(test_dat) <= 2) {warning("less than two cells are not allowed for further split")
      if(is.null(dim(test_dat))) {
        Env$sigclust_p[, j] <- sigclust_obj[[1]]
        Env$sigclust_z[, j] <- sigclust_obj[[2]]
        test_dat <- matrix(test_dat, nrow = 1)
        test_dat <- test_dat[, which(colSums(test_dat) > 0)]
        test_dat <- matrix(test_dat, nrow = 1)
      }
      test_dat <- test_dat[, which(colSums(test_dat) > 0)]
      Env$i <- Env$i + 1
    }

    else {
      test_dat <- test_dat[, which(colSums(test_dat) > 0)]
      dat <- adj_CDF_logit(test_dat)
      hc <- hclust(dist(dat, method = "euclidean"), method = "ward.D2")
      clust <- cutree(hc, k = 2)
      clust1 <- names(clust)[which(clust == 1)]
      clust2 <- names(clust)[which(clust == 2)]
      dat1 <- test_dat[which(rownames(dat) %in% clust1),]
      dat1 <- dat1[, which(colSums(dat1) > 0)]
      dat2 <- test_dat[which(rownames(dat) %in% clust2),]
      dat2 <- dat2[, which(colSums(dat2) > 0)]

      Env$res_split[which(rownames(Env$res_split) %in% clust1), j] <- 1
      Env$res_split[which(rownames(Env$res_split) %in% clust2), j] <- 2

      if(is.null(Env$clust1all[[j]]) & is.null(Env$clust2all[[j]])){
        Env$sigclust_p[, j] <- sigclust_obj[[1]]
        Env$sigclust_z[, j] <- sigclust_obj[[2]]
      }
      if(!is.null(Env$clust1all[[j]]) & all(rownames(test_dat) %in% Env$clust1all[[j]])){
        Env$sigclust_p[which(rownames(Env$res_split) %in% clust1), j] <- sigclust_obj[[1]]
        Env$sigclust_z[which(rownames(Env$res_split) %in% clust1), j] <- sigclust_obj[[2]]
      }

      if(!is.null(Env$clust2all[[j]]) & all(rownames(test_dat) %in% Env$clust2all[[j]])){
        Env$sigclust_p[which(rownames(Env$res_split) %in% clust2), j] <- sigclust_obj[[1]]
        Env$sigclust_z[which(rownames(Env$res_split) %in% clust2), j] <- sigclust_obj[[2]]
      }

      Env$clust1all[[j+1]] <- clust1
      Env$clust2all[[j+1]] <- clust2

      clustonly(dat1, j = j+1)
      clustonly(dat2, j = j+1)
    }
  }

  j <- 1
  Env$res_split <- matrix(NA, nrow = nrow(test_dat), ncol = J)
  Env$sigclust_p <- matrix(NA, nrow = nrow(test_dat), ncol = J)
  Env$sigclust_z <- matrix(NA, nrow = nrow(test_dat), ncol = J)
  rownames(Env$res_split) <- rownames(test_dat)
  Env$i <- 1
  Env$clust1all <- vector(mode = "list", length = J)
  Env$clust2all <- vector(mode = "list", length = J)

  clustonly(test_dat, j = 1)

  res_split_now <- Env$res_split
  sigclust_p <- Env$sigclust_p
  sigclust_z <- Env$sigclust_z

  res2 <- as.data.frame(res_split_now)
  res2$clust <- apply(res2[, 1:ncol(res2)], 1 , paste , collapse = "-" )
  res2 <- data.frame(names = rownames(res2), cluster = clust_clean(res2$clust))

  res_all <- list(res2, sigclust_p, sigclust_z)
  return(res_all)
}
