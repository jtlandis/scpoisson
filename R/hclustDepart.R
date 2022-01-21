#!/usr/bin/env Rscript

library(here)
library(SingleCellExperiment)
library(stringr)
library(purrr)

sample <-  "plate56" 
choose_plate <- c("Plate05A", "Plate06A")

load(here("data/sce_clean.rda"))
counts <- as.data.frame(assay(sce_clean))
test <- counts[rowSums(counts) > 0, ]
if(!is.null(choose_plate)) {test_plate <- test[, str_extract(colnames(test), "(?<=).+(?=_...)") %in% choose_plate]}
if(is.null(choose_plate)) {test_plate <- test}
df <- t(test_plate)
df <- df[,!grepl("MT-|RPL|RPS",colnames(df))]
test_plate <- t(df)
test_plate <- test_plate[which(rowSums(test_plate) > 0), ]
test_dat <- t(test_plate)



# Below from sigclust2
# Change code line move through nodes of dendrogram [for (k in 1:1)] to avoid calculation for each node

# shc test
shc_test <- function(x, metric = "euclidean", vecmet = NULL, matmet = NULL,
                     linkage = "ward.D2", l = 2,
                     alpha = 1, icovest = 1, bkgd_pca = FALSE, n_sim = 100,
                     n_min = 10, rcpp = FALSE, ci = "2CI", null_alg = "hclust",
                     ci_idx = 1, ci_emp = FALSE) {  
  
  n <- nrow(x)
  p <- ncol(x)
  
  if (n < 3) {
    stop("n must be >= 3")
  }
  
  n_ci <- length(ci)
  if (length(null_alg) != n_ci) {
    stop("ci and null_alg must be of same length")
  }
  
  for (ii in 1:n_ci) {
    if (ci[ii] == "linkage" && null_alg[ii] == "2means")
      stop("ci = 'linkage', null_alg = '2means' cannot be specified")
  }
  
  if (ci_idx > n_ci) {
    stop("invalid choice for ci_idx; ci_idx must be < length(ci)")
  }
  
  if (alpha > 1 || alpha < 0) {
    stop("invalid choice for alpha; alpha must be 0 < alpha < 1")
  }
  
  if (!is.matrix(x)) {
    stop("x must be a matrix; use as.matrix if necessary")
  }
  
  if (n_min < 3) {
    stop("n_min must be >= 3")
  }
  
  if (n_min > n) {
    stop("n_min must be <= n")
  }
  
  if (!is.null(vecmet) && !is.null(matmet)) {
    stop("only one of vecmet and matmet can be specified")
  }
  
  if (!is.null(vecmet)) {
    if (!is.function(vecmet)) {
      stop(paste("vecmet must be a function taking two vectors as input",
                 "and returning a real-valued dissimilarity"))
    }
    metric <- NULL
  }
  
  if (!is.null(matmet)) {
    if (!is.function(matmet)) {
      stop(paste("matmet must be a function taking a data matrix as input",
                 "and returning an object of class dist"))
    }
    metric <- NULL
  }
  
  if (rcpp && !requireNamespace("Rclusterpp", quietly = TRUE)) {
    stop("'Rclusterpp' package is not available.\n",
         "Either specify 'rcpp = FALSE' or install 'Rclusterpp' from GitHub using:\n",
         "> devtools::install_github('nolanlab/Rclusterpp')")
  }
  
  ## test vecmet and assign matmet if vecmet specified
  if (!is.null(vecmet)) {
    tryCatch({
      tmp <- vecmet(x[1, ], x[2, ])
    }, warning = function(e) {
      stop(paste0("warning for vecmet specification: ", e))
    }, error = function(e) {
      stop(paste0("error with vecmet specification: ", e))
    })
    matmet <- function(x) {
      as.dist(outer(split(x, row(x)), split(x, row(x)),
                    Vectorize(vecmet)))
    }
  }
  
  ## rclusterpp doesn't recognize 'ward.D2', stop and let user know
  if ((linkage == "ward.D2") & rcpp) {
    stop("Use 'ward' (in place of 'ward.D2') for linkage when rcpp = TRUE.")
  }
  
  ## apply initial clustering
  x_clust <- .initcluster(x, n, p, metric, matmet, linkage, l, 
                          n_ci, ci, rcpp)
  ci_dat <- x_clust$ci_dat
  hc_dat <- x_clust$hc_dat
  idx_hc <- x_clust$idx_hc
  
  ## for plotting purposes, change heights of dendrogram
  if ((linkage == "ward") & rcpp) {
    hc_dat$height <- sqrt(2*hc_dat$height)
  }
  
  ## p-values for all <= (n-1) tests
  p_emp <- matrix(2, nrow=n-1, ncol=n_ci)
  p_norm <- matrix(2, nrow=n-1, ncol=n_ci)
  colnames(p_emp) <- paste(null_alg, ci, sep="_")
  colnames(p_norm) <- paste(null_alg, ci, sep="_")
  
  ## null covariance parameters for all <= (n-1) tests
  eigval_dat <- matrix(-1, nrow=n-1, ncol=p)
  eigval_sim <- matrix(-1, nrow=n-1, ncol=p)
  backvar <- rep(-1, n-1)
  ci_sim <- array(-1, dim=c(n-1, n_sim, n_ci))
  
  ## determine parent nodes for all nodes
  pd_map <- .pd_map(hc_dat, n)
  
  ## compute Meinshausen cutoffs for significance at alpha
  cutoff <- fwer_cutoff(idx_hc, alpha)
  
  ## keep track of each node was tested
  nd_type <- rep("", n-1)
  
  ## move through nodes of dendrogram
  for (k in 1:1) {
    
    ## indices for subtree
    idx_sub <- unlist(idx_hc[k, ])
    n_sub <- length(idx_sub)
    
    ## only calc p-values for branches w/ more than n_min
    if (n_sub < n_min) {
      nd_type[k] <- "n_small"
      next
    }
    
    ## if parent wasn't significant, skip
    ## - placed after n_min check on purpose
    if ((alpha < 1) && (k > 1) && (nd_type[pd_map[k]] != "sig")) {
      nd_type[k] <- "no_test"
      next
    }
    
    ## estimate null Gaussian
    xk_null <- null_eigval(x[idx_sub, ], n_sub, p, icovest, bkgd_pca)
    
    ## prevent messages from looped application of clustering
    ## messages will be thrown once at initial clustering
    suppressMessages(
      ## simulate null datasets
      for (i in 1:n_sim) {
        xsim <- .simnull(xk_null$eigval_sim, n_sub, p)
        ci_sim[k, i, ] <- .calcCI_shc(xsim, p, metric, matmet, linkage, l, 
                                      n_ci, ci, null_alg, rcpp)
      }
    )
    
    ## compute p-values
    m_idx <- colMeans(as.matrix(ci_sim[k, , ]))
    s_idx <- apply(as.matrix(ci_sim[k, , ]), 2, sd)
    p_norm[k, ] <- pnorm(ci_dat[k, ], m_idx, s_idx)
    p_emp[k, ] <- colMeans(as.matrix(ci_sim[k, , ]) <= 
                             matrix(ci_dat[k, ], nrow=n_sim, 
                                    ncol=n_ci, byrow=TRUE))
    
    ## flip p-values for linkage based testing
    p_norm[k, ci == "linkage"] <- 1-p_norm[k, ci == "linkage"]
    p_emp[k, ci == "linkage"] <- 1-p_emp[k, ci == "linkage"]
    
    ## keep everything
    eigval_dat[k, ] <- xk_null$eigval_dat
    eigval_sim[k, ] <- xk_null$eigval_sim
    backvar[k] <- xk_null$backvar
    
    ## update nd_type (node type)
    if (alpha < 1) {
      if (ci_emp) {
        nd_type[k] <- ifelse(p_emp[k, ci_idx] < cutoff[k],
                             "sig", "not_sig")
      } else {
        nd_type[k] <- ifelse(p_norm[k, ci_idx] < cutoff[k],
                             "sig", "not_sig")
      }
    } else {
      nd_type[k] <- "cutoff_skipped"
    }
  }
  
  ## return shc S3 object
  structure(
    list(in_mat = x,
         in_args = list(metric = metric, linkage = linkage, alpha = alpha,
                        l = l, bkgd_pca = bkgd_pca, n_sim = n_sim,
                        n_min = n_min, icovest = icovest, ci = ci,
                        null_alg = null_alg, ci_idx = ci_idx, ci_emp = ci_emp),
         eigval_dat = eigval_dat,
         eigval_sim = eigval_sim,
         backvar = backvar,
         nd_type = nd_type,
         ci_dat = ci_dat,
         ci_sim = ci_sim,
         p_emp = p_emp,
         p_norm = p_norm,
         idx_hc = idx_hc,
         hc_dat = hc_dat),
    class = "shc")
}



## #############################################################################
## #############################################################################
## helper functions

## identify parent node of each node in dendrogram
.pd_map <- function(hc, n) {
  ## determine parent branch node for all children nodes along dendrogram
  pd_pairs <- rbind(cbind(hc$merge[, 1], 1:(n-1)), 
                    cbind(hc$merge[, 2], 1:(n-1)))
  pd_map <- data.frame(pd_pairs[pd_pairs[, 1] > 0, ])
  names(pd_map) <- c("dtr", "prt")
  pd_map <- pd_map$prt[order(pd_map$dtr)] #the parent of each daughter
  pd_map <- c(pd_map, n) #add final node without a parent
  
  ## flip index, hclust and shc use reversed ordering
  n - rev(pd_map)
}


## determine obs indices at each node of the dendrogram
.idx_hc <- function(hc, n) {
  ## list array of cluster indices at each of the n-1 merges
  idx_hc <- array(list(), c(2*n-1, 2))
  idx_hc[1:n, 1] <- as.list(n:1)
  idx_hc[(n+1):(2*n-1), ] <- hc$merge + n + (hc$merge<0)
  
  ## complete idx_hc
  for (k in 1:(n-1)) {
    idx_hc[[n+k, 1]] <- unlist(idx_hc[idx_hc[[n+k, 1]], ])
    idx_hc[[n+k, 2]] <- unlist(idx_hc[idx_hc[[n+k, 2]], ])
  }
  
  ## flip index, hclust and shc use revered ordering
  idx_hc[(2*n-1):(n+1), ]
}


## calculate sum of squares
.sumsq <- function(x) { norm(sweep(x, 2, colMeans(x), "-"), "F")^2 }


## calculate 2-means cluster index (n x p matrices)
.calc2CI <- function(x1, x2) {
  if (is.matrix(x1) && is.matrix(x2) && ncol(x1) == ncol(x2)) {
    (.sumsq(x1) + .sumsq(x2)) / .sumsq(rbind(x1, x2))
  } else {
    stop(paste("x1, x2 must be matrices with same ncols",
               "for 2CI calculation"))
  }      
}


## parse clustering parameters to produce hclust object
.cluster_shc <- function(x, metric, matmet, linkage, l, rcpp) {
  if (!is.null(matmet)) {
    hc_dat <- hclust(matmet(x), method=linkage)
  } else if (metric == "cor") {
    dmat <- 1 - WGCNA::cor(t(x))
    hc_dat <- hclust(as.dist(dmat), method=linkage)
  } else if (rcpp) {
    hc_dat <- Rclusterpp::Rclusterpp.hclust(x, method=linkage, 
                                            distance=metric, p=l)
  } else {
    hc_dat <- hclust(dist(x, method=metric, p=l), method=linkage)
  }
  hc_dat
}


## perform hierarchical clustering on the original data and 
## compute the corresponding cluster indices for each merge
.initcluster <- function(x, n, p, metric, matmet, linkage, l, 
                         n_ci, ci, rcpp) {
  
  ## obtain clustering solution
  hc_dat <- .cluster_shc(x, metric, matmet, linkage, l, rcpp)
  
  ## list array of cluster indices at each of the n-1 nodes
  idx_hc <- .idx_hc(hc_dat, n)
  
  ## matrix containing cluster indices
  ci_dat <- matrix(-1, nrow=n-1, ncol=n_ci)
  
  ## calculate cluster index(ices) for merge k
  for (i_ci in 1:n_ci) {
    if (ci[i_ci] == "2CI") {
      for (k in 1:(n-1)) {
        ci_dat[k, i_ci] <- .calc2CI(x[idx_hc[[k, 1]], , drop=FALSE],
                                    x[idx_hc[[k, 2]], , drop=FALSE])
      }
    } else if (ci[i_ci] == "linkage") {
      ## flip index, hclust and shc use revered ordering
      ci_dat[, i_ci] <- rev(hc_dat$height)
    }
  }
  
  list(hc_dat = hc_dat, 
       idx_hc = idx_hc,
       ci_dat = ci_dat)
}


## given null eigenvalues, simulate Gaussian dataset
.simnull <- function(eigval_sim, n, p) {
  simnorm <- matrix(rnorm(n*p, sd=sqrt(eigval_sim)), n, p, byrow=TRUE)
}


## perform hierarchical clustering on a simulated dataset and
## compute the correspond cluster indices for only the final merge
.calcCI_shc <- function(x, p, metric, matmet, linkage, l, 
                        n_ci, ci, null_alg, rcpp) {
  
  ##obtain clustering solution
  hc_isim <- .cluster_shc(x, metric, matmet, linkage, l, rcpp)
  split <- cutree(hc_isim, k=2)
  
  ##row vector containing cluster indices
  ci_isim <- matrix(-1, nrow=1, ncol=n_ci)
  
  for (i_ci in 1:n_ci) {
    if (ci[i_ci] == "2CI") {
      if (null_alg[i_ci] == "hclust") {
        ci_isim[i_ci] <- .calc2CI(x[split==1, , drop=FALSE],
                                  x[split==2, , drop=FALSE])        
      } else if (null_alg[i_ci] == "2means") {
        kmsol <- kmeans(x, centers=2)
        ci_isim[i_ci] <- kmsol$tot.withinss/kmsol$totss
      }
    } else if (ci[i_ci] == "linkage") {
      ci_isim[i_ci] <- hc_isim$height[nrow(x)-1]
    }
  }
  
  ci_isim
}

null_eigval <- function(x, n, p, icovest = 1, bkgd_pca = FALSE) {
  
  if (!(icovest %in% 1:3)) {
    warning("icovest should be 1, 2 or 3. Using default value: 1.")
    icovest <- 1
  }
  
  if (nrow(x) != n | ncol(x) != p)
    stop("Wrong size of matrix x!")
  
  ## compute background based on raw data
  ## or min of raw data and pca scores
  mad1 <- mad(as.matrix(x))
  if (bkgd_pca) {
    mad1 <- min(mad1, mad(as.matrix(prcomp(x)$x)) / sqrt(p/(n-1)))
  }
  backvar <- mad1^2
  
  avgx <- t(t(x) - colMeans(x))
  dv <- svd(avgx)$d
  eigval_dat <- dv^2/(n-1)
  
  ##pad with 0s
  eigval_dat <- c(eigval_dat, rep(0, p-length(eigval_dat)))
  eigval_sim <- eigval_dat
  
  if (icovest == 1) { #use soft 
    taub <- 0
    tauu <- .soft_covest(eigval_dat, backvar)$tau
    etau <- (tauu-taub) / 100
    ids <- rep(0, 100)
    ## tune to determine whether hard/soft more appropriate
    for (i in 1:100) {
      taus <- taub + (i-1)*etau
      eigval_temp <- eigval_dat - taus
      eigval_temp[eigval_temp < backvar] <- backvar
      ids[i] <- eigval_temp[1] / sum(eigval_temp)
    }
    tau <- taub + (which.max(ids)-1)*etau
    eigval_sim <- eigval_dat - tau
    eigval_sim[eigval_sim < backvar] <- backvar
    
  } else if (icovest == 2) { #use sample eigenvalues
    eigval_sim[eigval_sim < 0] <- 0
    
  } else if (icovest == 3) { #use hard thresholding
    eigval_sim[eigval_dat < backvar] <- backvar
  } else {
    stop("covest must be 1, 2 or 3")
  }
  
  list(eigval_dat = eigval_dat,
       backvar = backvar,
       eigval_sim = eigval_sim)
}



## helper function for computing soft thresholding estimator
.soft_covest <- function(vsampeigv, sig2b) {
  
  p <- length(vsampeigv)
  vtaucand <- vsampeigv - sig2b
  
  ##if all eigenvals > sig2b, just use sample eigenvals
  if (vtaucand[p] > 0) {
    return(list(veigvest = vsampeigv,
                tau = 0))
  }  
  
  ##if not enough power, just use flat est as in Matlab impl
  if (sum(vsampeigv) <= p*sig2b) {
    return(list(veigvest = rep(sig2b, p),
                tau = 0))
  }
  
  ##find threshold to preserve power
  which <- which(vtaucand <= 0)
  icut <- which[1] - 1
  powertail <- sum(vsampeigv[(icut+1):p])
  power2shift <- sig2b*(p-icut) - powertail
  
  vi <- c(1:icut)
  vcumtaucand <- sort(cumsum(sort(vtaucand[vi])), decreasing=TRUE)
  
  vpowershifted <- (vi-1)*vtaucand[vi] + vcumtaucand
  
  flag <- (vpowershifted < power2shift)
  if (sum(flag) == 0) {
    ## means decreasing everything still not enough
    itau <- 0
  } else {
    ## means for some index, decrease is sufficient
    itau <- which(flag)[1]
  }
  
  if (itau == 1) {
    powerprop <- power2shift/vpowershifted[1] #originally no [1] idx, PKK
    tau <- powerprop*vtaucand[1]
  } else if (itau == 0) {
    powerprop <- power2shift/vpowershifted[icut] 
    tau <- powerprop*vtaucand[icut] 
  } else {
    powerprop <- (power2shift-vpowershifted[itau]) /
      (vpowershifted[itau-1]-vpowershifted[itau])
    tau <- vtaucand[itau] + powerprop*(vtaucand[itau-1] - vtaucand[itau])
  }
  
  
  veigvest <- vsampeigv - tau 
  flag <- (veigvest > sig2b) 
  veigvest <- flag*veigvest + (1-flag)*(sig2b*rep(1, p))
  
  ##return eigenvalue estimate and soft threshold parameter, tau
  list(veigvest = veigvest,
       tau = tau)
}


#' return Family-Wise Error Rate (FWER) cutoffs
#' 
#' @name fwer_cutoff-generic
#' @docType methods
#' @export
#' @keywords internal
fwer_cutoff <- function(obj, ...) {
  UseMethod("fwer_cutoff", obj)
}



#' get FWER cutoffs for shc oject
#'
#' @param obj \code{shc} object
#' @param alpha numeric value specifying level
#' @param ... other parameters to be used by the function
#'
#' @name fwer_cutoff-shc
#' @export
#' @method fwer_cutoff shc
#' @author Patrick Kimes
fwer_cutoff.shc <- function(obj, alpha, ...) {
  fwer_cutoff(obj$idx_hc, alpha)
}



#' get FWER from idx_hc attribute of shc object 
#'
#' @param obj \code{shc} object
#' @param alpha numeric value specifying level
#' @param ... other parameters to be used by the function
#' 
#' @name fwer_cutoff-matrix
#' @method fwer_cutoff matrix
#' @author Patrick Kimes
#' @keywords internal
fwer_cutoff.matrix <- function(obj, alpha, ...) {
  alpha/(nrow(obj)+1) *
    apply(obj, 1, function(x) { length(unlist(x)) })
}


##########################################################################################

# Below are our clustering algorithm

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

`%notin%` <- Negate(`%in%`)

clust_clean <- function(clust){
  if(grepl("NA", clust, fixed = TRUE)){
    gsub("(.*?)(-NA.*)", "\\1", clust)
  }
  else(gsub('.{2}$', '', clust))
}

cluster_size <- function(test_dat){
  if(is.null(dim(test_dat))) {
    return(0)
  } else{
    test_dat <- test_dat[, which(colSums(test_dat) > 0)]
    return(nrow(test_dat))}
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
  F_adj <- function(a, b, .version = version){
      Fadj <- (ppois(a,b) + ppois(a-1,b)) / 2 
    return(Fadj)
  }
  #compare with adjusted CDF
  l_adj_test <- list(a = as.list(as.vector(as.matrix(test_set))), b = as.list(Pe))
  cdf_adj_test <- pmap(l_adj_test, function(a,b) F_adj(a,b)) %>% 
    matrix(ncol = ncol(test_set))
  mcdf_adj_test <- unlist(cdf_adj_test) %>% matrix(ncol = ncol(test_set))
  rownames(mcdf_adj_test) <- rownames(test_set)
  colnames(mcdf_adj_test) <- colnames(test_set)
  mcdf_adj_test[mcdf_adj_test < change] <- change
  mcdf_adj_test[mcdf_adj_test > (1 - change)] <- 1 - change
  return(psych::logit(mcdf_adj_test))
}

sigp <- function(test_dat, minSize = 10){
  
  if(is.null(dim(test_dat))) {
    return(list(NA,NA, NA, NA))
  } else if(nrow(test_dat) <= minSize) {
    return(list(NA,NA, NA, NA))
  } else{
    test_dat <- test_dat[, which(colSums(test_dat) > 0)]
  }
    clust_dat <- adj_CDF_logit(test_dat)
  
    shc_result <- shc_test(as.matrix(clust_dat), metric="euclidean", linkage="ward", icovest = 2, rcpp = T, n_sim = 1000)
    p <- round(shc_result$p_norm[1], 7)
    m_idx <- colMeans(as.matrix(shc_result$ci_sim[1, , ]))
    s_idx <- apply(as.matrix(shc_result$ci_sim[1, , ]), 2, sd)
    z <- (shc_result$ci_dat[1, ] - m_idx) / s_idx
    
    clust2 <- cutree(shc_result$hc_dat, k = 2)
    clust2 <- clust2[match(rownames(test_dat), names(clust2))]
    return(list(p, z, shc_result,
                clust2, clust_dat))
}

split_output_withsig <- function(test_dat, maxSplit = 10, minSize = 10){
  
  S <- minSize
  J <- maxSplit
  
  clustonly <- function(test_dat, j = j){
    sigclust_obj <- sigp(test_dat, S)
    if(j > J | 
       cluster_size(test_dat) <= S |
       sigclust_obj[[1]] > 0.05) {
      warning("finish required split or cluster homogeneous enough")
      
      if(is.null(.GlobalEnv$clust1all[[j]]) & is.null(.GlobalEnv$clust2all[[j]])){
        .GlobalEnv$sigclust_p[, j] <- sigclust_obj[[1]]
        .GlobalEnv$sigclust_z[, j] <- sigclust_obj[[2]]
      } 
      if(!is.null(.GlobalEnv$clust1all[[j]]) & all(rownames(test_dat) %in% .GlobalEnv$clust1all[[j]])){
        .GlobalEnv$sigclust_p[which(rownames(.GlobalEnv$res_split) %in% .GlobalEnv$clust1all[[j]]), j] <- sigclust_obj[[1]]
        .GlobalEnv$sigclust_z[which(rownames(.GlobalEnv$res_split) %in% .GlobalEnv$clust1all[[j]]), j] <- sigclust_obj[[2]]
      } 
      
      if(!is.null(.GlobalEnv$clust2all[[j]]) & all(rownames(test_dat) %in% .GlobalEnv$clust2all[[j]])){
        .GlobalEnv$sigclust_p[which(rownames(.GlobalEnv$res_split) %in% .GlobalEnv$clust2all[[j]]), j] <- sigclust_obj[[1]]
        .GlobalEnv$sigclust_z[which(rownames(.GlobalEnv$res_split) %in% .GlobalEnv$clust2all[[j]]), j] <- sigclust_obj[[2]]
      }   
      
      if(is.null(dim(test_dat))) {
        test_dat <- matrix(test_dat, nrow = 1)
        test_dat <- test_dat[, which(colSums(test_dat) > 0)]
        test_dat <- matrix(test_dat, nrow = 1)
      }
      test_dat <- test_dat[, which(colSums(test_dat) > 0)]
      #.GlobalEnv$res[[i]] <- test_dat
      .GlobalEnv$i <- .GlobalEnv$i + 1
    }
    else if (is.null(dim(test_dat))){
      .GlobalEnv$sigclust_p[, j] <- sigclust_obj[[1]]
      .GlobalEnv$sigclust_z[, j] <- sigclust_obj[[2]]
      test_dat <- matrix(test_dat, nrow = 1)
      test_dat <- test_dat[, which(colSums(test_dat) > 0)]
      test_dat <- matrix(test_dat, nrow = 1)
      #.GlobalEnv$res[[i]] <- test_dat
      .GlobalEnv$i <- .GlobalEnv$i + 1
    }
    else if (nrow(test_dat) <= 2) {warning("less than two cells are not allowed for further split")
      if(is.null(dim(test_dat))) {
        .GlobalEnv$sigclust_p[, j] <- sigclust_obj[[1]]
        .GlobalEnv$sigclust_z[, j] <- sigclust_obj[[2]]
        test_dat <- matrix(test_dat, nrow = 1)
        test_dat <- test_dat[, which(colSums(test_dat) > 0)]
        test_dat <- matrix(test_dat, nrow = 1)
      }
      test_dat <- test_dat[, which(colSums(test_dat) > 0)]
      .GlobalEnv$i <- .GlobalEnv$i + 1
    }
    
    else {
      test_dat <- test_dat[, which(colSums(test_dat) > 0)]
      dat <- adj_CDF_logit(test_dat)
      hc <- hclust(dist(dat, method = "euclidean"), method = "ward.D2")
      clust <- cutree(hc, k = 2)
      clust1 <- names(clust)[which(clust == 1)]
      clust2 <- names(clust)[which(clust == 2)]
      dat1 <- test_dat[which(rownames(dat) %in% clust1), ]
      dat2 <- test_dat[which(rownames(dat) %in% clust2), ]
      
      .GlobalEnv$res_split[which(rownames(.GlobalEnv$res_split) %in% clust1), j] <- 1
      .GlobalEnv$res_split[which(rownames(.GlobalEnv$res_split) %in% clust2), j] <- 2
      
      if(is.null(.GlobalEnv$clust1all[[j]]) & is.null(.GlobalEnv$clust2all[[j]])){
        .GlobalEnv$sigclust_p[, j] <- sigclust_obj[[1]]
        .GlobalEnv$sigclust_z[, j] <- sigclust_obj[[2]]
      } 
      if(!is.null(.GlobalEnv$clust1all[[j]]) & all(rownames(test_dat) %in% .GlobalEnv$clust1all[[j]])){
        .GlobalEnv$sigclust_p[which(rownames(.GlobalEnv$res_split) %in% clust1), j] <- sigclust_obj[[1]]
        .GlobalEnv$sigclust_z[which(rownames(.GlobalEnv$res_split) %in% clust1), j] <- sigclust_obj[[2]]
      } 
      
      if(!is.null(.GlobalEnv$clust2all[[j]]) & all(rownames(test_dat) %in% .GlobalEnv$clust2all[[j]])){
        .GlobalEnv$sigclust_p[which(rownames(.GlobalEnv$res_split) %in% clust2), j] <- sigclust_obj[[1]]
        .GlobalEnv$sigclust_z[which(rownames(.GlobalEnv$res_split) %in% clust2), j] <- sigclust_obj[[2]]
      } 
      
      .GlobalEnv$clust1all[[j+1]] <- clust1
      .GlobalEnv$clust2all[[j+1]] <- clust2
      
      clustonly(dat1, j = j+1)
      clustonly(dat2, j = j+1)
    }
  }
  
  j <- 1
  .GlobalEnv$res_split <- matrix(NA, nrow = nrow(test_dat), ncol = J)
  .GlobalEnv$sigclust_p <- matrix(NA, nrow = nrow(test_dat), ncol = J)
  .GlobalEnv$sigclust_z <- matrix(NA, nrow = nrow(test_dat), ncol = J)
  rownames(.GlobalEnv$res_split) <- rownames(test_dat)
  .GlobalEnv$i <- 1
  .GlobalEnv$clust1all <- vector(mode = "list", length = J)
  .GlobalEnv$clust2all <- vector(mode = "list", length = J)
  
  clustonly(test_dat, j = 1)
  
  res_split_now <- .GlobalEnv$res_split
  sigclust_p <- .GlobalEnv$sigclust_p
  sigclust_z <- .GlobalEnv$sigclust_z
  res_all <- list(res_split_now, sigclust_p, sigclust_z)
  return(res_all)
}

# split_output_all <- split_output_withsig(test_dat)
# save(split_output_all, file = paste("~/single-cell-analysis_yue/single cell_cell clustering/clean data II/results/", sample, "_hclustDepart.Rdata", sep = ""))


