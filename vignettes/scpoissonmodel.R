## ----setup--------------------------------------------------------------------
library(scpoisson)
library(dplyr)
library(ggplot2)

## ----poissoneity, message = F-------------------------------------------------
# UMI count matrix
test_dat <- as.data.frame(get_example_data("p5"))
scppp_obj <- scppp(test_dat, "rows")

# Q-Q envelope plot
qqplot_env_pois(scppp_obj, L = 10, lambda = 5) 


## ----clustering, warning=FALSE------------------------------------------------
# UMI count matrix
test_dat <- as.data.frame(get_example_data("p56"))
scppp_obj <- scppp(test_dat, "rows")

scppp_obj <- HclustDepart(scppp_obj, maxSplit = 3)

# cluster label for each cell
clust_res <- scppp_obj[["clust_results"]]$Hclust[[1]]
head(clust_res)
table(clust_res[,2])


## ----example3, message=FALSE, warning=FALSE-----------------------------------
scppp_obj <- LouvainDepart(scppp_obj)

# cluster label for each cell
clust_res2 <- scppp_obj[["clust_results"]]$Lclust[[4]]
head(clust_res2)
table(clust_res2[,2])

