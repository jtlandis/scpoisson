# Tests differential_expression.R

set.seed(1234)
para1 <- matrix(c(rep(2L, 20), rep(1L, 40)), ncol = 1)
para2 <- matrix(c(rep(2L, 10), rep(1L, 20)), nrow = 1)
dat <- para1 %*% para2
dat[31:60, 1:10] <- dat[31:60, 1:10] + 10L
counts <- map_dbl(dat, ~rpois(1, .x)) %>%
  matrix(ncol = 30)
storage.mode(counts) <- "integer"
colnames(counts) <- paste0("c", 1:30)
rownames(counts) <- paste0("g", 1:60)
scppp_obj <- scppp(counts)
scppp_obj <- suppressWarnings(HclustDepart(scppp_obj, maxSplit = 3))
scppp_obj <- adj_CDF_logit(scppp_obj)


test_that("differential_expression works as expected", {
  # require clustering results
  expect_error(diff_gene_list(scppp(t(counts))))
  # cluster label should match the ones from clustering results
  expect_error(diff_gene_list(scppp_obj, clust1 = "1", clust2 = "2-1"))
  scppp_obj <- suppressWarnings(diff_gene_list(scppp_obj, clust1 = "1", clust2 = "2"))
  # columns as as expected
  expect_equal(colnames(scppp_obj[["de_results"]]$Hclust),
               c("variable", "clust1_mean", "clust2_mean", "clust1_n", "clust2_n",
                 "mean_diff", "statistic", "p.value", "padj", "abs_diff"))
})
