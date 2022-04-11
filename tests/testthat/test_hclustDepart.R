
set.seed(1234)
para1 <- matrix(c(rep(2L, 20), rep(1L, 40)), ncol = 1)
para2 <- matrix(c(rep(2L, 10), rep(1L, 20)), nrow = 1)
dat <- para1 %*% para2
dat[31:60, 1:10] <- dat[31:60, 1:10] + 10L
counts <- map_dbl(dat, ~rpois(1, .x)) %>%
  matrix(ncol = 30)
storage.mode(counts) <- "integer"
colnames(counts) <- paste0("cell", 1:30)
rownames(counts) <- paste0("gene", 1:60)
scppp_obj <- scppp(counts)
scppp_obj <- suppressWarnings(HclustDepart(scppp_obj, maxSplit = 3))

test_that("HclustDepart works as expected", {
  # cell 1~10 in one cluster, the rest in the other cluster
  expect_true(all(scppp_obj[["clust_results"]]$Hclust[[1]]$cluster[1:10]  == "1"))
  expect_true(all(scppp_obj[["clust_results"]]$Hclust[[1]]$cluster[11:30]  == "2"))
  # p-value for first split less than 0.05
  expect_true(scppp_obj[["clust_results"]]$Hclust[[2]][1:1] < 0.05)
})
