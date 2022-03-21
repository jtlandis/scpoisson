
set.seed(1234)
dat <- matrix(rpois(20, 0.2), 5, 4)
colnames(dat) <- paste0("gene", 1:4)
rownames(dat) <- paste0("cell", 1:5)

test_that("model departure requires removing rows/columns with only zero values", {
  expect_error(adj_CDF_logit(dat), " zero values")
  expect_error(adj_CDF_logit(dat[, which(colSums(dat) > 0)]), " zero values")
  expect_error(adj_CDF_logit(dat[which(rowSums(dat) > 0), ]), " zero values")
  expect_error(adj_CDF_logit(dat[which(rowSums(dat) > 0), which(colSums(dat) > 0)]), NA)
})
