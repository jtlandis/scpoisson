

test_that("scppp expects integers", {
  #expect an error
  expect_error(scppp(matrix(character(10))), " is not TRUE")
  #expect an error
  expect_error(scppp(matrix(rnorm(10))), " is not TRUE")
  #expect no error
  expect_error(scppp(matrix(integer(10))), NA)
})

test_that("scppp S3 dispatch works", {
  expect_equal(scppp(c(1:10))[["data"]], c(1:10))
  expect_equal(scppp(matrix(1:10))[["data"]], matrix(1:10))
})

set.seed(1234)
dat <- matrix(rpois(20, 0.2), 5, 4)
rownames(dat) <- paste0("gene", 1:5)
colnames(dat) <- paste0("cell", 1:4)

test_that("access data in scppp", {
  expect_equal(scppp(dat)[["data"]], dat)
})

test_that("scppp expects cells as columns and genes as rows by default", {
  expect_true(all(startsWith(rownames(scppp(dat)[["data"]]), "gene")))
  expect_true(all(startsWith(rownames(scppp(t(dat), sample = "column")[["data"]]), "cell")))
})


