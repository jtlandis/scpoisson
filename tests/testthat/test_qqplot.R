library(vdiffr)

test_that("plot default call works properly", {

  # check output ggplot
  # Integer vector as input
  set.seed(1234)
  scppp_obj <- scppp(rpois(200, 3), "rows")
  p <- qqplot_env_pois(rpois(200, 3), 3, 100)
  expect_s3_class(p, "ggplot")
  vdiffr::expect_doppelganger("qq-env-plot-integer", p)

  # Count matrix as input
  set.seed(1234)
  dat <- matrix(c(rpois(300, 5), rpois(200, 1)), ncol = 20)
  scppp_obj <- scppp(dat, "rows")
  set.seed(1234)
  p2 <- qqplot_env_pois(scppp_obj, L = 2, lambda = 5)
  expect_s3_class(p2, "ggplot")
  vdiffr::expect_doppelganger("qq-env-plot-matrix", p2)
})
