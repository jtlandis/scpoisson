
set.seed(1234)
scppp_obj <- scppp(rpois(200, 3))
qqplot_env_pois(scppp_obj, 3, 100)

# Count matrix as input
set.seed(1234)
dat <- matrix(c(rpois(300, 5), rpois(200, 1)), ncol = 20)
scppp_obj <- scppp(dat)
qqplot_env_pois(scppp_obj, L = 2, lambda = 5)

test_that("plot default call works properly", {
  set.seed(1234)
  ## check output ggplot
  expect_is(qqplot_env_pois(rpois(200, 3), 3, 100), "ggplot")
  expect_is(qqplot_env_pois(scppp(matrix(c(rpois(300, 5), rpois(200, 1)), ncol = 20)), L = 2, lambda = 5), "ggplot")
})
