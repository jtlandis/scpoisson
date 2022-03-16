

test_that("scppp expects integers", {
  #expect an error
  expect_error(scppp(matrix(character(10))), " is not TRUE")
  #expect no error
  expect_error(scppp(matrix(integer(10))), NA)
})

test_that("scppp S3 dispatch works", {

})
