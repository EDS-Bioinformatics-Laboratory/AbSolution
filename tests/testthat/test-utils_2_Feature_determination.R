## TESTS: translate_fun #####

test_that("Basic translation works and prioritizes first ORF, AAA-> K", {
  expect_equal(translate_fun("AAAT", ORF_begins=1)[[1]], "K")
  expect_equal(translate_fun("AAAT", ORF_begins=1)[[2]], 1)
})

test_that("Basic translation works using the ORF_begins options", {
  expect_equal(translate_fun("TAAA", ORF_begins=1)[[1]], "NO_GOOD_ORF")
  expect_equal(translate_fun("TAAA", ORF_begins=1)[[2]], 0)
  expect_equal(translate_fun("TAAA", ORF_begins=0)[[1]], "K")
  expect_equal(translate_fun("TAAA", ORF_begins=0)[[2]], 2)
  expect_equal(translate_fun("TATAA",ORF_begins=1)[[1]], "Y")
  expect_equal(translate_fun("TATAA",ORF_begins=1)[[2]], 1)
})

test_that("Gives back empty AAs, A-> ", {
  expect_equal(translate_fun("A")[[1]], "")
  expect_equal(translate_fun("A")[[2]], "")
})
