test_that("Romeo and Juliet anchor extract", {

  r_and_j_anchor_set <- solve_nmf(
    create_input(Romeo_and_Juliet_tdm,
                 vocab = rownames(Romeo_and_Juliet_tdm),
                 covariates = acts,
                 topics = 25))$anchors
  expected_list <- c("time", "god", "light", "montague", "wilt",
                    "fair", "paris", "blood", "juliet", "speak",
                    "heaven", "sweet", "stand", "banished", "nurse",
                    "lady", "night", "tybalt", "death", "sir",
                    "dead", "child", "love", "romeo", "day")

  testthat::expect_identical(object = r_and_j_anchor_set,
                             expected = expected_list)



})
