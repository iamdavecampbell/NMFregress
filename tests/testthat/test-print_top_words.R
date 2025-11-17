test_that("Romeo and Juliet top words check", {

  r_and_j_juliet_top_words <- print_top_words(solve_nmf(
    create_input(Romeo_and_Juliet_tdm,
                 vocab = rownames(Romeo_and_Juliet_tdm),
                 covariates = acts,
                 topics = 25)))$juliet
  expected_list <- c("juliet", "beauty",
                    "grave", "talk",
                    "county", "kiss",
                    "hand", "tomb",
                    "poison", "true")

  testthat::expect_identical(object = r_and_j_juliet_top_words,
                             expected = expected_list)

})
