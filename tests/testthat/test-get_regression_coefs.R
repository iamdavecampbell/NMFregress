test_that("multiplication works", {

  my_input <- create_input(Romeo_and_Juliet_tdm,
                          rownames(Romeo_and_Juliet_tdm),
                          topics = 25,
                          covariates = acts)
  my_output <- solve_nmf(my_input)

  # BETA
  coefs <- get_regression_coefs(my_output, model = "BETA", topics = c("romeo","juliet"))
  expect_equal(coefs["romeo","mean.Intercept"],
               -1.807218, tolerance = .0000001)
  expect_equal(coefs["juliet","mean.act3"],
               5.249710, tolerance = .0000001)
  expect_equal(coefs["juliet","precision.act5"],
               -6.1965656, tolerance = .0000001)
  expect_equal(coefs["romeo","precision.act2"],
               0.5200418, tolerance = .0000001)


  # OLS
  coefs <- get_regression_coefs(my_output, model = "OLS", topics = c("romeo","juliet"))
  expect_equal(coefs["romeo","Intercept"],
               0.04868248, tolerance = .0000001)
  expect_equal(coefs["juliet","act2"],
               0.05582528, tolerance = .0000001)

  #GAM
  my_output$covariates <- matrix( acts_continuous, ncol = 1, dimnames = list(NULL,"acts"))
  coefs <- get_regression_coefs(my_output,
                                model = "GAM",
                                topics = c("romeo"))
  expect_equal(coefs["romeo",50],
               -2.028353, tolerance = .00001)
  expect_equal(coefs["X_pred_vals",100],
               1.658228, tolerance = .00001)
  expect_equal(coefs["X_pred_vals",150],
               2.339623, tolerance = .00001)
  expect_equal(coefs["romeo",200],
               -1.910082, tolerance = .00001)


})
