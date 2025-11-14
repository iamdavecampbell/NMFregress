#' @title Extract the Regression Coefficients
#'
#' @description Compute OLS coefficients, fitting a linear model between a
#' user's specified covariates and topic weight. Estimates are produced
#' simultaneously for all topics.
#'
#' @param output An object of class nmf_output.
#'
#' @param obs_weights Weights for the documents, as used by weighted least
#' squares. Defaults to null.
#'
#' @param model choose c("BETA", "GAM", "OLS") for OLS, beta regression or a
#' Generalized Additive model with family beta.  Default is "BETA".   OLS is
#' only useful if the covariates are categorical, to be passed to
#' get_regression_coefs
#'
#' @param return_just_coefs is a logical, if TRUE then just return the
#' coefficients, if FALSE, then return the output from the lm or betareg
#' function.
#'
#' @param formula is the formula to be passed into betareg or GAM.  With
#' betaregression, formula has the form Y~ model+for+mean | model+for+dispersion
#'
#' @param link is the link function for the GLM for the mean when using
#' betaregression or GAM.
#'
#' @param link.phi is the link function for the GLM for the precision when using
#'  betaregression.
#'
#' @param type is one of ML, BC, BR for betaregression to use Maximum
#' Likelihood, Bias Corrected, or Bias Reduced respectively
#'
#' @param topics which topics should be used. If empty, then all of them are
#' used.
#'
#' @param na.rm remove the problematic 'no topic' documents now so that we don't
#' end up with a bootstrap sample of all NAs
#'
#' @return A matrix of regression coefficients (named if column names have been
#' specified for the design matrix).
#'
#' @import dplyr
#' @importFrom stats as.formula coefficients fitted  weights
#' @importFrom betareg betareg betareg.control
#' @importFrom mgcv gam betar predict.gam
#'
#' @examples
#' my_input <- create_input(Romeo_and_Juliet_tdm,
#'                         rownames(Romeo_and_Juliet_tdm),
#'                         topics = 10,
#'                         covariates = acts)
#' my_output <- solve_nmf(my_input)
#' get_regression_coefs(my_output, model = "BETA")
#'
#'
#' @export
get_regression_coefs <- function(output,
                                 obs_weights = NULL,
                                 model = c("BETA", "GAM", "OLS"),
                                 return_just_coefs = TRUE,
                                 formula = NULL,
                                 link = "logit",
                                 link.phi = "log",
                                 type = "ML",
                                 topics = NULL,
                                 na.rm = TRUE) {
  ##### Check input types/whether they include covariates
  if (!inherits(output, "nmf_output")) {
    stop("Output must be of class nmf_output.")
  }
  if (is.null(output$covariates)) {
    stop("No design matrix specified with output.")
  }
  if (!(is.null(obs_weights))) {
    if (!(is.numeric(obs_weights))) {
      stop("Document weights must be in the form of a numeric vector.")
    }
    if (sum(obs_weights >= 0) != length(obs_weights)) {
      stop("Document weights must be positive or zero.")
    }
    if (length(obs_weights) != ncol(output$theta)) {
      stop("Document weights and documents have differing lengths.")
    }
  }
  if (is.null(topics)) {
    topics <- output$anchors
  }
  if (length(model) == 3) {
    model <- "BETA"
  }
  ##### set up matrices for regression/return value
  # select the rows corresponding to topics selected:
  if (na.rm) {
    #remove the problematic 'no topic' documents now so that we don't end up
    # with a bootstrap sample of all NAs
    theta <- t(output$theta[which(output$anchors %in% topics),
                            which(!is.na(1 / output$sum_theta_over_docs) &
                                    !is.infinite(1 / output$sum_theta_over_docs)
                                  )])
    covariates         <- output$covariates[which(
      !is.na(1 / output$sum_theta_over_docs) &
        !is.infinite(1 / output$sum_theta_over_docs)), ] |>
      as.data.frame()
    colnames(covariates) <- colnames(output$covariates)
    sum_theta_over_docs <- output$sum_theta_over_docs[which(
      !is.na(1 / output$sum_theta_over_docs) &
        !is.infinite(1 / output$sum_theta_over_docs))]
  } else {
    # keep all docs.
    theta               <- t(output$theta)
    covariates          <- output$covariates
    sum_theta_over_docs <- output$sum_theta_over_docs
  }


  # handle spaces and odd characters in column names
  colnames(covariates) <-  make.names(colnames(covariates))


  # Assume that if there is an intercept that it is provided by the user.
  if (is.null(formula) && model == "BETA") {
    factors <- colnames(covariates)
    formula <- stats::as.formula(paste("y~", paste(
      paste(factors, collapse = "+"),
      "-1 |",
      paste(factors, collapse = "+"),
      "-1"
    )))
  }
  if (is.null(formula) && model == "GAM") {
    factors <- colnames(covariates)
    formula <- stats::as.formula(paste("y~", paste(paste0(
      "s(", factors, ")"
    ), collapse = "+")))
  }

  #  Normalization of theta by it's row sums
  denominator <- sum_theta_over_docs
  normalize <- function(x, denominator) {
    if (is.matrix(x) && nrow(x) == 1) {
      x <- t(x)
    }
    diag(1 / denominator) %*% x
  }

  ##### set up a data frame for regression
  ##### fit a linear model on all topics using specified covariates

  if ((min(theta) <= 0 ||
         max(theta) >= 1) &&
         model %in% c("BETA", "GAM")) {
    # assume that normalization has not yet occurred
    theta_nonzero <- normalize(theta, denominator = denominator)
  }
  # in the off chance that this is still a problem try this:
  if ((min(theta) <= 0 ||
        max(theta) >= 1) &&
        model %in% c("BETA", "GAM")) {
    # increase all values by a tenth of fractional occurrence of a word
    # within a topic:
    # fractional occurrence = minimum of 1/nrow or the smallest nonzezro entry
    # of the column / 1000.
    theta_nonzero <- normalize(theta + min(1 / ncol(theta),
                                           min(theta[theta > 0] / 1000)),
                               denominator = denominator +
                                 2 * min(1 / ncol(theta),
                                         min(theta[theta > 0] / 1000)))
    warning(
      paste(
        " Some topic probabilities are 0 or 1. Increasing all
                  probabilities by (1/ncol.\n) and then dividing  by
                  '(sum_theta_over_docs' + 2/ncol.\n)",
        " Minimum was ",
        round(min(theta), 7),
        " Minimum is now ",
        round(min(theta_nonzero), 7),
        " Maximum was ",
        round(max(theta), 7),
        " Maximum is now ",
        round(max(theta_nonzero), 7)
      )
    )

  } else {
    theta_nonzero <- normalize(theta, denominator = denominator)
  }
  if (is.null(obs_weights)) {
    if (model == "OLS") {
      if (return_just_coefs) {
        beta <- stats::coef(stats::lm.fit(x = as.matrix(covariates),
                                          y = theta_nonzero))
        beta <- t(beta)
        rownames(beta) <- topics
      } else {
        # return the model output
        beta <- stats::lm.fit(x = covariates, y = theta_nonzero)
        colnames(beta$coefficients) <- topics
      }#end of using bootstrap T/F
    } else {
      if (model == "BETA") {
        # use a Beta regression model
        if (return_just_coefs) {
          beta <- matrix(NA,
                         nrow = ncol(theta_nonzero),
                         ncol = 2 * ncol(covariates))
          colnames(beta) <- c(paste0("mean.", colnames(covariates)),
                              paste0("precision.", colnames(covariates)))
          rownames(beta) <- topics
          for (thetaindex in seq_len(ncol(theta_nonzero))) {
            fail <- 1
            while (fail != 0 && fail < 10) {
              tryCatch({
                if (fail == 1) {
                  data =  data.frame(theta_nonzero[,thetaindex]
                                     ,covariates)
                }else{
                  # reconstruct theta_nonzero but pull it inwards away from boundaries.
                  data =  data.frame(
                    normalize(theta[,thetaindex] +
                                (fail-1)*min(1/ncol(theta[,thetaindex]),
                          min(theta[theta[,thetaindex] > 0,thetaindex] / 1000)),
                          denominator = denominator +
                            (fail-1+1)*min(1/ncol(theta[, thetaindex]),
                                           min(theta[theta[,thetaindex] > 0,
                                                     thetaindex] / 1000))
                    ),covariates)
                }
                colnames(data)  = c("y",colnames(covariates))
                beta[thetaindex,] <-
                  betareg::betareg(formula, data = data,
                                   link = link,
                                   link.phi = link.phi,  # link.phi is for the dispersion in betaregression
                                   type = type,
                                   control = betareg::betareg.control(
                                     fsmaxit = 10000)) |>
                  coefficients()|>
                  unlist()
                fail <- 0 #if it works
              },# end trycatch expression
              error = function(e){
                fail <<-fail+1
                cat(fail);
                cat(" It'll be ok... Sometimes beta-type regressions
                    fail because the smallest value is too close to zero.
                    Let's increase the smallest value and try again. \n ")},
              finally= {# completed
                # this was more useful when troubleshooting.  Otherwise it's a lot of output when in bootstrap
                # cat(paste("\n Completed topic", thetaindex))
              }
              )
            }# end while
          }

        }else{# return the betareg model output and not just the coefficients
          beta = list()
          for(thetaindex in seq_len(ncol(theta_nonzero))){
            fail = 1
            while (fail != 0 && fail < 10) {
              tryCatch({
                if (fail == 1) {
                  data =  data.frame(theta_nonzero[, thetaindex]
                                     ,covariates)
                } else {
                  # reconstruct theta_nonzero but pull it inwards away
                  # from boundaries.
                  data =  data.frame(
                    normalize(theta[,thetaindex] +
                        (fail-1)*min(1/ncol(theta[, thetaindex]) ,
                        min(theta[theta[, thetaindex] > 0, thetaindex] / 1000)),
                          denominator = denominator +
                          (fail-1+1)*min(1/ncol(theta[,thetaindex]) ,
                          min(theta[theta[, thetaindex] > 0, thetaindex] / 1000))
                    ),covariates)
                }
                colnames(data)  = c("y",colnames(covariates))
                beta[[thetaindex]] <-
                  betareg::betareg(formula, data = data,
                                   link = link,
                                   link.phi = link.phi,
                                   # the dispersion in betaregression
                                   type = type,
                                   control = betareg::betareg.control(
                                     fsmaxit = 10000))
                fail <- 0 #if it works
              },# end trycatch expression
              error = function(e){
                fail <<- fail+1
                cat(fail)
                cat(" It'll be ok... Sometimes beta-type regressions fail
                    because the smallest value is too close to zero.
                    Let's increase the smallest value and try again. \n ")},
              finally= {# completed
                # this was more useful when troubleshooting.
                # Otherwise it's a lot of output when in bootstrap
                # cat(paste("\n Completed topic", thetaindex))
              }
              )
            }# end while
            names(beta)[thetaindex] = topics[thetaindex]
          }

        }#end of return_just_coefs or the full model output
        #(typically if not using bootstrap)
      } else {# model == GAM
        if (return_just_coefs) {
          pred_X_vals = unique(covariates)
          # inferring what is meant here,
          # let's assume that it means predicting the GAM at
          # the input covariate values.
          if (is.null(covariates |> dim())) {
            nrow_x = length(unique(covariates))
            ncol_x = 1
          } else {
            nrow_x = nrow(unique(covariates))
            ncol_x = ncol(unique(covariates))
          }
          # make a place to put the predicted values.
          beta = matrix(NA, nrow = ncol(theta_nonzero), ncol = nrow_x)
          if(ncol_x==1){
            colnames(beta) = apply(pred_X_vals, 1, function(x) {paste0("X.",x)})
          }else{
            colnames(beta) = apply(pred_X_vals |> matrix(ncol = 1),
                                   1,
                                   function(x) {paste0("X.", x, collapse=".")})
          }
          rownames(beta) = topics
          for(thetaindex in seq_len(ncol(theta_nonzero))){
            fail = 1
            while (fail !=0 && fail < 10) {
              tryCatch({
                if(fail == 1) {
                  data =  data.frame(theta_nonzero[,thetaindex]
                                     ,covariates)
                } else {
                  # reconstruct theta_nonzero
                  # but pull it inwards away from boundaries.
                  data =  data.frame(
                    normalize(theta[, thetaindex] +
                          (fail - 1) * min(1/ncol(theta[,thetaindex]) ,
                        min(theta[theta[,thetaindex] > 0, thetaindex] / 1000)),
                              denominator = denominator +
                          (fail - 1 + 1) * min(1 / ncol(theta[, thetaindex]) ,
                              min(theta[theta[,thetaindex]>0,thetaindex]/1000))
                    ),covariates)
                }
                colnames(data)  = c("y",colnames(covariates))
                beta[thetaindex,] <- mgcv::gam(formula,
                                               family = betar,
                                               link = link,
                                               data = data)|>
                  mgcv::predict.gam(newdata = pred_X_vals)
                fail <- 0 # if it works, then this will
                # break out of the while loop.
              },# end trycatch expression,
              # if fails, then push fail back up to 1 from zero
              error = function(e){fail <<-fail+1
              cat(fail)
              cat("\n It'll be ok... Sometimes beta-type regressions
                  fail because the smallest value is too close to zero.
                  Let's increase the smallest value and try again. \n ")},
              finally= {# completed
                # this was more useful when troubleshooting.
                # Otherwise it's a lot of output when in bootstrap
                # cat(paste("\n Completed topic", thetaindex))
              }
              )
            }#end while
          }# end for loop
          if(ncol_x==1){
            beta = rbind(beta,X_pred_vals = c(as.matrix(pred_X_vals)))
          }
        }else{#return the whole model
          beta = list()
          for(thetaindex in 1:ncol(theta_nonzero)){
            fail = 1
            while(fail!=0 & fail < 10){
              tryCatch({
                if(fail == 1){
                  data =  data.frame(theta_nonzero[,thetaindex]
                                     ,covariates)
                }else{
                  # reconstruct theta_nonzero but pull it inwards away from boundaries.
                  data =  data.frame(
                    normalize(theta[,thetaindex] +
                      (fail-1)*min(1/ncol(theta[, thetaindex]) ,
                      min(theta[theta[, thetaindex] > 0, thetaindex] / 1000)),
                              denominator = denominator +
                    (fail-1+1)*min(1/ncol(theta[, thetaindex]),
                                   min(theta[theta[,thetaindex] > 0,
                                             thetaindex] / 1000))
                    ),covariates)
                }
                colnames(data)  = c("y",colnames(covariates))
                beta[[thetaindex]] <- mgcv::gam(formula,
                                                family = betar,
                                                  link = link,
                                                data = data)
                fail <- 0 # if it works, then this
                # will break out of the while loop.
              },# end trycatch expression,
              # if fails, then push fail back up to 1 from zero
              error = function(e){
                fail <<-fail+1
                cat(fail)
                cat(" It'll be ok... Sometimes beta-type regressions fail
                    because the smallest value is too close to zero.
                    Let's increase the smallest value and try again. \n ")},
              finally= {# completed
                # this was more useful when troubleshooting.
                # Otherwise it's a lot of output when in bootstrap
                # cat(paste("\n Completed topic", thetaindex))
              }
              )
            }# end while
            names(beta)[thetaindex] = topics[thetaindex]
          }# end forloop

          # end of GAM with returning the whole model.
        }
      }# end of GAM
    }# end of model choice

  }else{ # with weights
    if(model == "OLS"){
      if(return_just_coefs){
        beta = stats::coef(stats::lm.fit(x = as.matrix(covariates),
                                         y = theta_nonzero,
                                         weights = obs_weights))
        beta = t(beta)
        rownames(beta) = topics
      }else{# return the lm model output
        beta = stats::lm.fit(x = covariates,
                             y = theta_nonzero,
                             weights = obs_weights)
        colnames(beta$coefficients) = topics
      }

    }else{
      if(model == "BETA"){
        # use a Beta regression model with weights
        zero_block = matrix(0, nrow(covariates) - nrow(theta), ncol(theta))
        theta_nonzero = rbind(theta_nonzero, zero_block)
        obs_weights = c(obs_weights, rep(1, nrow(zero_block)))
        if(return_just_coefs){#
          beta = matrix(NA,
                        nrow = ncol(theta_nonzero),
                        ncol = 2*length(covariates))
          colnames(beta) = c(paste0("mean.", colnames(covariates)),
                             paste0("precision.", colnames(covariates)))
          rownames(beta) = topics
          for(thetaindex in 1:ncol(theta_nonzero)){
            fail = 0
            while(fail!=0 & fail < 10){
              fail = 1
              tryCatch({
                cat(paste0("working on ", thetaindex))
                beta[thetaindex,] =
                  betareg::betareg(formula, data = data,
                                   link = link,
                                   weights = obs_weights,
                                   link.phi = link.phi,
                                   # for the dispersion in betaregression
                                   type = type,
                                   control = betareg::betareg.control(
                                     fsmaxit = 10000)) |>
                  coefficients() |>
                  unlist()

                fail = -1 #it works

              },error = function(e){print(fail)},
              finally= {
                if(all(is.na(beta[thetaindex,]))){
                  if(fail != -1){
                    cat(paste(fail, "fail for index ",
                              thetaindex, " epsilon = ",
                              min(theta_nonzero[,thetaindex])))}
                  fail <<- fail + 1;
                  #if it worked now fail = 0,
                  # if it didn't work then fail is growing 2+
                  theta_nonzero[, thetaindex] <<- theta_nonzero[,thetaindex] +
                    fail*min(theta_nonzero[, thetaindex])*.5;
                }
              }
              )
            }

          }
        }else{
          beta = list()
          for(thetaindex in 1:ncol(theta_nonzero)){
            fail = 0
            while(fail!=0 & fail < 10){
              fail = 1
              tryCatch({
                if(fail == 1){
                  data =  data.frame(theta_nonzero[,thetaindex]
                                     ,covariates)
                }else{
                  # reconstruct theta_nonzero but
                  # pull it inwards away from boundaries.
                  data =  data.frame(
                    normalize(theta[,thetaindex] +
                                (fail - 1) * min(1 / ncol(theta[,thetaindex]) ,
                        min(theta[theta[, thetaindex] > 0, thetaindex] / 1000)),
                              denominator = denominator +
                                (fail - 1 + 1 ) * min(
                                  1 / ncol(theta[, thetaindex]),
                                  min(theta[theta[, thetaindex] > 0,
                                            thetaindex] / 1000))
                    ),covariates)
                }
                colnames(data)  = c("y",colnames(covariates))
                cat(paste0("working on ", thetaindex))
                beta[thetaindex] =
                  betareg::betareg(formula, data = data,
                                   link = link,
                                   weights = obs_weights,
                                   link.phi = link.phi,
                                   # for the dispersion in betaregression
                                   type = type,
                                   control = betareg::betareg.control(
                                     fsmaxit = 10000))

                fail = -1 #it works

              },error = function(e){print(fail)},
              finally= {
                if(all(is.na(beta[[thetaindex]]))){
                  if(fail != -1){
                    cat(paste(fail,
                              "fail for index ",
                              thetaindex,
                              " epsilon = ",
                              min(theta_nonzero[, thetaindex])))}
                  fail <<- fail + 1; #if it worked now fail = 0,
                  # if it didn't work then fail is growing 2+
                  theta_nonzero[,thetaindex] <<- theta_nonzero[, thetaindex] +
                    fail * min(theta_nonzero[, thetaindex]) * .5
                }
              }
              )
            }
            names(beta)[thetaindex] = topics[thetaindex]
          }
          # return the whole regression model
        }
      }else{# model == GAM
        if(return_just_coefs){
          pred_X_vals = unique(covariates)
          # inferring what is meant here,
          # let's assume that it means predicting the GAM at
          # the input covariate values.
          if(is.null(covariates|> dim())){
            nrow_x <- length(unique(covariates))
            ncol_x <- 1
          }else{
            nrow_x <- nrow(unique(covariates))
            ncol_x <- ncol(unique(covariates))
          }
          # make a place to put the predicted values.
          beta = matrix(NA, nrow = ncol(theta_nonzero), ncol = nrow_x)
          if (ncol_x == 1) {
            colnames(beta) <- apply(pred_X_vals, 1,
                                   function(x) {
                                     paste0("X.",x)
                                     })
          }else{
            colnames(beta) = apply(pred_X_vals |> matrix(ncol = 1),
                                   1,
                                function(x){
                                  paste0("X.", x, collapse = ".")
                                  })
          }
          rownames(beta) = topics
          for (thetaindex in seq_len(theta_nonzero)) {
            fail = 1
            while(fail != 0 && fail < 10) {
              tryCatch({
                data <- cbind(theta_nonzero[, thetaindex] + (fail - 1) *
                               min(theta_nonzero[, thetaindex]) * .5,
                             covariates) |> as.data.frame()
                colnames(data)  <- c("y", colnames(covariates))
                beta[thetaindex, ] <- mgcv::gam(formula,
                                                family = betar,
                                                link = link,
                                                data = data,
                                                weights) |>
                  mgcv::predict.gam(newdata = pred_X_vals)
                fail <- 0 # if it works, then this will break out of
                # the while loop.
              }, # end trycatch expression,
              # if fails, then push fail back up to 1 from zero
              error = function(e) {
                fail <<- fail + 1
                cat(fail)
                cat("\n It'll be ok... Sometimes beta-type regressions fail
                    because the smallest value is too close to zero.  Let's
                    increase the smallest value and try again. \n ")},
              finally = {# completed
              }
              )
            }#end while
          }# end for loop
          if(ncol_x == 1) {
            beta = rbind(beta, X_pred_vals = c(as.matrix(pred_X_vals)))
          }
        }else{#return the whole model
          beta = list()
          for(thetaindex in 1:ncol(theta_nonzero)){
            fail <- 1
            while(fail != 0 && fail < 10) {
              tryCatch({
                data <- cbind(theta_nonzero[, thetaindex] +
                            (fail - 1) * min(theta_nonzero[, thetaindex]) * .5,
                            covariates) |>
                  as.data.frame()
                colnames(data)  <- c("y", colnames(covariates))
                beta[[thetaindex]] <- mgcv::gam(formula,
                                          family = betar,
                                          link = link,
                                          data = data,
                                          weights = obs_weights)
                fail <- 0 # if it works,
                # then this will break out of the while loop.
              }, # end trycatch expression,
              # if fails, then push fail back up to 1 from zero
              error = function(e) {
                fail <<- fail + 1
                cat(fail)
                cat(" It'll be ok... Sometimes beta-type regressions fail
                because the smallest value is too close to zero.
                Let's increase the smallest value and try again. \n ")},
              finally = {# completed
                cat(paste("\n Completed topic", thetaindex))
              }
              )
            }# end while
            names(beta)[thetaindex] = topics[thetaindex]
          }# end forloop
          # end of GAM with returning the whole model.
        }
      }# end of GAM
    }# end of model choice
  }# end of obs weights
  ##### return
  class(beta) <- c("BRETT_ouput", class(beta))
  return(beta)
}

#' @title Bootstrap Regression
#'
#' @description Bootstrap regression coefficients in order to
#' estimate their sampling distribution.
#'
#' @param output An object of class nmf_output
#'
#' @param samples The number of bootstrap samples to use. If set to 1
#' then return the full regression model output.  If >1, then collect the
#' coefficients for bootstrap.
#'
#' @param obs_weights weights for observations to be passed to
#' get_regression_coefs
#'
#' @param model choose c("BETA", "GAM", "OLS") for OLS, beta regression or a
#' Generalized Additive model with family beta.  Default is "BETA".
#' OLS is only useful if the covariates are categorical, to be
#' passed to get_regression_coefs
#'
#' @param return_just_coefs returns the coefficients as opposed to the full
#' betaregression output to be passed to get_regression_coefs
#'
#' @param topics is a vector of the subset of topics of interest.  Default is to
#' perform regression on all topics.
#'
#' @param formula of the form Y = X |Z  for the model E(Y) = XB for the mean and
#' var(Y)=Zß for the precision to be passed to get_regression_coefs
#'
#' @param link is the link function for the betaregression GLM to be passed to
#' get_regression_coefs or when using GAM
#'
#' @param link.phi is the link function for the precision to be passed to
#' get_regression_coefs, only for betaregression
#'
#' @param na.rm remove the problematic 'no topic' documents now so that we don't
#' end up with a bootstrap sample of all NAs
#'
#' @param type  is one of ML, BR, BC for maximum likelihood, bias reduces, or
#' bias corrected estimates to be passed to get_regression_coefs
#'
#'
#' @return Most of the original call values and a boot_reg list containing
#' matrices/vectors, each of which contains regression coefficients produced by
#' get_regression_coefs(). Each list element corresponds to a bootstrap sample.
#' Combining a particular element across bootstrap iterates estimates the
#' sampling distribution of the associated estimator; see brett_plot() and
#' create_error_bars().
#'
#' @import dplyr
#' @importFrom stats coefficients  weights
#' @importFrom utils head tail
#'
#' @examples
#' my_input <- create_input(Romeo_and_Juliet_tdm,
#'                         rownames(Romeo_and_Juliet_tdm),
#'                         topics = 10,
#'                         covariates = acts)
#' my_output <- solve_nmf(my_input)
#' boot_samples = boot_reg(my_output, samples = 100, model = "OLS")
#'
#' @seealso [boot_reg_stratified()]
#'
#' @export
boot_reg <- function(output,
                    samples,
                    obs_weights = NULL,
                    model = c("BETA", "GAM", "OLS"),
                    return_just_coefs = TRUE,
                    topics = NULL,
                    formula = NULL,
                    link = "logit",
                    link.phi = "log",
                    type = "ML",
                    na.rm = TRUE) {
  ##### check input types/whether covariates are specified
  if (!inherits(output, "nmf_output")) {
    stop("Output must be of class nmf_output.")
  }
  if (is.null(output$covariates)) {
    stop("No design matrix specified.")
  }
  samples <- as.integer(samples)
  if (samples <= 0) {
    stop("Samples must be a positive integer.")
  }

  ##### set up matrices for regression/list for return value
  if (is.null(topics)) {
    topics <- output$anchors
  }
  ##### set up matrices for regression/return value
  # select the rows corresponding to topics selected:
  if (na.rm) {
    #remove the problematic 'no topic' documents now so that we don't end up
    # with a bootstrap sample of all NAs
    theta              <- output$theta[which(output$anchors %in% topics), which(
      !is.na(1 / output$sum_theta_over_docs) &
        !is.infinite(1 / output$sum_theta_over_docs)
    )]
    covariates         <- output$covariates[which(
      !is.na(1 / output$sum_theta_over_docs) &
        !is.infinite(1 / output$sum_theta_over_docs)
    ), ] |> as.data.frame()
    colnames(covariates) <- colnames(output$covariates)
    sum_theta_over_docs <- output$sum_theta_over_docs[which(
      !is.na(1 / output$sum_theta_over_docs) &
        !is.infinite(1 / output$sum_theta_over_docs)
    )]
  } else {
    # keep all docs.
    theta               <- output$theta
    covariates          <- output$covariates
    sum_theta_over_docs <- output$sum_theta_over_docs
  }
  to_return <- list()

  ##### use a while loop; this is because bootstrap sample
  # may not include all factor levels
  ##### throw out these samples
  ##### is there a better way of dealing with this?

  for (i in 1:samples) {
    ##### counter for samples where not all factor levels are present
    ##### counter for list elements (exclude bad samples in output)
    bad_samples <- 0
    if (i == 1) {
      j <- 0
    }

    ##### produce bootstrap sample and form associated theta and covariate
    # this is written so that a sum-to-one constraint is automatically pushed
    # through if it exists
    # note that sum to one constraints are incompatible with beta regression
    # since it requires an observed
    constraint_block <- utils::tail(covariates, nrow(covariates) - ncol(theta))

    covariate_block <- utils::head(covariates, ncol(theta))
    boot_docs <- sample(seq_len(ncol(theta)), replace = TRUE)
    boot_theta <- theta[, boot_docs]
    # put the constraint back if it exists:
    zero_block <- matrix(0, nrow(theta), nrow(covariates) - ncol(theta))
    boot_theta <- cbind(boot_theta, zero_block)

    boot_covariates <- covariate_block[boot_docs, ] |> as.data.frame()
    colnames(boot_covariates) <- colnames(covariate_block)
    boot_covariates <- rbind(boot_covariates, constraint_block)

    ##### check if all factor levels are used for appropriate variables
    ##### jump to next iteration if not
    if (0 %in% colSums(boot_covariates)) {
      bad_samples <- bad_samples + 1
      next
    }
    j <- j + 1

    ##### create a new nmf_output object with bootstrapped theta and covariate
    boot_output                     <- output
    boot_output$theta               <- boot_theta
    boot_output$covariates          <- boot_covariates
    boot_output$anchors             <- topics
    boot_output$sum_theta_over_docs <- sum_theta_over_docs[boot_docs]
    ##### call get_regression_coefs and append list element
    boot_coefs <- get_regression_coefs(
      output = boot_output,
      obs_weights = obs_weights,
      model = model,
      return_just_coefs = return_just_coefs,
      formula = formula,
      link = link,
      link.phi = link.phi,
      type = type,
      topics = topics,
      na.rm = na.rm
    )
    to_return[[j]] <- boot_coefs

    ##### progress of iterations
    if (i %% 50 == 0) {
      cat(i, " of ", samples, " bootstrap samples complete.\n")
    }

  }

  ##### print warning message if bad_samples > 0
  if (bad_samples > 0) {
    cat(
      "Warning: not all factor levels present in ",
      bad_samples,
      " bootstrap samples.
        These samples have been excluded from the output --
    inferences may be affected as a result.\n"
    )
  }

  ##### return
  toreturn <- list(
    boot_reg = to_return,
    obs_weights = obs_weights,
    model = model,
    return_just_coefs = return_just_coefs,
    formula = formula,
    link     = ifelse(model == "OLS", NA, link),
    link.phi = ifelse(model == "OLS", NA, link.phi),
    type = ifelse(model == "BETA", type, NA),
    topics = topics,
    na.rm = na.rm
  )
  class(toreturn) <- c("BRETT_ouput", "list")
  return(toreturn)
}


#' @title Stratified Bootstrap Regression
#'
#' @description Bootstrap regression coefficients in order to
#' estimate their sampling distribution. Stratified bootstrap is used to
#' maintain a constant number of individuals in each group when groups are
#' categorical
#'
#' @param output An object of class nmf_output
#'
#' @param samples The number of bootstrap samples to use.
#'
#' @param parallel number of cores to use.  Skip the parallel overhead if
#' using 1 core.
#'
#' @param obs_weights weights for observations to be passed to
#' get_regression_coefs
#'
#' @param model choose c("BETA", "GAM", "OLS") for OLS, beta regression or a
#' Generalized Additive model with family beta.  Default is "BETA".
#' OLS is only useful if the covariates are categorical, to be passed
#' to get_regression_coefs
#'
#' @param return_just_coefs returns the coefficients as opposed to the full
#' betaregression output to be passed to get_regression_coefs
#'
#' @param topics is a vector of the subset of topics of interest.  Default is to
#' perform regression on all topics.
#'
#' @param formula of the form Y = X |Z  for the model E(Y) = XB for the mean and
#' var(Y)=Zß for the precision to be passed to get_regression_coefs
#'
#' @param link is the link function for the betaregression GLM to be passed to
#' get_regression_coefs
#'
#' @param link.phi is the link function for the precision to be passed to
#' get_regression_coefs
#'
#' @param type  is one of ML, BR, BC for maximum likelihood, bias reduces, or
#' bias corrected estimates to be passed to get_regression_coefs
#'
#' @param na.rm remove the problematic 'no topic' documents now so that we don't
#' end up with a bootstrap sample of all NAs
#'
#' @return Most of the original call values and a boot_reg list containing
#' matrices/vectors, each of which contains regression coefficients produced by
#' get_regression_coefs(). Each list element corresponds to a bootstrap sample.
#' Combining a particular element across bootstrap iterates estimates the
#' sampling distribution of the associated estimator; see brett_plot() and
#' create_error_bars().
#'
#' @import dplyr
#' @importFrom stats coefficients  weights
#' @importFrom utils head tail
#' @importFrom prodlim row.match
#' @import doParallel
#' @import foreach
#' @importFrom foreach %dopar%
#' @import parallel
#'
#' @examples
#' my_input <- create_input(Romeo_and_Juliet_tdm,
#'                       rownames(Romeo_and_Juliet_tdm),
#'                       topics = 10,
#'                        covariates = acts)
#' my_output <- solve_nmf(my_input)
#' boot_samples <- boot_reg_stratified(my_output, samples=100,
#'                       model = "OLS", parallel =1)
#'
#' @seealso [boot_reg()]
#'
#' @export

boot_reg_stratified <- function(output,
                               samples,
                               parallel = 4,
                               obs_weights = NULL,
                               model = c("BETA", "GAM", "OLS"),
                               return_just_coefs = TRUE,
                               topics = NULL,
                               formula = NULL,
                               link = "logit",
                               link.phi = "log",
                               type = "ML",
                               na.rm = TRUE) {
  ##### check input types/whether covariates are specified
  if (!inherits(output, "nmf_output")) {
    stop("Output must be of class nmf_output.")
  }
  if (is.null(output$covariates)) {
    stop("No design matrix specified.")
  }
  samples <- as.integer(samples)
  if (samples <= 0) {
    stop("Samples must be a positive integer.")
  }
  parallel <- as.integer(parallel)
  if (parallel <= 0) {
    stop("parallel must be a positive integer.")
  }

  ##### set up matrices for regression/list for return value
  if (is.null(topics)) {
    topics <- output$anchors
  }

  ##### set up matrices for regression/return value
  # select the rows corresponding to topics selected:
  if (na.rm) {
    #remove the problematic 'no topic' documents now so that we don't end up
    # with a bootstrap sample of all NAs
    theta              <- output$theta[which(output$anchors %in% topics), which(
      !is.na(1 / output$sum_theta_over_docs) &
        !is.infinite(1 / output$sum_theta_over_docs)
    )]
    covariates         <- output$covariates[which(
      !is.na(1 / output$sum_theta_over_docs) &
        !is.infinite(1 / output$sum_theta_over_docs)
    ), ] |> as.data.frame()
    colnames(covariates) <- colnames(output$covariates)
    sum_theta_over_docs <- output$sum_theta_over_docs[which(
      !is.na(1 / output$sum_theta_over_docs) &
        !is.infinite(1 / output$sum_theta_over_docs)
    )]
  } else {
    # keep all docs.
    theta               <- output$theta
    covariates          <- output$covariates
    sum_theta_over_docs <- output$sum_theta_over_docs
  }


  to_return <- list()


  # this is written so that the constraint is automatically
  # pushed through if it exists
  constraint_block <- utils::tail(covariates, nrow(covariates) - ncol(theta))
  zero_block <- matrix(0, nrow(theta), nrow(covariates) - ncol(theta))
  covariate_block <- utils::head(covariates, ncol(theta))
  # identify the categories (remove the intercept)
  categorical_groups <- unique(covariate_block)
  groups <- apply(covariate_block, 1, function(x) {
    paste0(names(x[x == 1]), collapse = "_")
  })

  group_levels <- unique(groups)
  rownames(categorical_groups) <- group_levels
  group_count <- NULL
  for (group in seq_along(categorical_groups)) {
    group_count[group] <- sum(prodlim::row.match(
      covariate_block,
      categorical_groups[group, ], 0))
  }
  names(group_count) <- group_levels

  if (parallel == 1) {
    for (i in 1:samples) {
      ##### produce bootstrap sample and form associated theta and covariate

      boot_docs <- rep(NA, ncol(theta))
      start <- 0
      sampled_inds <- NULL
      for (group_strat in group_levels) {
        ####NOTE THIS NEEDS RETHINKING WHEN THERE ARE MULTIPLE COVARIATES
        sampled_indices <- which(prodlim::row.match(
          covariate_block,
          categorical_groups[group_strat, ], 0) == 1) [
                          sample(1:group_count[[group_strat]], replace = TRUE)]
        sampled_inds <- c(sampled_inds, sampled_indices)
        boot_docs[start + (1:group_count[[group_strat]])] <- sampled_indices
        start <- start + group_count[[group_strat]]
      }

      boot_theta <- theta[, boot_docs]
      # put the constraint back if it exists:
      boot_theta <- cbind(boot_theta, zero_block)
      boot_covariates <- covariate_block[boot_docs, ] |> as.data.frame()
      colnames(boot_covariates) <- colnames(covariate_block)
      boot_covariates <- rbind(boot_covariates, constraint_block)

      ##### create a new nmf_output object but with
      # bootstrapped theta and covariate
      boot_output                     <- output
      boot_output$theta               <- boot_theta
      boot_output$covariates          <- boot_covariates
      boot_output$anchors             <- topics
      boot_output$sum_theta_over_docs <- sum_theta_over_docs[boot_docs]
      boot_coefs <- get_regression_coefs(
        boot_output,
        obs_weights = obs_weights,
        model = model,
        return_just_coefs = return_just_coefs,
        formula = formula,
        link = link,
        link.phi = link.phi,
        type = type,
        topics = topics,
        na.rm = na.rm
      )
      to_return[[i]] <- boot_coefs

      ##### progress of iterations
      if (i %% 50 == 0) {
        cat(i, " of ", samples, " bootstrap samples complete.\n")
      }
    }
  } else {
    #parallel setup chunk 1 start
    #create the cluster
    my_cluster <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(cl = my_cluster)
    # defaults to fork on linux and psock on windows
    #parallel setup chunk 1 end

    #parallel loop start
    to_return <- foreach::foreach(
      i = 1:samples,
      .export = c("get_regression_coefs",
                  "betareg::betareg",
                  "betareg::betareg.control")
    ) %dopar% {
      ##### produce bootstrap sample and form associated theta and covariate

      boot_docs <- rep(NA, ncol(theta))
      start <- 0
      sampled_inds <- NULL
      for (group_strat in group_levels) {
        sampled_indices <- which(groups == group_strat)[
               sample(1:group_count[[group_strat]], replace = TRUE)]
        sampled_inds <- c(sampled_inds, sampled_indices)
        boot_docs[start + (1:group_count[[group_strat]])] <- sampled_indices
        start <- start + group_count[[group_strat]]
      }

      boot_theta <- theta[, boot_docs]
      # put the constraint back if it exists:
      boot_theta <- cbind(boot_theta, zero_block)
      boot_covariates <- covariate_block[boot_docs, ] |> as.data.frame()
      colnames(boot_covariates) <- colnames(covariate_block)
      boot_covariates <- rbind(boot_covariates, constraint_block)

      ##### create a new nmf_output object but with
      # bootstrapped theta and covariate
      boot_output                     <- output
      boot_output$theta               <- boot_theta
      boot_output$covariates          <- boot_covariates
      boot_output$anchors             <- topics
      boot_output$sum_theta_over_docs <- output$sum_theta_over_docs[boot_docs]
      return(
        get_regression_coefs(
          boot_output,
          obs_weights = obs_weights,
          model = model,
          return_just_coefs = return_just_coefs,
          formula = formula,
          link = link,
          link.phi = link.phi,
          type = type,
          topics = topics,
          na.rm = na.rm
        )
      )
    }
  }
  ##### return
  toreturn <- list(
    boot_reg = to_return,
    obs_weights = obs_weights,
    model = model,
    return_just_coefs = return_just_coefs,
    formula = formula,
    link     = ifelse(model == "OLS", NA, link),
    link.phi = ifelse(model == "OLS", NA, link.phi),
    type = type,
    topics = topics,
    na.rm = na.rm,
    type = "stratified bootstrap"
  )
  class(toreturn) <- c("BRETT_ouput", "list")
  return(toreturn)
}




#' @title Plots of Output
#'
#' @description Plot smoothed histograms of regression effects for a topic.
#' Use case 1: producing density plots of bootstrap OLS samples with
#' respect to discrete coefficients.  That's the only time OLS will be valid
#' anyways.
#' Use case 2: producing quantile plots for the fitted beta distribution with
#' respect to covariate coefficients. Works best when the coefficients are
#' discrete or when "newdata" is provided.
#'
#' @param brett_object A list of bootstrap samples, as produced by boot_reg() or
#' model coefficients from get_regression_coefficients()
#'
#' @param topic The topic to use as the response variable,
#'  labeled by its anchor word.
#'
#' @param newdata provided if the bootsamples are regression coefficients and we
#' wish to plot the fitted values at 'newdata' points.  newdata should be a
#' data.frame where the N rows are N new observations with values across the
#' columns.  The rownames are used for plotting as appropriate.
#'
#' @param model Used if 'boot_samples' is not from boot_reg
#' and therefore doesn't contain the model information already.  Used if
#' newdata is supplied so that the  function knows how to convert the
#' bootstrap samples into predictions at 'newdata' points.
#' MAY BE REMOVED IN FUTURE
#'
#' @param link Used ONLY if 'boot_samples' is not from boot_reg and therefore
#'  doesn't contain the model information already.  Used if newdata is
#'  supplied so that the function knows how to convert the bootstrap
#'  samples into predictions at 'newdata' points.  MAY BE REMOVED IN FUTURE
#'
#' @param from_stm logical hack so that the same plots from
#' can be made from the output of STM
#'
#' @return A list as produced by ggplot2. Calling this function without
#' assignment to a variable will send a plot to the plotting window.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom betareg predict
#' @importFrom mgcv gam betar predict.gam
#' @importFrom stringr str_to_title
#'
#' @examples
#' my_input <- create_input(Romeo_and_Juliet_tdm,
#'                         rownames(Romeo_and_Juliet_tdm),
#'                         topics = 10,
#'                         covariates = acts)
#' my_output <- solve_nmf(my_input)
#' brett_object <- boot_reg(my_output, 100, model = "OLS",
#'                  topics = c("death", "love"))
#' brett_plot(brett_object, topic = "death")
#'
#' brett_object = get_regression_coefs(my_output, model = "BETA",
#'              return_just_coefs = FALSE, topics = c("death", "love"))
#' brett_plot(brett_object, topic="death") +
#'               ggplot2::coord_cartesian(ylim = c(0,5))
#'
#'
#'
#' @export
brett_plot <- function(brett_object,
                       topic,
                       newdata = NULL,
                       model = NULL,
                       link = NULL,
                       from_stm = FALSE) {
  ##### check input types
  if (!inherits(brett_object, "BRETT_ouput")) {
    stop("Input must be an object of class nmf_input.")
  }
  stopifnot("one topic at a time please" = length(topic) == 1)
  # eventually loop this.
  if (inherits(brett_object, "list")) {
    # could be OLS or BETA, with return_just_coefs = FALSE
    if (is.null(brett_object$boot_reg)) {
      # didn't just return the coefs, so plot quantiles,
      boot_reg_list <- brett_object
      bootstrapped <- from_stm
      model <- "BETA"
      if ("gam" %in% class(boot_reg_list[[topic]])) {
        # check if GAM or Beta regression was used.
        model <- "GAM"
      }
    } else {
      # if brett_object is a bootstrap object using OLS or BETA
      boot_reg_list <- brett_object$boot_reg
      model         <- brett_object$model
      link          <- brett_object$link
      bootstrapped  <- TRUE
    }
  }

  if (bootstrapped && !(topic %in% row.names(boot_reg_list[[1]]))) {
    stop("Topic not among the anchor words.")
  }
  if (!bootstrapped && !(topic %in%
                         c(names(boot_reg_list),
                           rownames(boot_reg_list[[1]])))) {
    stop("Topic not among the anchor words.")
  }

  if (is.null(newdata)) {
    # if   = Using betaregression at the observation points
    # else = Just plot the effects not fits.
    if (inherits(boot_reg_list[[topic]], "betareg")) {

      # want to plot data fit at "newdata"
      # Beta regression model was used
      # plotting the prediction and prediction interval
      newdata <- boot_reg_list[[topic]]$model
      newdata <- subset(newdata, select = -c(y)) |> unique()
      # in case you have repeated binary measures make the labels nicer:
      if (all(apply(newdata, 1, function(x) {
        all(x %in% c(0, 1))
      }))) {
        names_2_use <- apply(newdata, 1, function(x) {
          ind <- which(x != 0)

          paste(names(x)[ind], collapse = "+")
        })
        rownames(newdata) <- names_2_use
      }
      topic_plot <- beta_distribution_plot(betamodel = boot_reg_list[[topic]],
                                           topic,
                                           newdata)

    } else {
      # just plot the effects not fits.
      ##### create a data frame with covariate effects as
      # columns and bootstrap samples as rows
      num_covar <-  dim(boot_reg_list[[1]])[2]
      num_samples <-  length(boot_reg_list)
      boot_mat <-  matrix(rep(0, num_samples * num_covar), ncol = num_covar)
      for (i in seq_along(boot_reg_list)) {
        boot_mat[i, ] <-  boot_reg_list[[i]][topic, ]
      }
      boot_effect <-  as.data.frame(boot_mat)
      names(boot_effect) <-  colnames(boot_reg_list[[1]])

      ##### melt the data frame so it can be passed to ggplot
      boot_effect <-  reshape2::melt(boot_effect)
      names(boot_effect) <-  c("covar", "weight")

      ##### create and evaluate plot
      topic_plot <- ggplot2::ggplot(data = boot_effect,
                                    ggplot2::aes(
                                      x = .data$weight,
                                      y = .data$covar,
                                      height = ggplot2::after_stat(density)
                                    )) +
        ggridges::stat_density_ridges(
          fill = "lightblue",
          rel_min_height = 0.005,
          scale = 1.25
        ) +
        ggplot2::geom_vline(xintercept = 0,
                            color = "red",
                            linetype = "dashed") +
        ggplot2::labs(
          title = paste(
            "Bootstrap Distribution of Coefficients: Topic =",
            stringr::str_to_title(topic)
          ),
          x = "P(topic|X)",
          y = stringr::str_to_title(attributes(dimnames(
            brett_object[[1]]
          ))[[1]][2])
        )


    }
  } else {
    if (is.matrix(newdata)) {
      newdata <- data.frame(newdata)
      }
    if (bootstrapped) {
      if (model == "BETA") {
        # want to plot data fit at "newdata"
        # Beta regression model was used
        # plotting the bootstrap distribution
        boot_reg_colnames <- boot_reg_list[[1]] |>
          colnames()
        mean_cols         <- boot_reg_colnames |>
          grep(pattern = "mean", value = TRUE)
        # eventually do this with the precision as well
        precision_cols    <- boot_reg_colnames |>
          grep(pattern = "mean",
               value = TRUE,
               invert = TRUE)

        x_beta         <- lapply(boot_reg_list, function(x) {
          x[, mean_cols] %*% t(newdata)
        })
        fitted_values <- lapply(x_beta, function(xb) {
          exp(xb) / (1 - exp(xb))
        })

        num_covar <- dim(fitted_values[[1]])[2]
        num_samples <- length(fitted_values)
        boot_mat <- matrix(rep(0, num_samples * num_covar), ncol = num_covar)
        for (i in seq_along(fitted_values)) {
          boot_mat[i, ] <- fitted_values[[i]][topic, ]
        }
        boot_effect <- as.data.frame(boot_mat)
        names(boot_effect) <- rownames(newdata)

        ##### melt the data frame so it can be passed to ggplot
        boot_effect <- reshape2::melt(boot_effect)
        names(boot_effect) <- c("newdatapoint", "fitted")

        ##### create and evaluate plot
        topic_plot <- ggplot2::ggplot(data = boot_effect,
                                      ggplot2::aes(x = newdatapoint,
                                                   y = fitted)) +
          ggplot2::geom_boxplot(fill = "lightblue") +
          ggplot2::labs(
            title = paste(
              "Beta model fit for Topic =",
              stringr::str_to_title(topic)
            ),
            x = "Covariate value",
            y = "P(topic | new data)"
          )

      }
      if (model == "OLS") {
        # want to plot data fit at "newdata"
        # OLS model was used
        # plotting the bootstrap distribution
        fitted_values <- lapply(boot_reg_list, function(x) {
          x %*% t(newdata)
        })
        num_covar <- dim(fitted_values[[1]])[2]
        num_samples <- length(fitted_values)
        boot_mat <- matrix(rep(0, num_samples * num_covar), ncol = num_covar)
        for (i in seq_along(fitted_values)) {
          boot_mat[i, ] <- fitted_values[[i]][topic, ]
        }
        boot_effect <- as.data.frame(boot_mat)
        names(boot_effect) <- rownames(newdata)
        ##### melt the data frame so it can be passed to ggplot
        boot_effect <- reshape2::melt(boot_effect)
        names(boot_effect) <- c("newdatapoint", "fitted")
        ##### create and evaluate plot
        topic_plot <- ggplot2::ggplot(data = boot_effect,
                                      ggplot2::aes(x = newdatapoint,
                                                   y = fitted)) +
          ggplot2::geom_boxplot(fill = "lightblue") +
          ggplot2::labs(
            title = paste(
              "Bootstrap distribution for Topic =",
              stringr::str_to_title(topic)
            ),
            x = "Covariate value",
            y = "P(topic | new data)"
          )

      }
    } else {
      #not bootstrapped
      if (any(class(boot_reg_list[[1]]) == "betareg")) {
        topic_plot <- beta_distribution_plot(betamodel = boot_reg_list[[topic]],
                                             topic = topic,
                                             newdata)
      } else {
        if (model == "GAM") {
          # want to plot data fit at "newdata"
          # Beta regression model was used
          # plotting the prediction and prediction interval

          fitted_values <- mgcv::predict.gam(
            boot_reg_list[[topic]],
            newdata = newdata,
            type = "response",
            se.fit = TRUE
          )
          if (ncol(newdata) == 1) {
            fitted_values <- data.frame(
              newdata = newdata,
              fitted = fitted_values$fit,
              plus_1sefit  = fitted_values$fit +
                fitted_values$se.fit,
              minus_1sefit = fitted_values$fit -
                fitted_values$se.fit
            )
          } else {
            fitted_values <- data.frame(
              newdata = rownames(newdata),
              fitted = fitted_values$fit,
              plus_1sefit  = fitted_values$fit +
                fitted_values$se.fit,
              minus_1sefit = fitted_values$fit -
                fitted_values$se.fit
            )
          }
          colnames(fitted_values)[1] <- "newdata"
          ##### melt the data frame so it can be passed to ggplot
          fitted_values <- reshape2::melt(fitted_values, "newdata")
          names(fitted_values) <- c("newdata", "variable", "value")
          ##### create and evaluate plot
          topic_plot <- ggplot2::ggplot(data = fitted_values,
                                        ggplot2::aes(
                                          x = newdata,
                                          y = value,
                                          colour = variable
                                        )) +
            ggplot2::geom_line(lwd = 2, alpha = .5) +
            ggplot2::labs(
              title = paste(
                "GAM model fit +/- 1 SE,
                                    Topic =",
                stringr::str_to_title(topic)
              ),
              x = "Covariate value",
              y = "P(topic | new data)"
            )

        }
      }
    }
  }
  eval(topic_plot)
}



#' @title Create Error Bars for Bootstrap Inference
#'
#' @description Produce a data frame with two-sided, 95\% confidence intervals
#' based on bootstrap estimates of sampling distributions.
#'
#' @param brett_object A list of bootstrap samples, as produced by boot_reg().
#'
#' @param topicvector the selected topic(s) of interest for inference, default
#' is all topics, though this is often too much output to be useful.
#'
#' @param coverage the level of coverage for the interval, default is 95\%
#'
#' @return A data frame. Each row corresponds to a covariate and the columns
#' give the CI.
#'
#' @importFrom dplyr filter
#' @importFrom stats quantile
#' @importFrom reshape2 melt
#'
#' @examples
#' my_input = create_input(Romeo_and_Juliet_tdm,
#'                         rownames(Romeo_and_Juliet_tdm),
#'                         topics = 10,
#'                         covariates = acts)
#' my_output = solve_nmf(my_input)
#' brett_object = boot_reg(my_output, samples = 100, model = "OLS")
#' bootstrap_error_bars(brett_object, topicvector = "death")
#'
#' @export
bootstrap_error_bars <- function(brett_object,
                                 topicvector = NULL,
                                 coverage = .95) {
  ##### check inputs
  if (!inherits(brett_object, "BRETT_ouput")) {
    stop("brett_object must be a BRETT_ouput object, as outputted by
              boot_reg() or get_regression_coefs.")
  }
  if (!inherits(brett_object, "list")) {
    stop("to build the error bars we need the bootstrap model,
         rerun get_regression_coefs with the input flag
         return_just_coefs = FALSE")
  }
  if (is.null(brett_object$boot_reg)) {
    # if brett_object comes from get_regression_coefs
    boot_reg_list <- brett_object
  } else {
    # if brett_object is a list of bootstrap runs
    boot_reg_list <- brett_object$boot_reg
  }
  if (!(is.matrix(boot_reg_list[[1]]))) {
    stop("brett_object must be a list of matrices, as outputted by
              boot_reg()")
  }
  if (!(0 < coverage && coverage < 1)) {
    stop("coverage must be a numerical value in the interval (0,1).")
  }

  ##### melt and name a data frame of samples
  sample_frame <- reshape2::melt(boot_reg_list)
  colnames(sample_frame) <- c("topic", "coef", "estimate", "sample")
  if (!is.null(topicvector)) {
    sample_frame <- dplyr::filter(sample_frame, topic %in% topicvector)
  }


  ##### construct data frame w/ all possible covariate/topic combinations
  ##### fill with corresponding quantiles
  error_frame <- expand.grid(unique(sample_frame$topic),
                             unique(sample_frame$coef))
  names(error_frame) <- c("topic", "coef")
  for (i in seq_len(nrow(error_frame))) {
    subset_frame <- sample_frame[sample_frame$topic == error_frame$topic[i] &
                                   sample_frame$coef == error_frame$coef[i], ]
    error_frame$lower[i]  <- stats::quantile(subset_frame$estimate,
                                             (1 - coverage) / 2)
    error_frame$median[i] <- stats::quantile(subset_frame$estimate,
                                             .5)
    error_frame$upper[i]  <- stats::quantile(subset_frame$estimate,
                                             coverage + (1 - coverage) / 2)
  }
  colnames(error_frame)[c(3, 5)] <- paste0(c("lower", "upper"), coverage)
  ##### return
  return(error_frame)
}

#' @title Internally used plot for the predicted Beta distribution
#'
#' @description Show the predicted Beta distribution at the covariate values.
#' Intended for internal use.  Used in the function `brett_plot`
#'
#' @param betamodel the output from betaregression so that betareg::predict
#' can be used
#'
#' @param newdata Design matrix of the new data points at which to predict.
#'
#' @param topic Name of the topic to plot
#'
#' @return A plot as produced by ggplot2. Calling this function without
#' assignment to a variable will send a plot to the plotting window.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom ggplot2 geom_ribbon ggplot geom_line labs
#'
#'
#' @examples
#' my_input <- create_input(Romeo_and_Juliet_tdm,
#'                         rownames(Romeo_and_Juliet_tdm),
#'                         topics = 10,
#'                         covariates = acts)
#' my_output <- solve_nmf(my_input)
#' brett_object <- boot_reg(my_output, 100, model = "OLS",
#'                  topics = c("death", "love"))
#' brett_plot(brett_object, topic = "death")
#'
#' brett_object = get_regression_coefs(my_output, model = "BETA",
#'                        return_just_coefs = FALSE,
#'                        topics = c("death", "love"))
#' brett_plot(brett_object, topic="death") +
#'               ggplot2::coord_cartesian(ylim = c(0,5))
#'
beta_distribution_plot <- function(betamodel,
                                   topic,
                                   newdata) {
  # importing ggpubr is for the vignette
  # set up to plot the full Beta regression distribution
  increment <- .005
  xgrid <- seq(from = increment, to = 1 - increment, by = increment)
  # make predictions
  beta_fit <- betareg::predict(betamodel,
                               newdata = newdata,
                               type = "density", at = xgrid)
  toplot <- beta_fit |>
    tibble::as_tibble() |>
    dplyr::mutate(covariate = colnames(newdata)) |>
    tidyr::pivot_longer(cols = -c("covariate"),
                        names_to = "gridpoint",
                        values_to = "density") |>
    dplyr::mutate(gridpoint = as.numeric(
      gsub(.data$gridpoint, pattern = "d_", replacement = "")))

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = toplot,
                         ggplot2::aes(x = .data$gridpoint,
                                      ymax = .data$density,
                                      ymin = 0,
                                      fill = .data$covariate,
                                      colour = .data$covariate),
                         alpha = .2) +
    ggplot2::geom_line(data = toplot,
                       aes(x = .data$gridpoint,
                           y = .data$density,
                           colour = .data$covariate),
                       alpha = 1, show.legend = FALSE) +
    ggplot2::labs(title = paste("Topic =", topic),
                  x = "P(topic | new data)",
                  y = "Density")
}
