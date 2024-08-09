require(mgcv)
require(dplyr)
require(prodlim) # for the row match function in the stratified bootstrap version.
#' get_regression_coefs
#'
#' Compute OLS coefficients, fitting a linear model between a user's specified covariates and topic
#' weight. Estimates are produced simultaneously for all topics.
#'
#' @param output An object of class nmf_output.
#'
#' @param obs_weights Weights for the documents, as used by weighted least squares. Defaults to null.
#'
#' @param Model choose c("BETA", "GAM", "OLS") for OLS, beta regression or a Generalized Additive model with family beta.  Default is "BETA".   OLS is only useful if the covariates are categorical, to be passed to get_regression_coefs
#'
#' @param return_just_coefs is a logical, if TRUE then just return the coefficients, if FALSE, then return the output from the lm or betareg function.
#' 
#' @param formula is the formula to be passed into betareg or GAM.  With betaregression, formula has the form Y~ model+for+mean | model+for+dispersion 
#' 
#' @param link is the link function for the GLM for the mean when using betaregression or GAM.
#' 
#' @param link.phi is the link function for the GLM for the precision when using betaregression.
#' 
#' @param type is one of ML, BC, BR for betaregression to use Maximum Likelihood, Bias Corrected, or Bias Reduced respectively
#' 
#' @param theta_transformation transformation of theta matrix, currently only NULL or 'log' are allowed.  This will apply the "log(theta)" transformation before regression
#' 
#' @return A matrix of regression coefficients (named if column names have been specified
#' for the design matrix).
#'
#' @examples
#' neurips_input = create_input(neurips_tdm, neurips_words,
#'    topics = 10, project = TRUE, proj_dim = 500, covariates = year_bins)
#' neurips_output = solve_nmf(neurips_input)
#' get_regression_coefs(neurips_output)
#'
#' @export
get_regression_coefs = function(output, obs_weights = NULL, 
                                Model = c("BETA","GAM", "OLS"), 
                                return_just_coefs = TRUE, 
                                formula = NULL,
                                link = "logit",
                                link.phi = "log", 
                                type = "ML",
                                theta_transformation = NULL,
                                topics = NULL,
                                na.rm = TRUE){
  ##### Check input types/whether they include covariates
  if(class(output) != "nmf_output"){
    stop("Output must be of class nmf_output.")
  }
  if(is.null(output$covariates)){
    stop("No design matrix specified with output.")
  }
  if(!(is.null(obs_weights))){
    if(!(is.numeric(obs_weights))){
      stop("Document weights must be in the form of a numeric vector.")
    }
    if(sum(obs_weights >= 0) != length(obs_weights)){
      stop("Document weights must be positive or zero.")
    }
    if(length(obs_weights) != ncol(output$theta)){
      stop("Document weights and documents have differing lengths.")
    }
  }
  if(is.null(topics)){
    topics = output$anchors
  }
  if(length(Model)==3){Model = "BETA"}
  ##### set up matrices for regression/return value
  # select the rows corresponding to topics selected:
  if(na.rm){
    #remove the problematic 'no topic' documents now so that we don't end up with a bootstrap sample of all NAs
    theta              = t(output$theta[which(output$anchors%in% topics),
                                      which(!is.na(1/output$sum_theta_over_docs) & 
                                              !is.infinite(1/output$sum_theta_over_docs))])
    covariates         = output$covariates[which(!is.na(1/output$sum_theta_over_docs) & 
                                                   !is.infinite(1/output$sum_theta_over_docs)),
    ] |> as.data.frame()
    colnames(covariates) = colnames(output$covariates)
    sum_theta_over_docs = output$sum_theta_over_docs[which(!is.na(1/output$sum_theta_over_docs) & 
                                                             !is.infinite(1/output$sum_theta_over_docs))
    ]
  }else{
    # keep all docs.
    theta               = t(output$theta)
    covariates          = output$covariates
    sum_theta_over_docs = output$sum_theta_over_docs
  }  

  
  # handle spaces and odd characters in column names
  colnames(covariates) =  make.names(colnames(covariates))
  
  
  
  # Deal with formulas for betaregression, put all covariates into the mean and precision
  # Assume that if there is an intercept that it is provided by the user.
  if(is.null(formula) & Model == "BETA"){
    factors = colnames(covariates)
    formula = as.formula(paste("y~", paste(
      paste(factors, collapse="+"), "-1 |",
      paste(factors, collapse="+"),"-1"))
    )
  }
  if(is.null(formula) & Model == "GAM"){
    factors = colnames(covariates)
    formula = as.formula(
      paste("y~",
        paste(
          paste0("s(",factors,")"),
                collapse="+")
            )
      )
  }  

  #  Normalization of theta by it's row sums
  denominator = sum_theta_over_docs
  normalize = function(x,denominator){
     return( diag( 1/denominator ) %*% x)
   }
  
  

  
  ##### set up a data frame for regression
  ##### fit a linear model on all topics using specified covariates

  if(!is.null(theta_transformation)){
    if(theta_transformation == "log"){
      print("This is a work in progress")
      print("log(.) expands the values to the real line")
      print("but maps 0 --> infinity")
      print("So we're adding 1 to all counts...")
      theta_nonzero = log(theta+1)
    }else{warning("un-recognized theta_transformation; skipping it.")}
  }
  if( (min(theta)<=0 | max(theta)>=1) & Model %in% c("BETA", "GAM")){
    # increase all values by a tenth of fractional occurrence of a word within a topic:
    # fractional occurrence = minimum of 1/nrow or the smallest nonzezro entry of the column / 1000.
    theta_nonzero = normalize(theta+min(1/ncol(theta) ,min(theta[theta>0]/1000)), 
                              denominator = denominator +
                                2*min(1/ncol(theta) ,min(theta[theta>0]/1000))
                              )
    warning(paste(" beta regresssion won't work with observed {0,1} since this will result in infinite betas; rescaling valued by 'sum_theta_over_docs' + nugget of 1/ncol.\n",
                  " Minimum is now ", round(min(theta_nonzero),7),
                  " Maximum is now ", round(max(theta_nonzero),7)))
    
  }else{
    theta_nonzero = normalize(theta, denominator = denominator)
  }
  if(is.null(obs_weights)){
    
    if(Model == "OLS"){
      if(return_just_coefs){
        beta = stats::coef(stats::lm.fit(x = as.matrix(covariates), y = theta_nonzero))
        beta = t(beta)
        rownames(beta) = topics 
      }else{# return the model output
        beta = stats::lm.fit(x = covariates, y = theta_nonzero)
        colnames(beta$coefficients) = topics 
      }#end of using bootstrap T/F
    }else{
      if(Model == "BETA"){# use a Beta regression model
        if(return_just_coefs){
          beta = matrix(NA, nrow = ncol(theta_nonzero), ncol = 2*ncol(covariates))
          colnames(beta) = c(paste0("mean.", colnames(covariates)),
                           paste0("precision.", colnames(covariates)))
          rownames(beta) = topics 
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
                                  (fail-1)*min(1/ncol(theta[,thetaindex]) ,min(theta[theta[,thetaindex]>0,thetaindex]/1000)), 
                                denominator = denominator +
                                          (fail-1+1)*min(1/ncol(theta[,thetaindex]) ,min(theta[theta[,thetaindex]>0,thetaindex]/1000))
                      ),covariates)
                   }
                    colnames(data)  = c("y",colnames(covariates))
                    beta[thetaindex,] <-
                              betareg::betareg(formula, data = data,
                                           link = link,
                                           link.phi = link.phi,  # link.phi is for the dispersion in betaregression
                                           type = type,
                                           control = betareg::betareg.control(fsmaxit = 10000))|> coefficients()|> unlist()
                             fail <- 0 #if it works
                   },# end trycatch expression
                   error = function(e){fail <<-fail+1; cat(fail); cat(" It'll be ok... Sometimes beta-type regressions fail because the smallest value is too close to zero.  Let's increase the smallest value and try again. \n ")},
                   finally= {# completed
                     # this was more useful when troubleshooting.  Otherwise it's a lot of output when in bootstrap
                     # cat(paste("\n Completed topic", thetaindex))
                   }
                )
              }# end while
          }
        
      }else{# return the betareg model output and not just the coefficients
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
                              (fail-1)*min(1/ncol(theta[,thetaindex]) ,min(theta[theta[,thetaindex]>0,thetaindex]/1000)), 
                            denominator = denominator +
                              (fail-1+1)*min(1/ncol(theta[,thetaindex]) ,min(theta[theta[,thetaindex]>0,thetaindex]/1000))
                  ),covariates)
              }
                colnames(data)  = c("y",colnames(covariates))
                beta[[thetaindex]] <-
                  betareg::betareg(formula, data = data,
                                   link = link,
                                   link.phi = link.phi,  # link.phi is for the dispersion in betaregression
                                   type = type,control = betareg::betareg.control(fsmaxit = 10000))
                fail <- 0 #if it works
            },# end trycatch expression
            error = function(e){fail <<-fail+1; cat(fail); cat(" It'll be ok... Sometimes beta-type regressions fail because the smallest value is too close to zero.  Let's increase the smallest value and try again. \n ")},
            finally= {# completed
              # this was more useful when troubleshooting.  Otherwise it's a lot of output when in bootstrap
              # cat(paste("\n Completed topic", thetaindex))
            }
          )
          }# end while
          names(beta)[thetaindex] = topics[thetaindex]
        }
        
      }#end of return_just_coefs or the full model output (typically if not using bootstrap)
      }else{# Model == GAM
        if(return_just_coefs){
          pred_X_vals = unique(covariates)
          # inferring what is meant here, let's assume that it means predicting the GAM at 
          # the input covariate values.
          if(is.null(covariates|> dim())){
            nrow_X = length(unique(covariates))
            ncol_X = 1
          }else{
            nrow_X = nrow(unique(covariates))
            ncol_X = ncol(unique(covariates))
          }
          # make a place to put the predicted values.
          beta = matrix(NA, nrow = ncol(theta_nonzero), ncol = nrow_X)
          if(ncol_X==1){
            colnames(beta) = apply(pred_X_vals,1,function(x){paste0("X.",x)})
          }else{
            colnames(beta) = apply(pred_X_vals|> matrix(ncol = 1),
                                   1,
                                   function(x){paste0("X.",x, collapse=".")})
          }
          rownames(beta) = topics 
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
                                (fail-1)*min(1/ncol(theta[,thetaindex]) ,min(theta[theta[,thetaindex]>0,thetaindex]/1000)), 
                              denominator = denominator +
                                (fail-1+1)*min(1/ncol(theta[,thetaindex]) ,min(theta[theta[,thetaindex]>0,thetaindex]/1000))
                    ),covariates)
                }
                colnames(data)  = c("y",colnames(covariates))
                beta[thetaindex,] <- mgcv::gam(formula,
                                                family=mgcv::betar(link=link),
                                                data = data)|>
                  predict(newdata = pred_X_vals)
                fail <- 0 # if it works, then this will break out of the while loop.
              },# end trycatch expression,
              # if fails, then push fail back up to 1 from zero
              error = function(e){fail <<-fail+1; cat(fail); cat("\n It'll be ok... Sometimes beta-type regressions fail because the smallest value is too close to zero.  Let's increase the smallest value and try again. \n ")},
              finally= {# completed
                # this was more useful when troubleshooting.  Otherwise it's a lot of output when in bootstrap
                # cat(paste("\n Completed topic", thetaindex))
              }
            )
          }#end while
        }# end for loop
          if(ncol_X==1){
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
                                (fail-1)*min(1/ncol(theta[,thetaindex]) ,min(theta[theta[,thetaindex]>0,thetaindex]/1000)), 
                              denominator = denominator +
                                (fail-1+1)*min(1/ncol(theta[,thetaindex]) ,min(theta[theta[,thetaindex]>0,thetaindex]/1000))
                    ),covariates)
                }
                colnames(data)  = c("y",colnames(covariates))
                beta[[thetaindex]] <- mgcv::gam(formula,
                                                    family=mgcv::betar(link=link),
                                                     data = data)
                fail <- 0 # if it works, then this will break out of the while loop.
              },# end trycatch expression,
              # if fails, then push fail back up to 1 from zero
              error = function(e){fail <<-fail+1; cat(fail); cat(" It'll be ok... Sometimes beta-type regressions fail because the smallest value is too close to zero.  Let's increase the smallest value and try again. \n ")},
              finally= {# completed
                # this was more useful when troubleshooting.  Otherwise it's a lot of output when in bootstrap
                    # cat(paste("\n Completed topic", thetaindex))
              }
            )
          }# end while
            names(beta)[thetaindex] = topics[thetaindex]
        }# end forloop
        
        # end of GAM with returning the whole model.
    }
      }# end of GAM
    }# end of Model choice
    
    }else{ # with weights
    if(Model == "OLS"){
      if(return_just_coefs){
        beta = stats::coef(stats::lm.fit(x = as.matrix(covariates), y = theta_nonzero, weights = obs_weights))
        beta = t(beta)
        rownames(beta) = topics 
      }else{# return the lm model output
        beta = stats::lm.fit(x = covariates, y = theta_nonzero, weights = obs_weights)
        colnames(beta$coefficients) = topics 
      }
        
    }else{
      if(Model == "BETA"){
      # use a Beta regression model with weights
      zero_block = matrix(0, nrow(covariates) - nrow(theta), ncol(theta))
      theta_nonzero = rbind(theta_nonzero, zero_block)
      obs_weights = c(obs_weights, rep(1, nrow(zero_block)))
      if(return_just_coefs){# 
        beta = matrix(NA, nrow = ncol(theta_nonzero), ncol = 2*length(covariates))
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
                               link.phi = link.phi,  # link.phi is for the dispersion in betaregression
                               type = type,control = betareg::betareg.control(fsmaxit = 10000))|> coefficients()|> unlist()
            
            fail = -1 #it works
            
          },error = function(e){print(fail)},
            finally= {
              if(all(is.na(beta[thetaindex,]))){
                if(fail != -1){cat(paste(fail, "fail for index ", thetaindex, " epsilon = ", min(theta_nonzero[,thetaindex])))}
                fail <<- fail + 1; #if it worked now fail = 0,  if it didn't work then fail is growing 2+
                theta_nonzero[,thetaindex] <<- theta_nonzero[,thetaindex]+fail*min(theta_nonzero[,thetaindex])*.5; 
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
                     # reconstruct theta_nonzero but pull it inwards away from boundaries.
                    data =  data.frame(
                      normalize(theta[,thetaindex] +     (fail-1)*min(1/ncol(theta[,thetaindex]) ,min(theta[theta[,thetaindex]>0,thetaindex]/1000)), 
                                denominator = denominator +
                                          (fail-1+1)*min(1/ncol(theta[,thetaindex]) ,min(theta[theta[,thetaindex]>0,thetaindex]/1000))
                      ),covariates)
                   }
              colnames(data)  = c("y",colnames(covariates))
              cat(paste0("working on ", thetaindex))
              beta[thetaindex] = 
                betareg::betareg(formula, data = data,
                                 link = link,
                                 weights = obs_weights,
                                 link.phi = link.phi,  # link.phi is for the dispersion in betaregression
                                 type = type,control = betareg::betareg.control(fsmaxit = 10000))
              
              fail = -1 #it works
              
            },error = function(e){print(fail)},
            finally= {
              if(all(is.na(beta[[thetaindex]]))){
                if(fail != -1){cat(paste(fail, "fail for index ", thetaindex, " epsilon = ", min(theta_nonzero[,thetaindex])))}
                fail <<- fail + 1; #if it worked now fail = 0,  if it didn't work then fail is growing 2+
                theta_nonzero[,thetaindex] <<- theta_nonzero[,thetaindex]+fail*min(theta_nonzero[,thetaindex])*.5; 
              }
            }
            )  
          }
          names(beta)[thetaindex] = topics[thetaindex]
        }
        # return the whole regression model
      }
    }else{# Model == GAM
      if(return_just_coefs){
        pred_X_vals = unique(covariates)
        # inferring what is meant here, let's assume that it means predicting the GAM at 
        # the input covariate values.
        if(is.null(covariates|> dim())){
          nrow_X = length(unique(covariates))
          ncol_X = 1
        }else{
          nrow_X = nrow(unique(covariates))
          ncol_X = ncol(unique(covariates))
        }
        # make a place to put the predicted values.
        beta = matrix(NA, nrow = ncol(theta_nonzero), ncol = nrow_X)
        if(ncol_X==1){
          colnames(beta) = apply(pred_X_vals,1,function(x){paste0("X.",x)})
        }else{
          colnames(beta) = apply(pred_X_vals|> matrix(ncol = 1),
                                 1,
                                 function(x){paste0("X.",x, collapse=".")})
        }
        rownames(beta) = topics 
        for(thetaindex in 1:ncol(theta_nonzero)){
          fail = 1
          while(fail!=0 & fail < 10){
            tryCatch({  
              data = cbind(theta_nonzero[,thetaindex]+(fail-1)*min(theta_nonzero[,thetaindex])*.5,
                           covariates)|> as.data.frame()
              colnames(data)  = c("y",colnames(covariates))
              beta[thetaindex,] <- mgcv::gam(formula,
                                             family=mgcv::betar(link=link),
                                             data = data,
                                             weights)|>
                predict(newdata = pred_X_vals)
              fail <- 0 # if it works, then this will break out of the while loop.
            },# end trycatch expression,
            # if fails, then push fail back up to 1 from zero
            error = function(e){fail <<-fail+1; cat(fail); cat("\n It'll be ok... Sometimes beta-type regressions fail because the smallest value is too close to zero.  Let's increase the smallest value and try again. \n ")},
            finally= {# completed
              # this was more useful when troubleshooting.  Otherwise it's a lot of output when in bootstrap
              # cat(paste("\n Completed topic", thetaindex))
            }
            )
          }#end while
        }# end for loop
        if(ncol_X==1){
          beta = rbind(beta,X_pred_vals = c(as.matrix(pred_X_vals)))
        }
      }else{#return the whole model
        beta = list()
        for(thetaindex in 1:ncol(theta_nonzero)){
          fail = 1
          while(fail!=0 & fail < 10){
            tryCatch({
              data = cbind(theta_nonzero[,thetaindex]+(fail-1)*min(theta_nonzero[,thetaindex])*.5,
                           covariates)|> as.data.frame()
              colnames(data)  = c("y",colnames(covariates))
              beta[[thetaindex]] <- mgcv::gam(formula,
                                              family=mgcv::betar(link=link),
                                              data = data,
                                              weights = obs_weights)
              fail <- 0 # if it works, then this will break out of the while loop.
            },# end trycatch expression,
            # if fails, then push fail back up to 1 from zero
            error = function(e){fail <<-fail+1; cat(fail); cat(" It'll be ok... Sometimes beta-type regressions fail because the smallest value is too close to zero.  Let's increase the smallest value and try again. \n ")},
            finally= {# completed
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
  return(beta)
}
#' boot_reg
#'
#' Bootstrap regression coefficients in order to estimate their sampling distribution.
#'
#' @param output An object of class nmf_output
#'
#' @param samples The number of bootstrap samples to use. If set to 1 then return the full regression model output.  If >1, then collect the coefficients for bootstrap.
#'
#' @param obs_weights weights for observations to be passed to get_regression_coefs
#' 
#' @param Model choose c("BETA", "GAM", "OLS") for OLS, beta regression or a Generalized Additive model with family beta.  Default is "BETA".  OLS is only useful if the covariates are categorical, to be passed to get_regression_coefs
#' 
#' @param return_just_coefs returns the coefficients as opposed to the full betaregression output to be passed to get_regression_coefs
#'
#' @param topics is a vector of the subset of topics of interest.  Default is to perform regression on all topics. 
#' 
#' @param formula of the form Y = X |Z  for the model E(Y) = XB for the mean and var(Y)=Zß for the precision to be passed to get_regression_coefs
#' 
#' @param link is the link function for the betaregression GLM to be passed to get_regression_coefs or when using GAM
#' 
#' @param link.phi is the link function for the precision to be passed to get_regression_coefs, only for betaregression
#' 
#' @param type  is one of ML, BR, BC for maximum likelihood, bias reduces, or bias corrected estimates to be passed to get_regression_coefs
#'
#' 
#' @return Most of the original call values and a boot_reg list containing matrices/vectors, each of which contains regression coefficients produced by
#' get_regression_coefs(). Each list element corresponds to a bootstrap sample. Combining a
#' particular element across bootstrap iterates estimates the sampling distribution
#' of the associated estimator; see boot_plot() and create_error_bars().
#'
#' @examples
#' neurips_input = create_input(neurips_tdm, neurips_words,
#'    topics = 10, project = TRUE, proj_dim = 500, covariates = year_bins)
#' neurips_output = solve_nmf(neurips_input)
#' boot_samples = boot_reg(neurips_output, 1000)
#'
#' @export
boot_reg = function(output, samples, 
                    obs_weights = NULL, 
                    Model = c("BETA", "GAM", "OLS"), 
                    return_just_coefs = TRUE,
                    topics = NULL,
                    formulas = TRUE, 
                    formula = NULL,
                    link = "logit",
                    link.phi = "log",
                    type = "ML",
                    theta_transformation = NULL,
                    na.rm = TRUE){
  
  ##### check input types/whether covariates are specified
  if(class(output) != "nmf_output"){
    stop("Output must be of class nmf_output.")
  }
  if(is.null(output$covariates)){
    stop("No design matrix specified.")
  }
  samples = as.integer(samples)
  if(samples <= 0){
    stop("Samples must be a positive integer.")
  }
  
  ##### set up matrices for regression/list for return value
  if(is.null(topics)){
    topics = output$anchors
  }
  ##### set up matrices for regression/return value
  # select the rows corresponding to topics selected:
  if(na.rm){
    #remove the problematic 'no topic' documents now so that we don't end up with a bootstrap sample of all NAs
    theta              = output$theta[which(output$anchors%in% topics),
                                      which(!is.na(1/output$sum_theta_over_docs) & 
                                              !is.infinite(1/output$sum_theta_over_docs))]
    covariates         = output$covariates[which(!is.na(1/output$sum_theta_over_docs) & 
                                                   !is.infinite(1/output$sum_theta_over_docs)),
    ] |> as.data.frame()
    colnames(covariates) = colnames(output$covariates)
    sum_theta_over_docs = output$sum_theta_over_docs[which(!is.na(1/output$sum_theta_over_docs) & 
                                                             !is.infinite(1/output$sum_theta_over_docs))
    ]
  }else{
    # keep all docs.
    theta               = output$theta
    covariates          = output$covariates
    sum_theta_over_docs = output$sum_theta_over_docs
  }
  to_return = list()
  
  ##### use a while loop; this is because bootstrap sample may not include all factor levels
  ##### throw out these samples
  ##### is there a better way of dealing with this?
  
    for(i in 1:samples){
    
    ##### counter for samples where not all factor levels are present
    ##### counter for list elements (exclude bad samples in output)
    bad_samples = 0
    if(i == 1){
      j = 0
    }
    
    ##### produce bootstrap sample and form associated theta and covariate
    # this is written so that a sum-to-one constraint is automatically pushed through if it exists
    # note that sum to one constraints are incompatible with beta regression since it requires an observed
    constraint_block = tail(covariates, nrow(covariates) - ncol(theta))
    
    covariate_block = head(covariates, ncol(theta))
    boot_docs = sample(1:ncol(theta), replace = T)
    boot_theta = theta[,boot_docs]
    # put the constraint back if it exists:
    zero_block = matrix(0, nrow(theta),nrow(covariates) - ncol(theta))
    boot_theta = cbind(boot_theta, zero_block)
    
    boot_covariates = covariate_block[boot_docs,]|> as.data.frame()
    colnames(boot_covariates) = colnames(covariate_block)
    boot_covariates = rbind(boot_covariates, constraint_block)
    
    ##### check if all factor levels are used for appropriate variables
    ##### jump to next iteration if not
    if(0 %in% colSums(boot_covariates)){
      bad_samples = bad_samples + 1
      next
    }
    j = j + 1
    
    ##### create a new nmf_output object but with bootstrapped theta and covariate
    boot_output                     = output
    boot_output$theta               = boot_theta
    boot_output$covariates          = boot_covariates
    boot_output$anchors             = topics
    boot_output$sum_theta_over_docs = sum_theta_over_docs[boot_docs]
    ##### call get_regression_coefs and append list element
    boot_coefs = get_regression_coefs(output = boot_output, 
                                      obs_weights = obs_weights, 
                                      Model = Model, 
                                      return_just_coefs = return_just_coefs, 
                                      formula = formula,
                                      link = link,
                                      link.phi = link.phi, 
                                      type = type, 
                                      theta_transformation = theta_transformation,
                                      topics = topics,
                                      na.rm = na.rm)
    to_return[[j]] = boot_coefs
    
    ##### progress of iterations
    if(i %% 50 == 0){
      cat(i, " of ", samples, " bootstrap samples complete.\n")
    }
    
  }
 
  ##### print warning message if bad_samples > 0
  if(bad_samples > 0){
    cat("Warning: not all factor levels present in ", bad_samples, " bootstrap samples.
        These samples have been excluded from the output -- inferences may be affected as a result.\n")
  }
  
  ##### return

  return(list(boot_reg =to_return,
              obs_weights = obs_weights, 
              Model = Model, 
              return_just_coefs = return_just_coefs, 
              formula = formula,
              link     = ifelse(Model == "OLS", NA, link),
              link.phi = ifelse(Model == "OLS", NA, link.phi), 
              type = ifelse(Model == "BETA",type, NA),  
              theta_transformation = theta_transformation,
              topics = topics,
              na.rm = na.rm)
      ) 
}


#' boot_reg_stratified
#'
#' Bootstrap regression coefficients in order to estimate their sampling distribution.
#' Stratified bootstrap is used to maintain a constant number of individuals in each group when groups
#' are categorical
#'
#' @param output An object of class nmf_output
#'
#' @param samples The number of bootstrap samples to use.
#'
#' @param parallel number of cores to use.  Skip the parallel overhead if using 1 core.
#'
#' @param obs_weights weights for observations to be passed to get_regression_coefs
#' 
#' @param Model choose c("BETA", "GAM", "OLS") for OLS, beta regression or a Generalized Additive model with family beta.  Default is "BETA".  OLS is only useful if the covariates are categorical, to be passed to get_regression_coefs
#'
#' @param return_just_coefs returns the coefficients as opposed to the full betaregression output to be passed to get_regression_coefs
#' 
#' @param topics is a vector of the subset of topics of interest.  Default is to perform regression on all topics.
#' 
#' @param formula of the form Y = X |Z  for the model E(Y) = XB for the mean and var(Y)=Zß for the precision to be passed to get_regression_coefs
#' 
#' @param link is the link function for the betaregression GLM to be passed to get_regression_coefs
#' 
#' @param link.phi is the link function for the precision to be passed to get_regression_coefs
#' 
#' @param type  is one of ML, BR, BC for maximum likelihood, bias reduces, or bias corrected estimates to be passed to get_regression_coefs
#'
#' @return Most of the original call values and a boot_reg list containing matrices/vectors, each of which contains regression coefficients produced by
#' get_regression_coefs(). Each list element corresponds to a bootstrap sample. Combining a
#' particular element across bootstrap iterates estimates the sampling distribution
#' of the associated estimator; see boot_plot() and create_error_bars().
#'
#' @examples
#' neurips_input = create_input(neurips_tdm, neurips_words,
#'    topics = 10, project = TRUE, proj_dim = 500, covariates = year_bins)
#' neurips_output = solve_nmf(neurips_input)
#' boot_samples = boot_reg(neurips_output, 1000)
#'
#' @export
#' 
#' 
boot_reg_stratified = function(output, samples, parallel = 4,
                               obs_weights = NULL, 
                               Model = c("BETA", "GAM", "OLS"), 
                               return_just_coefs = TRUE, 
                               topics = NULL,
                               formula = NULL,
                               link = "logit",
                               link.phi = "log", 
                               type = "ML",
                               theta_transformation = NULL,
                               na.rm = TRUE){
  
  ##### check input types/whether covariates are specified
  if(class(output) != "nmf_output"){
    stop("Output must be of class nmf_output.")
  }
  if(is.null(output$covariates)){
    stop("No design matrix specified.")
  }
  samples = as.integer(samples)
  if(samples <= 0){
    stop("Samples must be a positive integer.")
  }
  parallel = as.integer(parallel)
  if(parallel <= 0){
    stop("parallel must be a positive integer.")
  }

  ##### set up matrices for regression/list for return value
  if(is.null(topics)){
    topics = output$anchors
  }
  
  ##### set up matrices for regression/return value
  # select the rows corresponding to topics selected:
  if(na.rm){
    #remove the problematic 'no topic' documents now so that we don't end up with a bootstrap sample of all NAs
    theta              = output$theta[which(output$anchors%in% topics),
                                       which(!is.na(1/output$sum_theta_over_docs) & 
                                             !is.infinite(1/output$sum_theta_over_docs))]
    covariates         = output$covariates[which(!is.na(1/output$sum_theta_over_docs) & 
                                                 !is.infinite(1/output$sum_theta_over_docs)),
    ] |> as.data.frame()
    colnames(covariates) = colnames(output$covariates)
    sum_theta_over_docs = output$sum_theta_over_docs[which(!is.na(1/output$sum_theta_over_docs) & 
                                                             !is.infinite(1/output$sum_theta_over_docs))
                                                     ]
      }else{
    # keep all docs.
    theta               = output$theta
    covariates          = output$covariates
    sum_theta_over_docs = output$sum_theta_over_docs
  }
  
  
  to_return = list()
  
  
  # this is written so that the constraint is automatically pushed through if it exists
    constraint_block = tail(covariates, nrow(covariates) - ncol(theta))
    zero_block = matrix(0, nrow(theta),nrow(covariates) - ncol(theta))
  covariate_block = head(covariates, ncol(theta))
  # identify the categories (remove the intercept)
  categorical_groups = unique(covariate_block)
  groups = apply(covariate_block,1,function(x){paste0(names(x[x==1]), collapse = "_")})
  
  group_levels = unique(groups)
  rownames(categorical_groups) = group_levels
  group_count = NULL
  for( group in 1:nrow(categorical_groups)){
    group_count[group] = sum(row.match(covariate_block,categorical_groups[group,],0))
  }
  names(group_count) = group_levels
  
  if(parallel ==1){
    for(i in 1:samples){
      
      ##### produce bootstrap sample and form associated theta and covariate
      
      boot_docs = rep(NA, ncol(theta))
      start = 0
      sampled_inds = NULL
      for(group_strat in group_levels){
        ####NOTE THIS NEEDS RETHINKING WHEN THERE ARE MULTIPLE COVARIATES
        sampled_indices = which(row.match(covariate_block,categorical_groups[group_strat,],0)==1)[
                                        sample(1:group_count[[group_strat]], replace = T)
                                        ]
        sampled_inds = c(sampled_inds,sampled_indices)
        boot_docs[start+(1:group_count[[group_strat]])] = sampled_indices
        start = start + group_count[[group_strat]]
      }
      
      boot_theta = theta[,boot_docs]
      # put the constraint back if it exists:
      boot_theta = cbind(boot_theta, zero_block)
      boot_covariates = covariate_block[boot_docs,]|> as.data.frame()
      colnames(boot_covariates) = colnames(covariate_block)
      boot_covariates = rbind(boot_covariates, constraint_block)
      
      ##### create a new nmf_output object but with bootstrapped theta and covariate
      boot_output                     = output
      boot_output$theta               = boot_theta
      boot_output$covariates          = boot_covariates
      boot_output$anchors             = topics
      boot_output$sum_theta_over_docs = sum_theta_over_docs[boot_docs]
      boot_coefs = get_regression_coefs(boot_output,
                                        obs_weights = obs_weights, 
                                        Model = Model, 
                                        return_just_coefs = return_just_coefs, 
                                        formula = formula,
                                        link = link,
                                        link.phi = link.phi, 
                                        type = type, 
                                        theta_transformation = theta_transformation,
                                        topics = topics,
                                        na.rm = na.rm)
      to_return[[i]] = boot_coefs
      
      ##### progress of iterations
      if(i %% 50 == 0){
        cat(i, " of ", samples, " bootstrap samples complete.\n")
      }
   }
  }else{
    #parallel setup chunk 1 start
    library(doParallel)
    #create the cluster
    my.cluster <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(cl = my.cluster)
    # defaults to fork on linux and psock on windows
    #parallel setup chunk 1 end
    
    #parallel loop start
    to_return<- foreach(i=1:samples, .export = c("get_regression_coefs"),
                        .packages = "betareg") %dopar% {
      
      ##### produce bootstrap sample and form associated theta and covariate
      
      boot_docs = rep(NA, ncol(theta))
      start = 0
      sampled_inds = NULL
      for(group_strat in group_levels){
        sampled_indices = which(groups == group_strat)[sample(1:group_count[[group_strat]], replace = T)]
        sampled_inds = c(sampled_inds,sampled_indices)
        boot_docs[start+(1:group_count[[group_strat]])] = sampled_indices
        start = start + group_count[[group_strat]]
      }
      
      boot_theta = theta[,boot_docs]
      # put the constraint back if it exists:
      boot_theta = cbind(boot_theta, zero_block)
      boot_covariates = covariate_block[boot_docs,]|> as.data.frame()
      colnames(boot_covariates) = colnames(covariate_block)
      boot_covariates = rbind(boot_covariates, constraint_block)
      
      ##### create a new nmf_output object but with bootstrapped theta and covariate
      boot_output                     = output
      boot_output$theta               = boot_theta
      boot_output$covariates          = boot_covariates
      boot_output$anchors             = topics
      boot_output$sum_theta_over_docs = output$sum_theta_over_docs[boot_docs]
      return(get_regression_coefs(boot_output,
                                  obs_weights = obs_weights, 
                                  Model = Model, 
                                  return_just_coefs = return_just_coefs, 
                                  formula = formula,
                                  link = link,
                                  link.phi = link.phi, 
                                  type = type, 
                                  theta_transformation = theta_transformation,
                                  topics = topics,
                                  na.rm = na.rm)
      )
    }
  }
  ##### return
  return(list(boot_reg =to_return,
                  obs_weights = obs_weights, 
                  Model = Model, 
                  return_just_coefs = return_just_coefs, 
                  formula = formula,
                  link     = ifelse(Model == "OLS", NA, link),
                  link.phi = ifelse(Model == "OLS", NA, link.phi), 
                  type = type, 
                  theta_transformation = theta_transformation,
                  topics = topics,
                  na.rm = na.rm,
                  type = "stratified bootstrap")
        )
  
   
}




#' boot_plot
#'
#' Plot smoothed histograms of regression effects for a given topic.
#'
#' @param boot_samples A list of bootstrap samples, as produced by boot_reg().
#'
#' @param topic The topic to use as the response variable, labeled by its anchor word.
#'
#' @param newdata provided if the bootsamples are regression coefficients and we wish to plot the fitted values at 'newdata' points.  newdata should be a data.frame where the N rows are N new observations with values across the columns.  The rownames are used for plotting as appropriate.
#'
#' @param Model Used if 'boot_samples' is not from boot_reg and therefore doesn't contain the model information already.  Used if newdata is supplied so that the function knows how to convert the bootstrap samples into predictions at 'newdata' points.  WILL BE REMOVED IN FUTURE
#' 
#' @param link Used if 'boot_samples' is not from boot_reg and therefore doesn't contain the model information already.  Used if newdata is supplied so that the function knows how to convert the bootstrap samples into predictions at 'newdata' points.  WILL BE REMOVED IN FUTURE
#'
#' @return A list as produced by ggplot2. Calling this function without assignment to a variable
#' will send a plot to the plotting window.
#'
#' @examples
#' neurips_input = create_input(neurips_tdm, neurips_words,
#'    topics = 10, project = TRUE, proj_dim = 500, covariates = year_bins)
#' neurips_output = solve_nmf(neurips_input)
#' boot_samples = boot_reg(neurips_output, 1000)
#' boot_plot(boot_samples, "model")
#'
#' @export
boot_plot = function(boot_samples, 
                     topic, 
                     newdata=NULL, 
                     Model = NULL, 
                     link = NULL,
                     theta_transformation = NULL){
  
  ##### check inputs
  if(is.list(boot_samples) & is.null(boot_samples$boot_reg)){
    # print("looks like output from get_regression_coefs")
    # print("might eventually make boot_reg a class to fix this ad hoc handling")
    # if boot_samples comes from get_regression_coefs 
    boot_reg_list = boot_samples
    bootstrapped = FALSE 
    if("gam" %in% class(boot_reg_list[[topic]])){
      Model = "GAM"
    }
  }else{
    # if boot_samples comes from boot_reg
    boot_reg_list = boot_samples$boot_reg
    Model         = boot_samples$Model
    link          = boot_samples$link
    theta_transformation         = boot_samples$theta_transformation
    bootstrapped = TRUE
  }
  
  
  if(!(is.list(boot_samples))){
    stop("boot_samples must be a list of vectors/matrices, as output by
              boot_reg() or get_regression_coefs().")
  }
  if(bootstrapped & !(topic %in% row.names(boot_reg_list[[1]]))){
    stop("Topic not among the anchor words.")
  }
  if(!bootstrapped & !(topic %in% c(names(boot_reg_list),rownames(boot_reg_list[[1]])))){
    stop("Topic not among the anchor words.")
  }
  
  if( is.null(newdata)){
    # if   = Using betaregression at the observation points
    # else = Just plot the effects not fits.
    if(any(class(boot_reg_list[[topic]])=="betareg")){
      
      # want to plot data fit at "newdata"
      # Beta regression model was used
      # plotting the prediction and prediction interval
      newdata = boot_reg_list[[topic]]$model
      newdata = subset(newdata, select=-c(y))|> unique()
      # in case you have repeated binary measures make the labels nicer:
      if(all(apply(newdata, 1, function(x){all(x%in%c(0,1))}))){
        names_2_use = apply(newdata,1,function(x){ind=which(x!=0);paste(names(x)[ind], collapse = "+")})
        rownames(newdata) = names_2_use
      }
      q.025 = predict(boot_reg_list[[topic]], 
                      newdata = newdata, 
                      type = "quantile", at = c(.025))
      q.5   = predict(boot_reg_list[[topic]], 
                      newdata = newdata, 
                      type = "quantile", at = c(.5))
      q.975 = predict(boot_reg_list[[topic]], 
                      newdata = newdata, 
                      type = "quantile", at = c(0.975))
      if(ncol(newdata)==1){
        fitted_values=data.frame(newdata = newdata,q.025,q.5,q.975)
      }else{
        fitted_values=data.frame(newdata = rownames(newdata),q.025,q.5,q.975)
      }
      ##### melt the data frame so it can be passed to ggplot
      fitted_values = reshape2::melt(fitted_values)
      names(fitted_values) = c("newdatapoint","quantile", "value")
      
      ##### create and evaluate plot
      topic_plot = ggplot2::ggplot(data = fitted_values, 
                                   ggplot2::aes(x=newdatapoint, y=value,colour = quantile, group = quantile))+
        geom_line(lwd = 2, alpha = .5)+
        ggplot2::labs(title = paste("Model fit with 95% interval, Topic =", stringr::str_to_title(topic)),
                      x = "Covariate value",
                      y = "P(topic | new data)")
      
      
    }else{
      # just plot the effects not fits.
      ##### create a data frame with covariate effects as columns and bootstrap samples as rows
      num_topics = dim(boot_reg_list[[1]])[1]
      num_covar = dim(boot_reg_list[[1]])[2]
      num_samples = length(boot_reg_list)
      boot_mat = matrix(rep(0, num_samples*num_covar), ncol = num_covar)
      for(i in 1:length(boot_reg_list)){
        boot_mat[i,] = boot_reg_list[[i]][topic,]
      }
      boot_effect = as.data.frame(boot_mat)
      names(boot_effect) = colnames(boot_reg_list[[1]])
      
      ##### melt the data frame so it can be passed to ggplot
      boot_effect = reshape2::melt(boot_effect)
      names(boot_effect) = c("covar", "weight")
      
      ##### create and evaluate plot
      topic_plot = ggplot2::ggplot(data = boot_effect, ggplot2::aes(x=weight, y=covar, height = after_stat(density)))+
        ggridges::stat_density_ridges(fill = "lightblue",rel_min_height = 0.005,scale = 1.25, stat = "density")+
        ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
        ggplot2::labs(title = paste("Bootstrap Distribution of Coefficients: Topic =", stringr::str_to_title(topic)),
                      x = "P(topic|X)",#expression(beta),
                      y = stringr::str_to_title(attributes(dimnames(boot_samples[[1]]))[[1]][2]))
      
      
    }
  }else{
    if( bootstrapped){
    
      if( Model == "BETA"){
      # want to plot data fit at "newdata"
      # Beta regression model was used
      # plotting the bootstrap distribution
      boot_reg_colnames = boot_reg_list[[1]]|> colnames() 
      mean_cols         = boot_reg_colnames |> grep(pattern = "mean", value = TRUE)
      precision_cols    = boot_reg_colnames |> grep(pattern = "mean", value = TRUE, invert = TRUE)
      
      xBeta         = lapply(boot_reg_list, function(x){x[,mean_cols]%*%t(newdata)})
      fitted_values = lapply(xBeta, function(xB){exp(xB)/(1-exp(xB))})
      
      num_covar = dim(fitted_values[[1]])[2]
      num_samples = length(fitted_values)
      boot_mat = matrix(rep(0, num_samples*num_covar), ncol = num_covar)
      for(i in 1:length(fitted_values)){
        boot_mat[i,] = fitted_values[[i]][topic,]
      }
      boot_effect = as.data.frame(boot_mat)
      names(boot_effect) = rownames(newdata)
      
      # if(ncol(newdata)==1){
      #     boot_effect=cbind(newdata = newdata,boot_mat)
      #   else{
      #     boot_effect=cbind(newdata = newdata,boot_mat)
      #   }
      
      ##### melt the data frame so it can be passed to ggplot
      boot_effect = reshape2::melt(boot_effect)
      names(boot_effect) = c("newdatapoint", "fitted")
      
      ##### create and evaluate plot
      topic_plot = ggplot2::ggplot(data = boot_effect, ggplot2::aes(x=newdatapoint, y=fitted))+
        geom_boxplot(fill = "lightblue")+
        ggplot2::labs(title = paste("Beta model fit for Topic =", stringr::str_to_title(topic)),
                      x = "Covariate value",
                      y = "P(topic | new data)")
      
    }
      if( Model == "OLS"){
        # want to plot data fit at "newdata"
        # OLS model was used
        # plotting the bootstrap distribution
        fitted_values = lapply(boot_reg_list, function(x){x%*%t(newdata)})
        if(!is.null(theta_transformation)){
          fitted_values = lapply(fitted_values,function(x){exp(x)-1})
        }
        num_covar = dim(fitted_values[[1]])[2]
        num_samples = length(fitted_values)
        boot_mat = matrix(rep(0, num_samples*num_covar), ncol = num_covar)
        for(i in 1:length(fitted_values)){
          boot_mat[i,] = fitted_values[[i]][topic,]
        }
        boot_effect = as.data.frame(boot_mat)
        names(boot_effect) = rownames(newdata)
        # if(ncol(newdata)==1){
        #     boot_effect=cbind(newdata = newdata,boot_mat)
        #   else{
        #     boot_effect=cbind(newdata = newdata,boot_mat)
        #   }
        ##### melt the data frame so it can be passed to ggplot
        boot_effect = reshape2::melt(boot_effect)
        names(boot_effect) = c("newdatapoint", "fitted")
        
        ##### create and evaluate plot
        topic_plot = ggplot2::ggplot(data = boot_effect, ggplot2::aes(x=newdatapoint, y=fitted))+
          geom_boxplot(fill = "lightblue")+
          ggplot2::labs(title = paste("Bootstrap distribution for Topic =", stringr::str_to_title(topic)),
                        x = "Covariate value",
                        y = "P(topic | new data)")
        
    }
    }else{
      #not bootstrapped
      if( any(class(boot_reg_list[[1]]) == "betareg")){
      # want to plot data fit at "newdata"
      # Beta regression model was used
      # plotting the prediction and prediction interval
      
      q.025 = predict(boot_reg_list[[topic]], 
                      newdata = data.frame(newdata), 
                      type = "quantile", at = c(.025))
      q.5   = predict(boot_reg_list[[topic]], 
                      newdata = data.frame(newdata), 
                      type = "quantile", at = c(0.5))
      q.975 = predict(boot_reg_list[[topic]], 
                      newdata = data.frame(newdata), 
                      type = "quantile", at = c(0.975))
      if(ncol(newdata)==1){
        fitted_values=data.frame(newdata = newdata,q.025,q.5,q.975)
      }else{
        fitted_values=data.frame(newdata = rownames(newdata),q.025,q.5,q.975)
      }
      ##### melt the data frame so it can be passed to ggplot
      fitted_values = reshape2::melt(fitted_values)
      names(fitted_values) = c("newdatapoint","quantile", "value")
      
      ##### create and evaluate plot
      topic_plot = ggplot2::ggplot(data = fitted_values, 
                                   ggplot2::aes(x=newdatapoint, y=value,colour = quantile, group = quantile))+
        geom_line(lwd = 2, alpha = .5)+
        ggplot2::labs(title = paste("Model fit with 95% interval, Topic =", stringr::str_to_title(topic)),
                      x = "Covariate value",
                      y = "P(topic | new data)")
      
    }else{
      if( Model == "GAM"){
      # want to plot data fit at "newdata"
      # Beta regression model was used
      # plotting the prediction and prediction interval
      
      fitted_values = predict(boot_reg_list[[topic]], 
                              newdata = newdata, 
                              type = "response",se.fit=TRUE)
      if(ncol(newdata)==1){
        fitted_values=data.frame(newdata = newdata,
                                 fitted = fitted_values$fit,
                                 plus_1sefit  = fitted_values$fit + fitted_values$se.fit,
                                 minus_1sefit = fitted_values$fit - fitted_values$se.fit)
      }else{
        fitted_values=data.frame(newdata = rownames(newdata),
                                 fitted = fitted_values$fit,
                                 plus_1sefit  = fitted_values$fit + fitted_values$se.fit,
                                 minus_1sefit = fitted_values$fit - fitted_values$se.fit)
      }
      colnames(fitted_values)[1] = "newdata"
      ##### melt the data frame so it can be passed to ggplot
      fitted_values = reshape2::melt(fitted_values,"newdata")
      names(fitted_values) = c("newdata","variable", "value")
      
      ##### create and evaluate plot
      topic_plot = ggplot2::ggplot(data = fitted_values, 
                                   ggplot2::aes(x=newdata, y=value,colour = variable))+
        geom_line(lwd = 2, alpha = .5)+
        ggplot2::labs(title = paste("GAM Model fit +/- 1 SE, Topic =", stringr::str_to_title(topic)),
                      x = "Covariate value",
                      y = "P(topic | new data)")
      
      }
    }
    }
    
  }
  

eval(topic_plot)
}
#' create_error_bars
#'
#' Output a data frame with two-sided, 95% confidence intervals based on bootstrap estimates of
#' sampling distributions.
#'
#' @param boot_samples A list of bootstrap samples, as produced by boot_reg().
#'
#' @param topic the selected topic(s) of interest for inference, default is all topics, though this is often too much output to be useful.
#' 
#'  @param coverage the level of coverage for the interval, default is 95%
#'
#' @return A data frame. Each row corresponds to a covariate and the columns give the CI.
#'
#' @examples
#' neurips_input = create_input(neurips_tdm, neurips_words,
#'    topics = 10, project = TRUE, proj_dim = 500, covariates = year_bins)
#' neurips_output = solve_nmf(neurips_input)
#' boot_samples = boot_reg(neurips_output, 1000)
#' create_error_bars(boot_samples)
#'
#' @export
create_error_bars = function(boot_samples, topic = NULL, coverage = .95){
  
  ##### check inputs
  if(!(is.list(boot_samples))){
    stop("boot_samples must be a list of matrices, as outputted by
              boot_reg() or get_regression_coefs.")
  } 
  if(is.list(boot_samples) & is.null(boot_samples$boot_reg)){
    # print("looks like output from get_regression_coefs")
    # print("might eventually make boot_reg a class to fix this ad hoc handling")
    # if boot_samples comes from get_regression_coefs 
    boot_reg_list = boot_samples
  }else{
    # if boot_samples comes from boot_reg
    boot_reg_list = boot_samples$boot_reg
  }

  if(!(is.matrix(boot_reg_list[[1]]))){
    stop("boot_samples must be a list of matrices, as outputted by
              boot_reg(). This function does not yet support large design matrices.")
  }
  if(!(0 < coverage & coverage < 1)){
    stop("coverage must be a numerical value in the interval (0,1).")
  }
  
  ##### melt and name a data frame of samples
  sample_frame = reshape2::melt(boot_reg_list)
  if(!is.null(topic)){
    sample_frame = dplyr::filter(sample_frame, Var1 %in% topic)
  }
  names(sample_frame) = c("topic", "coef", "estimate", "sample")
  
  ##### construct data frame w/ all possible covariate/topic combinations
  ##### fill with corresponding quantiles
  error_frame = expand.grid(unique(sample_frame$topic), unique(sample_frame$coef))
  names(error_frame) = c("topic", "coef")
  for(i in 1:nrow(error_frame)){
    subset_frame = sample_frame[sample_frame$topic == error_frame$topic[i] &
                                  sample_frame$coef == error_frame$coef[i],]
    error_frame$lower[i]  = stats::quantile(subset_frame$estimate, (1-coverage)/2)
    error_frame$median[i] = stats::quantile(subset_frame$estimate, .5)
    error_frame$upper[i]  = stats::quantile(subset_frame$estimate, coverage + (1 - coverage)/2)
  }
  colnames(error_frame)[c(3,5)] = paste0(c("lower","upper"),coverage)
  
  
  ##### return
  return(error_frame)
  
}
