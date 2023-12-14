#' get_regression_coefs
#'
#' Compute OLS coefficients, fitting a linear model between a user's specified covariates and topic
#' weight. Estimates are produced simultaneously for all topics.
#'
#' @param output An object of class nmf_output.
#'
#' @param obs_weights Weights for the documents, as used by weighted least squares. Defaults to null.
#'
#' @param OLS Logical for using OLS.  Alternative is to use beta regression
#'
#' @param return_just_coefs is a logical, if TRUE then just return the coefficients, if FALSE, then return the output from the lm or betareg function.
#' 
#' @param formula is the formula to be passed into betareg.  Of the form Y~ model+for+mean | model+for+dispersion
#' 
#' @param link is the link function for the GLM for the mean when using betaregression
#' 
#' @param link.phi is the link function for the GLM for the precision when using betaregression
#' 
#' @param type is one of ML, BC, BR for betaregression to use Maximum Likelihood, Bias Corrected, or Bias Reduced respectively
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
get_regression_coefs = function(output, obs_weights = NULL, OLS = FALSE, return_just_coefs = TRUE, formula = NULL, link = "logit",
                                link.phi = "log", type = "ML"){
  
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
  
  ##### set up matrices for regression/return value
  theta = t(output$theta)
  covariates = output$covariates
  
  # handle spaces and odd characters in column names
  colnames(covariates) =  make.names(colnames(covariates))
  
  
  
  # Deal with formulas for betaregression, put all covariates into the mean and precision
  if(is.null(formula) & OLS != TRUE){
    factors = colnames(covariates)
    formula = as.formula(paste("y~", paste(
      paste(factors, collapse="+"), "|",
      paste(factors, collapse="+")
    )
    )
    )
  }
  
  
  ##### simple function to normalize rows of theta
  normalize = function(x){
    return(x/sum(x))
  }
  
  ##### set up a data frame for regression
  ##### fit a linear model on all topics using specified covariates
  theta_nonzero = apply(theta, FUN = normalize, MARGIN = 2)
  
  if(min(theta)==0 & OLS != TRUE){
    epsilon = max(theta_nonzero)*1e-10
    theta_nonzero = theta_nonzero + epsilon
    warning(paste("beta regresssion won't work with observed {0,1} since this will result in infinite betas; increasing all values by", epsilon))
  }
  if(is.null(obs_weights)){

    if(OLS == TRUE){
      if(return_just_coefs){
        beta = stats::coef(stats::lm.fit(x = covariates, y = theta))
        beta = t(beta)
        rownames(beta) = output$anchors 
      }else{# return the model output
        beta = stats::lm.fit(x = covariates, y = theta)
        colnames(beta$coefficients) = output$anchors 
      }#end of using bootstrap T/F
      
    }else{# use a Beta regression model

      if(return_just_coefs){
          beta = matrix(NA, nrow = ncol(theta_nonzero), ncol = 2*ncol(covariates)+1)
          colnames(beta) = c(paste0("mean.", colnames(covariates)),
                           paste0("precision.", colnames(covariates)), "epsilon")
          rownames(beta) = output$anchors 
        for(thetaindex in 1:ncol(theta_nonzero)){
          fail = 1
          while(fail==1){

          tryCatch({
            data = cbind(theta[,thetaindex]+(fail-1)*min(theta_nonzero[,thetaindex])*.5,
                         covariates)
            colnames(data)  = c("y",colnames(covariates))
                 beta[thetaindex,] <-
                  betareg::betareg(formula, data = data,
                               link = link,
                               link.phi = link.phi,  # link.phi is for the dispersion,
                               type = type,control = betareg::betareg.control(fsmaxit = 10000))|> coefficients()|> unlist()

                 fail <- 0 #if it works
          },error = function(e){fail <<-fail+1; cat(fail); cat(" It'll be ok... Sometimes beta regression fails because the smallest value is too close to zero.  Let's increase the smallest value and try again.")},
        finally= {
          if(all(is.na(beta[thetaindex,]))){
            if(fail != 1){cat(paste(fail, "fail for index ", thetaindex, " epsilon = ", min(theta_nonzero[,thetaindex])))}
            fail <<- fail + 1; #if it worked now fail = 0,  if it didn't work then fail is growing 0->1->2->...
            cat(fail)
          }
        }
         )
          }
        }
          

      }else{# return the betareg model output and not just the coefficients
        beta = list()
        for(thetaindex in 1:ncol(theta_nonzero)){
          fail = 1
          while(fail==1){
            tryCatch({
              data = cbind(theta[,thetaindex]+(fail-1)*min(theta_nonzero[,thetaindex])*.5,
                           covariates)
              colnames(data)  = c("y",colnames(covariates))
              beta[thetaindex,] <-
                betareg::betareg(formula, data = data,
                                 link = link,
                                 link.phi = link.phi,  # link.phi is for the dispersion,
                                 type = type,control = betareg::betareg.control(fsmaxit = 10000))|> coefficients()|> unlist()

              fail <- 0 #if it works
            },error = function(e){fail <<-fail+1; cat(fail); cat(" It'll be ok... Sometimes beta regression fails because the smallest value is too close to zero.  Let's increase the smallest value and try again.")},
            finally= {
              if(all(is.na(beta[thetaindex,]))){
                if(fail != 1){cat(paste(fail, "fail for index ", thetaindex, " epsilon = ", min(theta_nonzero[,thetaindex])))}
                fail <<- fail + 1; #if it worked now fail = 0,  if it didn't work then fail is growing 0->1->2->...
                cat(fail)
              }
            }
            )
          }
          
        }
        names(beta) = output$anchors
         
      }#end of return_just_coefs or the full model output (typically if not using bootstrap)
    }# end of OLS = TRUE / FALSE
    }else{ # with weights
    if(OLS == TRUE){
      if(return_just_coefs){
        beta = stats::coef(stats::lm.fit(x = covariates, y = theta, weights = obs_weights))
        beta = t(beta)
        rownames(beta) = output$anchors 
      }else{# return the lm model output
        beta = stats::lm.fit(x = covariates, y = theta, weights = obs_weights)
        colnames(beta$coefficients) = output$anchors 
      }
        
    }else{# use a Beta regression model with weights
      zero_block = matrix(0, nrow(covariates) - nrow(theta), ncol(theta))
      theta = rbind(theta, zero_block)
      obs_weights = c(obs_weights, rep(1, nrow(zero_block)))
      if(return_just_coefs){# 
        beta = matrix(NA, nrow = ncol(theta_nonzero), ncol = 2*length(covariates))
        colnames(beta) = c(paste0("mean.", colnames(covariates)),
                           paste0("precision.", colnames(covariates)))
        rownames(beta) = output$anchors
        for(thetaindex in 1:ncol(theta_nonzero)){
          fail = 0
          while(fail==0){
            fail = 1
            tryCatch({
            cat(paste0("working on ", thetaindex))
            beta[thetaindex,] = 
              c(betareg::betareg.fit(y=theta_nonzero[,thetaindex],
                                     x=covariates,z=covariates,
                                     weights = obs_weights,
                                     link = link, 
                                     link.phi = link.phi,  # link.phi is for the dispersion,
                                     type = type,control = betareg::betareg.control(fsmaxit = 10000))|> coefficients()|> unlist(), 
                min(theta_nonzero[,thetaindex]))
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
          while(fail==0){
            fail = 1
            tryCatch({
              cat(paste0("working on ", thetaindex))
              beta[[thetaindex]] = 
                betareg::betareg(y=theta_nonzero[,thetaindex],
                                       x=covariates,z=covariates,
                                       weights = obs_weights,
                                       link = link, 
                                       link.phi = link.phi,  # link.phi is for the dispersion,
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
          
        }
        names(beta) = output$anchors 
        
        # return the whole regression model
      }
    }
    }
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
#' @param ... additional inputs to be passed to get_regression_coefs, for now only available option is OLS = TRUE/FALSE or (in future) obs_weights for regression though that is not yet active
#'
#' @return A list containing matrices/vectors, each of which contains regression coefficients produced by
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
boot_reg = function(output, samples, ...){
  
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
  theta = output$theta
  covariates = output$covariates
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
    boot_theta = cbind(theta, zero_block)
    
    boot_covariates = covariate_block[boot_docs,]
    boot_covariates = rbind(boot_covariates, constraint_block)
    
    ##### check if all factor levels are used for appropriate variables
    ##### jump to next iteration if not
    if(0 %in% colSums(boot_covariates)){
      bad_samples = bad_samples + 1
      next
    }
    j = j + 1
    
    ##### create a new nmf_output object but with bootstrapped theta and covariate
    boot_output = output
    boot_output$theta = boot_theta
    boot_output$covariates = boot_covariates
    
    ##### call get_regression_coefs and append list element
    boot_coefs = get_regression_coefs(boot_output, ...)
    to_return[[j]] = boot_coefs
    
    ##### progress of iterations
    if(i %% 10 == 0){
      cat(i, " of ", samples, " bootstrap samples complete.\n")
    }
    
  }
 
  ##### print warning message if bad_samples > 0
  if(bad_samples > 0){
    cat("Warning: not all factor levels present in ", bad_samples, " bootstrap samples.
        These samples have been excluded from the output -- inferences may be affected as a result.\n")
  }
  
  ##### return
  return(to_return)
  
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
#' @param ... additional inputs to be passed to get_regression_coefs, for now only available option is OLS = TRUE/FALSE
#'
#' @return A list containing matrices/vectors, each of which contains regression coefficients produced by
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
boot_reg_stratified = function(output, samples, parallel = 4,...){
  # eventually deprecate the "..." options since the "..." for now is just to allow the OLS option
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
  ##### simple function to normalize rows of theta
  normalize = function(x){
    return(x/sum(x))
  }
  
  ##### set up matrices for regression/list for return value
  theta = output$theta
  covariates = output$covariates
  to_return = list()
  
  # this is written so that the constraint is automatically pushed through if it exists
    constraint_block = tail(covariates, nrow(covariates) - ncol(theta))
  
  covariate_block = head(covariates, ncol(theta))
  # identify the categories (remove the intercept)
  categorical_groups = covariate_block[, apply(covariate_block,2,function(x){length(unique(x))>1})]
  groups = apply(categorical_groups,1,function(x){names(x[x==1])})
  group_levels = unique(groups)
  group_count = table(groups)
  
  
  if(parallel ==1){
    for(i in 1:samples){
      
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
      boot_covariates = covariate_block[boot_docs,]
      boot_covariates = rbind(boot_covariates, constraint_block)
      
      ##### set up a data frame for regression
      ##### fit a linear model on all topics using specified covariates
      boot_theta = apply(boot_theta, FUN = normalize, MARGIN = 2)
      
      ##### create a new nmf_output object but with bootstrapped theta and covariate
      boot_output = output
      boot_output$theta = boot_theta
      boot_output$covariates = boot_covariates

      boot_coefs = get_regression_coefs(boot_output,...)
      to_return[[i]] = boot_coefs
      
      ##### progress of iterations
      if(i %% 10 == 0){
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
      boot_covariates = covariate_block[boot_docs,]
      boot_covariates = rbind(boot_covariates, constraint_block)
      
      ##### set up a data frame for regression
      ##### fit a linear model on all topics using specified covariates
      boot_theta = apply(boot_theta, FUN = normalize, MARGIN = 2)
      
      ##### create a new nmf_output object but with bootstrapped theta and covariate
      boot_output = output
      boot_output$theta = boot_theta
      boot_output$covariates = boot_covariates
      
      return(get_regression_coefs(boot_output,...))
    }
  }
  ##### return
  return(to_return)
}




#' boot_plot
#'
#' Plot smoothed histograms of regression effects for a given topic.
#'
#' @param boot_samples A list of bootstrap samples, as produced by boot_reg().
#'
#' @param topic The topic to use as the response variable, labeled by its anchor word.
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
boot_plot = function(boot_samples, topic){
  
  ##### check inputs
  if(!(is.list(boot_samples))){
    stop("boot_samples must be a list of vectors/matrices, as output by
              boot_reg().")
  }
  if(!(is.matrix(boot_samples[[1]]))){
    stop("boot_samples must be a list of matrices, as output by
              boot_reg().")
  }
  if(!(topic %in% row.names(boot_samples[[1]]))){
    stop("Topic not among the anchor words.")
  }
  
  ##### create a data frame with covariate effects as columns and bootstrap samples as rows
  num_topics = dim(boot_samples[[1]])[1]
  num_covar = dim(boot_samples[[1]])[2]
  num_samples = length(boot_samples)
  boot_mat = matrix(rep(0, num_samples*num_covar), ncol = num_covar)
  for(i in 1:length(boot_samples)){
    boot_mat[i,] = boot_samples[[i]][topic,]
  }
  boot_effect = as.data.frame(boot_mat)
  names(boot_effect) = colnames(boot_samples[[1]])
  
  ##### melt the data frame so it can be passed to ggplot
  boot_effect = reshape2::melt(boot_effect)
  names(boot_effect) = c("covar", "weight")
  
  ##### create and evaluate plot
  topic_plot = ggplot2::ggplot(data = boot_effect, ggplot2::aes(x=weight, y=covar))+
    ggridges::stat_density_ridges(fill = "lightblue")+
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed")+
    ggplot2::labs(title = paste("Topic =", stringr::str_to_title(topic)),
                  x = expression(beta),
                  y = stringr::str_to_title(attributes(dimnames(boot_samples[[1]]))[[1]][2]))
  eval(topic_plot)
  
}
#' create_error_bars
#'
#' Output a data frame with two-sided, 95% confidence intervals based on bootstrap estimates of
#' sampling distributions.
#'
#' @param boot_samples A list of bootstrap samples, as produced by boot_reg().
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
create_error_bars = function(boot_samples, coverage = .95){
  
  ##### check inputs
  if(!(is.list(boot_samples))){
    stop("boot_samples must be a list of matrices, as outputted by
              boot_reg().")
  }
  if(!(is.matrix(boot_samples[[1]]))){
    stop("boot_samples must be a list of matrices, as outputted by
              boot_reg(). This function does not yet support large design matrices.")
  }
  if(!(0 < coverage & coverage < 1)){
    stop("coverage must be a numerical value in the interval (0,1).")
  }
  
  ##### melt and name a data frame of samples
  sample_frame = reshape2::melt(boot_samples)
  names(sample_frame) = c("topic", "coef", "estimate", "sample")
  
  ##### construct data frame w/ all possible covariate/topic combinations
  ##### fill with corresponding quantiles
  error_frame = expand.grid(unique(sample_frame$topic), unique(sample_frame$coef))
  names(error_frame) = c("topic", "coef")
  for(i in 1:nrow(error_frame)){
    subset_frame = sample_frame[sample_frame$topic == error_frame$topic[i] &
                                  sample_frame$coef == error_frame$coef[i],]
    error_frame$lower[i] = stats::quantile(subset_frame$estimate, (1-coverage)/2)
    error_frame$upper[i] = stats::quantile(subset_frame$estimate, coverage + (1 - coverage)/2)
  }
  
  ##### return
  return(error_frame)
  
}
