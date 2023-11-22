#' get_regression_coefs
#'
#' Compute OLS coefficients, fitting a linear model between a user's specified covariates and topic
#' weight. Estimates are produced simultaneously for all topics.
#'
#' @param output An object of class nmf_output.
#'
#' @param obs_weights Weights for the documents, as used by weighted least squares. Defaults to null.
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
get_regression_coefs = function(output, obs_weights = NULL, OLS = TRUE){
  
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
  
  ##### set up matrices for OLS/return value
  theta = t(output$theta)
  covariates = output$covariates
  
  ##### simple function to normalize rows of theta
  normalize = function(x){
    return(x/sum(x))
  }
  
  ##### set up a data frame for regression
  ##### fit a linear model on all topics using specified covariates
  theta = apply(theta, FUN = normalize, MARGIN = 2)
  if(is.null(obs_weights)){
    zero_block = matrix(0, nrow(covariates) - nrow(theta), ncol(theta))
    theta = rbind(theta, zero_block)
    if(min(theta)==0){
      warning("beta regresssion won't work with obvserved {0,1} since this will result in infinite betas, setting smallest values to 1e-16")
      theta_nonzero = theta
      theta_nonzero[theta_nonzero<10^-16] = 10^-16
    }else{
      theta_nonzero = theta
    }
    if(OLS == TRUE){
      beta = stats::coef(stats::lm.fit(x = covariates, y = theta))
      
    }else{# use a Beta regression model
      ###IN BOOT_REG AND BOOT_REG_STRATIFIED USE THE INPUT constraint = FALSE
      # index = 4
      beta = betareg::betareg(theta_nonzero[,index]~covariates|covariates, # after "|" covariates included are for the dispersion.
                              link = "logit", 
                              link.phi = "log",  # link.phi is for the dispersion,
                              type = "BC") |> #type is for ML, BiasCorrected, BiasReduced
        coefficients()
      ### OBSERVED MUST BE IN (0,1), not [0,1] <-- set them to 10^-6 or maybe max(10^-10, {1st quartile ÷ # actual zero values}
    }
  }else{
    zero_block = matrix(0, nrow(covariates) - nrow(theta), ncol(theta))
    theta = rbind(theta, zero_block)
    obs_weights = c(obs_weights, rep(1, nrow(zero_block)))
    beta = stats::coef(stats::lm.wfit(x = covariates, y = theta, w = obs_weights))
  }
  
  ##### return
  beta = t(beta)
  rownames(beta) = output$anchors
  return(beta)
  
}
#' boot_reg
#'
#' Bootstrap OLS coefficients in order to estimate their sampling distribution.
#'
#' @param output An object of class nmf_output
#'
#' @param samples The number of bootstrap samples to use.
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
boot_reg = function(output, samples, constraint = TRUE, ...){
  
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
  
  ##### set up matrices for OLS/list for return value
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
    if(constraint){
      constraint_block = tail(covariates, nrow(covariates) - ncol(theta))
    }else{
      constraint_block = NULL
    }
    covariate_block = head(covariates, ncol(theta))
    boot_docs = sample(1:ncol(theta), replace = T)
    boot_theta = theta[,boot_docs]
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
#' Bootstrap OLS coefficients in order to estimate their sampling distribution.
#' Stratified bootstrap is used to maintain a constant number of individuals in each group when groups
#' are categorical
#'
#' @param output An object of class nmf_output
#'
#' @param samples The number of bootstrap samples to use.
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
boot_reg_stratified = function(output, samples, constraint = FALSE...){
  
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
  ##### simple function to normalize rows of theta
  normalize = function(x){
    return(x/sum(x))
  }
  
  ##### set up matrices for OLS/list for return value
  theta = output$theta
  covariates = output$covariates
  to_return = list()
  
  if(constraint){
    constraint_block = tail(covariates, nrow(covariates) - ncol(theta))
  }else{
    constraint_block = NULL
  }
  covariate_block = head(covariates, ncol(theta))
  # identify the categories (remove the intercept)
  categorical_groups = covariate_block[,
                                       apply(output$covariates,2,sum) != dim(covariate_block)[1]
  ]
  groups = apply(categorical_groups,1,function(x){names(x[x==1])})
  group_levels = unique(groups)
  group_count = table(groups)
  
  
  
  for(i in 1:samples){
    
    ##### produce bootstrap sample and form associated theta and covariate
    
    boot_docs = rep(NA, ncol(theta))
    start = 0
    sampled_inds = NULL
    for(group_strat in group_levels){
      sampled_indices = which(groups == group_strat)[sample(1:group_count[[group_strat]], replace = T)]
      sampled_inds = c(sampled_inds,sampled_indices)
      boot_docs[start+(1:group_count[[group_strat]])] = sampled_indices
      # if(any(is.na(boot_docs[start+(1:group_count[[group_strat]])]))){
      #   print("error")
      #   print(sampled_indices)
      #   print("start")
      #   print(start)
      #   print("inds")
      #
      # }
      start = start + group_count[[group_strat]]
    }
    
    boot_theta = theta[,boot_docs]
    boot_covariates = covariate_block[boot_docs,]
    boot_covariates = rbind(boot_covariates, constraint_block)
    
    ##### set up a data frame for regression
    ##### fit a linear model on all topics using specified covariates
    boot_theta = apply(boot_theta, FUN = normalize, MARGIN = 2)
    # if(0 %in% colSums(boot_covariates) | is.infinite(sum(boot_theta)) | is.na(sum(boot_theta))){
    #   bad_samples = boot_docs
    #   print(boot_docs)
    #   print(i)
    #   break
    #
    # }
    
    ##### create a new nmf_output object but with bootstrapped theta and covariate
    boot_output = output
    boot_output$theta = boot_theta
    boot_output$covariates = boot_covariates
    
    # note that theta should be normalized and sum to one.
    # however numerical precision isn't perfect and can lead to problems in the lm.fit
    # error handling for when numerical precision breaks things:
    if(all(boot_theta |> colSums()  == 1)){
      ##### call get_regression_coefs and append list element
      boot_coefs = get_regression_coefs(boot_output, ...)
      to_return[[i]] = boot_coefs
    }else{
      #potential problem cols where numerical precision was not perfect:
      nugget = 10* max(abs(1-(colSums(boot_theta) )))
      col_indices = which(colSums(boot_theta)!=1)
      boot_output$theta[,col_indices] = boot_theta[,col_indices] + nugget
      ##### call get_regression_coefs and append list element
      boot_coefs = get_regression_coefs(boot_output, ...)
      to_return[[i]] = boot_coefs
    }
    
    
    
    ##### progress of iterations
    if(i %% 10 == 0){
      cat(i, " of ", samples, " bootstrap samples complete.\n")
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
