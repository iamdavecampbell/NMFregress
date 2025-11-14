
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
#' @return A data frame. Each row corresponds to a covariate and the columns give the CI.
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
#' brett_object = boot_reg(my_output, samples = 100, Model = "OLS")
#' create_error_bars(brett_object, topicvector = "death")
#'
#' brett_object = get_regression_coefs(my_output, Model = "BETA", return_just_coefs = TRUE)
#' create_error_bars(brett_object, topicvector = "death")
#'
#' @export
create_error_bars = function(brett_object,
                             topicvector = NULL,
                             coverage = .95){
  ##### check inputs
  if(!inherits(brett_object, "BRETT_ouput")){
    stop("brett_object must be a BRETT_ouput object, as outputted by
              boot_reg() or get_regression_coefs.")
  }
  if(!inherits(brett_object, "list")){
    stop("to build the error bars we need the bootstrap model, rerun get_regression_coefs with the input flag return_just_coefs = FALSE")
  }
  if(is.null(brett_object$boot_reg)){
    # if brett_object comes from get_regression_coefs
    boot_reg_list <- brett_object
  }else{
    # if brett_object is a list of bootstrap runs
    boot_reg_list <- brett_object$boot_reg
  }
  if(!(is.matrix(boot_reg_list[[1]]))){
    stop("brett_object must be a list of matrices, as outputted by
              boot_reg()")
  }
  if(!(0 < coverage & coverage < 1)){
    stop("coverage must be a numerical value in the interval (0,1).")
  }

  ##### melt and name a data frame of samples
  sample_frame <- reshape2::melt(boot_reg_list)
  colnames(sample_frame) <- c("topic", "coef", "estimate", "sample")
  if(!is.null(topicvector)){
    sample_frame <- dplyr::filter(sample_frame, topic %in% topicvector)
  }


  ##### construct data frame w/ all possible covariate/topic combinations
  ##### fill with corresponding quantiles
  error_frame <- expand.grid(unique(sample_frame$topic), unique(sample_frame$coef))
  names(error_frame) <- c("topic", "coef")
  for(i in 1:nrow(error_frame)){
    subset_frame <- sample_frame[sample_frame$topic == error_frame$topic[i] &
                                   sample_frame$coef == error_frame$coef[i],]
    error_frame$lower[i]  <- stats::quantile(subset_frame$estimate, (1-coverage)/2)
    error_frame$median[i] <- stats::quantile(subset_frame$estimate, .5)
    error_frame$upper[i]  <- stats::quantile(subset_frame$estimate, coverage + (1 - coverage)/2)
  }
  colnames(error_frame)[c(3,5)] <- paste0(c("lower","upper"),coverage)
  ##### return
  return(error_frame)
}































#
#
#
#
#
#   old quantile plot for beta regression fit
#
#
# q_025 <- betareg::predict(
#   betamodel,
#   newdata = newdata,
#   type = "quantile",
#   at = c(.025)
# )
# q_5   <- betareg::predict(
#   boot_reg_list[[topic]],
#   newdata = newdata,
#   type = "quantile",
#   at = c(.5)
# )
# q_975 <- betareg::predict(
#   boot_reg_list[[topic]],
#   newdata = newdata,
#   type = "quantile",
#   at = c(0.975)
# )
# if (ncol(newdata) == 1) {
#   fitted_values <- data.frame(newdata = newdata,
#                               q_025,
#                               q_5,
#                               q_975)
# } else {
#   fitted_values <- data.frame(newdata = rownames(newdata),
#                               q_025,
#                               q_5,
#                               q_975)
# }
# ##### melt the data frame so it can be passed to ggplot
# fitted_values <-  reshape2::melt(fitted_values)
# names(fitted_values) <-  c("newdatapoint", "quantile", "value")
#
# ##### create and evaluate plot
# topic_plot <-  ggplot2::ggplot(
#   data = fitted_values,
#   ggplot2::aes(
#     x = .data$newdatapoint,
#     y = .data$value,
#     colour = .data$quantile,
#     group = .data$quantile
#   )
# ) +
#   ggplot2::geom_line(lwd = 2, alpha = .5) +
#   ggplot2::labs(
#     title = paste(
#       "model fit with 95% interval, Topic =",
#       stringr::str_to_title(topic)
#     ),
#     x = "Covariate value",
#     y = "P(topic | new data)"
#   )
# }
