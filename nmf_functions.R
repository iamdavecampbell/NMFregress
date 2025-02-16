#' solve_nmf
#'
#' Solve the non-negative matrix factorization problem with specified inputs.
#'
#' @param input An object of class nmf_input.
#'
#' @param user_anchors User-specified anchors. These can be specified if the user wishes to find anchors via their
#' own method. Defaults to null.
#'
#' @return An object of class nmf_output. This contains the word-topic matrix,
#' the block of the word-topic matrix associated with anchor words, the anchor words,
#' the order the anchor words were found, the entire vocabulary (sorted to match the word-topic
#' matrix), and the number of topics. These are all needed to print a summary of the topics.
#'
#' @examples
#' neurips_input = create_input(neurips_tdm, neurips_words,
#'    topics = 10, project = TRUE, proj_dim = 500)
#' neurips_output = solve_nmf(neurips_input)
#'
#' @export
solve_nmf = function(input, user_anchors = NULL){

  ##### check/convert input types
  if(class(input) != "nmf_input"){
    stop("Input must be an object of class nmf_input.")
  }
  if(!(is.null(user_anchors))){
    if(!(sum(user_anchors %in% input$vocab) == length(user_anchors))){
      stop("User_anchors must be a subset of vocab.")
    }
    if(length(user_anchors) > input$topics){
      stop("Number of user anchors should not exceed number of topics.")
    }
  }

  ##### simple function to get norm of vector
  #get_norm = function(x){                       ############ Replaced by sum (l1 normalization)
  #  return(sqrt(sum(x^2)))                      ############ More interpratable (probabilities)
  #}

  ##### find anchors w/ qr decomposition (using Gaussian random projection if specified)
  ##### store anchors
  if(!(is.null(user_anchors))){
    qr_tdm = qr(t(input$tdm), LAPACK = TRUE)
    extract_order_anchors = anchors = union(user_anchors, input$vocab[qr_tdm$pivot])[1:input$topics]
    anchors = sort(anchors)
    cat("Anchors recovered -- solving nmf.\n")
  }else if(input$project == TRUE){
    proj_mat = matrix(stats::rnorm(input$proj_dim*ncol(input$tdm), 0, 1), nrow = input$proj_dim)
    qr_tdm = qr(((1/sqrt(input$proj_dim))*proj_mat) %*% t(input$tdm), LAPACK = TRUE)
    extract_order_anchors = anchors = input$vocab[qr_tdm$pivot][1:input$topics]
    anchors = sort(anchors)
    cat("Anchors recovered -- solving nmf.\n")
  }else{
    qr_tdm = qr(t(input$tdm), LAPACK = TRUE)
    extract_order_anchors = anchors = input$vocab[qr_tdm$pivot][1:input$topics]
    anchors = sort(anchors)
    cat("Anchors recovered -- solving nmf.\n")
  }

  ##### set up necessary matrices for nmf
  anchor_rows = match(anchors, input$vocab)
  non_anchor_rows = setdiff(1:length(input$vocab), anchor_rows)
  anchor_block = input$tdm[anchor_rows,]
  other_block = input$tdm[non_anchor_rows,]
  Y = matrix(rep(0, length(anchor_rows)*length(non_anchor_rows)), ncol = length(anchor_rows))

  ##### prepare matrices for input to nnls (using Haddamard random projection if specified)
  if(input$project == TRUE){

      ##### set up parameters for random projection
      power_of_two = 2^ceiling(log(ncol(anchor_block), base = 2))
      num_zeros = power_of_two - ncol(anchor_block)
      d = sample(c(1,-1), power_of_two, replace = TRUE)

      ##### apply random projection
      ##### details are messy, refer to thesis
      A = d*rbind(t(anchor_block), matrix(rep(0, num_zeros*nrow(anchor_block)), nrow = num_zeros))
      A = as.matrix(apply(A, FUN = phangorn::fhm, MARGIN = 2))
      B = d*rbind(t(other_block), matrix(rep(0, num_zeros*nrow(other_block)), nrow = num_zeros))
      B = as.matrix(apply(B, FUN = phangorn::fhm, MARGIN = 2))
      selected = sample(c(1,0), power_of_two,
                  prob = c(input$proj_dim/power_of_two,1-input$proj_dim/power_of_two), replace = TRUE)
      selected = which(selected == 1)
      A = rep(sqrt(power_of_two/input$proj_dim), length(selected))*A[selected,]
      B = rep(sqrt(power_of_two/input$proj_dim), length(selected))*B[selected,]

  }else{

      ##### transpose matrices for input to nnls
      A = t(anchor_block)
      B = t(other_block)

  }

  ##### solve non-negative least squares problem for each word
  for(i in 1:nrow(Y)){
    Y[i,] = as.vector(nnls::nnls(A = A, b = B[,i])$x)
    if(i %% 500 == 0){cat(i, " rows of word-topic matrix recovered.\n")}
  }

  ##### solve for diagonal scaling matrix (lambda), phi, and theta
  Y_with_ones = rbind(diag(length(anchor_rows)), Y)
  colnames(Y_with_ones) = input$vocab[anchor_rows]
  lambdas = 1/apply(Y_with_ones, FUN = sum, MARGIN = 2)
  anchor_order = order(lambdas)
  phi = Y_with_ones %*% diag(lambdas)
  colnames(phi) = colnames(Y_with_ones)
  theta = diag(1/lambdas) %*% anchor_block
  phi = phi[,anchor_order]
  theta = theta[anchor_order,]
  rownames(theta) = colnames(phi)
  ##### return object of class nmf_output
  to_return = list(phi = phi,
                   theta = theta,
                   anchors = anchors[anchor_order],
                   extract_order_anchors = extract_order_anchors,
                   lambdas = lambdas,
                   vocab = input$vocab[c(anchor_rows, non_anchor_rows)],
                   topics = input$topics,
                   covariates = input$covariate,
                   sum_theta_over_docs =  apply(theta,2,sum))
  class(to_return) = "nmf_output"
  cat("Complete -- outputting object of class nmf_output.\n")
  return(to_return)

}
#' print_top_words
#'
#' Print the most important words associated with each topic.
#'
#' @param output An object of class nmf_output.
#'
#' @param n The number of words to print.
#'
#' @return A list with one element for each topic. These are labeled by anchor words. Each
#' list element contains the top words in that topic, printed in order of importance (largest
#' pseudo-count in the word-topic matrix.)
#'
#' @examples
#' neurips_input = create_input(neurips_tdm, neurips_words,
#'    topics = 10, project = TRUE, proj_dim = 500)
#' neurips_output = solve_nmf(neurips_input)
#' print_top_words(neurips_output, n = 25)
#'
#' @export
print_top_words = function(output, n = 10){

  ##### check/convert input types
  if(class(output) != "nmf_output"){
    stop("Output must be an object of class nmf_output.")
  }
  n = as.integer(n)

  ##### prepare, fill, and return matrix of top words
  word_list = list()
  for(i in 1:output$topics){
    word_list[[i]]= output$vocab[order(output$phi[,i], decreasing = TRUE)][1:n]
    names(word_list)[i] = output$anchors[i]
  }
  return(word_list)

}




#' get_lambda
#'
#' Extracts the inverse lambdas and matches them to the anchor words.
#'
#' @param output An object of class nmf_output.
#'
#' @return A data frame of the anchor words and the inverse lambda (lambdacrit)
#'
#' @examples
#' neurips_input = create_input(neurips_tdm, neurips_words,
#'    topics = 10, project = FALSE)
#' neurips_output = solve_nmf(neurips_input)
#' get_lambda(neurips_output)
#'
#' @export
get_lambda = function(output){
   lambdas = output$phi[1:output$topics,
                       1:output$topics]
   words_and_lambda = data.frame(anchors = output$vocab[1:output$topics],
          lambdacrit = apply(lambdas,1,function(x){1/x[x!=0]})
   )
   return(words_and_lambda[ order(words_and_lambda$lambdacrit),])
}






#' get_reconstruction_error.
#'
#' Calculates the goodness of fit measure for matrix reconstruction
#' The Frobenius norm is computed based on a TDM reconstructed with incrementally increasing topics
#'
#' @param output An object of class nmf_output from solve_nmf
#'
#' @param input all of the information passed into the original solve_nmf function
#'
#' @param Ntopics optional, the maximum number of topics to consider
#'
#' @return a data frame with columns 'topics' taking values from 0: Ntopics.  Note that
#' zero is actually the Frobenius norm of the original TDM.  When 'topics'>0 the
#' norm is applied to "reconstruction(based on # topics) - original"
#'
#' @examples
#' neurips_input = create_input(neurips_tdm, neurips_words,
#'    topics = 10, project = FALSE)
#' neurips_output = solve_nmf(neurips_input)
#' get_reconstruction_error(neurips_output, neurips_input)
#'
#' @export
get_reconstruction_error = function(output, input, Ntopics = ncol(output$phi)){
  # find the indices of the anchors from the original data vocabulary


  frobeniusnorm = data.frame(topics = 0:Ntopics,Frobenius_norm = rep(NA, Ntopics+1))
  docs_2_keep = which(output$sum_theta_over_docs>0)
  rownames(input$tdm) = input$vocab
  # handle origin; no topics case
  frobeniusnorm[1,"Frobenius_norm"] = norm(tdm_hat[,docs_2_keep],
                                           type = "f")
  rownames(output$phi) = output$vocab
  # Ntopics= 1. handle differently since indexing breaks the matrix and requires a transpose
  tdm_hat = output$phi[output$vocab,1] %*% t(output$theta[1,]) ##<-- there has been a shuffling of rows of phi
  frobeniusnorm[2,"Frobenius_norm"] = norm(tdm_hat[,docs_2_keep] -
                                           input$tdm[output$vocab,docs_2_keep],
                                           type = "f")

    for (topic_index in 2:ncol(output$phi)){
              tdm_hat = output$phi[output$vocab,1:topic_index] %*% output$theta[1:topic_index,]
              frobeniusnorm[topic_index+1,"Frobenius_norm"] = norm(tdm_hat[,docs_2_keep] -
                                                                   input$tdm[output$vocab,docs_2_keep],
                                                                   type = "f")
  }

  p = frobeniusnorm|> ggplot(aes(x = topics, y = Frobenius_norm))+geom_point()
  print(p)
  return(list(frobeniusnorm = frobeniusnorm,
              plot = p))
}

