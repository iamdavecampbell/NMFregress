#' create_input 
#' 
#' This function creates an object of class nmf_input, later passed to the solve_nmf
#' function.
#' 
#' @param tdm A term-document matrix. This should be a standard R matrix, not a
#' term-document matrix as produced by some text mining packages. 
#' 
#' @param vocab A character vector of all words appearing in the corpus. In other words, 
#' the entire vocabulary. Make sure the words appear in the same order as the corresponding
#' rows in the term-document matrix. 
#' 
#' @param topics The number of topics to fit in the topic model. 
#' 
#' @param project Should a random projection be performed? If yes, random projections are used in 
#' two places. First, they reduce the dimensionality of the term-doducment matrix to 
#' assist in the finding of anchor words. Second, they are used to speed up the non-negative
#' least squares problems that are needed to solve for the word-topic matrix. Defaults to false. 
#' 
#' @param proj_dim What dimension should the random projections take? Both instances effectively
#' reduce the number of documents in the corpus. Despite formal guarantees,
#' experience suggests projections with as few as 100 dimensions can give good results
#' for collections on the order of 10000 documents. Defaults to null. 
#' 
#' @param covariates A design matrix with at least as many rows as the tdm has columns (one for each
#' document, plus any rows encoding constraints on categorical variables). Any trailing rows beyond the 
#' number of doucments will be interpreted as constriants. Must be of the standard R matrix type. 
#' Matrix columns should be named if the user wants regression coefficients named accordingly. 
#' Intercept column must be explicitly included if desired. 
#' 
#' @return Object of class nmf_input. This contains all the information needed to solve the 
#' NMF problem as provided by the function's arguments 
#' arguments 
#' 
#' @examples neurips_input = create_input(neurips_tdm, neurips_words, 
#'    topics = 10, project = TRUE, proj_dim = 100) 
#' 
#' @export 
create_input = function(tdm, vocab, topics, project = FALSE, proj_dim = NULL, covariates = NULL){
  
  ##### check basic input types 
  stopifnot(is.matrix(tdm)) 
  stopifnot(is.character(vocab)) 
  stopifnot(is.numeric(topics)) 
  topics = as.integer(topics) 
  stopifnot(is.logical(project)) 
  
  ##### make sure the user doesn't specify any contradictory inputs 
  if(project == TRUE & is.null(proj_dim)){
    stop("No projection dimensions specified.")
  }
  if(project == TRUE & !(is.null(proj_dim))){
    proj_dim = as.integer(proj_dim)
    if(!(proj_dim > topics & proj_dim < ncol(tdm))){
      stop("Projection dimension must be between the number of topics and the number of documents.")
    }
  }
  if(project == FALSE){
    proj_dim = NULL
  }
  if(!(is.null(covariates))){
    if(!(is.matrix(covariates))){
      stop("Covariates must be in the form of a matrix.")
    }
    if(nrow(covariates) < ncol(tdm)){
      stop("Design matrix must have at least as many rows as there are documents.")
    }
  }
  if(!(topics < ncol(tdm) & topics < nrow(tdm))){
    stop("Number of topics must be less than the both the number of rows and columns in the TDM.")
  } 
  if(length(vocab) != nrow(tdm)){
    stop("Vocab must be the same length as the number of rows in the TDM.")
  }
  
  ##### return object of class nmf_input 
  to_return = list(tdm, vocab, topics, project, proj_dim, covariates)
  names(to_return) = c("tdm", "vocab", "topics", "project", "proj_dim", "covariates")
  class(to_return) = "nmf_input"
  return(to_return)
  
}

