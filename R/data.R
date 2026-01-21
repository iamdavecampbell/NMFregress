#' THE TRAGEDY OF ROMEO AND JULIET by William Shakespeare
#' The text as a term document matrix with covariates
#'
#' @format ## `matrix`
#' A matrix with 788 rows and 5 columns
#' \describe{
#'   \item{(Intercept)}{Act 1 is the baseline for the design matrix}
#'   \item{actstwo}{indicator for Act 2}
#'   \item{actsthree}{indicator for Act 3}
#'   \item{actsfour}{indicator for Act 4}
#'   \item{actsfive}{indicator for Act 5}
#'   ...
#' }
#'
#' @source <https://www.gutenberg.org/files/1513/1513-h/1513-h.htm>
"acts"


#' THE TRAGEDY OF ROMEO AND JULIET by William Shakespeare
#' The text as a term document matrix with covariates
#'
#' @format ## `vector`
#' A vector with 788 elements
#' \describe{
#'   \item{values}{Continuous conversion of line number as a fraction of the
#'         way through each act.}
#'   ...
#' }
#'
#' @source <https://www.gutenberg.org/files/1513/1513-h/1513-h.htm>
"acts_continuous"




#' THE TRAGEDY OF ROMEO AND JULIET by William Shakespeare
#' The text as a term document matrix with covariates
#'
#' @format ## `matrix`
#' A matrix with 2997 rows and  788 columns:
#' \describe{
#'   \item{rows}{unique terms, sorted alphabetically, as defined by the row
#'   labels. Terms have been filtered to remove stopwords}
#'   \item{columns}{Each column is one line of the play}
#'   \item{elements}{Counts of occurence of each term within each line
#'   of the play.}
#' }
#'
#' @source <https://www.gutenberg.org/files/1513/1513-h/1513-h.htm>
"Romeo_and_Juliet_tdm"
