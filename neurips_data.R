#' Neurips TDM
#' 
#' Term-document matrix compiled from papers at the Neurips conference, 1987 - 2015.  
#' 
#' @docType data 
#' 
#' @usage data(neurips_tdm) 
#' 
#' @format An integer matrix. 
#' 
#' @references Poisson Random Fields for Dynamic Feature Models'. Perrone V., Jenkins P. A., Spano D., Teh Y. W. (2016)
#' 
#' @source \href{https://archive.ics.uci.edu/ml/datasets/NIPS+Conference+Papers+1987-2015}{UCI Repository}
"neurips_tdm"



#' Neurips Vocab 
#' 
#' #' The entire vocabulary of words used in papers at the Neurips conference, 1987 - 2015.  
#' 
#' @docType data 
#' 
#' @usage data(neurips_words) 
#' 
#' @format A character vector. 
#' 
#' @references Poisson Random Fields for Dynamic Feature Models'. Perrone V., Jenkins P. A., Spano D., Teh Y. W. (2016)
#' 
#' @source \href{https://archive.ics.uci.edu/ml/datasets/NIPS+Conference+Papers+1987-2015}{UCI Repository}
"neurips_words"

#' Neurips Covariate
#' 
#' The year each neurips paper was published, binned into five year intervals (inclusive). 
#' 
#' 
#' @docType data 
#' 
#' @usage data(neurips_years) 
#' 
#' @format A data frame of 1 varaible: the year each paper was published, in five year intervals (i.e., 
#' this is a discrete covariate).  
#' 
#' @references Poisson Random Fields for Dynamic Feature Models'. Perrone V., Jenkins P. A., Spano D., Teh Y. W. (2016)
#' 
#' @source \href{https://archive.ics.uci.edu/ml/datasets/NIPS+Conference+Papers+1987-2015}{UCI Repository}
"neurips_years"