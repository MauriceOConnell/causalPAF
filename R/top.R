#' @title Show the top part of an object
#'
#' @description This combines head with a number of columns
#' @param mat a matrix or dataframe
#' @param n Number of rows/columns to display
#' @export
#' @keywords models Regression
#' @seealso \code{\link[utils]{head}}
#' @return NULL
#' @aliases causal
#' @examples \dontrun{
#' # I don't want you to run this
#' }
#' x = matrix(rnorm(100), ncol = 10)
#' top(x)
top <- function(mat, n = 5){
  print(mat[1:n, 1:n])
  return(NULL)
}
