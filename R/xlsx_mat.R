#' @title Blah
#'
#' @description Blah
#' @param n Number of rows/columns to display
#' @export
#' @import readxl
#' @importFrom utils head
#' @return NULL
xlsx_mat <- function(n = 2){
  x = read_excel(path = system.file("ExcelData.xlsx",
                                    package = "causalPAF"))
  x = head(x, n = n)
  # print(x)
  x
}
