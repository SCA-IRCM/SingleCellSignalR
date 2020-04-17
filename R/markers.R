#' @title Markers
#' @description Provides a table of cell type specific markers.
#'
#' @details To use this function the user must have some knowledge about the
#' cell composition of her data set. For instance, if the dataset comes from a
#' breast cancer tumor, the user may select "bc", but also "immune" and "tme"
#' for cells of the microenvironment (see **Examples of use** in the User's
#' Guide vignette).
#'
#' @param category one or several of the following "immune", "tme", "melanoma",
#' "bc"
#'
#' @return The function returns a cell type specific markers table.
#' @export
#'
#' @examples
#' markers(c("immune","bc"))
#'
markers <- function(category=c("immune","tme","melanoma","bc")){
  a <- NULL
  category <- match.arg(category, several.ok=TRUE)
  for (i in category){
    if (i=="immune"){
      a <- c(a,c(4:6,10:15))
    }
    if (i=="tme"){
      a <- c(a,c(7,8))
    }
    if (i=="melanoma"){
      a <- c(a,9)
    }
    if (i=="bc"){
      a <- c(a,c(seq_len(3)))
    }
  }
  m <- markers_default[,a]
  return(m)
}
