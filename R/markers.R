#' @title Markers
#' @description Provides a table of cell type specific markers.
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
markers = function(category=c("immune","tme","melanoma","bc")){
  a=NULL
  category = match.arg(category, several.ok = TRUE)
  for (i in category){
    if (i=="immune"){
      a=c(a,c(4:6,10:15))
    }
    if (i=="tme"){
      a=c(a,c(7,8))
    }
    if (i=="melanoma"){
      a=c(a,9)
    }
    if (i=="bc"){
      a=c(a,c(seq_len(3)))
    }
  }
  m=markers_default[,a]
  return(m)
}
