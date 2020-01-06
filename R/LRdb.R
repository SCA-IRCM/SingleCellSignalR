#' Ligand/Receptor interactions data table
#'
#' @format A data frame with 3369 rows of 13 variables:
#' \describe{
#'   \item{ligand}{ligand gene symbol}
#'   \item{receptor}{receptor gene symbol}
#'   \item{ligand.name}{receptor gene name}
#'   \item{receptor.name}{ligand gene name}
#'   \item{ligand.synonyms}{ligand synonymes}
#'   \item{receptor.synonyms}{receptor synonymes}
#'   \item{ligand.altern.names}{ligand alternative names}
#'   \item{receptor.altern.names}{receptor alternative names}
#'   \item{source}{provenance of the interaction}
#'   \item{PMIDs}{PubmedID of the publication reporting the interaction}
#'   \item{cells.L}{cell types from which the ligand is exclusive}
#'   \item{cells.R}{cell types from which the receptor is exclusive}
#'   \item{remarks}{remarks}
#'   ...
#' }
#' @source \url{EDF R\&D}
"LRdb"
