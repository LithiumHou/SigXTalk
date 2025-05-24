#' Filter the prior database using genes expressed in the dataset
#'
#' @param DB The raw prior database, whose first two columns contain the "from" and "to" genes, respectively
#' @param allgenes All the expressed genes in the dataset.
#' @return The filtered database.
#' @export
#'
Filter_DB <- function(DB, allgenes) {

  DB_new <- DB[, 1:2]
  colnames(DB_new) <- c("From", "To")
  selfregu <- which(DB[, 1] == DB[, 2])
  if (!is_empty(selfregu)) {
    DB_new <- DB_new[-selfregu, ]
  }
  DB_new <- DB_new[DB_new$From %in% allgenes, ]
  DB_new <- DB_new[DB_new$To %in% allgenes, ]

  return(DB_new)
}


#' Human Receptor-TF interactions
#'
#' A database for the biological interactions from receptors to TFs in Human
#' @references Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nature Methods (2019)
#' @format A dataframe containing receptors, TFs, sources and databases
#' @source \url{https://github.com/LithiumHou/SigXTalk}
"RTF_human"


#' Mouse Receptor-TF interactions
#'
#' A database for the biological interactions from receptors to TFs in Mouse
#' @references Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nature Methods (2019)
#' @format A dataframe containing receptors, TFs, sources and databases
#' @source \url{https://github.com/LithiumHou/SigXTalk}
"RTF_mouse"


#' Human TF-Target gene interactions
#'
#' A database for the biological interactions from TFs to Targets in Human
#' @references Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nature Methods (2019)
#' @format A dataframe containing TFs, target genes, sources and databases
#' @source \url{https://github.com/LithiumHou/SigXTalk}
"TFT_human"


#' Mouse TF-Target gene interactions
#'
#' A database for the biological interactions from TFs to Targets in Mouse
#' @references Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nature Methods (2019)
#' @format A dataframe containing TFs, target genes, sources and databases
#' @source \url{https://github.com/LithiumHou/SigXTalk}
"TFT_mouse"
