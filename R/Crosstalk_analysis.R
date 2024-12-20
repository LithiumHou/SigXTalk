#' Count the crosstalk pathways for each gene.
#'
#' @param CC_results The dataframe of PRS values.
#' @param KeyGenes The genes to be considered. By default, will consider all the genes in the pathways
#' @param data_type The type of the KeyGene (SSC or Target).
#' @param verbose Whether show notifications
#' @return The dataframe containing the number of crosstalk pathways of each gene.
#' @export
#' @examples
#' \donttest{
#' counts <- Count_Crosstalk(CC_results, data_type = "Target")
#' }
Count_Crosstalk <- function(CC_results, KeyGenes = NULL, data_type = "SSC", verbose = T) {
  if (data_type == "SSC") {
    if (is.null(KeyGenes)) {
      KeyGenes <- CC_results$SSC %>% unique()
    }
    KeySSC <- KeyGenes
    SSC_used <- intersect(KeySSC, CC_results$SSC) %>% unique()
    if (!length(SSC_used)) {
      stop("The input SSCs have no crosstalk!")
    }
    Count_Xmodule <- array(0, length(SSC_used))
    for (i in 1:length(SSC_used)) {
      SSC <- SSC_used[i]
      temp_df <- CC_results[CC_results$SSC == SSC, ]
      nXTalk <- nrow(temp_df)
      if (nXTalk > 1) {
        if (verbose) {
          message(paste0("There are ", nXTalk, "crosstalk pathways for SSC:", SSC, "!\n"))
        }
        Count_Xmodule[i] <- nXTalk
      } else {
        if (verbose) {
          message(paste0("There is no crosstalk for SSC:", SSC, "!\n"))
        }
        Count_Xmodule[i] <- 0
      }
    }
    names(Count_Xmodule) <- SSC_used
  } else if (data_type == "Target") {
    if (is.null(KeyGenes)) KeyGenes <- CC_results$Target %>% unique()
    KeyTG <- KeyGenes
    TG_used <- intersect(KeyTG, CC_results$Target) %>% unique()
    if (!length(TG_used)) {
      stop("The input TGs have no crosstalk!")
    }
    Count_Xmodule <- array(0, length(TG_used))
    for (i in 1:length(TG_used)) {
      tg <- TG_used[i]
      temp_df <- CC_results[CC_results$Target == tg, ]
      nXTalk <- nrow(temp_df)
      if (nXTalk > 1) {
        if (verbose) {
          message(paste0("There are", nXTalk, " crosstalk pathways for Target:", tg, "!\n"))
        }
        Count_Xmodule[i] <- nXTalk
      } else {
        if (verbose) {
          message(paste0("There is no crosstalk for Target:", tg, "!\n"))
        }
        Count_Xmodule[i] <- 0
      }
    }
    names(Count_Xmodule) <- TG_used
  }else{
    stop("Must specify the type of gene to be counted (SSC or Target)!")
  }

  return(Count_Xmodule)
}


#' Calculate the TRS values
#'
#' @param CC_results The dataframe of PRS values.
#' @param data_type The type of the KeyGene (SSC or target).
#' @return The dataframe of TRS.
#' @export
#' @examples
#' \donttest{
#' CC_pair_results <- Aggregate_Causality(CC_results, data_type = "Target")
#' }
Aggregate_Causality <- function(CC_results, data_type = "Target") {
  if (data_type == "Target") {
    CC_df <- CC_results[, c("Receptor", "Target", "Weight")]
    CC_df$Weight <- 0.999 * CC_df$Weight / max(CC_df$Weight) # Normalize
    CC_aggregated <- aggregate(CC_df$Weight, list(CC_df$Rec, CC_df$Target), function(x) Nonlinear_Sum(x, type = "sum"))
    colnames(CC_aggregated) <- c("Receptor", "Target", "Weight")
    CC_sorted <- CC_aggregated[order(CC_aggregated$Rec, -CC_aggregated$Weight), ]
  } else if (data_type == "SSC") {
    CC_df <- CC_results[, c("Receptor", "SSC", "Weight")]
    CC_df$Weight <- 0.999 * CC_df$Weight / max(CC_df$Weight) # Normalize
    CC_aggregated <- aggregate(CC_df$Weight, list(CC_df$Receptor, CC_df$TF), function(x) Nonlinear_Sum(x, type = "sum"))
    colnames(CC_aggregated) <- c("Receptor", "TF", "Weight")
    CC_sorted <- CC_aggregated[order(CC_aggregated$Receptor, -CC_aggregated$Weight), ]
  }

  return(CC_sorted)
}


#' Obtain a FIDELITY matrix for a fixed target gene.
#'
#' @param CC_results The dataframe of PRS values.
#' @param KeyTG The target to be focused on.
#' @return The list containing the fidelity matrix within 3 types crosstalk modules.
#' @export
#' @examples
#' \donttest{
#' Fid_mat <- Calculate_Fidelity_Matrix(CC_results, KeyTG)
#' }
Calculate_Fidelity_Matrix <- function(CC_results, KeyTG) {

  colnames(CC_results) <- c("Receptor","SSC","Target","Weight")

  if (!KeyTG %in% CC_results$Target) {
    stop("The input gene must be a target gene!")
    return(NULL)
  }
  if (length(KeyTG) != 1) {
    stop("The input gene must be a target gene!")
    return(NULL)
  }
  CC_used <- CC_results %>% dplyr::filter(Target == KeyTG)
  CC_mat <- df_mat(CC_used, row = Receptor, col = SSC, value = Weight)
  CC_mat[is.na(CC_mat)] <- 0
  Fid_SSC <- sweep(CC_mat, 1, rowSums(CC_mat), FUN = "/") # for a given Receptor-Target, the fidelity for each SSC
  Fid_Receptor <- sweep(CC_mat, 2, colSums(CC_mat), FUN = "/") # for a given SSC-Target, the fidelity for each rec
  Fid_all <- CC_mat / sum(CC_mat) # fidelity for each pathway Receptor-SSC-Target

  Fid_results <- list(
    SSC = Fid_SSC,
    Receptor = Fid_Receptor,
    all = Fid_all
  )
  return(Fid_results)
}

#' Obtain a SPECIFICITY matrix for a fixed signal(receptor).
#'
#' @param CC_results The dataframe of PRS values.
#' @param KeyRec The signal(receptor) to be focused on.
#' @return The list containing the specificity matrix within 3 types crosstalk modules.
#' @export
#' @examples
#' \donttest{
#' Spe_mat <- Calculate_Specificity_Matrix(CC_results, KeyRec)
#' }
Calculate_Specificity_Matrix <- function(CC_results, KeyRec) {

  if (!KeyRec %in% CC_results$Receptor) {
    message("The input gene must be a receptor!")
    return(NULL)
  }
  if (length(KeyRec) != 1) {
    stop("The input gene must be a receptor!")
    return(NULL)
  }
  CC_used <- CC_results %>% dplyr::filter(Receptor == KeyRec)
  CC_used <- CC_used[!duplicated(CC_used[,1:3]),]
  CC_mat <- tidyfst::df_mat(CC_used, row = Target, col = SSC, value = Weight)
  CC_mat[is.na(CC_mat)] <- 0

  Spe_SSC <- sweep(CC_mat, 1, rowSums(CC_mat), FUN = "/") # for a given Receptor-Target, the specificity for each SSC
  Spe_Target <- sweep(CC_mat, 2, colSums(CC_mat), FUN = "/") # for a given Receptor-SSC, the specificity for each Target
  Spe_all <- CC_mat / sum(CC_mat) # specificity for each pathway Receptor-SSC-Target

  Spe_results <- list(
    SSC = Spe_SSC,
    Target = Spe_Target,
    all = Spe_all
  )
  return(Spe_results)
}

#' Calculate the L-R-T pathway fidelity out of all the pathways with the same SSC and Target
#'
#' @param CC_results The dataframe of PRS values.
#' @param KeyTG The target to be focused on.
#' @param KeySSC The SSC to be focused on
#' @return The array containing the fidelity of the signal, out of all the pathways with the same SSC and Target.
#' @export
#'
Calculate_Pathway_Fidelity <- function(CC_results, KeyTG, KeySSC) {

  if (length(KeyRec) + length(KeySSC) != 2) stop("Must input one SSC and one Target! \n")
  CC_used <- CC_results %>%
    dplyr::filter(Target == KeyTG) %>%
    dplyr::filter(SSC == KeySSC)
  Fid_pathway <- CC_used$Weight / sum(CC_used$Weight)
  names(Fid_pathway) <- CC_used$Target
  return(Fid_pathway)
}

#' Calculate the L-R-T pathway specificity out of all the pathways with the same Signal and SSC
#'
#' @param CC_results The dataframe of PRS values.
#' @param KeyRec The signal (receptor) to be focused on.
#' @param KeySSC The SSC to be focused on
#' @return The array containing the fidelity of the signal, out of all the pathways with the same Signal and SSC.
#' @export
#'
Calculate_Pathway_Specificity <- function(CC_results, KeyRec, KeySSC) {

  if (length(KeyRec) + length(KeySSC) != 2) stop("Must input one Receptor and one SSC! \n")
  CC_used <- CC_results %>%
    dplyr::filter(Receptor == KeyRec) %>%
    dplyr::filter(SSC == KeySSC)
  Spe_pathway <- CC_used$Weight / sum(CC_used$Weight)
  names(Spe_pathway) <- CC_used$Target
  return(Spe_pathway)
}


#' Calculate the pairwise pathway fidelity/specificity
#'
#' @param CC_pair_results The dataframe of TRS values.
#' @param KeyGene The Signal/Target to be focused on
#' @param type Whether to calculate fidelity (Fid) or specificity (Spe)
#' @return The array containing the fidelity of the signal, out of all the pathways with the same Signal and SSC.
#' @export
#' @examples
#' \donttest{
#' df <- Calculate_Pairwise(CC_pair_results, KeyGene = NULL, type = "Fid")
#' }
Calculate_Pairwise <- function(CC_pair_results, KeyGene = NULL, type = "Fid"){

  if(type == "Fid" || type == "Fidelity"){
    if(is.null(KeyGene)){
      KeyGene <- CC_pair_results$Target
    }else{
      KeyGene <- intersect(KeyGene, CC_pair_results$Target)
    }
    if(length(KeyGene) == 0){
      stop("The input gene is not a target gene!\n")
    }else{
      CC_used <- dplyr::filter(CC_pair_results, Target %in% KeyGene)
      temp <- CC_used %>% group_by(Target) %>% mutate(Fidelity = Weight / sum(Weight))
      temp <- temp[,c("Receptor","Target","Fidelity")]
    }
  }
  if(type == "Spe" || type == "Specificity"){
    if(is.null(KeyGene)){
      KeyGene <- CC_pair_results$Receptor
    }else{
      KeyGene <- intersect(KeyGene, CC_pair_results$Receptor)
    }
    if(length(KeyGene) == 0){
      stop("The input gene is not a Receptor!\n")
    }else{
      CC_used <- dplyr::filter(CC_pair_results, Receptor %in% KeyGene)
      temp <- CC_used %>% group_by(Receptor) %>% mutate(Specificity = Weight / sum(Weight))
      temp <- temp[,c("Receptor","Target","Specificity")]
    }
  }
  return(temp)
}
