#' Infer the Cell-cell communication pattern in a given dataset.
#'
#' @import dplyr
#' @import Seurat
#' @import CellChat
#' @param SeuratObj The Seurat object of the dataset.
#' @param use_spatial Whether use spatial information. Set to True only when using spatial transcriptomics data.
#' @param db_use The species of the object. Set to human or mouse.
#' @param scale_factors The spatial scaling factors of the dataset. Use only when use_spatial = True.
#' @return The CellChat Object containing Cell-cell communication information.
#' @export
#' @examples
#' \donttest{
#' Infer_CCI(SeuratObj, use_spatial = F, db_use = "human", scale_factors = NULL)
#' }
Infer_CCI <- function(SeuratObj, use_spatial = F, db_use = "human", scale_factors = NULL) {

  if ("RNA" %in% names(SeuratObj@assays)) {
    data <- SeuratObj@assays$RNA$data
  } else if ("SCT" %in% names(SeuratObj@assays)) {
    data <- SeuratObj@assays$SCT$data
  } else {
    data <- SeuratObj@assays$Spatial$data
  }
  meta <- data.frame(labels = Idents(SeuratObj), row.names = names(Idents(SeuratObj)))

  if (use_spatial) {
    meta$slices <- "slice1" %>% as.factor()
    spatial_locs <- Seurat::GetTissueCoordinates(SeuratObj, scale = NULL, cols = c("imagerow", "imagecol"))
    if (is.null(scale_factors)) {
      stop("No spatial factors input!\n")
    } else {
      spot_size <- 65 # the theoretical spot size (um) in 10X Visium
      conversion_factor <- spot_size / scale_factors$spot_diameter_fullres
      spatial_factors <- data.frame(ratio = conversion_factor, tol = spot_size / 2)
    }
    MLcellchat <- createCellChat(
      object = data, meta = meta, group.by = "labels",
      datatype = "spatial", coordinates = spatial_locs, spatial.factors = spatial_factors
    )
  } else {
    MLcellchat <- createCellChat(object = data, meta = meta, group.by = "labels", datatype = "RNA")
  }
  # Input the LR interaction database from CellChat
  if (db_use == "human") {
    cellchatdb_use <- CellChatDB.human
  } else if (db_use == "mouse") {
    cellchatdb_use <- CellChatDB.mouse
    # cellchatdb_use <- subsetDB(cellchatdb_use, search = "Secreted Signaling") # use Secreted Signaling
  } else {
    stop("Must choose a validated CCI database!")
  }

  MLcellchat@DB <- cellchatdb_use
  MLcellchat <- subsetData(MLcellchat) # This step is necessary even if using the whole database
  MLcellchat <- identifyOverExpressedGenes(MLcellchat)
  MLcellchat <- identifyOverExpressedInteractions(MLcellchat)
  if (use_spatial) {
    MLcellchat <- computeCommunProb(MLcellchat,
      type = "truncatedMean", trim = 0.1,
      distance.use = TRUE, interaction.range = 250, scale.distance = 0.01,
      contact.dependent = TRUE, contact.range = 100
    )
  } else {
    MLcellchat <- computeCommunProb(MLcellchat)
  }

  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups

  MLcellchat <- filterCommunication(MLcellchat, min.cells = 10)
  MLcellchat <- computeCommunProbPathway(MLcellchat)
  MLcellchat <- aggregateNet(MLcellchat)
  Weight <- MLcellchat@net$prob

  return(MLcellchat)
}

#' Infer regulator and their corresponding TGs.
#'
#' @param DB The gene-gene interaction database, normally between TF and TG.
#' @param exp The expression matrix.
#' @param selected_genes Select specific TGs. If set to NULL, will by default use all TGs in the database.
#' @param pv_threshold Threshold for the p-value of the Fisher test. Set to >= 1 if without extra filtering.
#' @import rlang
#' @import parallel
#' @return The filtered database
FindRegulator <- function(DB, exp, selected_genes = NULL,pv_threshold = 1) {

  # Load necessary packages

  message("Filtering regulator ... \n")
  colnames(DB) <- c("From","To")

  # Clarify the genes used and update the prior dataset
  all_genes <- rownames(exp)
  Regulator_All <- intersect(DB$From, all_genes) %>% unique()
  Target_All <- intersect(DB$To, all_genes) %>% unique()

  if (!is.null(selected_genes)) {
    selected_genes <- intersect(Target_All, selected_genes)
  } else {
    stop("Must specify target genes!")
  }

  DB_filtered <- DB %>%
    dplyr::filter(From %in% Regulator_All) %>%
    dplyr::filter(To %in% Target_All)

  # if(pv_threshold >= 1){
  #   colnames(DB_filtered) <- c("TF","TG")
  #   return(DB_filtered)
  # }
  Regulator_All <- DB_filtered$From %>% unique()
  Target_All <- DB_filtered$To %>% unique()

  # Run the fisher's test
  Regulator_of_Target <- lapply(Regulator_All, function(x) {
    yofx <- DB_filtered[DB_filtered$From == x, 2] %>%
      unlist() %>% as.vector()
    intersect(yofx, all_genes)
  })
  names(Regulator_of_Target) <- Regulator_All

  Regulator_pval <- lapply(Regulator_of_Target, function(x){
    Fisher_Test(subset1 = x, subset2 = selected_genes, background = Target_All)
  })

  Regulator_pval <- unlist(Regulator_pval)
  names(Regulator_pval) <- Regulator_All
  Regulator_filtered <- Regulator_pval[as.numeric(Regulator_pval) < pv_threshold] %>% names()

  Regulator_of_Target <- Regulator_of_Target[Regulator_filtered]
  Regulator_of_Target <- lapply(Regulator_of_Target, function(x){intersect(x, selected_genes)})

  # Calculate the co-expression
  DB_new <- data.frame(TF = rep(names(Regulator_of_Target), times = sapply(Regulator_of_Target, length)),TG = unlist(Regulator_of_Target))

  return(DB_new)
}


#' Prepare all the input to the HGNN module.
#'
#' @param SeuratObj The Seurat object of the dataset.
#' @param target_type The receiver cell type.
#' @param TGs The target genes used.
#' @param CCC_results The Cellchat object.
#' @param RTFDB The receptor-TF database.
#' @param TFTGDB The TF-TG database.
#' @param data_dir The output directory of this function
#' @param assay The assay of the Seurat Object to use. By default use the "RNA" assay.
#' @param datatype The datatype of the assay. By default use the "scale.data".
#' @param exp_threshold The threshold for filtering out the low-expression genes.
#' @param CCC_threshold The threshold for filtering out lowly activated Ligand-Receptor pairs.
#' @param Fisher_threshold The threshold for the p-value of the Fisher test. Set to >= 1 if without extra filtering
#' @param species The species of the dataset. If set to "mouse" the gene names will be converted.
#' @return The prepared databases, and the expression matrix that is stored in the data_dir for HGNN module.
#' @import tidyfst
#' @import rlang
#' @import stringr
#' @import igraph
#' @import stats
#' @import utils
#' @export
#' @examples
#' \donttest{
#' Prepare_Input_New <- function(SeuratObj, target_type, TGs, CCC_results, RTFDB, TFTGDB, data_dir)
#' }
Prepare_Input_New <- function(
    SeuratObj, target_type, TGs, CCC_results, RTFDB, TFTGDB, data_dir,
    assay = "RNA", datatype = "scale.data",
    exp_threshold = 0.1, CCC_threshold = 0.05, Fisher_threshold = 1, species = "human") {

  # Pre-process the expression matrix
  message("Getting the expression matrix ... \n")
  Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = assay, datatype = datatype, cutoff = exp_threshold)

  gc()
  # Filter the Lig-Rec database
  message("Filter the Ligand-Receptor results ... \n")
  all_genes <- rownames(Exp_clu)
  LR_Pairprob <- Extract_LR_Prob(CCC_results, target_type = target_type)
  if(species == "mouse"){
    LR_Pairprob[,1:2] <- data.frame(lapply(LR_Pairprob[,1:2], function(x) stringr::str_to_title(tolower(x))))
  }
  LR_Pairprob <- LR_Pairprob[LR_Pairprob$From %in% all_genes, ]
  LR_Pairprob <<- LR_Pairprob[LR_Pairprob$To %in% all_genes, ]
  LR_Pairprob <- LR_Pairprob[which(LR_Pairprob$Weight >= CCC_threshold * max(LR_Pairprob$Weight)), ]


  # Filter the TF-TG database
  message("Filter the TF-TG database ... \n")
  Ncores <- parallel::detectCores()
  TFTGDB <- Filter_DB(TFTGDB, all_genes)
  TFTG_filtered <- FindRegulator(TFTGDB, exp = Exp_clu, selected_genes = TGs,
    pv_threshold = Fisher_threshold)

  # Filter the Rec-TF database
  message("Filter the Receptor-TF database ... \n")

    if(is.data.frame(RTFDB)){

    RTFDB <- Filter_DB(RTFDB, all_genes)
    Recs <- intersect(LR_Pairprob$To %>% unique(), RTFDB$From)
    TFs <- TFTG_filtered$TF %>% unique()
    genes <- unlist(RTFDB) %>% unique()
    gg <- graph_from_data_frame(RTFDB, directed = T, vertices = genes)

    RTF_filtered <- myPageRank(gg, Recs, TFs)

  }else{ # The database is from KEGG

    Recs <- LR_Pairprob$To %>% unique()
    TFs <- TFTG_filtered$TF %>% unique()
    graph_metas <- names(RTFDB)
    RTF_list <- lapply(graph_metas, function(x){
      gg <- RTFDB[[x]]
      temp <- myPageRank(gg, Recs, TFs)
      temp <- temp[,c("regulator","target","PPR")]
      colnames(temp) <- c("Receptor","TF","PPR")
      if(nrow(temp) > 0){
        temp$pathway <- x
        temp[temp$PPR>0,]
      }else{
        temp$pathway <- character(0)
      }

    })
    RTF_all <- do.call(rbind, RTF_list)
    RTF_filtered <- aggregate(RTF_all$PPR, list(RTF_all$Receptor, RTF_all$TF), sum)

  }

  colnames(RTF_filtered) <- c("Receptor","TF","PPR")
  RTF_filtered <- na.omit(RTF_filtered)
  RTF_filtered <<- RTF_filtered[RTF_filtered$PPR>0,]
  TFTG_filtered <<- TFTG_filtered %>% filter(TF %in% TFs)

  filen <- paste0(data_dir,"/RecTFDB.csv")
  write.table(RTF_filtered, file = filen, quote = F, sep = " ")
  filen <- paste0(data_dir,"/TFTGDB.csv")
  write.table(TFTG_filtered, file = filen, quote = F, sep = " ")

  gene_used <- c(RTF_filtered$Receptor, RTF_filtered$TF, TFTG_filtered$TF, TFTG_filtered$TG) %>% unlist() %>% unique()
  Exp_clu <<- Exp_clu[gene_used, ]

  filen <- paste0(data_dir,"/ExpressionCount.csv")
  write.table(Exp_clu, file = filen, quote = F, sep = " ")

}

#' Calculate the PRS values.
#'
#' @param Exp_clu The expression matrix.
#' @param RTT_results The dataframe of activated pathways, generated from the HGNN module.
#' @param cutoff The threshold for filtering out lowly co-expresssed TF-TG pairs.
#' @return The dataframe of PRS values.
#' @import ranger
#' @import tidymodels
#' @import foreach
#' @export
#' @examples
#' \donttest{
#' RecTFTG <- PRS_calc(Exp_clu, ress, cutoff = 0.15)
#' }
PRS_calc <- function(Exp_clu, RTT_results, cutoff = 0.15) {
  method <- tolower(method)

  # make sure gene lists are in DB
  message("Preprocessing data...\n")
  colnames(RTT_results)[1:3] <- c("Receptor", "TF", "TG")
  RecAll <- RTT_results$Receptor %>% unique()
  TFAll <- RTT_results$TF %>% unique()
  TGAll <- RTT_results$TG %>% unique()

  RTT_results$TFTG <- paste0(RTT_results$TF, "_", RTT_results$TG)

  RecofTFTG <- split(RTT_results, RTT_results$TFTG)
  TFTG_All <- names(RecofTFTG)

  # Use parallel computing
  Ncores <- parallel::detectCores()
  cl <- snow::makeSOCKcluster(min(24, Ncores))
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(min = 0, max = length(TFTG_All), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n);
  opts <- list(progress = progress)
  message(paste0("\n Estimating the nonlinear model, including ", length(TFTG_All), " TF-TG pairs. \n"))

  ModelList <- foreach::foreach(
    i = 1:length(TFTG_All), .packages = c(
      "dplyr", "ranger",
      "tidymodels", "rlang", "parsnip",
      "recipes"
    ),
    .export = c(
      "TFTG_All", "Exp_clu", "cutoff", "Run_Regression"
    ),
    .options.snow = opts, .errorhandling = "pass"
  ) %dopar% {

    tftg_name <- TFTG_All[i]
    tf_name <- strsplit(tftg_name, "_", fixed = T)[[1]][1]
    tg_name <- strsplit(tftg_name, "_", fixed = T)[[1]][2]
    message("Current TF-TG pair: ", tftg_name, "\n")

    Exp_tftg <- Exp_clu[tf_name, ] * Exp_clu[tg_name, ]
    Exp_tftg <- as.data.frame(Exp_tftg)
    tar_allrec <- RecofTFTG[[tftg_name]]$Receptor %>% unique()
    tar_allrec <- setdiff(tar_allrec,c(tf_name,tg_name))
    ncells <- sum(Exp_clu[tg_name, ]>0)

    if (is_empty(tar_allrec)) {
      cat("No corresponding Rec, skip this TF... \n")
      next
    }
    if (sum(Exp_tftg > 1e-05) < cutoff * ncells) {
      cat("No enough expression, skip this TF... \n")
      next
    }

    Exp_rec_used <- Exp_clu[tar_allrec, ] %>% as.matrix()

    if (min(dim(Exp_rec_used)) > 1) {
      Exp_rec_used <- t(Exp_rec_used)
      Exp_rec_used <- scale(Exp_rec_used)
    }
    if (!max(Exp_rec_used)) {
      cat("No expression of receptors, skip this TF... \n")
      next
    }

    data_use <- cbind(Exp_tftg, Exp_rec_used)
    colnames(data_use)[1:2] <- c("y", tar_allrec[1])
    colnames(data_use) <- gsub("-", "_", colnames(data_use))
    Run_Regression(data_use = data_use)

  }
  names(ModelList) <- TFTG_All

  close(pb)
  parallel::stopCluster(cl)
  rm(cl)
  gc()

  message("Estimating importance scores of Recs...\n")

  coefList <- Calculate_Importance(ModelList)
  RecTFTG <- do.call(rbind, coefList) %>% as.data.frame()
  RecTFTG_final <- RecTFTG %>% tidyr::separate(To, into = c("SSC", "Target"), sep = "_")

  RecTFTG_final <- RecTFTG_final[,c("From", "SSC", "Target", "Weight")]
  colnames(RecTFTG_final) <- c("Receptor", "SSC", "Target", "Weight")

  return(RecTFTG_final)
}

#' Run the single expression.
#'
#' @param data_use The expression of genes.
#' @return The regression model object.
#' @import recipes
#' @import parsnip
#'
Run_Regression <- function(data_use) {

      data_rec <-
        recipe(y ~ ., data = data_use)
      data_juiced <- data_rec %>%
        prep() %>%
        juice()
      rf_model <- rand_forest(mtry = 3, trees = 500, mode = "regression") %>%
        set_engine("ranger", importance = "permutation")
      rf_fit <- fit(rf_model, y ~ ., data_juiced)
      Model_fit <- rf_fit$fit

  return(Model_fit)
}

#' Extracting the PRS values from the regression models.
#'
#' @param ModelList The list of regression model objects.
#' @return The dataframe of PRS values.
#'
Calculate_Importance <- function(ModelList) {
  coefList <- list()

  for (x in 1:length(ModelList)) {
    cur_model <- ModelList[[x]]
    cur_gene <- names(ModelList)[x]
    coefs <- cur_model$variable.importance
    if (!is.null(coefs)) {
      temp_df <- cbind(names(coefs), rep(cur_gene, length(coefs)), coefs %>% as.numeric()) %>% as.data.frame()
      colnames(temp_df) <- c("From", "To", "Weight")
      temp_df$From <- gsub("_", "-", temp_df$From)
      temp_df$Weight <- as.numeric(temp_df$Weight)
      temp_df <- temp_df[temp_df$Weight > 1e-05, ] %>% na.omit()
    } else {
      temp_df <- NULL
    }
    coefList[[cur_gene]] <- temp_df
  }
  return(coefList)
}
