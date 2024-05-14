Infer_CCI <- function(SeuratObj, LRDB = NULL, cellchat_output = T, use_spatial = 0,db_use = "human", scale_factors = NULL, Kh = 0.5, nHill = 1, zero.omit = F) {
  
  if ("RNA" %in% names(SeuratObj@assays)) {
    data <- LayerData(SeuratObj, assay = "RNA", layer = "data")
  } else if ("SCT" %in% names(SeuratObj@assays)) {
    data <- LayerData(SeuratObj, assay = "SCT", layer = "data")
  } else {
    data <- LayerData(SeuratObj, assay = "Spatial", layer = "counts")
  }
  meta <- data.frame(labels = Idents(SeuratObj), row.names = names(Idents(SeuratObj)))

  if (cellchat_output) {
    require(CellChat)
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
    } else {
      stop("Must choose a valida CCI database!")
    }
    # cellchatdb_use <- subsetDB(cellchatdb_use, search = "Secreted Signaling") # use Secreted Signaling

    MLcellchat@DB <- cellchatdb_use
    MLcellchat <- subsetData(MLcellchat) # This step is necessary even if using the whole database
    MLcellchat <- identifyOverExpressedGenes(MLcellchat)
    MLcellchat <- identifyOverExpressedInteractions(MLcellchat)
    if(use_spatial){
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

  } else {

    pairLRsig <- LRDB
    allLigand <- as.character(pairLRsig$From)
    allRec <- as.character(pairLRsig$To)
    nPairs <- nrow(pairLRsig)
    nCluster <- nlevels(meta$labels)

    # Subset the data
    data_use <- data / max(data)
    gene_use <- unique(c(allLigand, allRec))
    data_use <- data_use[row.names(data_use) %in% gene_use, ] %>% as.matrix()

    rm(data)

    data_use_avg <- aggregate(t(data_use), list(meta$labels), FUN = mean)
    data_use_avg <- t(data_use_avg[, -1])
    colnames(data_use_avg) <- levels(meta$labels)

    rm(data_use)
    data_use_avg[is.na(data_use_avg)] <- 0

    dataLavg <- data_use_avg[allLigand, ]
    dataRavg <- data_use_avg[allRec, ]

    Weight <- array(0, dim = c(nCluster, nCluster, nPairs))

    for (i in 1:nPairs) {
      # ligand/receptor
      dataLR <- Matrix::crossprod(matrix(dataLavg[i, ], nrow = 1), matrix(dataRavg[i, ], nrow = 1))
      P1 <- dataLR^nHill / (Kh^nHill + dataLR^nHill)
      Weight[, , i] <- P1 %>% as.numeric()
    }
    dimnames(Weight)[[1]] <- levels(meta$labels)
    dimnames(Weight)[[2]] <- levels(meta$labels)
    dimnames(Weight)[[3]] <- lapply(1:nrow(LRDB), function(x) {
      paste0(LRDB[x, 1], "_", LRDB[x, 2])
    })
  }
  if (cellchat_output) {
    return(MLcellchat)
  } else {
    return(Weight)
  }
}
