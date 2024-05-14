# Count the number of crosstalk pathways

Count_Crosstalk <- function(CC_results, KeyGenes = NULL, verbose = T, type = "SSC") {
  if (type == "SSC") {
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
  } else if (type == "TG") {
    if (is.null(KeyGenes)) KeyGenes <- CC_results$TG %>% unique()
    KeyTG <- KeyGenes
    TG_used <- intersect(KeyTG, CC_results$TG) %>% unique()
    if (!length(TG_used)) {
      stop("The input TGs have no crosstalk!")
    }
    Count_Xmodule <- array(0, length(TG_used))
    for (i in 1:length(TG_used)) {
      tg <- TG_used[i]
      temp_df <- CC_results[CC_results$TG == tg, ]
      nXTalk <- nrow(temp_df)
      if (nXTalk > 1) {
        if (verbose) {
          message(paste0("There are", nXTalk, " crosstalk pathways for TG:", tg, "!\n"))
        }
        Count_Xmodule[i] <- nXTalk
      } else {
        if (verbose) {
          message(paste0("There is no crosstalk for TG:", tg, "!\n"))
        }
        Count_Xmodule[i] <- 0
      }
    }
    names(Count_Xmodule) <- TG_used
  }

  return(Count_Xmodule)
}


# Calculate the signal->target Causality using the copula sum

Aggregate_Causality <- function(CC_results, sum_type, data_type = "TG") {
  if (data_type == "TG") {
    CC_df <- CC_results[, c("signal", "TG", "Weight_all")]
    CC_df$Weight_all <- 0.999 * CC_df$Weight_all / max(CC_df$Weight_all) # Normalize
    CC_aggregated <- aggregate(CC_df$Weight_all, list(CC_df$signal, CC_df$TG), function(x) Nonlinear_Sum(x, type = sum_type))
    colnames(CC_aggregated) <- c("signal", "TG", "Weight")
    CC_sorted <- CC_aggregated[order(CC_aggregated$signal, -CC_aggregated$Weight), ]
  } else if (data_type == "SSC") {
    CC_df <- CC_results[, c("signal", "SSC", "Weight_all")]
    CC_df$Weight_all <- 0.999 * CC_df$Weight_all / max(CC_df$Weight_all) # Normalize
    CC_aggregated <- aggregate(CC_df$Weight_all, list(CC_df$signal, CC_df$SSC), function(x) Nonlinear_Sum(x, type = sum_type))
    colnames(CC_aggregated) <- c("signal", "SSC", "Weight")
    CC_sorted <- CC_aggregated[order(CC_aggregated$signal, -CC_aggregated$Weight), ]
  }


  return(CC_sorted)
}

Aggregate_SSCTG_Causality <- function(CC_results, sum_type) {
  CC_df <- CC_results[, c("SSC", "TG", "Weight_all")]
  CC_df$Weight_all <- 0.999 * CC_df$Weight_all / max(CC_df$Weight_all) # Normalize
  CC_aggregated <- aggregate(CC_df$Weight_all, list(CC_df$SSC, CC_df$TG), function(x) Nonlinear_Sum(x, type = "Copula"))
  colnames(CC_aggregated) <- c("SSC", "TG", "Weight")
  return(CC_aggregated)
}

# Obtain a FIDELITY matrix for a fixed target gene.

Calculate_Fidelity <- function(CC_results, KeyTG) {
  if (!KeyTG %in% CC_results$TG) {
    stop("The input gene must be a target gene!")
    return(NULL)
  }
  if (length(KeyTG) != 1) {
    stop("The input gene must be a target gene!")
    return(NULL)
  }
  CC_used <- CC_results[CC_results$TG %in% KeyTG, ]
  CC_mat <- df_mat(CC_used, row = signal, col = SSC, value = Weight_all)
  CC_mat[is.na(CC_mat)] <- 0
  Fid_SSC <- CC_mat / rowSums(CC_mat) # for a given signal-TG, the fidelity for each SSC
  Fid_signal <- CC_mat / colSums(CC_mat) # for a given SSC-TG, the fidelity for each signal
  Fid_all <- CC_mat / sum(CC_mat) # fidelity for each pathway signal-SSC-TG

  Fid_results <- list(
    SSC = Fid_SSC,
    signal = Fid_signal,
    all = Fid_all
  )
  return(Fid_results)
}

Calculate_Specificity <- function(CC_results, Keysignal) {
  Keysignal <- unique(Keysignal)
  if (!Keysignal %in% CC_results$signal) {
    message("The input gene must be a signal!")
    return(NULL)
  }
  CC_used <- CC_results[CC_results$signal %in% Keysignal, ]
  CC_mat <- df_mat(CC_used, row = TG, col = SSC, value = Weight_all)
  CC_mat[is.na(CC_mat)] <- 0

  Spe_SSC <- CC_mat / rowSums(CC_mat) # for a given signal-TG, the specificity for each SSC
  Spe_tg <- CC_mat / colSums(CC_mat) # for a given signal-SSC, the specificity for each tg
  Spe_all <- CC_mat / sum(CC_mat) # specificity for each pathway signal-SSC-TG

  Spe_results <- list(
    SSC = Spe_SSC,
    tg = Spe_tg,
    all = Spe_all
  )
  return(Spe_results)
}

# Calculate the L-R-T pathway specificity

Calculate_Pathway_Specificity <- function(CC_results, Keysignal, KeySSC) {
  if (length(Keysignal) + length(KeySSC) != 2) stop("Must input one signal and one SSC! \n")
  CC_used <- CC_results %>%
    filter(signal == Keysignal) %>%
    filter(SSC == KeySSC)
  Spe_pathway <- CC_used$Weight_all / sum(CC_used$Weight_all)
  names(Spe_pathway) <- CC_used$TG
  return(Spe_pathway)
}

# Compare the crosstalk patterns of multiple crosstalk results
Compare_Multitype_CC <- function(CC_pair_list, type = "Fid") {
  alltypes <- names(CC_pair_list)
  signal_list <- lapply(CC_pair_list, function(x) {
    x$signal %>% unique()
  })
  names(signal_list) <- alltypes

  TG_list <- lapply(CC_pair_list, function(x) {
    x$TG %>% unique()
  })
  names(TG_list) <- alltypes

  common_signals <- Reduce(intersect, signal_list)
  common_TGs <- Reduce(intersect, TG_list)

  if (type == "Fid") {
    Fid_all <- lapply(alltypes, function(x) {
      CC_pair_res <- CC_pair_list[[x]]
      CC_pair_used <- CC_pair_res %>% filter(TG %in% common_TGs)
      CC_pair_mat <- df_mat(CC_pair_used, row = signal, col = TG, value = Weight)
      CC_pair_mat[is.na(CC_pair_mat)] <- 0
      Fid_pair_mat <- CC_pair_mat / colSums(CC_pair_mat)
      Fid_pair_df <- mat_df(Fid_pair_mat)
      colnames(Fid_pair_df) <- c("signal", "TG", "Fidelity")
      Fid_pair_df$celltype <- x
      Fid_pair_df
    })
    Fid_df <- do.call(rbind, Fid_all)

    return(Fid_df)
  }
  if (type == "Spe") {
    Spe_all <- lapply(alltypes, function(x) {
      CC_pair_res <- CC_pair_list[[x]]
      CC_pair_used <- CC_pair_res %>% filter(TG %in% common_TGs)
      CC_pair_mat <- df_mat(CC_pair_used, row = signal, col = TG, value = Weight)
      CC_pair_mat[is.na(CC_pair_mat)] <- 0
      Spe_pair_mat <- CC_pair_mat / rowSums(CC_pair_mat)
      Spe_pair_df <- mat_df(Spe_pair_mat)
      colnames(Spe_pair_df) <- c("signal", "TG", "Specificity")
      Spe_pair_df$celltype <- x
      Spe_pair_df
    })
    Spe_df <- do.call(rbind, Spe_all)

    return(Spe_df)
  }
}

# Subcluster the cells using the expression of highly specific genes
XT_Subclustering <- function(SeuratObj, pair_results, Keysignals, target_type, topk = 10) {
  
  require(flexclust)

  Tar_weight <- pair_results %>% filter(signal %in% Keysignals)
  Tar_weight <- aggregate(Tar_weight$Weight, list(Tar_weight[,2]), sum)
  colnames(Tar_weight) <- c("Tar", "Weight")
  Tar_weight <- Tar_weight[order(Tar_weight$Weight, decreasing = T), ]
  topk <- min(topk, nrow(Tar_weight))
  Tar_weight <- Tar_weight[1:topk,]
  SeuratObj$celltype <- Idents(SeuratObj)
  subobject <- subset(SeuratObj, celltype == target_type)
  Exp_Tars <- subobject@assays$SCT$counts[Tar_weight$Tar, ]
  ws <- Tar_weight$Weight
  ws <- ws/sum(ws)

  cluster_results <- cclust(Exp_Tars %>% t(), k = 2, dist = "manhattan", method = "hardcl", weights = ws)
  cluster_center <- cluster_results@centers
  cluster_df <- data.frame(cell = colnames(Exp_Tars), clusterID = cluster_results@second, interacted = NA)
  norm1 <- sum(cluster_center[1, ])
  norm2 <- sum(cluster_center[2, ])

  if (norm1 < norm2) {
    cluster_df[which(cluster_df$clusterID == 1), 3] <- FALSE
    cluster_df[which(cluster_df$clusterID == 2), 3] <- TRUE
  } else {
    cluster_df[which(cluster_df$clusterID == 1), 3] <- TRUE
    cluster_df[which(cluster_df$clusterID == 2), 3] <- FALSE
  }

  return(cluster_df)
}
