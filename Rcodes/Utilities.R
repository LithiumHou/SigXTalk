# Get the expression matrix from a seurat object
Get_Exp_Clu <- function(SeuratObj, clusterID = NULL, assay = "SCT", datatype = "data", aggregated = F, cutoff = 0,
                        k_neigh = 50, atacbinary = TRUE, max_overlap = 0.8, imputation = F) {
  if (imputation) {
    require(SeuratWrappers)
    SeuratObj <- RunALRA(object = SeuratObj)
  }
  if (is.null(clusterID)) {
    subobject <- SeuratObj
  } else {
    subobject <- subset(SeuratObj, idents = clusterID)
  }
  if("SCT" %in% names(SeuratObj@assays)){
    Count_clu <- subobject@assays$SCT$counts
  }else{
    Count_clu <- subobject@assays$RNA$counts
  }
  
  allgenes <- rownames(Count_clu)
  ncells <- ncol(Count_clu)
  expressedgenes <- rownames(Count_clu)[which(rowSums(Count_clu > 0) > cutoff * ncells)]
  expressedgenes <- intersect(expressedgenes, allgenes)

  if (aggregated) {
    if ("wnn.umap" %in% names(subobject@reductions)) {
      cell_coord_i <- subobject@reductions$wnn.umap@cell.embeddings
    } else {
      cell_coord_i <- subobject@reductions$umap.rna@cell.embeddings
    }
    sub_aggregated_data <- Aggregation_Single(subobject, cell_coord_i, k_neigh, atacbinary, max_overlap)
    Exp_clu <- sub_aggregated_data$rna
  } else {
    if (imputation) {
      Exp_clu <- subobject@assays$alra@data
    } else {
      if (datatype == "counts"){
        Exp_clu <- subobject@assays$SCT$counts
      }else {
        if (assay == "SCT"){
          Exp_clu <- subobject@assays$SCT$data
        }else{
          Exp_clu <- subobject@assays$RNA$data
        }
        
      }
    }
  }

  # filter low-expressed genes
  expressedgenes <- intersect(expressedgenes,rownames(Exp_clu))
  Exp_clu <- Exp_clu[expressedgenes, ]
  Exp_clu <- as(Exp_clu, "sparseMatrix")
  gc()

  return(Exp_clu)
}

# Convert interaction matrices to From-To lists
Convert_Mat_To_Relationship <- function(input_mat, rm_zero = F) {
  tryCatch(
    {
      output_df <- tidyfst::mat_df(input_mat)
      colnames(output_df) <- c("From", "To", "Weight")
      output_df$Weight <- as.numeric(output_df$Weight)
      return(output_df)
    },
    error = function(e) {
      message("Warning: package 'tidyfst' is not installed, using other methods...\n")
      output_df <- c()
      from_size <- nrow(input_mat)
      to_size <- ncol(input_mat)
      for (i in 1:from_size) {
        tempfrom <- rownames(input_mat)[i]
        tempdata <- cbind(rep(tempfrom, to_size), colnames(input_mat), input_mat[i, ])
        tempdata <- as.data.frame(tempdata)
        colnames(tempdata) <- c("From", "To", "Weight")
        if (rm_zero) {
          tempdata <- tempdata[tempdata$Weight != 0, ]
        }
        output_df <- rbind(output_df, tempdata)
      }
      output_df$Weight <- as.numeric(output_df$Weight)
      rownames(output_df) <- NULL
      return(output_df)
    }
  )
}

# Convert from a cellchat object or network prob 3-dim tensor to an LR prob matrix
Extract_LR_Prob <- function(result, source_type = NULL, target_type = NULL, cellchat_use) {
  if(is.null(source_type) & is.null(target_type)){
    stop("Must input the source or the target!\n")
  }
  if (cellchat_use) {
    probs <- result@net$prob
    pvals <- result@net$pval
    if (is.null(source_type)) { # consider all sender types
      LR_vector <- c()
      for (type in dimnames(probs)[[1]]) {
        LR_df <- data.frame(interaction = names(probs[type, target_type, ]),Weight = probs[type, target_type, ],pval = pvals[type, target_type, ])
        LR_df <- LR_df %>%
          dplyr::filter(Weight > 0) %>%
          filter(pval < 0.1)
        if (!rlang::is_empty(LR_df)) LR_vector <- rbind(LR_vector, LR_df)
      }
    } else {
      if (!is.null(target_type)) {
        LR_df <- data.frame(interaction = names(probs[source_type, target_type, ]),Weight = probs[source_type, target_type, ],pval = pvals[source_type, target_type, ])
        LR_df <- LR_df %>%
          dplyr::filter(Weight > 0) %>%
          filter(pval < 0.1)
        LR_vector <- LR_df
      } else {
        LR_vector <- c()
        for (type in dimnames(probs)[[2]]) {
          LR_df <- data.frame(interaction = names(probs[source_type, type, ]),Weight = probs[source_type, type, ],pval = pvals[source_type, type, ])
          LR_df <- LR_df %>%
            dplyr::filter(Weight > 0) %>%
            filter(pval < 0.1)
          if (!rlang::is_empty(LR_df)) LR_vector <- rbind(LR_vector, LR_df)
        }
      }
    }
  } else {
    if (is.null(source_type)) { # consider all sender types
      LR_vector <- c()
      probs <- result
      for (type in dimnames(probs)[[1]]) {
        LR_df <- data.frame(interaction = names(probs[type, target_type, ]),Weight = probs[type, target_type, ])
        if (!rlang::is_empty(LR_df)) LR_vector <- rbind(LR_vector, LR_df)
      }
    } else { #consider all receiver types
      LR_vector <- c()
      probs <- result
      if(!is.null(target_type)){
        LR_vector <- data.frame(interaction = names(probs[source_type, target_type, ]),Weight = probs[source_type, target_type, ])
        LR_vector <- LR_vector %>%
            dplyr::filter(Weight > 0)
      }else{
        for (type in dimnames(probs)[[2]]) {
          LR_df <- data.frame(interaction = names(probs[source_type, type, ]),Weight = probs[source_type, type, ],pval = pvals[source_type, type, ])
          LR_df <- LR_df %>%
            dplyr::filter(Weight > 0) %>%
            filter(pval < 0.1)
          if (!rlang::is_empty(LR_df)) LR_vector <- rbind(LR_vector, LR_df)
        }
      }
      
    }
  }
  LR_vector$Weight <- as.numeric(LR_vector$Weight)
  Pairprob_final <- c()

  for (i in 1:nrow(LR_vector)) {
    tempname <- strsplit(LR_vector$interaction[i], "_") %>% unlist()
    if (length(tempname) == 3) {
      temp_frame1 <- data.frame("From" = tempname[1], "To" = tempname[2], "Weight" = LR_vector$Weight[i] / 2)
      temp_frame2 <- data.frame("From" = tempname[1], "To" = tempname[3], "Weight" = LR_vector$Weight[i] / 2)
      Pairprob_final <- rbind(Pairprob_final, temp_frame1, temp_frame2)
    } else {
      temp_frame1 <- data.frame("From" = tempname[1], "To" = tempname[2], "Weight" = LR_vector$Weight[i])
      Pairprob_final <- rbind(Pairprob_final, temp_frame1)
    }
  }
  Pairprob_final <- as.data.frame(Pairprob_final)
  Pairprob_final$Weight <- as.numeric(Pairprob_final$Weight)
  Pairprob_sum <- aggregate(Pairprob_final$Weight, list(Pairprob_final$From, Pairprob_final$To), sum)
  colnames(Pairprob_sum) <- c("From", "To", "Weight")
  return(Pairprob_sum)
}

# Calculate the distance matrix between cells
Calculate_Distance <- function(cellnames, cell_coords, p = 1){
  cellid <- cellnames[which(cellnames %in% rownames(cell_coords))]
  cell_coords_use <- cell_coords[cellid,]
  # Use the Manhattan distance
  distance_cells <- dist(cell_coords_use, method = "euclidean", diag = 1, upper = 1, p) %>% as.matrix()
  rownames(distance_cells) <- cellnames
  colnames(distance_cells) <- cellnames
  return(distance_cells)
}

# Calculate the copula sum
Nonlinear_Sum <- function(x, type = "Copula"){
  # input should be a numerical vector
  # sum type: copula, prob, sum
  if(!is.numeric(x)){
    stop("The input should be a numerical vector!\n")
  }
  x_len <- length(x)
  if(x_len < 2){
    if(x_len == 1){
      return(x)
    }else{
      stop("The input should at least contain 1 number!\n ")
    }
  }

  if(x_len == 2){
    if(type == "Copula"){
      return((x[1]+x[2]-2*x[1]*x[2])/(1-x[1]*x[2]))
    }
    else{
      return(1-(1-x[1])*(1-x[2]))
    }
  }else{
    if(type == "Copula"){
      Csum <- (x[1]+x[2]-2*x[1]*x[2])/(1-x[1]*x[2])
      for(i in 2:(x_len-1)){
        temp <- (Csum+x[i+1]-2*Csum*x[i+1])/(1-Csum*x[i+1])
        Csum <- temp
      }
    }else if(type == "prob"){
      Csum <- (1-(1-x[1])*(1-x[2]))
      for(i in 2:(x_len-1)){
        temp <- (1-(1-Csum)*(1-x[i+1]))
        Csum <- temp
      }
    }else if(type == "sum"){
      Csum <- sum(x)
    }

    return(Csum)
  }
}
