#' Get the expression matrix from a seurat object
#'
#' @param SeuratObj The Seurat object of the dataset.
#' @param clusterID The cell type used (i.e. target_type).
#' @param assay The assay of the Seurat Object to use. By default use the "RNA" assay.
#' @param datatype The datatype of the assay. By default use the "scale.data".
#' @param cutoff The threshold for filtering out the low-expression genes.
#' @return The expression matrix of the certain cell type.
#' @import SeuratObject
#' @export
#'
Get_Exp_Clu <- function(SeuratObj, clusterID = NULL, assay = "RNA", datatype = "data", cutoff = 0.05) {

  if (is.null(clusterID)) {
    subobject <- SeuratObj
  } else {
    subobject <- subset(SeuratObj, idents = clusterID)
  }
  if("SCT" %in% names(SeuratObj@assays)){
    Count_clu <- subobject@assays$SCT$counts
  }else if("RNA" %in% names(SeuratObj@assays)){
    Count_clu <- subobject@assays$RNA$counts
  }else{
    Count_clu <- subobject@assays$Spatial$counts
  }

  allgenes <- rownames(Count_clu)
  ncells <- ncol(Count_clu)
  expressedgenes <- rownames(Count_clu)[which(rowSums(Count_clu > 0) > cutoff * ncells)]
  if(assay %in% names(subobject@assays)){
    if(datatype == "data"){
      Exp_clu <- subobject@assays[[assay]]$data
    }else if(datatype == "counts"){
      Exp_clu <- subobject@assays[[assay]]$counts
    }else{
      Exp_clu <- subobject@assays[[assay]]$scale.data
    }
  }else if ("RNA" %in% names(subobject@assays)){
    if(datatype == "data"){
      Exp_clu <- subobject@assays$RNA$data
    }else if(datatype == "counts"){
      Exp_clu <- subobject@assays$RNA$counts
    }else{
      Exp_clu <- subobject@assays$RNA$scale.data
    }
  }else if ("SCT" %in% names(subobject@assays)){
    if(datatype == "data"){
      Exp_clu <- subobject@assays$SCT$data
    }else if(datatype == "counts"){
      Exp_clu <- subobject@assays$SCT$counts
    }else{
      Exp_clu <- subobject@assays$SCT$scale.data
    }
  }else if ("Spatial" %in% names(subobject@assays)){
    if(datatype == "data"){
      Exp_clu <- subobject@assays$Spatial$data
    }else if(datatype == "counts"){
      Exp_clu <- subobject@assays$Spatial$counts
    }else{
      Exp_clu <- subobject@assays$Spatial$scale.data
    }
  }

  # filter low-expressed genes
  expressedgenes <- intersect(expressedgenes,rownames(Exp_clu))
  Exp_clu <- Exp_clu[expressedgenes, ]
  Exp_clu <- as(Exp_clu, "sparseMatrix")
  gc()

  return(Exp_clu)
}

#' Convert from a cellchat object an Ligand-Receptor probability dataframe
#'
#' @param result The CellChat object
#' @param source_type The sender cell type.
#' @param target_type The receiver cell type.
#' @param pv_threshold The p-value threshold for filtering out the insignificant LR pairs.
#' @return The expression matrix of the certain cell type.
#' @export
Extract_LR_Prob <- function(result, source_type = NULL, target_type = NULL, pv_threshold = 0.05) {
  if(is.null(source_type) & is.null(target_type)){
    stop("Must specify the source or the target cell type!\n")
  }

  probs <- result@net$prob
  pvals <- result@net$pval
  if (is.null(source_type)) { # consider all sender types
    LR_vector <- c()
    for (type in dimnames(probs)[[1]]) {
      LR_df <- data.frame(interaction = names(probs[type, target_type, ]),Weight = probs[type, target_type, ],pval = pvals[type, target_type, ])
      LR_df <- LR_df %>%
        dplyr::filter(Weight > 0) %>%
        dplyr::filter(pval < pv_threshold)
      if (!rlang::is_empty(LR_df)) LR_vector <- rbind(LR_vector, LR_df)
    }
  } else {
    if (!is.null(target_type)) {
      LR_df <- data.frame(interaction = names(probs[source_type, target_type, ]),Weight = probs[source_type, target_type, ],pval = pvals[source_type, target_type, ])
      LR_df <- LR_df %>%
        dplyr::filter(Weight > 0) %>%
        dplyr::filter(pval < pv_threshold)
      LR_vector <- LR_df
    } else {
      LR_vector <- c()
      for (type in dimnames(probs)[[2]]) {
        LR_df <- data.frame(interaction = names(probs[source_type, type, ]),Weight = probs[source_type, type, ],pval = pvals[source_type, type, ])
        LR_df <- LR_df %>%
          dplyr::filter(Weight > 0) %>%
          dplyr::filter(pval < pv_threshold)
        if (!rlang::is_empty(LR_df)) LR_vector <- rbind(LR_vector, LR_df)
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

#' Run the fisher's exact test
#'
#' @param subset1 The first set of genes.
#' @param subset2 The second set of genes.
#' @param background All the expressed genes
#' @return The p-value of the fisher test
#'
Fisher_Test <- function(subset1,subset2,background){
  subset1 <- intersect(subset1, background)
  subset2 <- intersect(subset2, background)
  a=intersect(subset1,subset2)
  la <- length(a)
  b=setdiff(subset2,a)
  lb <- length(b)
  c=setdiff(subset1,a)
  lc <- length(c)
  abc <- union(a,b)
  abc <- union(abc,c)
  d <- setdiff(background,abc)
  ld <- length(d)
  matrix=matrix(c(la,lc,lb,ld),nrow=2)
  fisher.test(matrix,alternative="greater")$p.value
}

#' Run the Chi-square's test
#' @param subset1 The first set of genes.
#' @param subset2 The second set of genes.
#' @param background All the expressed genes
#' @return The p-value of the test
Chisq_Test <- function(subset1,subset2,background){
  subset1 <- intersect(subset1, background)
  subset2 <- intersect(subset2, background)
  a=intersect(subset1,subset2)
  la <- length(a)
  b=setdiff(subset2,a)
  lb <- length(b)
  c=setdiff(subset1,a)
  lc <- length(c)
  abc <- union(a,b)
  abc <- union(abc,c)
  d <- setdiff(background,abc)
  ld <- length(d)
  matrix=matrix(c(la,lc,lb,ld),nrow=2)
  stats::chisq.test(matrix)$p.value
}

#' Perform the pagerank algorithm
#'
#' @param gg The graph object of the database
#' @param regulators The receptors, which are also the seed of the algorithm.
#' @param targets The destinations of the paths.
#' @references Identify, quantify and characterize cellular communication from single-cell RNA sequencing data with scSeqComm. Bioinformatics (2022)
#' @references Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nature Methods (2019)
#' @return The PPR values of the certain regulator-target pair.
#'
myPageRank <- function(gg, regulators, targets) {

  # Modified from scSeqComm and NicheNet
  # Convert dataframe to igraph object

  association <- data.frame(regulator = character(0), target = character(0), PPR = numeric(0))
  nod <- sapply(igraph::V(gg)$name, function(z) strsplit(strsplit(z,":")[[1]][2],",")) #nodes of graph and their gene components
  if(is.na(nod[[1]])[1]){
      nod <- sapply(igraph::V(gg)$name, function(z) strsplit(strsplit(z,":")[[1]][1],",")) #nodes of graph and their gene components
  }
  r <- nodes_find(regulators, gg, nod)
  tf <- nodes_find(targets, gg, nod)

  if (length(r) == 0 | length(tf) == 0) return(association) # pathway does not contains at least one R and one TF

  # for each receptor
  for (j in 1:length(r)){
    # node reachable by all nodes associated to the current receptor
    tf_connected <- lapply(seq_along(r[[j]]), function(i) {
      return(is_connected_to(r[[j]][i], gg))
    })

    E <- rep(0, times = length(igraph::V(gg))) # preference vector for pagerank algorithm

    # computation of PERSONALIZED PAGERANK setting receptor node as seed node
    ppr <- single_receptor_ppr_wrapper(r[[j]], E, gg, delta = 0.85)
    colnames(ppr) <- igraph::V(gg)$name

    # for each transcription factor
    for (k in 1:length(tf)) {

      # a TF can be contained in more than one node
      for (t in 1:length(tf[[k]])) {
        ans <- c()

        # if multiple receptor node (a receptorr R1 contained in more than on one): the Persornalized PageRank algoritm is computed taking all receptor nodes as seed nodes
        if (length(r[[j]]) > 1) {
          # which receptor node is connected to the current tf?
          recp <- unlist(lapply(tf_connected, function(i) tf[[k]][t] %in% i))

          if (any(recp)) {
            # recomputation of pageranks setting the receptor nodes connected to tf as seed nodes
            ppr <- single_receptor_ppr_wrapper(r[[j]][recp], E, gg, delta = 0.85)
            # if(any(ppr<0)) print(gg$name)

            colnames(ppr) <- igraph::V(gg)$name

            # result
            ans <- c(ans, ppr[1, tf[[k]]][t])
          }
        } else {
          # the receptor is contained in only one node

          if (tf[[k]][t] %in% unlist(tf_connected)) # tf node is reachable by receptor node?
            {
              # result
              ans <- c(ans, ppr[1, tf[[k]]][t])
            }
        }
      } # for multiple tf

      # save results

      if (length(ans) != 0) {
        association <- rbind(association, data.frame(regulator = names(r)[j], target = names(tf)[k], PPR = mean(ans)))
      }
    } # for tf
  } # for r

  association <- association[order(-association$PPR), ]

  return(association)
}


#' Do single time PPR algorithm
#' @param r The seed node.
#' @param E The edges of the graph.
#' @param graph The graph object.
#' @param delta The dampling factor.
#' @return The PPR matrix.
single_receptor_ppr_wrapper <- function(r,E,graph,delta){

  # set receptor as seed in preference vector E
  E[which(igraph::V(graph) %in% r,arr.ind = T)] <- 1

  #Personalized PageRank algorithm
  ppr <- igraph::page_rank(graph, algo = c("prpack"), vids = igraph::V(graph),directed = TRUE, damping = delta, personalized = E, weights = NULL) %>% .$vector
  ppr_matrix = matrix(unlist(ppr), ncol = length(E), byrow = TRUE)
  return(ppr_matrix)
}

#' Find the nodes that contains the certain gene
#' @param x The target gene.
#' @param gg The graph object.
#' @param nod The nodes.
#' @return All the nodes that contain gene x.
nodes_find <- function(x,gg,nod){

  ans <- list()
  for(j in 1:length(x)) {
    rec <- strsplit(x[[j]],",")[[1]]
    X <- c()

    #multi-subunit regulators
    for (r in rec)
    {
      #find nodes that contain current gene
      redok <- vapply(seq_along(nod),function(i) any(nod[[i]]==r), logical(1))

      if (any(redok))
      {
        X <- unique(c(X,igraph::V(gg)[igraph::V(gg)$name %in% names(nod[redok])]))  #ID node
      }
    }
    #save results
    if (length(X)>0)
    {
      t <- list(X)
      names(t) <- paste(x[j])
      ans <- c(ans,t)
    }
  }
  return(ans)
}

#' Find all nodes reachable by the given node.
#' @param r The seed node.
#' @param gg The graph object.
#' @return All the nodes that the seed node can reach via the graph.
is_connected_to <- function(r,gg)
{
  R_connection <- igraph::subcomponent(graph = gg, v = r, mode = "out") # NOTE: node "r" is reachable by itself, so it is included in the output
  R_connection <- sort(R_connection[R_connection!=r]) # remove node "r" and sort the vertex
  return(R_connection)
}


#' Calculate the different types of "sum"
#'
#' @param x A numerical array
#' @param type The type of "sum": normal sum, copula sum, probability sum
#' @return The "sum" of the data
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

#' Calculate the average expression of a gene
#'
#' @param targene The gene to be calculated
#' @param Expmat The expression matrix
#' @param nonzero Whether consider cells with 0 targene expression
#' @return The "sum" of the data
Calculate_Average_Expression <- function(targene, Expmat, nonzero = T){

  Ave_exp <- lapply(targene, function(gene) {
    if(!gene %in% rownames(Expmat)){
      0
    }else{
      Exp_gene <- Expmat[gene,]
      Trimean(Exp_gene, nonzero)
    }

  }) %>% as.data.frame()

  rownames(Ave_exp) <- "Exp"
  colnames(Ave_exp) <- targene
  return(Ave_exp)
}

#' Calculate the average value of a vector
#' @param x The vector.
#' @param nonzero Whether filter out non-zero values in the vector.
#' @return The tri-mean of the vector.
Trimean <- function(x,nonzero = T) {
  if(nonzero){
    x[x == 0] <- NA
  }
  mean(stats::quantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = T))
}

#' Filter the Signal-SSC-Target pathways with low PRS results
#'
#' @param results The raw PRS results generated by PRS_calc() function
#' @param PRS_thres The threshold for the PRS value. Pathways with PRS lower than PRS_thres*max(PRS) will be removed
#' @param remove_genes Boolean variables on whether remove mitochondrial and ribosomal genes
#' @return The filtered dataframe containing the PRS values of pathways
#' @export
#' @examples
#' \donttest{
#' RecTFTG_filtered <- Filter_results(ress, PRS_thres = 0.05)
#' }
Filter_results <- function(results, PRS_thres = 0.05, remove_genes = T){

  if(remove_genes){ # for human and mouse, separately
    results <- results[!grepl("^MT-", results$Target), ]
    results <- results[!grepl("^RPL", results$Target), ]
    results <- results[!grepl("^RPS", results$Target), ]
    results <- results[!grepl("^mt-", results$Target), ]
    results <- results[!grepl("^Rpl", results$Target), ]
    results <- results[!grepl("^Rps", results$Target), ]
  }
  # Filter out low-PRS pathways
  results_filtered <- results[results$Weight > PRS_thres * max(results$Weight), ]

  return(results_filtered)
}

#' Run a Python script with given conda environment and arguments
#'
#' @param script_path Character. Path to the Python script.
#' @param conda_env Character. The name of the conda environment.
#' @param args Character vector. Arguments to pass to the script (e.g., \code{c("--arg1", "value1","--arg2", "value2")}).
#'
#' @return No return value.
#' @export
#'
#' @examples
#' \dontrun{
#' Run_py_script(
#'   script_path = "./main.py",
#'   conda_env = "SigXTalk_py",
#'   args = c("--project",args.project,"--target_type",target_type)
#' )
#' }
Run_py_script <- function(script_path, conda_env, args = character()) {

  # Activate the conda environment
  reticulate::use_condaenv(conda_env, required = TRUE)

  # Construct full sys.argv with script name as the first element
  full_args <- c(basename(script_path), args)

  # Format as Python sys.argv assignment string
  sys_argv_code <- sprintf(
    "import sys; sys.argv = %s",
    paste0("['", paste(full_args, collapse = "', '"), "']")
  )

  # Set sys.argv and run the Python script
  reticulate::py_run_string(sys_argv_code)
  reticulate::py_run_file(script_path)

}
