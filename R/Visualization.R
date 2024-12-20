#' Visualize the LRI using the chord diagram
#'
#' @param result The LR pair probability extracted from the Extract_LR_Prob() function
#' @param topk The number of LR pairs that are visualized
#' @import ggplot2
#' @import circlize
#' @return A chord diagram of cell-cell communication
#' @export
PlotCCI_ChordPlot <- function(result, topk = 10) {
  requireNamespace("circlize")
  requireNamespace("grDevices")
  result.sorted <- result[order(result$Weight, decreasing = TRUE), ]
  topk <- min(topk, dim(result.sorted)[1])
  result.sorted <- result.sorted[1:topk, ]
  colnames(result.sorted)[3] <- "value"
  circos.clear()
  circos.par(start.degree = 85)
  chordp <- chordDiagram(result.sorted,
    directional = 1, direction.type = c("arrows"),
    big.gap = 30, order = c(rev(result.sorted$From), rev(result.sorted$To))
  )
  return(chordp)
}

#' Visualize the LRI using the circlized diagram
#'
#' @param result The LR pair probability extracted from the Extract_LR_Prob() function
#' @param topk The number of LR pairs that are visualized
#' @return A circle diagram of cell-cell communication
#' @export
PlotCCI_CirclePlot <- function(result, topk = 10) {
  vertex.label.color <- "black"
  vertex.weight <- 20
  vertex.label.cex <- 1
  arrow.width <- 2
  arrow.size <- 0.4
  alpha.edge <- 0.6

  result.sorted <- result[order(result$Weight, decreasing = TRUE), ]
  topk <- min(topk, dim(result)[1])
  result.sorted <- result.sorted[1:topk, ]
  result.sorted$Weight <- result.sorted$Weight / max(result.sorted$Weight)

  g <- graph_from_data_frame(result.sorted)
  igraph::E(g)$weight <- 2*result.sorted$Weight
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  layout <- in_circle()
  coords <- layout_(g, layout)
  coords_scale <- scale(coords)

  vertex.weight <- 5
  edge.alpha <- 0.6
  color.use <- grDevices::rainbow(length(igraph::V(g)))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  edge.weight.max <- 4 * max(result.sorted$Weight)
  igraph::E(g)$width <- 0.3 + edge.weight.max * igraph::E(g)$weight
  igraph::E(g)$label.cex <- 0.8
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::V(g)$color[edge.start[, 1]], edge.alpha)
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))

  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)), direction = -1, start = 0)
  label.dist <- 2
  edge.curved <- 0.2
  plot(g,
    edge.curved = edge.curved, vertex.shape = "circle", layout = coords_scale, margin = 0.2, vertex.label.dist = label.dist,
    vertex.label.degree = label.locs, vertex.label.family = "Times", edge.label.family = "Times"
  )
  gg <- grDevices::recordPlot()
  return(gg)
}

#' Show how many xtalk modules in the dataset
#'
#' @param CC_results The dataframe of PRS values.
#' @param KeyGenes The genes to be considered. By default, will consider all the genes in the pathways
#' @param data_type The type of the KeyGene (SSC or target).
#' @param top_percent In the inner panel, will show the the genes with most crosstalk pathways
#' @return Two diagrams. Outer: the distribution of number of crosstalk pathways; Inner: the top-k targets
#' @export
#'
PlotXT_Counts <- function(CC_results, KeyGenes = NULL, data_type = "Target", top_percent = 10) {

  Counts_pathway <- Count_Crosstalk(CC_results, KeyGenes, verbose = F, data_type = data_type)
  No_SSCs <- length(Counts_pathway)
  Counts_sorted <- sort(Counts_pathway)
  results <- data.frame(orders = 1:No_SSCs, gene = names(Counts_sorted), pathways = Counts_sorted)
  results$pathways <- as.numeric(results$pathways)
  results$pathways[results$pathways == 0] <- 1
  title_name <- paste0("Number of crosstalk pathways for ", data_type, "s")
  # windowsFonts(A = windowsFont("Arial"),T = windowsFont("Times New Roman"))

  #  p <- ggplot(results, aes(x=orders, y=pathways)) +
  #    geom_point(size = 2, aes(colour = factor(ifhub))) +
  #    ggrepel::geom_text_repel(
  #      data=results %>% filter(pathways >= Counts_sorted[No_SSCs-topk]), # Filter data first
  #      aes(label=gene),box.padding = 0.5,min.segment.length = 0.5,
  #      point.padding = 0.8,hjust = 1,vjust = 1,
  #      size = 6
  #    )+
  #    labs(x = data_type, y = "Pathways") +
  #    theme(axis.title = element_text(size = 24))+
  #    theme(axis.text = element_text(size = 18))+
  #    theme(legend.position = "none")

  p1 <- ggplot(results, aes(x = pathways)) +
    geom_histogram(binwidth = 1, fill = "#69b3a2", color = "black", alpha = 0.9) +
    theme_bw() +
    labs(x = "Number of crosstalk pathways", y = "Frequency") +
    theme(text = element_text(family = "Arial")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.position = "none")

  threshold <- quantile(Counts_sorted, 1 - top_percent / 100)
  if (threshold == max(Counts_sorted)) {
    results_topk <- dplyr::filter(results, pathways >= threshold)[, 2:3]
  } else {
    results_topk <- dplyr::filter(results, pathways > threshold)[, 2:3]
  }
  results_topk$colororder <- match(results_topk$pathways, results_topk$pathways %>% unique())
  allcolors <- grDevices::rainbow(results_topk$pathways %>% unique() %>% length())
  results_topk$colors <- allcolors[results_topk$colororder]
  p2 <- results_topk %>%
    mutate(gene = forcats::fct_reorder(gene, pathways)) %>%
    ggplot(aes(x = gene, y = pathways, fill = colors)) +
    geom_bar(stat = "identity", alpha = .6, width = .4) +
    coord_flip() +
    theme(text = element_text(family = "Arial")) +
    xlab("Target gene") +
    ylab("Number of pathways") +
    theme_bw() +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.position = "none")
  return(list(outer = p1, inner = p2))
}

#' The heatmap showing the TRS of signal-target pairs
#'
#' @param CC_pair_results The dataframe of TRS values.
#' @param Exp_clu The expression matrix of the target type.
#' @param KeyTG The genes to be considered.
#' @param topk The number of signal-target pairs to be plotted.
#' @return Two diagrams. Upper: the contribution of each signal on the expression of each target; Lower: the heatmap of the TRS values.
#' @export
#'
PlotXT_RecTGHeatmap <- function(CC_pair_results, Exp_clu, KeyTG, topk = 25) {

  CC_used <- CC_pair_results %>% dplyr::filter(Target %in% KeyTG)
  CC_used <- CC_used[order(CC_used$Weight,decreasing = T),]

  topk <- min(topk, nrow(CC_used))
  CC_used <- CC_used[1:topk, ]
  KeyTG <- CC_used$Target %>% unique() %>% sort()
  CC_mat <- df_mat(CC_used, row = Receptor, col = Target, value = Weight)
  CC_mat[is.na(CC_mat)] <- 0

  CC_df <- mat_df(CC_mat)
  colnames(CC_df) <- c("Receptor", "Target", "Weight")

  Ave_exp <- Calculate_Average_Expression(targene = KeyTG, Expmat = Exp_clu, nonzero = F)
  CC_normal <- apply(CC_mat, 2, function(x) x / sum(x)) %>% mat_df()
  colnames(CC_normal) <- c("Receptor", "Target", "Weight")
  CC_normal$Exp <- (CC_normal$Weight * (Ave_exp[CC_normal$Target])) %>% as.numeric()
  # windowsFonts(A = windowsFont("Arial"), T = windowsFont("Times New Roman"))
  requireNamespace("extrafont")
  requireNamespace("patchwork")

  p1 <- ggplot(CC_normal, aes(x = Target, weight = Exp, fill = Receptor)) +
    geom_bar(position = "stack") +
    labs(x = NULL, y = "Expression", fill = "Signal") +
    theme_minimal() +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_text(size = 14, face = "bold")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = 20, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 18, face = "bold")) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.position = "top")

  p2 <- ggplot(CC_df, aes(x = Target, y = Receptor, fill = Weight)) +
    geom_tile(color = "black", size = 1.5) +
    coord_equal() +
    scale_fill_gradient(low = "whitesmoke", high = "red") +
    labs(x = "Target gene", y = "Signal", fill = "Activity") +
    theme_minimal() +
    theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent")) +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 14, angle = 90, face = "bold", vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 14, face = "bold")) +
    theme(axis.title.x = element_text(size = 20, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 20, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 18, face = "bold")) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.position = "right")

  return(p1 / p2)
}

#' The river plot showing how the regulation flows
#'
#' @param CC_results The dataframe of PRS values.
#' @param KeyTG The target genes to be visualized.
#' @param min_weight The lower bound for the PRS of the pathways to be visualized
#' @return An alluvial diagram of crosstalk
#' @export
#'
PlotXT_Alluvial <- function(CC_results, KeyTG, min_weight = 0.45) {
  requireNamespace("ggalluvial")

  CC_used <- CC_results %>% dplyr::filter(Target %in% KeyTG)
  threshold <- quantile(CC_used$Weight, min_weight)
  CC_used2 <- CC_used %>% dplyr::filter(Weight > threshold)
  # windowsFonts(A = windowsFont("Arial"),T = windowsFont("Times New Roman"))

  p <- ggplot(
    data = CC_used2,
    aes(
      axis1 = Target, # First variable on the X-axis
      axis2 = SSC, # Second variable on the X-axis
      axis3 = Receptor, # Third variable on the X-axis
      y = Weight
    )
  ) +
    coord_flip() +
    ggalluvial::geom_alluvium(aes(fill = Target), width = 0.2, reverse = T) +
    ggalluvial::geom_stratum(alpha = .3, width = 0.2) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)), cex = 6
    ) +
    theme_void() +
    theme(text = element_text(family = "Arial")) +
    theme(legend.position = "none")
  return(p)
}

#' The heatmap of fidelity and specificity
#' @param CC_results The dataframe of PRS values.
#' @param KeyTG The target gene to be visualized.
#' @param threshold The lower bound for the fidelity of the pathways to be visualized
#' @return Two diagrams. Left: the fidelity matrix; Right: the specificity matrix.
#' @export
#'
PlotXT_FidSpe <- function(CC_results, KeyTG, threshold = 0.15) {

  if (!KeyTG %in% CC_results$Target) stop("Must input a target gene! \n")
  if (length(KeyTG) != 1) stop("Error: input more than 1 genes. \n")
  Fid_mat <- Calculate_Fidelity_Matrix(CC_results, KeyTG)
  Fid_mat <- Fid_mat$all
  if(nrow(Fid_mat) == 1){ # only one receptor
    Fid_df <- data.frame(From = rownames(Fid_mat), To = colnames(Fid_mat), Weight = Fid_mat %>% as.vector())
  }else if(ncol(Fid_mat) == 1){ # only one SSC
    Fid_df <- data.frame(From = rownames(Fid_mat), To = colnames(Fid_mat), Weight = Fid_mat %>% as.vector())
  }else{
    Fid_filtered <- Fid_mat[, colSums(Fid_mat) >= threshold] %>% as.matrix()
    Fid_df <- tidyfst::mat_df(Fid_filtered)
  }

  colnames(Fid_df) <- c("Receptor", "SSC", "Fidelity")

  Fid_df$Specificity <- apply(Fid_df, 1, function(x) {
    Spe_x <- Calculate_Pathway_Specificity(CC_results, KeyRec = x[1], KeySSC = x[2])
    Spe_x[KeyTG]
  })
  Fid_df[is.na(Fid_df$Specificity), "Specificity"] <- 0

  titlename <- paste0("Crosstalk for Target ", KeyTG)

  # windowsFonts(A = windowsFont("Arial"),T = windowsFont("Times New Roman"))
  p1 <- ggplot(Fid_df, aes(x = SSC, y = Receptor, fill = Fidelity)) +
    geom_tile(color = "black", size = 2) +
    coord_equal() +
    scale_fill_continuous(low = "whitesmoke", high = "green4", breaks = c(0, max(Fid_df$Fidelity)), labels = c("0.0", as.character(signif(max(Fid_df$Fidelity), digits = 2)))) +
    labs(x = "SSC", y = "Signal", fill = "Fidelity") +
    theme_minimal() +
    theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent")) +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 20, angle = 90, face = "bold", vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 20, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 18, face = "bold", vjust = 1)) +
    theme(legend.text = element_text(size = 16)) +
    theme(legend.key.width = unit(2, "cm")) +
    theme(legend.position = "top")
  p2 <- ggplot(Fid_df, aes(x = SSC, y = Receptor, fill = Specificity)) +
    geom_tile(color = "black", size = 2) +
    coord_equal() +
    scale_fill_continuous(low = "whitesmoke", high = "mediumorchid", breaks = c(0, max(Fid_df$Specificity)), labels = c("0.0", as.character(signif(max(Fid_df$Specificity), digits = 2)))) +
    labs(x = "SSC", y = "Signal", fill = "Specificity") +
    theme_minimal() +
    theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent")) +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 20, angle = 90, face = "bold", vjust = 0.5)) +
    # theme(axis.text.y = element_text(size = 20,face = "bold")) +
    theme(axis.text.y = element_blank()) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    # theme(axis.title.y = element_text(size = 24, face = "bold",vjust = 1) ) +
    theme(axis.title.y = element_blank()) +
    theme(legend.title = element_text(size = 18, face = "bold", vjust = 1)) +
    theme(legend.text = element_text(size = 16)) +
    theme(legend.key.width = unit(2, "cm")) +
    theme(legend.position = "top")
  p <- p1 + p2
  return(p)
}

#' The heatmap of PRS values for a given gene that is shared
#' @param CC_results The dataframe of PRS values.
#' @param gene_used The gene that is shared by pathways.
#' @param genetype The type of the shared gene: Receptor, SSC or Target
#' @param topk The maximum number of genes to be plotted.
#' @references Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics 2016.
#' @return A heatmap of the PRS or the normalized PRS.Side bars contain the distribution of fidelity or specificity values.
#' @export
#'
PlotXT_HeatMap <- function(CC_results, gene_used, genetype, topk = 25){
  requireNamespace("ComplexHeatmap")
  requireNamespace("grid")
  if(genetype == "Target" | genetype == "TG"){
    results_TG <- dplyr::filter(CC_results, Target == gene_used)
    results_TG <- results_TG[,c('Receptor','SSC','Weight')]
    temp_mat <- df_mat(results_TG, row = Receptor, col = SSC, value = Weight)
    legend_name <- "Fidelity"
    temp_mat[is.na(temp_mat)] <- 0
    temp_mat <- temp_mat/sum(temp_mat)
    temp_mat2 <- temp_mat[order(apply(temp_mat/sum(temp_mat),1,mean),decreasing = T),order(apply((temp_mat)/sum(temp_mat),2,mean),decreasing = T)]
    topk_used <- min(topk,ncol(temp_mat2))
    temp_mat2 <- temp_mat2[,1:topk_used]
    topk_used <- min(topk,nrow(temp_mat2))
    temp_mat2 <- temp_mat2[1:topk_used, ]
    # row (receptor) annotation
    ra =  ComplexHeatmap::rowAnnotation(Fid = ComplexHeatmap::anno_boxplot((temp_mat2/sum(temp_mat2)), height = unit(4, "cm")))
    # column (SSC) annotation
    ha =  ComplexHeatmap::HeatmapAnnotation(Fid = ComplexHeatmap::anno_boxplot((temp_mat2/sum(temp_mat2)), which='column',height = unit(4, "cm")))

    p <- ComplexHeatmap::Heatmap(temp_mat2,right_annotation = ra, top_annotation = ha,
                          cluster_rows = F, cluster_columns = F,name =legend_name,
                          row_names_gp = grid::gpar(fontsize = 18),
                          column_names_gp = grid::gpar(fontsize = 18),
                          width = unit(18, "cm"), height = unit(18, "cm"))
  }else if(genetype == "SSC" | genetype == "TF"){
    results_TG <- dplyr::filter(CC_results, SSC == gene_used)
    results_TG <- results_TG[,c('Receptor','Target','Weight')]
    temp_mat <- df_mat(results_TG, row = Receptor, col = Target, value = Weight)
    legend_name <- "PRS"
    temp_mat[is.na(temp_mat)] <- 0
    # temp_mat <- temp_mat/sum(temp_mat) # Do not normalize as we plot the PRS
    temp_mat2 <- temp_mat[order(apply(temp_mat/sum(temp_mat),1,mean),decreasing = T),order(apply((temp_mat)/sum(temp_mat),2,mean),decreasing = T)]
    topk_used <- min(topk,ncol(temp_mat2))
    temp_mat2 <- temp_mat2[,1:topk_used]
    topk_used <- min(topk,nrow(temp_mat2))
    temp_mat2 <- temp_mat2[1:topk_used, ]
    # row (receptor) annotation
    ra = ComplexHeatmap::rowAnnotation(Fid = ComplexHeatmap::anno_boxplot((temp_mat2/sum(temp_mat2)), height = unit(4, "cm")))
    # column (target) annotation
    ha = ComplexHeatmap::HeatmapAnnotation(Spe = ComplexHeatmap::anno_boxplot((temp_mat2/sum(temp_mat2)), which='column',height = unit(4, "cm")))

    p <- ComplexHeatmap::Heatmap(temp_mat2,right_annotation = ra, top_annotation = ha,
                          cluster_rows = F, cluster_columns = F,name =legend_name,
                          row_names_gp = grid::gpar(fontsize = 18),
                          column_names_gp = grid::gpar(fontsize = 18),
                          width = unit(18, "cm"), height = unit(18, "cm"))
  }else if(genetype == "Receptor" | genetype == "Rec"){
    results_TG <- dplyr::filter(CC_results, Receptor == gene_used)
    results_TG <- results_TG[,c('SSC','Target','Weight')]
    temp_mat <- df_mat(results_TG, row = SSC, col = Target, value = Weight)
    legend_name <- "Specificity"
    temp_mat[is.na(temp_mat)] <- 0
    temp_mat <- temp_mat/sum(temp_mat)
    temp_mat2 <- temp_mat[order(apply(temp_mat/sum(temp_mat),1,mean),decreasing = T),order(apply((temp_mat)/sum(temp_mat),2,mean),decreasing = T)]
    topk_used <- min(topk,ncol(temp_mat2))
    temp_mat2 <- temp_mat2[,1:topk_used]
    topk_used <- min(topk,nrow(temp_mat2))
    temp_mat2 <- temp_mat2[1:topk_used, ]
    # row (SSC) annotation
    ra = ComplexHeatmap::rowAnnotation(Spe = ComplexHeatmap::anno_boxplot((temp_mat2/sum(temp_mat2)), height = unit(4, "cm")))
    # column (SSC) annotation
    ha = ComplexHeatmap::HeatmapAnnotation(Spe = ComplexHeatmap::anno_boxplot((temp_mat2/sum(temp_mat2)), which='column',height = unit(4, "cm")))

    p <- ComplexHeatmap::Heatmap(temp_mat2,right_annotation = ra, top_annotation = ha,
                          cluster_rows = F, cluster_columns = F,name =legend_name,
                          row_names_gp = grid::gpar(fontsize = 18),
                          column_names_gp = grid::gpar(fontsize = 18),
                          width = unit(18, "cm"), height = unit(18, "cm"))
  }else{
    stop("Must specify a type of the key gene! \n")
  }

  return(p)
}

#' The ridgeline plot of the distribution of fidelity
#' @param CC_results The dataframe of PRS values.
#' @param KeyTG The target gene to be visualized.
#' @return A ridgeline plot of the distribution of fidelity for the signals that regulate the target via all the possible SSCs.
#' @export
#'
PlotXT_Ridgeline <- function(CC_results, KeyTG){

  requireNamespace("ggridges")
  Fid_mat <- Calculate_Fidelity_Matrix(CC_results, KeyTG)$all
  df <- mat_df(Fid_mat)
  colnames(df) <- c("Receptor","SSC","value")

  p <- df %>%
    mutate(Receptor = forcats::fct_reorder(Receptor, value)) %>%
    ggplot(aes(y=Receptor, x=value,  fill=Receptor)) +
      ggridges::geom_density_ridges(alpha=0.6) +
      ggridges::theme_ridges() +
      theme(
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
      ) +
      xlab("Fidelity") +
      ylab("Signal") +
      theme(text = element_text(family = "Arial")) +
      theme(axis.text.x = element_text(size = 20, angle = 90, face = "bold", vjust = 0.5)) +
      theme(axis.text.y = element_text(size = 20, face = "bold")) +
      theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
      theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
      theme(legend.title = element_text(size = 18, face = "bold", vjust = 1)) +
      theme(legend.text = element_text(size = 16))

  return(p)
}

#' The circulized bar charts of the specificity of Top-k specific targets for different conditions or signals
#' @param df The specificity dataframe
#' @param KeyFactors The conditions or signals. By default will use the second column of df.
#' @param topk The maximum number of targets to be plotted.
#' @param label_max The maximum specificity values shown in the charts.
#' @return A circulized bar charts of the specificity of Top-k specific targets for different conditions or signals.
#' @export
#'
PlotXT_MultiCircularBar <- function(df, KeyFactors = NULL, topk = 5, label_max = NULL) {

  colnames(df) <- c("individual", "group", "Specificity")
  if(is.null(KeyFactors)){
    KeyFactors <- df$group %>% as.array() %>% unique()
  }else KeyFactors <- intersect(KeyFactors,df$group)
  df <- dplyr::filter(df, group %in% KeyFactors)

  data <- c()
  for (type in KeyFactors) {
    temp_df <- dplyr::filter(df, group == type)
    temp_df <- temp_df[order(temp_df$Specificity, decreasing = T),]
    data <- rbind(data, temp_df[1:topk, ])
  }
  # Create dataset
  colnames(data) <- c("individual", "group", "value")
  data$group <- as.factor(data$group)
  rownames(data) <- NULL

  empty_bar <- 3
  to_add <- data.frame(matrix(NA, empty_bar * nlevels(data$group), ncol(data)))
  colnames(to_add) <- colnames(data)
  to_add$group <- rep(levels(data$group), each = empty_bar)
  data <- rbind(data, to_add)
  data <- data %>% arrange(group)
  data$id <- seq(1, nrow(data))

  label_data <- data
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle + 180, angle)

  base_data <- data %>%
    group_by(group) %>%
    summarize(start = min(id), end = max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title = mean(c(start, end)))

  grid_data <- base_data
  grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data) - 1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1, ]
  # windowsFonts(A = windowsFont("Arial"), T = windowsFont("Times New Roman"))

  if(is.null(label_max)){
    label_max <- max(data$value, na.rm = T)*1.25
    label_max <- signif(label_max,2)
  }

  p <- ggplot(data, aes(x = as.factor(id), y = value, fill = group)) + # Note that id is a factor. If x is numeric, there is some space between the first bar

    geom_bar(aes(x = as.factor(id), y = value, fill = group), stat = "identity", alpha = 0.5) +
    geom_segment(data = base_data, aes(x = start, y = label_max * 0.8, xend = end, yend = label_max * 0.8), colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE) +
    geom_segment(data = base_data, aes(x = start, y = label_max * 0.6, xend = end, yend = label_max * 0.6), colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE) +
    geom_segment(data = base_data, aes(x = start, y = label_max * 0.4, xend = end, yend = label_max * 0.4), colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE) +
    geom_segment(data = base_data, aes(x = start, y = label_max * 0.2, xend = end, yend = label_max * 0.2), colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE) +
    # Add text showing the value of each 100/75/50/25 lines
    annotate("text",
      x = rep(max(data$id), 4), y = c(label_max * 0.2, label_max * 0.4, label_max * 0.6, label_max * 0.8),
      label = c(paste0("Spe=", label_max * 0.2), paste0("Spe=", label_max * 0.4), paste0("Spe=", label_max * 0.6), paste0("Spe=", label_max * 0.8)),
      color = "black", size = 7.5, angle = 0, fontface = "bold", hjust = 1
    ) +
    geom_bar(aes(x = as.factor(id), y = value, fill = group), stat = "identity", alpha = 0.5) +
    ylim(-0.8 * label_max, 1.1 * label_max) +
    theme_void() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1, 6), "cm")
    ) +
    coord_polar() +
    geom_text(data = label_data, aes(x = id, y = label_max, label = individual, hjust = hjust), color = "black", fontface = "bold", alpha = 0.6, size = 7.5, angle = label_data$angle, inherit.aes = FALSE) +
    # Add base line information
    geom_segment(data = base_data, aes(x = start, y = -0.15 * label_max, xend = end, yend = -0.15 * label_max), colour = "black", alpha = 0.8, size = 0.6, inherit.aes = FALSE) +
    geom_text(data = base_data, aes(x = title, y = -0.5 * label_max, label = group), hjust = 0.5, colour = "black", alpha = 1, size = 7.5, fontface = "bold", inherit.aes = FALSE)
  return(p)
}

#' The chord diagram of PRS/fidelity/specificity values for a given gene that is shared
#' @param mat The matrix of PRS/fidelity/specificity values. Direction: columns to rows
#' @param orders The order of the genes. By default, orders <- c(rownames(mat), colnames(mat))
#' @param edge_colors The colors of the edges.
#' @return A chord diagram of PRS/fidelity/specificity values for a given gene that is shared.
#' @export
#'
PlotXT_Chord <- function(mat, orders = NULL,edge_colors = NULL){
  if(is.null(orders)){
    orders <- c(rownames(mat), colnames(mat))
  }
  circos.clear()
  if(is.null(edge_colors)){
    chordDiagram(t(mat),
      transparency = 0.25,
      order = orders,
      big.gap = 30,
      annotationTrack = "grid",  # Show sector labels
      annotationTrackHeight = 0.05,
      preAllocateTracks = 1    # Reserve space for labels
    )
  }else{
    chordDiagram(t(mat),
      col = edge_colors,
      transparency = 0.25,
      order = orders,
      big.gap = 30,
      annotationTrack = "grid",  # Show sector labels
      annotationTrackHeight = 0.05,
      preAllocateTracks = 1    # Reserve space for labels
    )
  }

  circos.track(track.index = 1, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      xplot = get.cell.meta.data("xplot")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")

      if(abs(xplot[2] - xplot[1]) < 15) {
          circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), col = "blue",cex = 2.5)
      } else {
          circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
              niceFacing = TRUE, adj = c(0.5, 0), col= "blue",cex = 2.5)
      }
  }, bg.border = NA)
}
