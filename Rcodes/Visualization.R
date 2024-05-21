###### Visualization

# Prepare the color sets
scPalette <- function(n) {
  colorSpace <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", "#222F75", "#1B9E77", "#B2DF8A",
    "#E3BE00", "#FB9A99", "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", "#8DD3C7", "#999999"
  )
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}

# Visualize the LRI
PlotCCI_ChordPlot <- function(result, topk = 10) {
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
  color.use <- scPalette(length(igraph::V(g)))
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
    vertex.label.degree = label.locs, vertex.label.family = "Helvetica", edge.label.family = "Arial"
  )
  gg <- recordPlot()
  return(gg)
}

# Show how many xtalk modules in the dataset

PlotXT_Counts <- function(CC_results, KeyGenes = NULL, datatype = "TG", top_percent = 10) {
  Counts_pathway <- Count_Crosstalk(CC_results, KeyGenes, verbose = F, type = datatype)
  No_SSCs <- length(Counts_pathway)
  Counts_sorted <- sort(Counts_pathway)
  results <- data.frame(orders = 1:No_SSCs, gene = names(Counts_sorted), pathways = Counts_sorted)
  results$pathways <- as.numeric(results$pathways)
  results$pathways[results$pathways == 0] <- 1
  title_name <- paste0("Number of crosstalk pathways for ", datatype, "s")
  # windowsFonts(A = windowsFont("Arial"),T = windowsFont("Times New Roman"))

  #  p <- ggplot(results, aes(x=orders, y=pathways)) +
  #    geom_point(size = 2, aes(colour = factor(ifhub))) +
  #    ggrepel::geom_text_repel(
  #      data=results %>% filter(pathways >= Counts_sorted[No_SSCs-topk]), # Filter data first
  #      aes(label=gene),box.padding = 0.5,min.segment.length = 0.5,
  #      point.padding = 0.8,hjust = 1,vjust = 1,
  #      size = 6
  #    )+
  #    labs(x = datatype, y = "Pathways") +
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
    results_topk <- filter(results, pathways >= threshold)[, 2:3]
  } else {
    results_topk <- filter(results, pathways > threshold)[, 2:3]
  }
  results_topk$colororder <- match(results_topk$pathways, results_topk$pathways %>% unique())
  results_topk$colors <- brewer.pal(results_topk$pathways %>% unique() %>% length(), "Paired")[results_topk$colororder]
  p2 <- results_topk %>%
    mutate(gene = fct_reorder(gene, pathways)) %>%
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

PlotXT_RecTGHeatmap <- function(CC_pair_results, Exp_clu, topk = 25, ranktype = "expression") {
  if (ranktype == "weight") {
    CC_TGrank <- aggregate(CC_pair_results$Weight, list(CC_pair_results$TG), FUN = function(x) Nonlinear_Sum(x, type = "Copula"))
    colnames(CC_TGrank) <- c("TG", "Weight")
    CC_TGrank <- CC_TGrank[order(-CC_TGrank$Weight), ]
    TG_used <- CC_TGrank[1:topk, ]$TG
  } else if (ranktype == "expression") {
    allTGs <- CC_pair_results$TG %>% unique()
    Exp_TGs <- Exp_clu[allTGs, ] %>% rowSums()
    TG_used <- sort(Exp_TGs, decreasing = T) %>% names()
  } else if (ranktype == "number") {
    Counts_pathway <- Count_Crosstalk(CC_results, KeyGenes = NULL, verbose = F, type = "TG")
    TG_used <- sort(Counts_pathway, decreasing = T) %>% names()
  }

  topk <- min(topk, length(TG_used))
  TG_used <- TG_used[1:topk]
  CC_used <- CC_pair_results %>% filter(TG %in% TG_used)
  CC_mat <- df_mat(CC_used, row = Signal, col = TG, value = Weight)
  CC_mat[is.na(CC_mat)] <- 0

  CC_df <- mat_df(CC_mat)
  colnames(CC_df) <- c("Signal", "TG", "Weight")

  Ave_exp <- Calculate_Average_Expression(targene = sort(TG_used), Expmat = Exp_clu, nonzero = F)
  CC_normal <- apply(CC_mat, 2, function(x) x / sum(x)) %>% mat_df()
  colnames(CC_normal) <- c("Signal", "TG", "Weight")
  CC_normal$Exp <- (CC_normal$Weight * (Ave_exp[CC_normal$TG])) %>% as.numeric()
  # windowsFonts(A = windowsFont("Arial"), T = windowsFont("Times New Roman"))
  require(extrafont)

  p1 <- ggplot(CC_normal, aes(x = TG, weight = Exp, fill = Signal)) +
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

  p2 <- ggplot(CC_df, aes(x = TG, y = Signal, fill = Weight)) +
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

  require(patchwork)
  return(p1 / p2)
}

# The river plot shows how the causality flows
PlotXT_Alluvial <- function(CC_results, KeyTG, topk = 8, ranktype = "weight", min_weight = 0.45) {
  require(ggalluvial)
  if (!is.null(KeyTG)) {
    TG_used <- KeyTG
  } else {
    if (ranktype == "weight") {
      CC_pair_results <- Aggregate_Causality(CC_results, sum_type = "Copula")
      CC_TGrank <- aggregate(CC_pair_results$Weight, list(CC_pair_results$TG), FUN = function(x) Nonlinear_Sum(x, type = "Copula"))
      colnames(CC_TGrank) <- c("TG", "Weight")
      CC_TGrank <- CC_TGrank[order(-CC_TGrank$Weight), ]
      TG_used <- CC_TGrank$TG
    } else if (ranktype == "expression") {
      CC_pair_results <- Aggregate_Causality(CC_results, sum_type = "Copula")
      allTGs <- CC_pair_results$TG %>% unique()
      Exp_TGs <- Exp_clu[allTGs, ] %>% rowSums()
      TG_used <- sort(Exp_TGs, decreasing = T) %>% names()
    } else if (ranktype == "number") {
      Counts_pathway <- Count_Crosstalk(CC_results, KeyGenes = NULL, verbose = F, type = "TG")
      TG_used <- sort(Counts_pathway, decreasing = T) %>% names()
    }

    topk <- min(topk, length(TG_used))
    TG_used <- TG_used[1:topk]
  }

  CC_used <- CC_results %>% filter(TG %in% TG_used)
  threshold <- quantile(CC_used$Weight_all, min_weight)
  CC_used2 <- CC_used %>% filter(Weight_all > threshold)
  # windowsFonts(A = windowsFont("Arial"),T = windowsFont("Times New Roman"))

  p <- ggplot(
    data = CC_used2,
    aes(
      axis1 = TG, # First variable on the X-axis
      axis2 = SSC, # Second variable on the X-axis
      axis3 = Signal, # Third variable on the X-axis
      y = Weight_all
    )
  ) +
    coord_flip() +
    geom_alluvium(aes(fill = TG), width = 0.2, reverse = T) +
    geom_stratum(alpha = .3, width = 0.2) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)), cex = 6
    ) +
    theme_void() +
    theme(text = element_text(family = "Arial")) +
    theme(legend.position = "none")
  return(p)
}

# The heatmap of fidelity and specificity
PlotXT_FidSpe <- function(CC_results, KeyTG,threshold = 0.15) {
  if (!KeyTG %in% CC_results$TG) stop("Must input a target gene! \n")
  if (length(KeyTG) != 1) stop("Error: input more than 1 genes. \n")
  Fid_mat <- Calculate_Fidelity(CC_results, KeyTG)
  Fid_mat <- Fid_mat$SSC
  if(nrow(Fid_mat) == 1){ # only one receptor
    Fid_df <- data.frame(From = rownames(Fid_mat), To = colnames(Fid_mat), Weight = Fid_mat %>% as.vector())
  }else if(ncol(Fid_mat) == 1){ # only one SSC
    Fid_df <- data.frame(From = rownames(Fid_mat), To = colnames(Fid_mat), Weight = Fid_mat %>% as.vector())
  }else{
    Fid_filtered <- Fid_mat[, colSums(Fid_mat) >= threshold] %>% as.matrix()
    Fid_df <- Convert_Mat_To_Relationship(Fid_filtered)
  }

  colnames(Fid_df) <- c("Signal", "SSC", "Fidelity")

  Fid_df$Specificity <- apply(Fid_df, 1, function(x) {
    Spe_x <- Calculate_Pathway_Specificity(CC_results, KeySignal = x[1], KeySSC = x[2])
    Spe_x[KeyTG]
  })
  Fid_df[is.na(Fid_df$Specificity), "Specificity"] <- 0

  titlename <- paste0("Crosstalk for TG ", KeyTG)

  # windowsFonts(A = windowsFont("Arial"),T = windowsFont("Times New Roman"))
  p1 <- ggplot(Fid_df, aes(x = SSC, y = Signal, fill = Fidelity)) +
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
  p2 <- ggplot(Fid_df, aes(x = SSC, y = Signal, fill = Specificity)) +
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

PlotXT_SingleCircularBar <- function(CC_pair_results, LR_Pairprob, topk = 5, label_max = NULL) {
  require(tidyverse)
  LR_Pairprob_used <- LR_Pairprob %>% filter(To %in% CC_pair_results$Signal)
  Signal_top <- aggregate(LR_Pairprob_used$Weight, list(LR_Pairprob_used$To), sum)

  SignalAll <- Signal_top$Group.1[order(Signal_top$x, decreasing = T)][1:topk]
  data <- c()
  for (signal in SignalAll) {
    temp_df <- CC_pair_results %>% filter(Signal == signal)
    temp_df$Specificity <- temp_df$Weight
    temp_df <- temp_df[order(temp_df$Specificity, decreasing = T), c(2, 1, 4)]
    data <- rbind(data, temp_df[1:5, ])
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
  # windowsFonts(A = windowsFont("Arial"),T = windowsFont("Times New Roman"))

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


PlotXT_MultiCircularBar <- function(Spe_all_df, Keyfactor, topk = 5, label_max = NULL) {
  require(tidyverse)

  if (!Keyfactor %in% colnames(Spe_all_df)) stop("Please select a key factor (celltype, SSC, etc)... \n")
  if(colnames(Spe_all_df)[1] == Keyfactor){
    Spe_all_df[,c(1,2)] <- Spe_all_df[,c(2,1)] # Make sure key factor is on the 2nd column
    colnames(Spe_all_df)[c(1,2)] <- colnames(Spe_all_df)[c(2,1)]
  } 
  colnames(Spe_all_df)[3] <- "Specificity"
  allfactors <- Spe_all_df[,Keyfactor] %>% unique() %>% as.character()
  data <- c()
  for (type in allfactors) {
    temp_df <- Spe_all_df[which(Spe_all_df[,Keyfactor] == type),]
    temp_df$Specificity <- temp_df$Specificity
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

PlotXT_ComparisonBubble <- function(Fid_all_df, SeuratObj, KeyTG) {
  Fid_all_df$Signal <- as.character(Fid_all_df$Signal)
  Fid_all_df$TG <- as.character(Fid_all_df$TG)
  Fid_all_mat <- df_mat(Fid_all_df %>% as.data.frame(), row = Signal, col = celltype, value = Fidelity)
  Fid_all_mat[Fid_all_mat %>% is.na()] <- 0
  Fid_all_mat <- Fid_all_mat / colSums(Fid_all_mat)
  Fid_df_used <- mat_df(Fid_all_mat)
  colnames(Fid_df_used) <- c("Signal", "celltype", "value")

  # calculate the average expressions
  alltypes <- Fid_df_used$celltype %>% unique()
  Exp_clus <- lapply(alltypes, function(x) {
    subobject <- subset(SeuratObj, idents = x)
    subobject@assays$RNA$counts
  })
  names(Exp_clus) <- alltypes
  Exps <- apply(Fid_df_used, 1, function(x) {
    Calculate_Average_Expression(
      targene = x[1][[1]] %>% as.character(),
      Expmat = Exp_clus[[x[2][[1]] %>% as.character()]],
      nonzero = F
    )[[1]] %>% as.numeric()
  })
  Fid_df_used$Exp <- Exps

  # windowsFonts(A = windowsFont("Arial"), T = windowsFont("Times New Roman"))

  p <- ggplot(Fid_df_used, aes(x = celltype, y = Signal)) +
    geom_point(aes(size = Exp, color = value)) +
    theme_bw() +
    scale_color_gradient(low = "lightgray", high = "blue") +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    labs(x = "", y = "Signal", color = "Fidelity", size = "Average \nExpression") +
    # theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.1)) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 18, face = "bold", vjust = 1)) +
    theme(legend.text = element_text(size = 16))
  return(p)
}

PlotXT_KEGG <- function(genes, type = "human") {
  if (type == "human") {
    id_list <- mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL") %>% na.omit()
    kegg <- enrichKEGG(
      gene = id_list,
      organism = "hsa", keyType = "kegg",
      pAdjustMethod = "BH", pvalueCutoff = 0.5, qvalueCutoff = 0.8,
      minGSSize = 10, maxGSSize = 500, use_internal_data = F
    )
    kegg <- data.frame(kegg)
    kegg_used <- data.frame(fun = kegg$Description, value = kegg$Count)
    kegg_used <- kegg_used %>%
      mutate(fun = fct_reorder(fun, value))

    ylabel <- paste0("Number of genes (total=", length(genes), ")")
    p <- ggplot(kegg_used, aes(x = fun, y = value, fill = "red")) +
      geom_bar(stat = "identity", width = 0.9) +
      coord_flip() +
      theme_bw() +
      ylim(0, 1.1 * max(kegg_used$value %>% ceiling())) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
      labs(x = "KEGG terms", y = ylabel) +
      theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.1)) +
      theme(axis.text.y = element_text(size = 16, face = "bold")) +
      theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
      theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
      theme(legend.position = "none") +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  }
  return(p)
}

PlotXT_Subcluster <- function(SeuratObj, Subcluster_df, source_type, target_type) {

  SeuratObj$newID <- SeuratObj$celltype %>% as.character()
  true_cells <- Subcluster_df[which(Subcluster_df$interacted), "cell"] %>% as.vector()
  false_cells <- Subcluster_df[which(!Subcluster_df$interacted), "cell"] %>% as.vector()
  senders <- names(SeuratObj$celltype[which(SeuratObj$celltype == source_type)])
  others <- setdiff(names(SeuratObj$newID), c(senders, true_cells, false_cells))
  SeuratObj$newID[true_cells] <- paste0(target_type, "-interacted")
  SeuratObj$newID[false_cells] <- paste0(target_type, "-not interacted")
  SeuratObj$newID[senders] <- source_type
  SeuratObj$newID[others] <- "Other cells"
  object_used <- SeuratObj
  Idents(object_used) <- object_used$newID
  colors <- c("green","blue","red","gray")
  names(colors) <- c(source_type,paste0(target_type, "-interacted"),paste0(target_type, "-not interacted"),"Other cells")
  p <- SpatialDimPlot(object_used, label = F, cols = colors)
  return(p)
}

PlotXT_distances <- function(SeuratObj, Subcluster_df, source_type, target_type) {
  true_cells <- Subcluster_df[which(Subcluster_df$interacted), "cell"] %>% as.vector()
  false_cells <- Subcluster_df[which(!Subcluster_df$interacted), "cell"] %>% as.vector()
  senders <- names(SeuratObj$celltype[which(SeuratObj$celltype == source_type)])
  cell_coords <- Seurat::GetTissueCoordinates(SeuratObj, scale = NULL, cols = c("imagerow", "imagecol"))
  dist_mat <- Calculate_Distance(c(true_cells, false_cells, senders), cell_coords, p = 1)

  true_min <- apply(dist_mat[true_cells, senders], 1, min)
  true_mean <- rowMeans(dist_mat[true_cells, senders])
  false_min <- apply(dist_mat[false_cells, senders], 1, min)
  false_mean <- rowMeans(dist_mat[false_cells, senders])

  df_min <- rbind(
    data.frame(type = "interacted", distance = true_min),
    data.frame(type = "non-interacted", distance = false_min)
  )

  p <- ggplot(df_min, aes(x = type, y = distance, fill = type)) +
    geom_boxplot(alpha = 0.3) +
    theme(legend.position = "none") +
    labs(x = "", y = "Distance to senders (units)") +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.1)) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1))

  return(p)
}

### Plot for supp figures
PlotXT_Counts_SSC <- function(CC_results, KeyGenes = NULL, datatype = "SSC", top_percent = 10) {
  Counts_df <- aggregate(CC_results$Weight_all, list(CC_results$SSC), length)
  colnames(Counts_df) <- c("SSC","Pathways")
  Counts_df <- mutate(Counts_df, SSC = fct_reorder(SSC, Pathways))
  ggplot(Counts_df, aes(x=SSC, y=Pathways)) +
    geom_point(size=6, color="#69b3a2") + # Show dots
    geom_text(
      label=Counts_df$SSC, 
      nudge_x = .25, nudge_y = 15.75, 
      size = 4, color = "red",
      check_overlap = T
    )+
    theme_bw()+
    labs(x = "SSC", y = "Number of Pathways") +
    scale_x_discrete(label = element_blank())+
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.1)) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1))
}
