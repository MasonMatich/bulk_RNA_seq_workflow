# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #                                                                   #   #
#   #             Visualize PCA by Experimental Condition               #   #
#   #                                                                   #   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

visualizePCsByCondition <- function(rlog, metadata, title) {

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #               Get PCs for Scree Plot                # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Compute PCA
  pca <- pca(assay(rlog), metadata = metadata)
  
  # Compute Total Variance and Remove 0s For Log Scale Plotting
  pca$variance <- ifelse(pca$variance == 0, NA, pca$variance)
  
  # Get Number of PCs Calculated for qplot()
  numPCs <- length(pca$variance)
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #             Organize Principle Components Into Df             # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Create DF With All Experiments, Strains, Replicates, and Associated PCs
  metaData <- pca$metadata %>%
    rownames_to_column("sample_id")
  
  metaData <- pivot_longer(
    metaData, 
    names(pca$metadata), 
    names_to = "grouping", 
    values_to = "ident"
  )
  
  # Find Elbow
  elb <- findElbowPoint(pca$variance)
  
  relevantPCs <- pca$rotated[,c(1:elb)] %>%
    rownames_to_column("sample_id")
  
  PCswMetadata <- left_join(metaData, relevantPCs)
  
  PCswMetadata <- pivot_longer(
    PCswMetadata, 
    starts_with("PC"), 
    names_to = "PCs"
  )
  
  PCswMetadata <- PCswMetadata %>%
    mutate(strain = word(sample_id, 1, sep = "_"))
  
  # Visualize PCs
  selectedPCs <- ggplot(PCswMetadata, aes(x = ident, y = value, colour = strain)) +
    geom_point() +
    stat_summary(
      fun = mean,
      geom = "line",
      aes(group = 1),
      linetype = "solid",
      show.legend = FALSE,
      colour = "grey40"
    ) +
    facet_grid(PCs ~ grouping, scales = "free") +
    ggtitle(sprintf("PCs %s", title))+
    theme(axis.text.x = element_text(angle = -90))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #                          Scree Plots                          # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  variance <- pca$variance
  numPCs <- length(pca$components)
  
  logScaleScree <- qplot(c(1:(numPCs -1)), variance[1:length(variance)-1]) + 
    geom_line() + 
    geom_point(size=2)+
    xlab("Principal Component") + 
    ylab("Variance Explained (Log10 Scale)") +
    ggtitle(sprintf("Log10 Scree Plot %s", title)) +
    scale_y_log10()
  
  linearScaleScree <- qplot(c(1:(numPCs -1)), variance[1:length(variance)-1]) + 
    geom_line() + 
    geom_point(size=2)+
    xlab("Principal Component") + 
    ylab("Variance Explained") +
    ggtitle(sprintf("Linear Scree Plot %s", title))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #                        Return All Plots                       # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  list(
    logScaleScree = logScaleScree, 
    linearScaleScree = linearScaleScree, 
    selectedPCs = selectedPCs
  )
}