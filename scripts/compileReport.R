# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #                                                                     #   #
#   #         Calls translateBioIDs.R and geneOntologyAnalysis.R          #   #
#   #                                                                     #   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

compileReport <- function(result.obj, contrast, sample_name = NULL, l2fc_filter, padj_filter, semantic_data, OrgDb){
  
  # Create comment structure for bookdown render
  cat(sprintf("## %s {.contrast}\n\n", contrast))
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #   #                                                                   #   #
  #   #                            DEG Report                             #   # 
  #   #                                                                   #   #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  cat("\n\n### DEG Report")
  cat("\n\n")
  knitDataTable(
    df = list(result.obj$DataFrame),
    tableName = "Differentially Expressed Genes",
    sample_name = sample_name,
    contrast = contrast,
    fileExtension = "DEGs",
    pageLength = 10
  )
  
  # # Volcano Plot
  volcano_file <- sprintf("%s_%s_volcano.png", sample_name, gsub(" ", "_", tolower(contrast)))
  volcano_path <- file.path(project_dirs$figures, volcano_file)

  # Relative path for HTML display
  volcano_rel_path <- file.path("figures", volcano_file)

  volcano <- EnhancedVolcano(
    result.obj$DataFrame,
    lab = result.obj$DataFrame$gene_name,
    x = 'log2FoldChange',
    y = 'padj',
    title = contrast,
    subtitle = sprintf('padj < %s, |log2FC| > %s', padj_filter, l2fc_filter),
    pCutoff = padj_filter,
    FCcutoff = l2fc_filter
  )

  ggsave(volcano_path, volcano, height = 12, width = 10, units = "in", dpi = 300)
  cat(sprintf('<img src="%s" width="80%%"/>\n\n', volcano_rel_path))
  
  # # MA Plot
  plotMA(
    result.obj$DESeqDataSet,
    ylim=c(-3,3), 
    main = sprintf("%s  MA plot", contrast)
    )
  cat("\n\n")
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #   #                                                                   #   #
  #   #                         Ontology Report                           #   #
  #   #                                                                   #   #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  cat("\n\n### Ontology Report {.tabset .analysisreport}\n\n")
  
  # Translate Wormbase IDs to ENTREZID
  result.obj.list <- translateBioIDs(
    DESeqResults.obj = result.obj$DESeqDataSet,
    bioID = "WORMBASE",
    l2fc_filter = l2fc_filter,
    padj_filter = padj_filter
  )
  
  # Run gruopGO(), enrichGO(), gseGO(), and goplot()
  go.result.list <- callClusterProfilerFunc(
    result.obj.list = result.obj.list,
    OrgDb = OrgDb,
    semantic_data = semantic_data)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # #             Dysregulated Reactome Pathway GSE, GO GSE             # # # 
  cat(paste0("\n\n#### Dysregulated\n\n"))
  
  # GO GSE
  # Open the GO GSEA div container (for styling)
  cat('<div class="analysis">\n\n')
  cat("\n\n##### Gene Ontology GSEA\n\n")
  
  if (any(!is.na(go.result.list$go_gse_red[[1]])) && nrow(go.result.list$go_gse_red[[1]]) > 0) {
    treemapPlot(go.result.list$go_gse_red[[1]])
  }
  
  knitDataTable(
    df = go.result.list$go_gse,
    tableName = " GSE Of All Differentially Expressed Genes",
    sample_name = sample_name,
    contrast = contrast,
    fileExtension = "GO_gsea"
  )

  # # # GSEA plot condiitonal rendering
  if (nrow(go.result.list$go_gse[[1]]) > 0){
    # Tabset for GSEA plots
    cat("##### {.tabset .tabset-dropdown}\n\n")

    # Print GSEA plot for each set identified
    for (i in 1:nrow(data.frame(go.result.list$go_gse[[1]]))) {
      if (i > 25) {break} # Dont print over 30 plots
      cat(sprintf("###### %s \n\n", go.result.list$go_gse[[1]]$Description[i]))

      print(gseaplot(
        go.result.list$go_gse[[1]],
        by = "all",
        title = go.result.list$go_gse[[1]]$Description[i],
        geneSetID = go.result.list$go_gse[[1]]$ID[i]))

      cat("\n\n")
    }
  }
  cat('\n\n</div>\n\n')
  
  
  # Reactome Pathway GSEA
  # Open the Reactome Pathway GSEA div container (for styling)
  cat('<div class="analysis">\n\n')
  
  cat("\n\n##### Reactome Pathway GSEA\n\n")
  
  knitDataTable(
    df = go.result.list$react_gse,
    tableName = " Reactome Pathway GSE Of All Differentially Expressed Genes",
    sample_name = sample_name,
    contrast = contrast,
    fileExtension = "Reactome_gsea"
  )
  
  # # # GSEA plot condiitonal rendering of RP GSEA plots
  if (nrow(go.result.list$react_gse[[1]]) > 0){
    # Tabset for Reactome GSEA plots
    cat("##### {.tabset .tabset-dropdown}\n\n")

    # Print GSEA plot for each react pathway set identified
    for (i in 1:nrow(data.frame(go.result.list$react_gse[[1]]))) {
      if (i > 25) {break} # Dont print over 30 plots
      cat(sprintf("###### %s \n\n", go.result.list$react_gse[[1]]$Description[i]))

      print(gseaplot(
        go.result.list$react_gse[[1]],
        by = "all",
        title = go.result.list$react_gse[[1]]$Description[i],
        geneSetID = go.result.list$react_gse[[1]]$ID[i]))

      cat("\n\n")
    }
  }
  
  # Close the Reactome Pathway GSEA div container (for styling)
  cat('\n\n</div>\n\n')
  

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # #    Upregulated GO classification, GO over-representation, RP over-representation    # # #
  cat(paste0("\n\n#### Upregulated\n\n"))
  
  # GO classification
  cat("\n\n##### GO Classification {.analysis}\n\n")
  knitDataTable(
    df = go.result.list$go_class_upreg,
    tableName = " GO Classification of Upregulated Genes",
    sample_name = sample_name,
    contrast = contrast,
    fileExtension = "GO_classification_upreg"
  )
  
  # GO over-representation
  cat("\n\n##### GO Over-representation {.analysis}\n\n")
  # # Treemap Plot of simplified terms
  if (any(!is.na(go.result.list$go_overrep_upreg_red[[1]])) && nrow(go.result.list$go_overrep_upreg_red[[1]]) > 0) {
    treemapPlot(go.result.list$go_overrep_upreg_red[[1]])
  }
  
  # # Table of simplified terms
  knitDataTable(
    df = go.result.list$go_overrep_upreg,
    tableName = " GO Over-representation Analysis of Upregulated Genes",
    sample_name = sample_name,
    contrast = contrast,
    fileExtension = "GO_overrepresented_upreg"
  )
  
  # RP over-representation
  cat("\n\n##### Reactome Pathway Over-representation {.analysis}\n\n")
  
  knitDataTable(
    df = go.result.list$react_overrep_upreg,
    tableName = " Reactome Pathway Over-representation Analysis of Upregulated Genes",
    sample_name = sample_name,
    contrast = contrast,
    fileExtension = "Reactome_overrepresented_upreg"
  )
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # #    Downregulated GO classification, GO over-representation, RP over-representation    # # #
  cat(paste0("\n\n#### Downregulated\n\n"))
  
  # GO classification
  cat("\n\n##### GO Classification {.analysis}\n\n")
  
  knitDataTable(
    df = go.result.list$go_class_downreg,
    tableName = " GO Classification of Downregulated Genes",
    sample_name = sample_name,
    contrast = contrast,
    fileExtension = "GO_classification_downreg"
  )

  # GO over-representation
  cat("\n\n##### GO Over-representation {.analysis}\n\n")
  
  # # Treemap Plot of simplified terms
  if (any(!is.na(go.result.list$go_overrep_downreg_red[[1]])) && nrow(go.result.list$go_overrep_downreg_red[[1]]) > 0) {
    treemapPlot(go.result.list$go_overrep_downreg_red[[1]])
    cat("\n\n")
  }
  
  # # Table of simplified terms
  knitDataTable(
    df = go.result.list$go_overrep_downreg,
    tableName = " GO Over-representation Analysis of Downregulated Genes",
    sample_name = sample_name,
    contrast = contrast,
    fileExtension = "GO_overrepresented_downreg"
  )
  
  # RP over-representation
  cat("\n\n##### Reactome Pathway Over-representation {.analysis}\n\n")
  
  knitDataTable(
    df = go.result.list$react_overrep_downreg,
    tableName = " Reactome Pathway Over-representation Analysis of Downregulated Genes",
    sample_name = sample_name,
    contrast = contrast,
    fileExtension = "Reactome_overrepresented_downreg"
  )

}