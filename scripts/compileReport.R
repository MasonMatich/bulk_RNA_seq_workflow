# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# #         Calls translateBioIDs.R and geneOntologyAnalysis.R        # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

compileReport <- function(result.obj, contrast, l2fc_filter, padj_filter){
  
  # Create comment structure for bookdown render
  cat(sprintf("# %s {.title}\n\n", contrast))

  # Show exact call which invoked function
  call_txt <- paste(deparse(match.call()), collapse = "\n")
  cat("**Inputs for function call:**\n\n")
  cat("```r\n")
  cat(call_txt, "\n\n")
  cat("```\n\n")
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #   #                                                                   #   #
  #   #                            DEG Report                             #   # 
  #   #                                                                   #   #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  cat("\n\n## DEG Report {.subtitle}")
  cat("\n\n")
  knitDataTable(
    df = list(result.obj$DataFrame),
    tableName = "Differentially Expressed Genes",
    contrast = contrast,
    fileExtension = "DEGs",
    pageLength = 10
  )
  
  # # Volcano Plot
  volcano <- EnhancedVolcano(
    result.obj$DataFrame,
    lab = result.obj$DataFrame[,1],
    x = 'log2FoldChange',
    y = 'pvalue'
  )
  print(volcano)
  cat("\n\n")
  
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
  
  cat("\n\n## Ontology Report {.tabset .subtitle}")
  
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
    OrgDb = org.Ce.eg.db)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # #             Dysregulated Reactome Pathway GSE, GO GSE             # # # 
  cat(paste0("\n\n### Dysregulated\n\n"))
  
  # Reactome Pathway GSE
  cat("**Reactome Pathway GSEA**")
  knitDataTable(
    df = go.result.list$react_gse,
    tableName = " Reactome Gene Set Enrichment Analysis",
    contrast = contrast,
    fileExtension = "Reactome_gsea"
  )
  
  # # # GSEA plot condiitonal rendering of RP GSEA plots
  # Tabset for Reactome GSEA plots
  cat("#### {.tabset .tabset-dropdown}\n\n")

  # Print GSEA plot for each react pathway set identified 
  for (i in 1:nrow(data.frame(go.result.list$react_gse[[1]]))) {
    if (i > 30) {break} # Dont print over 30 plots
    cat(sprintf("##### %s \n\n", go.result.list$react_gse[[1]]$Description[i]))
    print(gseaplot(
      go.result.list$react_gse[[1]], 
      by = "all", 
      title = go.result.list$react_gse[[1]]$Description[i], 
      geneSetID = go.result.list$react_gse[[1]]$ID[i]))
    cat("\n\n")
  }
  
  # GO GSE
  cat("#### \n\n") # Break out of previous tabset
  cat("**GO GSEA**")
  knitDataTable(
    df = go.result.list$go_gse,
    tableName = " GSE Of All Differentially Expressed Genes",
    contrast = contrast,
    fileExtension = "GO_gsea"
  )
  
  # # # GSEA plot condiitonal rendering
  # Tabset for GSEA plots
  cat("#### {.tabset .tabset-dropdown}\n\n")

  # Print GSEA plot for each set identified 
  for (i in 1:nrow(data.frame(go.result.list$go_gse[[1]]))) {
    if (i > 30) {break} # Dont print over 30 plots
    cat(sprintf("##### %s \n\n", go.result.list$go_gse[[1]]$Description[i]))
    print(gseaplot(
      go.result.list$go_gse[[1]], 
      by = "all", 
      title = go.result.list$go_gse[[1]]$Description[i], 
      geneSetID = go.result.list$go_gse[[1]]$ID[i]))
    cat("\n\n")
  }

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # #    Upregulated GO classification, GO over-representation, RP over-representation    # # #
  cat(paste0("\n\n### Upregulated\n\n"))
  
  # GO classification
  knitDataTable(
    df = go.result.list$go_class_upreg,
    tableName = " GO Classification of Upregulated Genes",
    contrast = contrast,
    fileExtension = "GO_classification_upreg"
  )
  
  # GO over-representation
  knitDataTable(
    df = go.result.list$go_overrep_upreg,
    tableName = " GO Over-representation Analysis of Upregulated Genes",
    contrast = contrast,
    fileExtension = "GO_overrepresented_upreg"
  )
  
  # RP over-representation
  knitDataTable(
    df = go.result.list$react_overrep_upreg,
    tableName = " Reactome Pathway Over-representation Analysis of Upregulated Genes",
    contrast = contrast,
    fileExtension = "Reactome_overrepresented_upreg"
  )
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # #    Downregulated GO classification, GO over-representation, RP over-representation    # # #
  cat(paste0("\n\n### Downregulated\n\n"))
  
  # GO classification
  knitDataTable(
    df = go.result.list$go_class_downreg,
    tableName = " GO Classification of Downregulated Genes",
    contrast = contrast,
    fileExtension = "GO_classification_downreg"
  )

  # GO over-representation
  knitDataTable(
    df = go.result.list$go_overrep_downreg,
    tableName = " GO Over-representation Analysis of Downregulated Genes",
    contrast = contrast,
    fileExtension = "GO_overrepresented_downreg"
  )
  
  # RP over-representation
  knitDataTable(
    df = go.result.list$react_overrep_downreg,
    tableName = " Reactome Pathway Over-representation Analysis of Downregulated Genes",
    contrast = contrast,
    fileExtension = "Reactome_overrepresented_downreg"
  )

}