# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# #         Calls translateBioIDs.R and geneOntologyAnalysis.R        # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

compileOntologyReport <- function(result.obj, strain, fraction, l2fc_filter, padj_filter){
  
  # Create comment structure for bookdown render
  cat(sprintf("# %s %s {.title .tabset}\n\n", strain, fraction))
  #cat(paste0("# ", strain," ", fraction, " {.title .tabset}\n\n"))

  # Show exact call which invoked function
  call_txt <- paste(deparse(match.call()), collapse = "\n")
  cat("**Inputs for function call:**\n\n")
  cat("```r\n")
  cat(call_txt, "\n\n")
  cat("```\n\n")
  
  # Translate Wormbase IDs to ENTREZID
  result.obj.list <- translateBioIDs(
    DESeqResults.obj = result.obj,
    bioID = "WORMBASE",
    l2fc_filter = l2fc_filter,
    padj_filter = padj_filter
  )
  
  # Run gruopGO(), enrichGO(), gseGO(), and goplot()
  go.result.list <- callClusterProfilerFunc(
    result.obj.list = result.obj.list,
    OrgDb = org.Ce.eg.db)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #     Print text and tables to console for bookdown rendering     # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  knitDataTable <- function(df, tableName, fileExtension){
    if(nrow(data.frame(df[[1]])) != 0){
      table <- list(datatable(data.frame(df[[1]]),
                            extensions = 'Buttons',
                            rownames = FALSE,
                            caption = tableName,
                            options = list(
                              dom = "Bfrtip",
                              buttons = list('copy',
                                 list(extend = 'csv', filename = sprintf("%s_%s_%s", gsub(" ", "_", tolower(strain)), tolower(fraction), fileExtension)),
                                 list(extend = 'excel', filename = sprintf("%s_%s_%s", gsub(" ", "_", tolower(strain)), tolower(fraction), fileExtension))
                              ),
                              pageLength = 5,
                              scrollX = TRUE)))
      print(htmltools::tagList(table[[1]]))
      cat("\n\n")
      write.csv(df[[1]], file = sprintf("./results/%s_%s_%s.csv", gsub(" ", "_", tolower(strain)), tolower(fraction), fileExtension))
    }else{emptyTableMessage(tableName)}
    
    # clean up environment
    rm(table)
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # #          Dysregulated Reactome Pathway GSE, GO GSE          # # # 
  cat(paste0("\n\n## Dysregulated\n\n"))
  
  # Reactome Pathway GSE
  cat("**Reactome Pathway GSEA**")
  knitDataTable(
    df = go.result.list$react_gse,
    tableName = " Reactome Gene Set Enrichment Analysis",
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
  cat("### \n\n") # Break out of previous tabset
  cat("**GO GSEA**")
  knitDataTable(
    df = go.result.list$go_gse,
    tableName = " GSE Of All Differentially Expressed Genes",
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
  cat(paste0("\n\n## Upregulated\n\n"))
  
  # GO classification
  knitDataTable(
    df = go.result.list$go_class_upreg,
    tableName = " GO Classification of Upregulated Genes",
    fileExtension = "GO_classification_upreg"
  )
  
  # GO over-representation
  knitDataTable(
    df = go.result.list$go_overrep_upreg,
    tableName = " GO Over-representation Analysis of Upregulated Genes",
    fileExtension = "GO_overrepresented_upreg"
  )
  
  # RP over-representation
  knitDataTable(
    df = go.result.list$react_overrep_upreg,
    tableName = " Reactome Pathway Over-representation Analysis of Upregulated Genes",
    fileExtension = "Reactome_overrepresented_upreg"
  )
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # #    Downregulated GO classification, GO over-representation, RP over-representation    # # #
  cat(paste0("\n\n## Downregulated\n\n"))
  
  # GO classification
  knitDataTable(
    df = go.result.list$go_class_downreg,
    tableName = " GO Classification of Downregulated Genes",
    fileExtension = "GO_classification_downreg"
  )

  # GO over-representation
  knitDataTable(
    df = go.result.list$go_overrep_downreg,
    tableName = " GO Over-representation Analysis of Downregulated Genes",
    fileExtension = "GO_overrepresented_downreg"
  )
  
  # RP over-representation
  knitDataTable(
    df = go.result.list$react_overrep_downreg,
    tableName = " Reactome Pathway Over-representation Analysis of Downregulated Genes",
    fileExtension = "Reactome_overrepresented_downreg"
  )

}