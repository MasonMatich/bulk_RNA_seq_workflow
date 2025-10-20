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
  
  knitDataTable <- function(df, )
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # #          Dysregulated Reactome Pathway GSE, GO GSE          # # # 
  cat(paste0("\n\n## Dysregulated\n\n"))
  
  # Reactome Pathway GSE
  cat("**Reactome Pathway GSEA**")
  if(nrow(data.frame(go.result.list$react_gse[[1]])) != 0){
    td1 <- list(datatable(data.frame(go.result.list$react_gse[[1]]),
                         extensions = 'Buttons',
                         rownames = FALSE,
                         caption = " Reactome Gene Set Enrichment Analysis.",
                         options = list(
                           dom = "Bfrtip",
                           buttons = list('copy',
                                       list(extend = 'csv', filename = sprintf("%s_%s_Reactome_gsea",gsub(" ", "_", tolower(strain)),tolower(fraction))),
                                       list(extend = 'excel', filename = sprintf("%s_%s_Reactome_gsea.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
                                       ),
                           pageLength = 5,
                           scrollX = TRUE)))
    print(htmltools::tagList(td1[[1]]))
    #htmltools::browsable(td1)
    cat("\n\n")
    write.csv(go.result.list$react_gse[[1]], file = sprintf("./results/%s_%s_Reactome.gsea.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{emptyTableMessage(' Reactome Gene Set Enrichment Analysis')}
  
  # # # GSEA plot condiitonal rendering of RP GSEA plots
  # Tabset for Reactome GSEA plots
  cat("#### {.tabset .tabset-dropdown}\n\n")

  # Print GSEA plot for each set identified 
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
  if(nrow(data.frame(go.result.list$go_gse[[1]])) != 0){
    td2 <- list(datatable(data.frame(go.result.list$go_gse[[1]]),
                          extensions = 'Buttons',
                          rownames = FALSE,
                          caption = " GSE Of All Differentially Expressed Genes",
                          options = list(
                            dom = "Bfrtip",
                            buttons = c('copy', 'csv', 'excel'),
                            pageLength = 5,
                            scrollX = TRUE)))
    print(htmltools::tagList(td2[[1]]))
    cat("\n\n")
    write.csv(go.result.list$go_gse[[1]], file = sprintf("./results/%s_%s_GO.gsea.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
    }else{emptyTableMessage(' GSE Of All Differentially Expressed Genes')}
  
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
  if(nrow(data.frame(go.result.list$go_class_upreg[[1]])) != 0){
    tu1 <- list(datatable(data.frame(go.result.list$go_class_upreg[[1]]),
                    extensions = 'Buttons',
                    rownames = FALSE,
                    caption = " GO Classification of Upregulated Genes.",
                    options = list(
                      dom = "Bfrtip",
                      buttons = c('copy', 'csv', 'excel'),
                      pageLength = 5,
                      scrollX = TRUE)))
    print(htmltools::tagList(tu1[[1]]))
    cat("\n\n")
    write.csv(go.result.list$go_class_upreg[[1]], file = sprintf("./results/%s_%s_GO.classification.upreg.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{emptyTableMessage(' GO Classification of Upregulated Genes')}
  
  # GO over-representation
  if(nrow(data.frame(go.result.list$go_overrep_upreg[[1]])) != 0){
    tu2 <- list(datatable(data.frame(go.result.list$go_overrep_upreg[[1]]),
                    extensions = "Buttons",
                    rownames = FALSE,
                    caption = " GO Over-representation Analysis of Upregulated Genes.",
                    options = list(
                      dom = "Bfrtip",
                      buttons = c('copy', 'csv', 'excel'),
                      pageLength = 5,
                      scrollX = TRUE)))
    print(htmltools::tagList(tu2[[1]]))
    cat("\n\n")
    write.csv(go.result.list$go_overrep_upreg[[1]], file = sprintf("./results/%s_%s_GO.overrepresented.upreg.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{emptyTableMessage(' GO Over-representation Analysis of Upregulated Genes')}
  
  # RP over-representation
  if(nrow(data.frame(go.result.list$react_overrep_upreg[[1]])) != 0){
    tu3 <- list(datatable(data.frame(go.result.list$react_overrep_upreg[[1]]),
                         extensions = 'Buttons',
                         rownames = FALSE,
                         caption = " Reactome Pathway Over-representation Analysis.",
                         options = list(
                           dom = "Bfrtip",
                           buttons = c('copy', 'csv', 'excel'),
                           pageLength = 5,
                           scrollX = TRUE)))
    print(htmltools::tagList(tu3[[1]]))
    cat("\n\n")
    write.csv(go.result.list$react_overrep_upreg[[1]], file = sprintf("./results/%s_%s_Reactome.overrepresented.upreg.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{emptyTableMessage(' Reactome Pathway Over-representation Analysis')}
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # # #    Downregulated GO classification, GO over-representation, RP over-representation    # # #
  cat(paste0("\n\n## Downregulated\n\n"))
  
  # GO classification
  if(nrow(data.frame(go.result.list$go_class_downreg[[1]])) != 0){
    tdo1 <- list(datatable(data.frame(go.result.list$go_class_downreg[[1]]),
                    extensions = "Buttons",
                    rownames = FALSE,
                    caption = " GO Classification of Downregulated Genes.",
                    options = list(
                      dom = "Bfrtip",
                      buttons = c('copy', 'csv', 'excel'),
                      pageLength = 5,
                      scrollX = TRUE)))
    print(htmltools::tagList(tdo1[[1]]))
    cat("\n\n")
    write.csv(go.result.list$go_class_downreg[[1]], file = sprintf("./results/%s_%s_GO.classification.downreg.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{emptyTableMessage(' GO Classification of Downregulated Genes')}

  # GO over-representation
  if(nrow(data.frame(go.result.list$go_overrep_downreg[[1]])) != 0){
    tdo2 <- list(datatable(data.frame(go.result.list$go_overrep_downreg[[1]]),
                    extensions = "Buttons",
                    rownames = FALSE,
                    caption = " GO Over-representation Analysis of Downregulated Genes.",
                    options = list(
                      dom = "Bfrtip",
                      buttons = c('copy', 'csv', 'excel'),
                      pageLength = 5,
                      scrollX = TRUE)))
    print(htmltools::tagList(tdo2[[1]]))
    cat("\n\n")
    write.csv(go.result.list$go_overrep_downreg[[1]], file = sprintf("./results/%s_%s_GO.overrepresented.downreg.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{emptyTableMessage(' GO Over-representation Analysis of Downregulated Genes')}
  
  # RP over-representation
  if(nrow(data.frame(go.result.list$react_overrep_downreg[[1]])) != 0){
    tdo3 <- list(datatable(data.frame(go.result.list$react_overrep_downreg[[1]]),
                         extensions = 'Buttons',
                         rownames = FALSE,
                         caption = " Reactome Pathway Over-representation Analysis.",
                         options = list(
                           dom = "Bfrtip",
                           buttons = c('copy', 'csv', 'excel'),
                           pageLength = 5,
                           scrollX = TRUE)))
    print(htmltools::tagList(tdo3[[1]]))
    cat("\n\n")
    write.csv(go.result.list$react_overrep_downreg[[1]], file = sprintf("./results/%s_%s_Reactome.overrepresented.downreg.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{emptyTableMessage(' Reactome Pathway Over-representation Analysis')}
  
  # Cleanup
  rm(td1, td2, tu1, tu2, tu3, tdo1, tdo2, tdo3)
}