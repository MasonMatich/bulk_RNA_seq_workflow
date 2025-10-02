# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# #         Calls translateBioIDs.R and geneOntologyAnalysis.R        # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

strainFractionGO <- function(result.obj, strain, fraction, l2fc_filter, padj_filter){
  
  # Create comment structure for bookdown render
  cat(paste0("# ", strain," ", fraction, " {.tabset .tabset}\n\n"))

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
  go.result.list <- runGO(
    geneNames.up = result.obj.list$upregulated_genes[,1],
    geneNames.down = result.obj.list$downregulated_genes[,1],
    geneList = result.obj.list$vectorized_all_DE,
    universe = result.obj.list$all_genes[,1],
    OrgDb = org.Ce.eg.db)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #     Print text and tables to console for bookdown rendering     # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  
  # # #                             GSEA                            # # # 
  cat(paste0("## GSEA \n\n"))
  # Print GSEA table
  t0 <- list(datatable(data.frame(go.result.list$go_gse[[1]]),
                       extensions = 'Buttons',
                       caption = " GSE Of All DEGs.",
                       options = list(
                         dom = "Bfrtip",
                         buttons = c('csvHtml5'),
                         pageLength = 5)))
  
  print(htmltools::tagList(t0[[1]]))
  cat("\n\n")
  
  # Tabset for GSEA plots
  cat(sprintf("#### %s %s GSEA Plots {.tabset .tabset-pills}\n\n", strain, fraction))
  
  
  # Print GSEA plot for each set identified 
  for (i in 1:nrow(data.frame(go.result.list$go_gse[[1]]))) {
    cat(sprintf("##### %s \n\n", go.result.list$go_gse[[1]]$Description[i]))
    print(gseaplot(go.result.list$go_gse[[1]], by = "all", title = go.result.list$go_gse[[1]]$Description[i], geneSetID = go.result.list$go_gse[[1]]$ID[i]))
    cat("\n\n")
  }

  
  # # #    Upregulated GO classification and GO over-representation   # # #
  cat(paste0("\n\n## Upregulated\n\n"))
  
  if(nrow(data.frame(go.result.list$go_class_upreg[[1]])) != 0){
    t1 <- list(datatable(data.frame(go.result.list$go_class_upreg[[1]]),
                    extensions = 'Buttons',
                    caption = " GO Classification of Upregulated Genes.",
                    options = list(
                      dom = "Bfrtip",
                      buttons = c('csvHtml5'),
                      pageLength = 5)))
    print(htmltools::tagList(t1[[1]]))
    cat("\n\n")
    write.csv(go.result.list$go_class_upreg[[1]], file = sprintf("./results/%s_%s_GO.classification.upreg.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{
    emptyTableMessage(' GO Classification of Upregulated Genes')
  }
  
  if(nrow(data.frame(go.result.list$go_overrep_upreg[[1]])) != 0){
    t2 <- list(datatable(data.frame(go.result.list$go_overrep_upreg[[1]]),
                    extensions = "Buttons",
                    caption = " GO Over-representation Analysis of Upregulated Genes.",
                    options = list(
                      dom = "Bfrtip",
                      buttons = c('csvHtml5'),
                      pageLength = 5)))
    print(htmltools::tagList(t2[[1]]))
    cat("\n\n")
    write.csv(go.result.list$go_overrep_upreg[[1]], file = sprintf("./results/%s_%s_GO.overrepresented.upreg.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{
    emptyTableMessage(' GO Over-representation Analysis of Upregulated Genes')
  }
  
  
  # # #    Downregulated GO classification and GO over-representation   # # #
  cat(paste0("\n\n## Downregulated\n\n"))
  
  if(nrow(data.frame(go.result.list$go_class_downreg[[1]])) != 0){
    t3 <- list(datatable(data.frame(go.result.list$go_class_downreg[[1]]),
                    extensions = "Buttons",
                    caption = " GO Classification of Downregulated Genes.",
                    options = list(
                      dom = "Bfrtip",
                      buttons = c('csvHtml5'),
                      pageLength = 5)))
    print(htmltools::tagList(t3[[1]]))
    cat("\n\n")
    write.csv(go.result.list$go_class_downreg[[1]], file = sprintf("./results/%s_%s_GO.classification.downreg.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{
    emptyTableMessage(' GO Classification of Downregulated Genes')
  }

  if(nrow(data.frame(go.result.list$go_overrep_downreg[[1]])) != 0){
    t4 <- list(datatable(data.frame(go.result.list$go_overrep_downreg[[1]]),
                    extensions = "Buttons",
                    caption = " GO Over-representation Analysis of Downregulated Genes.",
                    options = list(
                      dom = "Bfrtip",
                      buttons = c('csvHtml5'),
                      pageLength = 5)))
    print(htmltools::tagList(t4[[1]]))
    cat("\n\n")
    write.csv(go.result.list$go_overrep_downreg[[1]], file = sprintf("./results/%s_%s_GO.overrepresented.downreg.csv",gsub(" ", "_", tolower(strain)),tolower(fraction)))
  }else{
    emptyTableMessage(' GO Over-representation Analysis of Downregulated Genes')
  }
  rm(t0, t1, t2, t3, t4)
}