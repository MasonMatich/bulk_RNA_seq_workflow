# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# #         Calls translateBioIDs.R and geneOntologyAnalysis.R        # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

strainFractionGO <- function(result.obj, strain, fraction){
  
  # Create comment structure for bookdown render
  cat(paste0("# ", strain," ", fraction, " {.tabset .tabset-pills}\n\n"))

  # Show exact call which invoked function
  call_txt <- paste(deparse(match.call()), collapse = "\n")
  cat("Inputs for function call:\n\n")
  cat("```r\n")
  cat(call_txt, "\n\n")
  cat("```\n\n")
  
  # Translate Wormbase IDs to ENTREZID
  result.obj.list <- translateBioIDs(
    DESeqResults.obj = result.obj,
    bioID = "WORMBASE"
  )
  
  # Run gruopGO(), enrichGO(), gseGO(), and goplot()
  go.result.list <- geneOntologyAnalysis(
    geneNames.up = result.obj.list$upregulated_genes[,1],
    geneNames.down = result.obj.list$downregulated_genes[,1],
    geneList = result.obj.list$vectorized_all_DE,
    universe = result.obj.list$all_genes[,1],
    OrgDb = org.Ce.eg.db)
  
  # Print out Gene set enrichment analysis plot
  cat(paste0("## GSEA \n\n"))
  print(go.result.list$go_plot[[1]])
  
  # Print out GO classification and GO over-representation of unregulated genes
  cat(paste0("\n\n## Upregulated\n\n"))
  
  print(knitr::kable(head(go.result.list$go_class_upreg[[1]], 10),
    booktabs = TRUE,
    caption = ' GO Classification of Upregulated Genes.'))
  cat("\n\n")
  
  print(knitr::kable(
    head(go.result.list$go_overrep_upreg[[1]], 10),
    booktabs = TRUE,
    caption = ' GO Over-representation Analysis of Upregulated Genes.'))
  
  # Print out GO classification and GO over-representation of downregulated genes
  cat(paste0("\n\n## Downregulated\n\n"))
  
  print(knitr::kable(
    head(go.result.list$go_class_downreg[[1]], 10),
    booktabs = TRUE,
    caption = ' GO Classification of Downregulated Genes.'))
  
  cat("\n\n")
  
  print(knitr::kable(
    head(go.result.list$go_overrep_downreg[[1]], 10),
    booktabs = TRUE,
    caption = ' GO Over-representation Analysis of Downregulated Genes.'))
}