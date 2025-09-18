# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #         Gene Ontology Enrichment Analysis         # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
geneOntologyAnalysis <- function(geneNames.up, geneNames.down, geneList, universe, OrgDb){
  
  # # GO Classification # #
  go_classification <- groupGO(gene = geneNames,
          OrgDb = OrgDb,
          ont = "CC",
          level = 3,
          readable = TRUE)
  
  # # GO Over-Representation Analysis # #
  go_overrepresentation <- enrichGO(gene = geneNames,
           universe = universe,
           OrgDb = OrgDb,
           ont = "CC",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.01,
           qvalueCutoff = 0.05,
           readable = TRUE)
  
  # # GO Gene Set Enrichment Analysis # #
  go_gene_set_enrichment <- gseGO(geneList = geneList,
        OrgDb = OrgDb,
        ont = "CC",
        minGSSize = 100,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose = FALSE)
  
  go_plot <- goplot(go_gene_set_enrichment)
  
  go_results <- list(go_classification, go_overrepresentation, go_gene_set_enrichment, go_plot)
  names(go_results) <- c("go_class", "go_overrep", "go_gse", "go_plot")
  
  return(go_results)
}