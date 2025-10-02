# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #         Gene Ontology Enrichment Analysis         # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

runGO <- function(geneNames.up, geneNames.down, geneList, universe, OrgDb){
  
  
  # # # # # # # # # # # # # # # # # # # # # # #
  # # #   Analysis of Upregulated Genes   # # #
  # # # # # # # # # # # # # # # # # # # # # # #
  
  # # GO Classification # #
  go_class_upreg <- groupGO(gene = geneNames.up,
          OrgDb = OrgDb,
          ont = "CC",
          level = 3,
          readable = TRUE) %>%
    arrange(desc(Count))
  
  # # GO Over-Representation Analysis # #
  go_overrep_upreg <- enrichGO(gene = geneNames.up,
           universe = universe,
           OrgDb = OrgDb,
           ont = "CC",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.01,
           qvalueCutoff = 0.05,
           readable = TRUE) %>%
    arrange(desc(Count))
  
  
  # # # # # # # # # # # # # # # # # # # # # # # #
  # # #   Analysis of Downregulated Genes   # # #
  # # # # # # # # # # # # # # # # # # # # # # # #
  
  # # GO Classification # #
  go_class_downreg <- groupGO(gene = geneNames.down,
                               OrgDb = OrgDb,
                               ont = "CC",
                               level = 3,
                               readable = TRUE) %>%
    arrange(desc(Count))
  
  # # GO Over-Representation Analysis # #
  go_overrep_downreg <- enrichGO(gene = geneNames.down,
                                    universe = universe,
                                    OrgDb = OrgDb,
                                    ont = "CC",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.01,
                                    qvalueCutoff = 0.05,
                                    readable = TRUE) %>%
    arrange(desc(Count))
  
  
  
  # # # # # # # # # # # # # # # # # # # # # # # #
  # #     GO Gene Set Enrichment Analysis     # #
  # # # # # # # # # # # # # # # # # # # # # # # #
  go_gene_set_enrichment <- gseGO(geneList = geneList,
        OrgDb = OrgDb,
        ont = "CC",
        minGSSize = 100,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose = FALSE)
  
  go_results <- list(
    list(go_class_upreg),
    list(go_overrep_upreg), 
    list(go_class_downreg),
    list(go_overrep_downreg), 
    list(go_gene_set_enrichment)
    )
  
  names(go_results) <- c(
    "go_class_upreg", 
    "go_overrep_upreg", 
    "go_class_downreg", 
    "go_overrep_downreg", 
    "go_gse"
    )
  
  return(go_results)
}