# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #         Gene Ontology Enrichment Analysis         # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

callClusterProfilerFunc <- function(result.obj.list, OrgDb){
  
  # # Extract filtered lists from list container
  geneNames.variable <- result.obj.list$variable_genes[,1]
  geneNames.up <- result.obj.list$upregulated_genes[,1]
  geneNames.down <- result.obj.list$downregulated_genes[,1]
  geneList <- result.obj.list$vectorized_all_DE
  universe <- result.obj.list$all_genes[,1]
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #                       Analysis of Upregulated Genes                   # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  ## RUN WITH BP ONTOLOGY
  # # GO Classification # #
  go_class_upreg <- groupGO(
      gene = geneNames.up,
      OrgDb = OrgDb,
      ont = "BP",
      level = 3,
      readable = TRUE) %>%
    arrange(desc(Count))
  
  # # GO Over-Representation Analysis # #
  go_overrep_upreg <- enrichGO(
      gene = geneNames.up,
      universe = universe,
      OrgDb = OrgDb,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.1,
      readable = TRUE) %>%
    arrange(desc(Count))
  
  # # Reactome Pathway Over-Representation Analysis # # 
  reactome_overrep_upreg <- enrichPathway(
      gene = geneNames.up, 
      universe = universe,
      organism = 'celegans',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.1,
      minGSSize = 20,
      maxGSSize = 500,
      readable = TRUE) %>%
    arrange(desc(Count))
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #                      Analysis of Downregulated Genes                  # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # # GO Classification # #
  go_class_downreg <- groupGO(
      gene = geneNames.down,
      OrgDb = OrgDb,
      ont = "BP",
      level = 3,
      readable = TRUE) %>%
    arrange(desc(Count))
  
  # # GO Over-Representation Analysis # #
  go_overrep_downreg <- enrichGO(
      gene = geneNames.down,
      universe = universe,
      OrgDb = OrgDb,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.1,
      readable = TRUE) %>%
    arrange(desc(Count))
  
  # # Reactome Pathway Over-Representation Analysis # # 
  reactome_overrep_downreg <- enrichPathway(
      gene = geneNames.down, 
      universe = universe,
      organism = 'celegans',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.1,
      minGSSize = 20,
      maxGSSize = 500,
      readable = TRUE) %>%
    arrange(desc(Count))
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #                    GO Gene Set Enrichment Analysis                    # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  go_gene_set_enrichment <- gseGO(
    geneList = geneList,
        OrgDb = OrgDb,
        ont = "BP",
        minGSSize = 20,
        maxGSSize = 500,
        pvalueCutoff = 0.1,
        verbose = FALSE
        )
  
  # change gene IDs from ENTREZID to symbol
  go_gene_set_enrichment <- setReadable(go_gene_set_enrichment, OrgDb = OrgDb, keyType="ENTREZID")
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #                       Reactome Enrichment Analysis                    # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  reactome_gene_set_enrichment <- gsePathway(
    geneList = geneList, 
    organism = 'celegans',
    pvalueCutoff = 0.2,
    pAdjustMethod = "BH", 
    minGSSize = 20,
    maxGSSize = 500,
    verbose = FALSE
    )
  
  # change Reactome GSE gene IDs from ENTREZID to symbol
  reactome_gene_set_enrichment <- setReadable(reactome_gene_set_enrichment, OrgDb = OrgDb, keyType="ENTREZID")
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #                             Return Results                            # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  go_results <- list(
    list(go_class_upreg),
    list(go_overrep_upreg), 
    list(reactome_overrep_upreg),
    list(go_class_downreg),
    list(go_overrep_downreg), 
    list(reactome_overrep_downreg),
    list(go_gene_set_enrichment),
    list(reactome_gene_set_enrichment)
    )
  
  names(go_results) <- c(
    "go_class_upreg", 
    "go_overrep_upreg", 
    "react_overrep_upreg",
    "go_class_downreg", 
    "go_overrep_downreg", 
    "react_overrep_downreg",
    "go_gse",
    "react_gse"
    )
  
  return(go_results)
}