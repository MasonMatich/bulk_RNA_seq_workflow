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
  
  # # # Reduce GO Over-representation Redundancy
  if(length(go_overrep_upreg$ID) > 2){
    go_overrep_upreg_reduced <- reduceGORedundancy(
      GOResult = go_overrep_upreg,
      OrgDb = OrgDb,
      semdata = semantic_data,
      threshold = .7
    )
  }else{go_overrep_upreg_reduced <- NA}
  
  # # # Simplify GO Over-representation Redundancy
  go_overrep_upreg <- clusterProfiler::simplify(
    go_overrep_upreg,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
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
  
  # # # Reduce GO Over-representation Redundancy
  if(length(go_overrep_downreg$ID) > 2){
    go_overrep_downreg_reduced <- reduceGORedundancy(
      GOResult = go_overrep_downreg,
      OrgDb = OrgDb,
      semdata = semantic_data,
      threshold = .7
    )
  }else{go_overrep_downreg_reduced <- NA}
  
  # # # Simplify GO Over-representation Redundancy
  go_overrep_downreg <- clusterProfiler::simplify(
    go_overrep_downreg,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
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
  
  # # # Reduce GO Redundancy
  if(length(go_gene_set_enrichment$ID) > 2){
    gsea_reduced <- reduceGORedundancy(
      GOResult = go_gene_set_enrichment,
      OrgDb = OrgDb,
      semdata = semantic_data,
      threshold = .7
    )
  }else{gsea_reduced <- NA}
  
  # # # Simplify GO Over-representation Redundancy
  go_gene_set_enrichment <- clusterProfiler::simplify(
    go_gene_set_enrichment,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
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
    list(go_overrep_upreg_reduced),
    list(reactome_overrep_upreg),
    list(go_class_downreg),
    list(go_overrep_downreg), 
    list(go_overrep_downreg_reduced),
    list(reactome_overrep_downreg),
    list(go_gene_set_enrichment),
    list(gsea_reduced),
    list(reactome_gene_set_enrichment)
  )
  
  names(go_results) <- c(
    "go_class_upreg", 
    "go_overrep_upreg", 
    "go_overrep_upreg_red",
    "react_overrep_upreg",
    "go_class_downreg", 
    "go_overrep_downreg", 
    "go_overrep_downreg_red", 
    "react_overrep_downreg",
    "go_gse",
    "go_gse_red",
    "react_gse"
  )
  
  return(go_results)
}