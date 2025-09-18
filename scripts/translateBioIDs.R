# # # # # # # # # # # # # # # # # # # # # # # # #
# #         Translating Biological IDS        # #
# # # # # # # # # # # # # # # # # # # # # # # # #
# # # From WormBase to ENTREZID

translateBioIDs <- function (DESeqResults.obj, bioID, return.success = TRUE) {
  
  if (bioID != "WORMBASE"){
    print("⛔️ ERROR: Only WORMBASE biological ID supported as input.")
    return(0)
  }
  
  # Convert DESeq2 results object to DataFrame, and edit structure
  gene.df <- DESeqResults.obj %>%
    data.frame() %>%
    rownames_to_column("WORMBASE") %>%
    dplyr::select(WORMBASE, log2FoldChange)
  
  # Biological ID TranslatoR
  gene_ids <- bitr(gene.df[,'WORMBASE'],
                   fromType = "WORMBASE",
                   toType = "ENTREZID",
                   OrgDb = org.Ce.eg.db)
  cat("\n")
  
  # Record transfer success
  if (return.success == TRUE){
    t1 <- table(gene.df$WORMBASE %in% gene_ids$WORMBASE)
    cat(paste0(t1[["FALSE"]]," gene IDs were not transfered from WORMBASE to ENTREZID \n\n"))
    cat(paste0(t1[["TRUE"]]," gene IDs were successfully transfered from WORMBASE to ENTREZID \n\n"))
    cat("")
  }
  
  # Clear duplicate genes following ID transfer
  gene_ids <- gene_ids %>%
    distinct(WORMBASE, .keep_all = T) %>%
    column_to_rownames("WORMBASE")
  
  # Add ENTREZID to L2FC carrying DataFrame
  gene.df <- gene.df %>%
    mutate(ENTREZID = gene_ids[WORMBASE, 1]) %>%
    dplyr::select(ENTREZID, log2FoldChange) %>%
    subset(!is.na(ENTREZID)) # Remove NA ENTREZ ID values
  rownames(gene.df) <- NULL
  
  # #    Subset L2FC Data    # # 
  
  # Subset for highly varirable genes
  variable.genes <- gene.df %>%
    subset(abs(log2FoldChange) > 1.5)
  rownames(variable.genes) <- NULL
  
  # Subset for highly upregulated genes
  upregulated.genes <- gene.df %>%
    subset(log2FoldChange > 1.5)
  rownames(upregulated.genes) <- NULL
  
  # Subset for highly downregulated genes
  downregulated.genes <- gene.df %>%
    subset(log2FoldChange < -1.5)
  rownames(downregulated.genes) <- NULL
  
  # Subset for all genes (universe DataFrame)
  all.genes <- gene.df
  
  # Vectorized differentially expressed genes
  vectorized <- all.genes[,2]
  names(vectorized) <- all.genes[,1]
  vectorized <- sort(vectorized, decreasing = TRUE)
  
  # Save all data to list
  gene.list <- list(variable.genes, upregulated.genes, downregulated.genes, all.genes, vectorized)
  names(gene.list) <- c("variable_genes", "upregulated_genes", "downregulated_genes", "all_genes", "vectorized_all_DE")
  
  return(gene.list)
}