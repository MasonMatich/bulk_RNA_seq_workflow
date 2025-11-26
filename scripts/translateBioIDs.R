# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #                                                                   #   #
#   #                      Translating Biological IDS                   #   #
#   #                                                                   #   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # From WormBase to ENTREZID

translateBioIDs <- function (DESeqResults.obj, bioID, return.success = TRUE, l2fc_filter, padj_filter) {
  
  if (bioID != "WORMBASE"){
    stop("⛔️ ERROR: Only WORMBASE biological ID supported as input.")
  }
  
  # Convert DESeq2 results object to DataFrame, and edit structure
  gene.df <- DESeqResults.obj %>%
    data.frame() %>%
    rownames_to_column("WORMBASE") %>%
    dplyr::select(WORMBASE, log2FoldChange, pvalue, padj)
  
  # Biological ID TranslatoR
  gene_ids <- bitr(gene.df[,'WORMBASE'],
                   fromType = "WORMBASE",
                   toType = "ENTREZID",
                   OrgDb = org.Ce.eg.db)
  
  # Record success of ID transfer
  t1 <- table(gene.df$WORMBASE %in% gene_ids$WORMBASE)
  
  # Clear duplicate genes following ID transfer
  gene_ids <- gene_ids %>%
    distinct(WORMBASE, .keep_all = T) %>%
    column_to_rownames("WORMBASE")
  
  # Add ENTREZID to L2FC carrying DataFrame
  gene.df <- gene.df %>%
    mutate(ENTREZID = gene_ids[WORMBASE, 1]) %>%
    dplyr::select(ENTREZID, log2FoldChange, pvalue, padj) %>%
    subset(!is.na(ENTREZID)) %>% # remove NA ENTREZ ID values
    arrange(desc(log2FoldChange)) # <--- arrange genes by L2FC descending
  rownames(gene.df) <- NULL
  
  # #    Subset L2FC Data    # # 
  
  # Subset for highly varirable genes
  variable.genes <- gene.df %>%
    subset(abs(log2FoldChange) > l2fc_filter & padj < padj_filter & !is.na(padj))
  rownames(variable.genes) <- NULL

  # Subset for highly upregulated genes
  upregulated.genes <- gene.df %>%
    subset(log2FoldChange > l2fc_filter & padj < padj_filter & !is.na(padj))
  rownames(upregulated.genes) <- NULL

  # Subset for highly downregulated genes
  downregulated.genes <- gene.df %>%
    subset(log2FoldChange < -l2fc_filter & padj < padj_filter & !is.na(padj))
  rownames(downregulated.genes) <- NULL

  # Subset for all genes (universe DataFrame)
  all.genes <- gene.df %>%
    subset(!is.na(padj)) %>%
    arrange(desc(log2FoldChange))
  
  # Vectorized differentially expressed genes
  vectorized <- all.genes[,2]
  names(vectorized) <- all.genes[,1]
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #                     Print Results Following Subset                    # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  # Record transfer success
  if (return.success == TRUE){
    cat("\n")
    cat("**Transfering WORMBASE gene IDs to ENTREZID:**\n\n")
    cat("```r\n")
    cat(paste0(t1[["FALSE"]]," gene IDs were not transfered from WORMBASE to ENTREZID \n"))
    cat(paste0(t1[["TRUE"]]," gene IDs were successfully transfered from WORMBASE to ENTREZID \n\n"))
    cat("```\n\n")
  }
  
  # report number of highly variable genes used for analysis
  cat(sprintf("**Number of highly variable genes identified (with | L2FC | > %s and padj < %s) is:**\n\n", toString(l2fc_filter) , toString(padj_filter)))
  cat("```r\n")
  cat(length(variable.genes[[1]]))
  cat("\n\n```\n\n")
  
  # report number of highly variable genes used for analysis
  cat(sprintf("**Number of highly upregulated genes identified (with L2FC > %s and padj < %s) is:**\n\n", toString(l2fc_filter) , toString(padj_filter)))
  cat("```r\n")
  cat(length(upregulated.genes[[1]]))
  cat("\n\n```\n\n")
  
  # report number of highly variable genes used for analysis
  cat(sprintf("**Number of highly downregulated genes identified (with L2FC < -%s and padj < %s) is:**\n\n", toString(l2fc_filter) , toString(padj_filter)))
  cat("```r\n")
  cat(length(downregulated.genes[[1]]))
  cat("\n\n```\n\n")
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #                               Return Data                             # # 
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Save all data to list
  gene.list <- list(
    variable.genes, 
    upregulated.genes, 
    downregulated.genes, 
    all.genes, 
    vectorized
    )
  
  names(gene.list) <- c(
    "variable_genes", 
    "upregulated_genes", 
    "downregulated_genes", 
    "all_genes", 
    "vectorized_all_DE"
    )
  
  return(gene.list)
}