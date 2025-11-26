# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #                                                                   #   #
#   #        Run DESeq2 Analysis with Optional Batch Correction         #   #
#   #                                                                   #   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

runDESeq2Analysis <- function(count.matrix, colData, design, ref_level_cond, ref_level_value, metadata) {

  # Create DESeq Data Set
  dds <- DESeqDataSetFromMatrix(
    countData = count.matrix,
    colData = colData,
    design = as.formula(design)
  )

  # Set Reference Level
  dds[[ref_level_cond]] <- relevel(dds[[ref_level_cond]], ref = ref_level_value)

  # Subset out genes with low overall counts
  smallest.group.size <- length(grep(sprintf("^%s$", ref_level_value), metadata[[ref_level_cond]]))
  keep <- rowSums(counts(dds) >= 10) >= smallest.group.size
  dds <- dds[keep,]

  # Run DESeq
  dds <- DESeq(dds)

  # Regularized Log Transform
  rld <- rlog(dds)

  # Return both objects
  list(dds = dds, rld = rld)
}
