# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #                                                                   #   #
#   #                           DEG Analysis                            #   # 
#   #                                                                   #   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script applies DESeq2 to input bulk RNA-seq data
DEAnalysis <- function (dataset){
  
  knitr::asis_output("\n\n# DESeq2 Differential Expression Analysis {.title}\n\n")
  cat("\n\n## Analysis Methods\n\nHere, count data from the RNA-seq experiment is read in the form of a counts matrix. Each column holds data from one sample, and each row represents a gene, such that the i-th row and n-th column tells how reads of gene i were measured in sample n.
  The values received should be un-normalized counts of sequencing reads or fragments.")
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # #                         Format Count Matrix                       # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # # Read Count Matrix
  count.matrix <- read_tsv(trimws(dataset$counts_matrix_file_path)) %>%
    mutate(across(3:ncol(.), ~ as.integer(.x))) %>%
    column_to_rownames("gene_id")
  
  gene.reference <- dplyr::select(count.matrix, 1)
  
  count.matrix <- count.matrix %>%
    dplyr::select(-1)
  
  # # Input Validation
  # Check for agreement in cond_col_name and num_cond
  if (nrow(metadata) != (ncol(count.matrix))){
    stop("number of rows entered as metadata doesn't match number of columns in count matrix.")
  }
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # #                           Format Metadata                         # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  cat("Information entered as metadata describes the samples (columns) of the count matrix. This metadata is combined with the sample/column names from the count matrix so the accuracy of the metadata can be reviewed. Take time to review the table below.\n\n")
  # Design ColData Matrix
  sample <- colnames(count.matrix)
  colData <- data.frame(sample)
  
  # Add Columns for Condition to colData
  for (i in colnames(metadata)){
    colData[i] <- metadata[i]
  }
  
  colData <- column_to_rownames(colData, "sample")
  
  # Report Column Data
  print(as.matrix(colData))
  
  # # Input Validation
  # Check matrix design against count matrix columns
  if(all(colnames(count.matrix) != rownames(colData))){
    stop(sprintf("The structure of the colData matrix doesn't match the structure of the count matrix.\n   The received colData is %s.\n   The received count matrix columns are %s", rownames(colData), colnames(count.matrix)))
  }
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # #                Run Differential Expression Test                   # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  cat("Using the count matrix, column metadata, and user-inputed design formula 
  (expressing the variables to be used in modeling), a DESeq data set is created.\n
  The experimental variable and the associated control/experimental conditions 
  are given by the user in 'ref_level_cond' and 'ref_level_value' in the params_file respectively.\n
  Pre-filtering is performed, keeping only genes that have a count of at least 10 in a minimum number of samples. 
  The minimum number of samples is decided by calculating the number of times the reference level value 
  (control condition) is listed in the metadata with the assumption that experimental conditions will be repeated 
  the same number of times.\n 
  Standard differential expression analysis is performed with the DESeq function, and regularized 
  log transformation (rlog) transforms count data to a log2 scale for PCA.")
  
  # Create DESeq Data Set
  dds <- DESeqDataSetFromMatrix(
    countData = count.matrix,
    colData = colData,
    design = as.formula(dataset$design)
  )
  
  # Set Reference Level
  reflevelcond <- trimws(dataset$ref_level_cond)
  reflevelvalue <- trimws(dataset$ref_level_value)
  dds[[reflevelcond]] <- relevel(dds[[reflevelcond]], ref = reflevelvalue)
  
  # Subset out genes with low overall counts
  smallest.group.size <- length(grep(sprintf("^%s$",reflevelvalue), metadata[[reflevelcond]]))
  keep <- rowSums(counts(dds) >= 10) >= smallest.group.size
  dds <- dds[keep,]
  
  # Run DESeq
  dds <- DESeq(dds)
  
  # Regularized Log Transform
  rld <- rlog(dds)

  cat("If input data is designated to be batch corrected, the ComBat_seq() function from the SVA Bioconductor package 
  is used. The batch (condition) to be regressed out is given in the params_file as 'batch_cond', and the sample 
  information for this condition is taken from the metadata. \n
  The group argument of ComBat_seq() specifies biological covariates, whos signals should be preserved in adjusted data. 
  All conditions (columns) from the metadata which are not to be regressed out are passed to the group argument.\n
  Differential expression analysis is otherwise performed as described above on batch-corrected counts. \n")
  
  if (as.logical(trimws(dataset$batch_correct))) {
    
    # Check that condition for batch correction matches column metadata
    if(all(trimws(dataset$batch_cond) != colnames(metadata))){
      stop(sprintf("The given condition for batch correction:%s, doesn't match any prior conditions: %s", dataset$batch_cond, colnames(metadata)))
    }
    
    # Replicate / Batch Correction
    batch <- trimws(dataset$batch_cond)
    
    batch.correct <- ComBat_seq(
      counts = as.matrix(count.matrix),
      batch = metadata[[batch]],
      group = metadata[[colnames(metadata)[colnames(metadata) != batch]]]) %>%
      as.data.frame()
    
    # Create DESeq Data Set
    batch.correct.dds <- DESeqDataSetFromMatrix(
      countData = batch.correct,
      colData = colData,
      design = as.formula(dataset$design))
    
    # Set Factor Level
    reflevelcond <- trimws(dataset$ref_level_cond)
    reflevelvalue <- trimws(dataset$ref_level_value)
    batch.correct.dds[[reflevelcond]] <- relevel(batch.correct.dds[[reflevelcond]], ref = reflevelvalue)
    
    # Subset out genes with low overall counts
    smallest.group.size <- length(grep(sprintf("^%s$",reflevelvalue), metadata[[reflevelcond]]))
    keep <- rowSums(counts(batch.correct.dds) >= 10) >= smallest.group.size
    batch.correct.dds <- batch.correct.dds[keep,]
    
    # Run DESeq
    batch.correct.dds <- DESeq(batch.correct.dds)
    
    # Regularized Log Transform
    batch.correct.rld <- rlog(batch.correct.dds)
  }

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # #                         Calculate Contrasts                       # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  cat("Reference level condition and value (ref_level_cond and ref_level_value respectively) inputs from params_file
    are used to create contrast terms to compare treatment samples with control samples. This relies on the design 
    formula being formatted correctly to produce comparisons relevant to the experimental condition.\n")
  
  # # Find relevant contrasts
  contrasts <- list()
  
  # Factor values of experimental condition
  pattern <- factor(metadata[[reflevelcond]])
  comp <- levels(pattern)
  
  # Remove control state from experimental condition vector
  comp <- comp[comp != reflevelvalue]
  
  # Create Contrasts
  for (i in 1:length(comp)){
    contrasts[[i]] <- c(reflevelcond, comp[i], reflevelvalue)
  }
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # #                         Select DDS Object                         # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Select non-batch corrected dds for further analysis
  dds.for.analysis <- dds
  
  # Select batch corrected dds for further analysis if sample_params[i]$batch_correct = TRUE
  if (as.logical(trimws(dataset$batch_correct))){
    dds.for.analysis <- batch.correct.dds
  }

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # #                          Retrieve Results                         # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  cat("Based on the number of experimental conditions (identified in the ref_level_cond column of the metadata) 
    and the contrast terms comparing them to the experimental control (ref_level_value), differential expression 
    results are extracted from the DESeq data set object.")

  # Find results from relevant contrasts
  results <- list()
  
  # Pull results from DESeq object with contrasts
  for (i in 1:length(contrasts)){
    results[[i]] <- list()
    results[[i]][["DataFrame"]] <- results(dds.for.analysis, contrast = contrasts[[i]]) %>%
      data.frame() %>%
      rownames_to_column("gene_id") %>%
      mutate(gene_name = gene.reference[gene_id, 1], .before=baseMean) %>% # Add gene names to DESeq results matrices
      column_to_rownames("gene_id")
    
    results[[i]][["DESeqDataSet"]] <- results(dds.for.analysis, contrast = contrasts[[i]])
    
    # name results
    names(results)[i] <- paste(contrasts[[i]][1], contrasts[[i]][2], contrasts[[i]][3], sep = '_')
  }
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # #                               Save Files                                # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Save DESeq Output
  #saveRDS(dds, file = sprintf("results/%s_dds.rds", dataset$project))
  #saveRDS(rld, file = sprintf("results/%s_rld.rds", dataset$project))
  
  # Save Batch Corrected DESeq Output
  if(as.logical(trimws(dataset$batch_correct))){
    #saveRDS(batch.correct.dds, file = sprintf("results/%s_batch_correct_dds.rds", dataset$project))
    #saveRDS(batch.correct.rld, file = sprintf("results/%s_batch_correct_rld.rds", dataset$project))
  }
  
  
  # Return Results
  return.results <- list(
    results,
    colData,
    rld
  )
  
  names(return.results) <- c(
    "results",
    "metadata",
    "rld"
  )
  if(as.logical(trimws(dataset$batch_correct))){
    return.results["batch correct rld"] <- batch.correct.rld
  }
  
  return(return.results)
}