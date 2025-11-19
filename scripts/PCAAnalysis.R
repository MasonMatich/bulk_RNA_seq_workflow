# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #                                                                   #   #
#   #                           PCA Analysis                            #   # 
#   #                                                                   #   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script performs a Principle Component Analysis on bulk RNA-seq 

PCAAnalysis <- function(dataset, rld, batch.rld, metadata){
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # #                           No Regression PCA                       # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  cat("\n\n# Principle Component Analysis {.title}\n\n")
  cat("The Bioconductor package PCAtools is used to perform principle component analysis on regularized log transformed count data. \n
    Scree plots show principle component numbers on the x-axis, and their respective eigenvalues on the y-axis. \n
    The third graph plots experimental metadata against a number of PCs and their gene loadings to visualize the agreement with 
    the gene expression pattern of each condition for each PC.")
  
  cat("\n\n## No Regression\n\n")
  # View non-batch corrected
  pca <- visualizePCsByCondition(
    rlog = rld,
    metadata = metadata,
    title = sprintf("Before Regressing out%s", str_to_title(dataset$batch_cond)))
  
  print(pca$logScaleScree)
  print(pca$linearScaleScree)
  print(pca$selectedPCs)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # #                            Regression PCA                         # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  if(as.logical(trimws(dataset$batch_correct))){
    cat(sprintf("\n\n##%s Regressed Out\n\n", str_to_title(dataset$batch_cond)))
  }
  
  # View Rep effect corrected
  if(as.logical(trimws(dataset$batch_correct))){
    pca.batch <- visualizePCsByCondition(
      rlog = batch.rld,
      metadata = metadata,
      title = sprintf("After Regressing out%s", str_to_title(dataset$batch_cond)))
    
    print(pca.batch$logScaleScree)
    print(pca.batch$linearScaleScree)
    print(pca.batch$selectedPCs)
  }
  
  # Save PCA to PDF
  pdf(sprintf("figures/%s_PCA.pdf", dataset$project))
  
  # View non-batch corrected
  pca$logScaleScree
  pca$linearScaleScree
  pca$selectedPCs
  
  if(as.logical(trimws(dataset$batch_correct))){
    # View Rep effect corrected
    pca.batch$logScaleScree
    pca.batch$linearScaleScree
    pca.batch$selectedPCs
  }
  dev.off()

}