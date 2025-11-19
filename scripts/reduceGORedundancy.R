# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #                                                                   #   #
#   #                       Reduce GO Redundancy                        #   #
#   #                                                                   #   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

reduceGORedundancy <- function(GOResult, OrgDb, semdata = GOSemSim::godata(orgdb, ont = ont, keytype = keytype), threshold = .7){
  # if semdata isn't provided, generate it
  if (is.null(semdata)) {
    semdata <- GOSemSim::godata(
      OrgDb,
      ont = "BP",
      keytype = "ENTREZID"  # or "WORMBASE" if your OrgDb uses that
    )
  }
  
  # calculate similarity matrix
  similarity_matrix <- calculateSimMatrix(
    GOResult$ID,
    orgdb = OrgDb,
    ont = "BP",
    semdata = semdata,
    method = "Rel")
  
  # associate scores with each GO term
  scores <- setNames(-log10(GOResult$qvalue), GOResult$ID)
  
  if(any(is.na(scores))){
    scores <- "uniqueness"
  }
  
  # group GO terms based on similarity and statistical significance
  reduced_terms <- reduceSimMatrix(
    simMatrix = similarity_matrix,
    scores = scores,
    threshold = threshold,
    orgdb = OrgDb)
  
  return(reduced_terms)
}