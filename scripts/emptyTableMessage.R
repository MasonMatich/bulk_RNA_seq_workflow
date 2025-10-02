# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #     Uses cat() to print that a table is empty to console    # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

emptyTableMessage <- function(name_of_table){
  cat(sprintf("**Table %s has no rows.**", name_of_table))
}