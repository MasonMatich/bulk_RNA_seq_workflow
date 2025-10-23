# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #        Print text and tables to console for bookdown rendering        # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

knitDataTable <- function(df, tableName, contrast, fileExtension, pageLength = 5){
  if(nrow(data.frame(df[[1]])) != 0){
    table <- list(datatable(data.frame(df[[1]]),
                            extensions = 'Buttons',
                            rownames = FALSE,
                            caption = tableName,
                            options = list(
                              dom = "Bfrtip",
                              buttons = list('copy',
                                             list(extend = 'csv', filename = sprintf("%s_%s", gsub(" ", "_", tolower(contrast)), fileExtension)),
                                             list(extend = 'excel', filename = sprintf("%s_%s", gsub(" ", "_", tolower(contrast)), fileExtension))
                              ),
                              pageLength = pageLength,
                              scrollX = TRUE)))
    print(htmltools::tagList(table[[1]]))
    cat("\n\n")
    write.csv(df[[1]], file = sprintf("./results/%s_%s.csv", gsub(" ", "_", tolower(contrast)), fileExtension))
  }else{emptyTableMessage(tableName)}
  
  # clean up environment
  rm(table)
}