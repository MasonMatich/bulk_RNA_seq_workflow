# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #                                                                     #   #
#   #       Print text and tables to console for bookdown rendering       #   #
#   #                                                                     #   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

knitDataTable <- function(df, tableName, sample_name, contrast, fileExtension, pageLength = 5){
  if(nrow(data.frame(df[[1]])) != 0){
    table <- list(datatable(data.frame(df[[1]]),
                            extensions = 'Buttons',
                            rownames = FALSE,
                            caption = tableName,
                            options = list(
                              dom = "frtipB",
                              buttons = list('copy',
                                             list(extend = 'csv', filename = sprintf("%s_%s_%s", sample_name, gsub(" ", "_", tolower(contrast)), fileExtension)),
                                             list(extend = 'excel', filename = sprintf("%s_%s_%s", sample_name, gsub(" ", "_", tolower(contrast)), fileExtension))
                              ),
                              pageLength = pageLength,
                              scrollX = TRUE
                              )
                            )
                  )
    # Print DataTable
    print(htmltools::tagList(table[[1]]))
    cat("\n\n")
    
    # Write csv to results with table contents 
    write.csv(df[[1]], file = sprintf("../results/%s_%s_%s.csv", sample_name, gsub(" ", "_", tolower(contrast)), fileExtension))
    }else{emptyTableMessage(tableName)}
  
  # clean up environment
  rm(table)
}