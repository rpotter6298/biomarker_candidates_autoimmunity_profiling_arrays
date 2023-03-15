# Function to subset the dataframe based on the lowest (n) values
subset_top_n <- function(df, n) {
  if ("Q.Value" %in% colnames(df)) {
    subset_df <- df %>% top_n(-n, Q.Value)
  } else if ("adj.P.Val" %in% colnames(df)) {
    subset_df <- df %>% top_n(-n, adj.P.Val)
  } else {
    return(df)
  }
  
  return(subset_df)
}

subset_export <- function(dflist, n=15){
  # Apply the function to the dflist
  subsetted_dflist <- lapply(dflist, subset_top_n, n)
  # Get the name of the dflist
  dflist_name <- deparse(substitute(dflist))
  # Write the list of subsetted dataframes to an Excel file
  excel_file_name <- paste0(dflist_name, "_Top_", n, ".xlsx")
  # Create a new workbook
  wb <- createWorkbook()
  # Add worksheets for each subsetted dataframe and write data
  for (i in seq_along(subsetted_dflist)) {
    sheet_name <- names(dflist)[i]
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, subsetted_dflist[[i]])
  }
  
  # Save the workbook
  saveWorkbook(wb, file = file.path("stats", excel_file_name), overwrite = TRUE)
}


export_excel <- function(dflist, name = deparse(substitute(dflist))) {
  file_path <- file.path(getwd(), "stats", paste0(name, ".xlsx"))
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Iterate through the list of dataframes and add worksheets
  for (setid in seq_along(dflist)) {
    sheet_name <- paste0("Dataset_", setid)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, dflist[[setid]])
  }
  
  # Save the workbook
  saveWorkbook(wb, file = file_path, overwrite = TRUE)
}
