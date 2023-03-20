df <- P03_Transformed$Set_3

analytes <- df[-1:-2]
spearman_correlation <- cor(analytes, method = "spearman") 

# Create a color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(25)

# Plot heatmap
pheatmap(spearman_correlation, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = color_palette, 
         show_rownames = TRUE, 
         show_colnames = TRUE)

spearman_genes <- replace_analyte_ID_rcnames(spearman_correlation, I00_Antigens)


strongest_correlation <- function(correlation_table, n, genes_of_interest=NULL){
  lowertri = correlation_table[lower.tri(correlation_table, diag = FALSE)]
  sorted_abs_corr <- sort(abs(lowertri), decreasing = TRUE)
  top_n_indices <- head(order(abs(lowertri), decreasing = TRUE), n)
  
  # Get row and column indices
  row_indices <- row(correlation_table)[top_n_indices]
  col_indices <- col(correlation_table)[top_n_indices]
  
  # Create a data frame of the top n correlated pairs
  top_n_correlations <- data.frame(
    row = rownames(correlation_table)[row_indices],
    col = colnames(correlation_table)[col_indices],
    correlation = lowertri[top_n_indices])
    return(top_n_correlations)
}

strongest_correlation <- function(cor_matrix, n, genes_of_interest = NULL) {
  # Get the lower triangle of the correlation matrix and set the diagonal to NA
  cor_matrix_lower <- cor_matrix
  cor_matrix_lower[upper.tri(cor_matrix_lower, diag = TRUE)] <- NA
  
  # Melt the correlation matrix into a long format dataframe
  cor_df <- reshape2::melt(cor_matrix_lower, na.rm = TRUE)
  
  # Calculate the absolute value of correlations
  cor_df$abs_correlation <- abs(cor_df$value)
  
  # Filter based on the genes_of_interest if provided
  if (!is.null(genes_of_interest)) {
    cor_df <- cor_df[(cor_df$Var1 %in% genes_of_interest) | 
                       (cor_df$Var2 %in% genes_of_interest), ]
  }
  
  # Sort the dataframe by absolute value of correlations and select the top n
  top_correlations <- cor_df[order(-cor_df$abs_correlation),][1:n,]
  
  return(top_correlations)
}


strong = strongest_correlation(spearman_genes, genes_of_interest = top_genes, n=15)

library(dplyr)

replace_analyte_ID_rcnames <- function(df, antigens_list, type = "Gene.name") {
  lookup_df <- do.call(rbind, antigens_list)
  lookup_df <- distinct(lookup_df, Antigen.name, .keep_all = TRUE)
  
  # Replace row names
  for (i in 1:nrow(df)) {
    antigen_name <- rownames(df)[i]
    gene_name <- lookup_df[lookup_df$Antigen.name == antigen_name, type]
    
    if (length(gene_name) > 0) {
      rownames(df)[i] <- gene_name
    }
  }
  
  # Replace column names
  for (i in 1:ncol(df)) {
    antigen_name <- colnames(df)[i]
    gene_name <- lookup_df[lookup_df$Antigen.name == antigen_name, type]
    
    if (length(gene_name) > 0) {
      colnames(df)[i] <- gene_name
    }
  }
  
  return(df)
}



# replace_analyte_ID <- function(df, cols, antigens_list, type = "Gene.name") {
#   lookup_df <- do.call(rbind, antigens_list)
#   lookup_df <- distinct(lookup_df, Antigen.name, .keep_all = TRUE)
#   
#   for (col in cols) {
#     print(col)
#     if (col %in% colnames(df)) {
#       for (i in 1:nrow(df)) {
#         antigen_name <- df[i, col]
#         gene_name <- lookup_df[lookup_df$Antigen.name == antigen_name, type]
#         
#         if (length(gene_name) > 0) {
#           df[i, col] <- gene_name
#         }
#       }
#     } else {
#       message(paste("Column", col, "not found in the dataframe. Skipping..."))
#     }
#   }
#   
#   return(df)
# }

strong2 = replace_analyte_ID(strong, list("row", "col"),I00_Antigens)


signif.genes <- reports_P03_Transformed_bcrsn$Set_3_limma
top_genes <- head(signif.genes[order(signif.genes$adj.P.Val),], 25)
top_genes <- fill_analyte_info(top_genes)
top_genes <- top_genes$Gene.name