pp_dflist_wrapper <- function(dflist, pp_function){
  dflist_name = deparse(substitute(dflist))
  plot_name = paste(strsplit(deparse(substitute(pp_function)), "_")[[1]][-1], collapse = "")
  dir_path = file.path(getwd(),"plots",dflist_name,plot_name)
  dir.create(dir_path, recursive=TRUE, showWarnings = FALSE)
  for (setid in seq_along(dflist)){
    set_name = paste("Set",setid,sep="_")
    name = (file.path(dir_path,set_name))
    pp_function(dflist[[setid]], name)
  }
}

pca_plot <- function(data, show_ellipse = TRUE) {
  library(ggplot2)
  library(ggfortify)
  
  data$group <- factor(data$group)
  pca_data <- prcomp(data[-1:-2], center = TRUE, scale = TRUE)
  
  # Extract PCA scores
  pca_scores <- as.data.frame(pca_data$x)
  pca_scores$group <- data$group
  
  plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
    geom_point() +
    theme_classic() +
    labs(x = "PC1", y = "PC2", title = "PCA Plot")
  
  if (show_ellipse) {
    plot <- plot + stat_ellipse(aes(fill = group), geom = "polygon", level = 0.95, alpha = 0.2) +
      labs(title = "PCA Plot with 95% Confidence Ellipses")
  }
  
  return(plot)
}



pp_box_plot_multi <- function(data, name){
  # Set up plot area and device
  dim = display_division(ncol(data)-2)
  if ((ncol(data)-2)<=(dim[1]*dim[2])){
    filename = name
    col_range = 3:ncol(data)
    png(filename = paste0(filename,".png"), width = 1200+300*dim[1], height = 900+100*dim[2], res=250)
    par(mfrow = c(dim[2],dim[1]))
    for (i in col_range) {
      plot = boxplot(data[,i] ~ data$group, main = colnames(data)[i], xlab = "Group")
    }
    dev.off()
  }
  else{
    modifier = dim[1]*dim[2]
    for (d in seq(dim[3])){
      filename = paste(name,d,sep="_")
      png(filename = paste0(filename,".png"), width = 1200+300*dim[1], height = 900+300*dim[2], res=250)
      col_range = 3:(modifier+2)+(modifier*(d-1))
      col_range = col_range[col_range<=ncol(data)]
      par(mfrow = c(dim[2],dim[1]))
      print(col_range)
      for (i in col_range) {
        plot = boxplot(data[,i] ~ data$group, main = colnames(data)[i], xlab = "Group")
      }
      dev.off()
    }
  }
}

pp_heatmap <- function(df, subset = "full", show_row_names = TRUE, show_col_names = TRUE){
  # Apply the appropriate subset function based on the 'subset' parameter
  if(subset == "limma"){
    df = limma_subset(df)
  } else if(subset == "compstat"){
    df = compstat_subset(df)
  }
  # Remove the first two columns to create the 'subset' dataframe
  subset = df[-1:-2]
  # Standardize the columns of the 'subset' dataframe
  subset = apply(subset, 2, function(x) (x - mean(x)) / sd(x))
  # Create row names for the 'subset' dataframe based on the 'group' and 'Internal.LIMS.ID' columns
  rownames(subset) = paste0(df$group, "-", df$Internal.LIMS.ID)
  # Order the columns of the 'subset' dataframe by decreasing column means
  subset = subset[, order(colMeans(subset), decreasing = TRUE)]
  # Transpose the 'subset' dataframe
  subset = t(subset)
  # # Create empty annotations for rows and columns if labels should be hidden
  # annotation_row <- if (!show_row_labels) data.frame(labels = factor(rep("", nrow(subset))), stringsAsFactors = FALSE) else NULL
  # annotation_col <- if (!show_col_labels) data.frame(labels = factor(rep("", ncol(subset))), stringsAsFactors = FALSE) else NULL
  
  pheatmap(subset, scale = "none", cluster_rows = FALSE, cluster_cols = TRUE, show_rownames = show_row_names, show_colnames = show_col_names)
}

pp_multipca <- function(df, name){
  datalist <- list(base_data=df,limma_data=limma_subset(df),compstat_data = compstat_subset(df))
  plots <- lapply(names(datalist), function(name) {
    plot <- pca_plot(datalist[[name]])
    plot + ggtitle(name)
  })
  combined_plot <- do.call(grid.arrange, c(plots, ncol = 3))
  ggsave(paste(name, "multipca.png", sep="_"), combined_plot, width = 12, height = 4, dpi = 300)
  }


  