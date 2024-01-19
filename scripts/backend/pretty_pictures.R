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
  
  data$group <- factor(data$group, levels = c(0, 1), labels = c("Healthy", "Diseased"))
  pca_data <- prcomp(data[-1:-2], center = TRUE, scale = TRUE)
  
  # Extract PCA scores
  pca_scores <- as.data.frame(pca_data$x)
  pca_scores$group <- data$group
  # Calculate percentage of variance explained by each PC
  var_exp <- round(pca_data$sdev^2 / sum(pca_data$sdev^2) * 100, 2)
  
  # Update axis labels with percentage of variance explained
  x_label <- paste0("PC1 (", var_exp[1], "%)")
  y_label <- paste0("PC2 (", var_exp[2], "%)")
  
  plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
    geom_point() +
    theme_classic() +
    labs(x = x_label, y = y_label, title = "PCA Plot") + 
    scale_color_manual(values = c("Healthy" = "blue", "Diseased"="red"))
  
  if (show_ellipse) {
    plot <- plot + stat_ellipse(aes(fill = group), geom = "polygon", level = 0.95, alpha = 0.2) +
      labs(title = "") + 
      scale_fill_manual(values = c("Healthy" = "palegreen", "Diseased"="palegoldenrod"))
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
      plot = boxplot(data[,i] ~ data$group, main = colnames(data)[i], xlab = "Group", ylab="")
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
        plot = boxplot(data[,i] ~ data$group, main = colnames(data)[i], xlab = "Group", ylab="")
      }
      dev.off()
    }
  }
}

library(ggplot2)
library(ggbreak)

pp_gg_box_plot_multi <- function(data, name){
  
  data <- data[,-1]
  # Transform data to long format
  data_long <- tidyr::pivot_longer(data, -group, names_to="Variable", values_to="Value")
  
  # Calculate upper limit for normal scale before break
  upper_limit <- quantile(data_long$Value, 0.95)
  
  # Create the plot
  p <- ggplot(data_long, aes(x=group, y=Value)) + 
    geom_boxplot(aes(group=group, fill=factor(group)), outlier.shape = NA, color="black") +  # Set the boxes to neutral color and outline them in black
    geom_point(aes(color=factor(group)), position = position_jitter(width = 0.3), alpha=0.7) +  # Keep the points colored
    facet_wrap(~Variable, scales="free_y") + 
    coord_trans(y="log10") +  # Implement log transform
    theme_bw() + 
    scale_fill_manual(values = c("white", "white"), guide=FALSE) +  # Use white color for box fill and disable its legend
    scale_color_manual(values = c("#E57373", "#4DB6AC"), 
                       labels = c("Healthy Controls", "XFG Patients"),
                       name = "") +  # Return to the preferred colors
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12, face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = c(0.9, 0.1),
          legend.justification = c(1, 0),
          legend.background = element_blank(),
          legend.key = element_blank())
  
  # Save the plot to a file
  ggsave(filename = paste0(name, ".png"), plot = p, width = 10, height = 6)
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
  rownames(subset) = paste0(ifelse(df$group == 0, "Healthy", "Diseased"), "-", substr(df$Internal.LIMS.ID, start = nchar(df$Internal.LIMS.ID)-3, stop = nchar(df$Internal.LIMS.ID)))
  # Order the columns of the 'subset' dataframe by decreasing column means
  subset = subset[, order(colMeans(subset), decreasing = TRUE)]
  # Transpose the 'subset' dataframe
  subset = t(subset)
  
  pheatmap(subset, 
           scale = "none", 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           show_rownames = show_row_names, 
           show_colnames = show_col_names, 
           treeheight_row = 0,
           clustering_distance_cols = "euclidean", 
           clustering_distance_rows = "euclidean",
           clustering_method = "complete")
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


  