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

pca_plot <- function(data){
  data$group <- factor(data$group)
  pca_data <- prcomp(data[-1:-2], center = TRUE, scale = TRUE)
  plot <- autoplot(pca_data, data=data, colour = "group")
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

pp_heatmap <- function(df, subset = "full"){
  if(subset == "limma"){
    df = limma_subset(df)
  } else if(subset == "compstat"){
    df = compstat_subset(df)
  }
  subset = df[-1:-2]
  subset = apply(subset,2, function(x) (x - mean(x))/sd(x))
  rownames(subset) = paste0(df$group,"-",df$Internal.LIMS.ID)
  # subset = as.matrix(comp_data[comp_data$Internal.LIMS.ID!="GLA_02-0016",])
  subset = subset[,order(colMeans(subset), decreasing = TRUE)]
  # subset2 = subset[,-1]
  # subset = log2(subset)
  subset=t(subset)
  # subset = scale(subset)
  # subset2 = t(subset2)
  pheatmap(subset, scale = "none", cluster_rows = FALSE, cluster_cols = TRUE)
  
  # pheatmap(subset, scale = "column", cluster_rows=FALSE, cluster_cols=TRUE, show_colnames = TRUE)
  # pheatmap(subset, scale = "row", cluster_cols = TRUE, cluster_rows=FALSE, show_colnames = TRUE)
  # pheatmap(subset, scale = "none", show_colnames = TRUE)
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


