test_calc = P02_Merged$Set_3
reporting = comparative_statistics(test_calc)
reporting = fill_analyte_info(reporting)

analyte ="HPRA006876"
X = mean(test_calc[analyte][test_calc$group==0,])
Y = mean(test_calc[analyte][test_calc$group==1,])


fold.change <- function(x,y){
  fc = y/x
  return(fc)
}

fold.change(X, Y)


limma_test = limma_funct(P02_Merged$Set_3)
comp_test = comparative_statistics((P02_Merged$Set_3))


subset_top_n(reports_P03_Transformed_bcrsn$Set_3_limma, 15)
plotset = limma_subset(P03_Transformed$Set_3)
plots = pp_box_plot_multi(plotset, file.path("plots", "checkplots")

                          
                          #visualizations of limma results

limma_results <- limma_funct(P02_Merged$Set_3)
limma_results <- limma_funct(P03_Transformed$Set_3)

#volcano
threshold <- 0.05
limma_results$significant <- limma_results$adj.P.Val < threshold
ggplot(limma_results, aes(x = logFC, y= -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(x = "log2 Fold Change", y="-log10(Adjusted P-Value")+
  scale_color_manual(values = c("black", "red"))

#heatmap of FC values
top_genes <- limma_results[order(abs(limma_results$logFC), decreasing=TRUE)[1:50],]
logFC_matrix <- matrix(top_genes$logFC, nrow=1, ncol = length(top_genes$logFC), byrow = TRUE)
rownames(logFC_matrix) <- "logFC"
colnames(logFC_matrix) <- row.names(top_genes)

pheatmap(logFC_matrix, cluster_rows = FALSE, cluster_cols = FALSE)


#heatmap of normalized expression
