source(file.path("scripts", "workflow_frontend.R"))

#All the top 15s
excel_subset_export(reports_P03_Transformed_bcrsn)
#Report with good names
report = reports_P03_Transformed_bcrsn$Set_3_limma
anti_names = row.names((report))
gene_names = unname(replace_antigen_names(anti_names))
report["Gene_Names"] = gene_names
#Set 3 Pictures
pca_plot(P03_Transformed$Set_3, show_ellipse = TRUE)
pp_heatmap(P03_Transformed$Set_3, show_row_names=FALSE)
lim = limma_subset(P03_Transformed$Set_3, mode="raw")
names = replace_antigen_names(colnames(lim[-1:-2]))
colnames(lim)[-1:-2]<- unname(names)
pp_heatmap(lim)
pp_ROC_curve_limma(P03_Transformed$Set_3)


#LR and ROC
logit = roc_curve(P03_Transformed$Set_3, report = reports_P03_Transformed_bcrsn$Set_3_limma, stepwise = FALSE, validation = "LOOCV")
coefficients <- coef(logit$logistic_regression$logistic_model)
equation <- paste("logit(p) = ", round(coefficients[1], 3), "+", paste(round(coefficients[-1], 3), "* x", names(coefficients)[-1], collapse = " + "))
logit_BR = roc_curve(P03_Transformed$Set_3, report = reports_P03_Transformed_bcrsn$Set_3_limma, stepwise = TRUE, validation = "LOOCV")
coefficients_BR <- coef(logit_BR$logistic_regression$logistic_model)
equation_BR <- paste("logit(p) = ", round(coefficients_BR[1], 3), "+", paste(round(coefficients_BR[-1], 3), "* x", names(coefficients_BR)[-1], collapse = " + "))
analyte_names <- names(coefficients)[-1]
unname(replace_antigen_names(analyte_names))

fancy_roc(P03_Transformed$Set_3, report = reports_P03_Transformed_bcrsn$Set_3_limma, validation = "LOOCV")


#ROC
source(file.path("scripts", "roc_rough.R"))
#GSEA
source(file.path("scripts", "gsea_rough.R"))
C2 = msigdb_workflow("C2")
C2$data
C2$plot
C3 = msigdb_workflow("C3")
C3$data
C3$plot
df = P03_Transformed$Set_3


##Table
  difftable <- reports_P03_Transformed_bcrsn$Set_3_limma
  difftable <- difftable [difftable$adj.P.Val <0.1,]
  difftable_compstats <- comparative_statistics(P02_Merged$Set_3)
  difftable_fc <- difftable_compstats[c('log2FC', 'FC')]
  difftable_merge = merge(difftable, difftable_fc, by="row.names")
  rownames(difftable_merge) <- difftable_merge$Row.names
  #difftable_merge$Row.names <- NULL
  difftable_merge$Gene.names <- replace_antigen_names(difftable_merge$Row.names)
  # Create a new workbook
  wb <- createWorkbook()
  addWorksheet(wb, "Difftable_Merge")
  writeData(wb, "Difftable_Merge", difftable_merge)
  saveWorkbook(wb, file="stats/difftable_merge.xlsx")
  
#Boxplots
  pp_box_plot_multi(lim, "testplots")
  
  
#Pearson Correlations:
diff_results_filtered <- subset(reports_P03_Transformed_bcrsn$Set_3_limma, P.Value<0.05)
expr_matrix <- diff_results_filtered %>%
  select(-P.Value, -logFC) %>%
  pivot_wider(names_from = group, values_from = normalized_expression)

cor_matrix <- cor(expr_matrix, method = "pearson")

top_genes <- diff_results
logFC <- diff_results$logFC
expr_matrix <- t(P03_Transformed$Set_3[-1])
expr_matrix <- as.matrix(expr_matrix[row.names(top_genes),])
expr_matrix <- apply(expr_matrix, 2, as.numeric)
# compute correlation coefficients for top genes
cor_mat <- cor(t(expr_matrix))
row.names(cor_mat) = replace_antigen_names(row.names(top_genes))
colnames(cor_mat) = replace_antigen_names(row.names(top_genes))

#Bar Chart
library(reshape2)
library(plyr)
limmasubset = limma_subset(P03_Transformed$Set_3)
crosslima_subset = P02_Merged$Set_3[colnames(limmasubset)]
names = colnames(crosslima_subset[-1:-2])
newnames = unname(replace_antigen_names(names))
colnames(crosslima_subset)[-1:-2] = newnames 
bar_plot(crosslima_subset)
bar_plot <- function(set){
  melted = melt(set)
  melted = na.omit(melted)
  means <- ddply(melted, c("group", "variable"), summarise,
                 mean=mean(value))
  colnames(melted)[1] = "ID"
  melted$ID = as.factor(melted$ID)
  ggplot(data=means,aes(x=variable, y=mean, fill=factor(group, levels = c(0,1), labels = c("Healthy", "Diseased"))))+
    geom_col(position=position_dodge()) + theme(axis.text.x=element_text(angle = 270, vjust = 0.5, hjust = 0)) + 
    ggtitle( "Significantly Differentially Expressed Antigens") +
    scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
    xlab("Gene Name") + ylab("Adjusted Expression Level") +
    labs(fill = "Group")
}
