source(file.path("scripts", "workflow_frontend.R"))

#All the top 15s
excel_subset_export(reports_P03_Transformed_bcrsn)
#Set 3 Pictures
pca_plot(P03_Transformed$Set_3)
pp_heatmap(P03_Transformed$Set_3, show_row_names=FALSE)
pp_heatmap(limma_subset(P03_Transformed$Set_3, mode = "top"))
pp_ROC_curve_limma(P03_Transformed$Set_3)
#ROC
source(file.path("scripts", "roc_rough.R"))
#GSEA
source(file.path("scripts", "gsea_rough.R"))
C2 = msigdb_workflow("C2")
C2$plot
df = P03_Transformed$Set_3





