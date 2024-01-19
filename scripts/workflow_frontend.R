dataset <- "GLA02"
controls <- c("Anti-human IgG", "EBNA1", "Bare-bead", "His6ABP")
for (file in list.files(file.path("scripts", "backend"))) {
  print(file)
  source(file.path("scripts", "backend", file))
}
source(file.path("scripts", "importer.R"))


### Pipeline
import(dataset)
stage_1(I01_Import, "P01_Preprocessed")
stage_2(P01_Preprocessed, "P02_Merged")
stage_3(P02_Merged, "P03_Transformed")


# export_excel(P03_Transformed)
# Report with good names
report <- reports_P03_Transformed_bcrsn$Set_3_limma
anti_names <- row.names((report))
gene_names <- unname(replace_antigen_names(anti_names))
report["Gene_Names"] <- gene_names

# ROC with P.Value 0.1
fancy_roc(P03_Transformed$Set_3, report = reports_P03_Transformed_bcrsn$Set_3_limma, validation = "LOOCV")
roc_breakdown <- roc_curve(P03_Transformed$Set_3, report = reports_P03_Transformed_bcrsn$Set_3_limma, validation = "LOOCV", stepwise = TRUE)
roc_breakdown$roc
roc_breakdown$validation$accuracy
roc_breakdown$validation


# Heatmap
lim <- limma_subset(P03_Transformed$Set_3, mode = "raw")
names <- replace_antigen_names(colnames(lim[-1:-2]))
colnames(lim)[-1:-2] <- unname(names)
pp_heatmap(lim)
# Boxplots
pp_box_plot_multi(lim, "Boxplots")

#PCA
pca_plot(P03_Transformed$Set_3)

# GSEA
sig_adj <- reports_P03_Transformed_bcrsn$Set_3_limma %>%
  filter(adj.P.Val < 0.05)
sig_raw <- reports_P03_Transformed_bcrsn$Set_3_limma %>%
  filter(P.Value < 0.05)
# C2 = msigdb_workflow(reports_P03_Transformed_bcrsn$Set_3_limma, category = "C2")
C2 <- msigdb_workflow(sig_adj, category = "C2")
C2 <- msigdb_workflow(sig_raw, category = "C2")
C3 <- msigdb_workflow(sig_adj, category = "C3")
C3 <- msigdb_workflow(sig_raw, category = "C3")
C2$plot
C3$plot

# OTHER
#pp_box_plot_multi(limma_subset(P03_Transformed$Set_3), "Limma_Subset")
untrans <- untransform_subset(P03_Transformed$Set_3, P02_Merged$Set_3, subset_method = limma_subset, P = 0.05, mode = "default")
#pp_box_plot_multi(untrans, "untrans_limma")
pp_gg_box_plot_multi(untrans, "untrans_limma")
