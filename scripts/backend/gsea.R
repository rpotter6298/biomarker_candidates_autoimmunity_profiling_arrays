msigdb_workflow <- function(category = "C2") {
  library(msigdbr)
  library(ggplot2)
  
  msigdb_data <- msigdbr(species = "Homo sapiens", category = category)
  signif.genes <- reports_P03_Transformed_bcrsn$Set_3_limma
  top_genes <- head(signif.genes[order(signif.genes$adj.P.Val),], 25)
  top_genes <- fill_analyte_info(top_genes)
  top_genes <- top_genes$Gene.name
  
  msigdb_genes = select(msigdb_data, gs_name, gene_symbol)
  enrich_msigdb <- enricher(top_genes, TERM2GENE = msigdb_genes)
  
  enrich_msigdb_df <- enrich_msigdb@result %>% 
    separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
    separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"), sep="/") %>% 
    mutate_at(vars("size.term","size.category", "size.overlap.term","size.overlap.category"), as.numeric) %>% 
    mutate("k.K"=size.overlap.term/size.term)
  
  enrich_plot <- enrich_msigdb_df %>% 
    filter(p.adjust <= 0.05) %>% 
    ggplot(aes(x=reorder(Description, k.K), y=k.K)) +
    geom_col() +
    theme_classic() +
    coord_flip() +
    labs(y="Significant genes in set / Total genes in set \nk/K", x="Gene set",
         title = paste("Differentially expressed genes enriched in", category, "Gene sets (P<0.05)"))
  
  return(list(data = enrich_msigdb_df, plot = enrich_plot))
}

