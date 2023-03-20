# library(msigdbr)
# C2 <- msigdbr(species <- "Homo sapiens", category = "C2")
# signif.genes <- reports_P03_Transformed_bcrsn$Set_3_limma
# top_genes <- head(signif.genes[order(signif.genes$adj.P.Val),], 25)
# top_genes <- fill_analyte_info(top_genes)
# top_genes <- top_genes$Gene.name
# 
# 
# C2.genes = select(C2, gs_name, gene_symbol)
# enrich.C2 <- enricher(top_genes, TERM2GENE = C2.genes)
# head(enrich.C2)
# 
# enrich.C2.df <- enrich.C2@result %>% 
#   #separate ratios into 2 columns of data
#   separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
#   separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
#            sep="/") %>% 
#   #convert to numeric
#   mutate_at(vars("size.term","size.category",
#                  "size.overlap.term","size.overlap.category"),
#             as.numeric) %>% 
#   #Calculate k/K
#   mutate("k.K"=size.overlap.term/size.term)
# 
# enrich.C2.df %>% 
#   filter(p.adjust <= 0.05) %>% 
#   #Beautify descriptions by removing _ and HALLMARK
#   #mutate(Description = gsub("HALLMARK_","", Description),
#   #       Description = gsub("_"," ", Description)) %>% 
#   
#   ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
#              y=k.K)) +
#   geom_col() +
#   theme_classic() +
#   #Some more customization to pretty it up
#   #Flip x and y so long labels can be read
#   coord_flip() +
#   #fix labels
#   labs(y="Significant genes in set / Total genes in set \nk/K",
#        x="Gene set",
#        title = "Differentially expressed genes enriched in C2 Gene sets (P<0.05)")
# 


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

