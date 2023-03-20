library(tidyverse)
library(msigdbr)
library(clusterProfiler)

C2 <- msigdbr(species = "Homo sapiens", category = "C2")
class(H)

signif.genes = unique(sumdoc$Set_3$Gene)
signif.genes = unique(settrim_highlightsq$Gene.name)
signif.genes = c("DGCR2", "LGSN", "FUT2", "LOXL3", "CLTCL1", "SCG5","CDH5")
C2.genes = select(C2, gs_name, gene_symbol)

enrich.C2 <- enricher(signif.genes, TERM2GENE = C2.genes)

head(enrich.C2@result)

enrich.C2.df <- enrich.C2@result %>% 
  #separate ratios into 2 columns of data
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 
  #Calculate k/K
  mutate("k.K"=size.overlap.term/size.term)


enrich.C2.df %>% 
  filter(p.adjust <= 0.05) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  #mutate(Description = gsub("HALLMARK_","", Description),
  #       Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = "Differentially expressed genes enriched in C2 Gene sets (P<0.05)")





C7 <- msigdbr(species = "Homo sapiens", category = "C7")
C7.genes = select(C7, gs_name, gene_symbol)
enrich.C7 <- enricher(signif.genes, TERM2GENE = C7.genes)
enrich.C7.df <- enrich.C7@result %>% 
  #separate ratios into 2 columns of data
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 
  #Calculate k/K
  mutate("k.K"=size.overlap.term/size.term)
head(enrich.C7@result)

enrich.C7.df %>% 
  filter(p.adjust <= 0.05) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  #mutate(Description = gsub("HALLMARK_","", Description),
  #       Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = "Differentially expressed genes enriched in C7 Gene sets (P<0.05)")


C3 <- msigdbr(species = "Homo sapiens", category = "C3")
C3.genes = select(C3, gs_name, gene_symbol)
enrich.C3 <- enricher(signif.genes, TERM2GENE = C3.genes)
enrich.C3.df <- enrich.C3@result %>% 
  #separate ratios into 2 columns of data
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 
  #Calculate k/K
  mutate("k.K"=size.overlap.term/size.term)
head(enrich.C3@result)

enrich.C3.df %>% 
  filter(p.adjust <= 0.05) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  #mutate(Description = gsub("HALLMARK_","", Description),
  #       Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = "Differentially expressed genes enriched in C3 Gene sets (P<0.05)")
