dataset = "GLA02"
controls = c("Anti-human IgG", "EBNA1", "Bare-bead","His6ABP")
source(paste(getwd(),"/scripts/workflow_backend_1.R", sep=""))

### Workflow
import(dataset)

# Adjust for empty values, remove some noise
D01_Dflist = lapply(I01_Import, initial_clean)
# Generate Merged Set
D02_Dflist = generate_combined_set(D01_Dflist)
### Generate Set 3 From combined Set 1 & 2
#D01_Dflist = generate_combined_set(I01_Import)

export_excel(D02_Dflist)
### Decide cutoff value
cutoff_val = 0

### Remove empty, remove variables with no reading above cutoff
#D02_dflist = lapply(D01_Dflist, initial_clean)

pp_boxplot_dflist(D02_Dflist)

#Apply Fisher's Test
cutoff_val = 250
FishReport = lapply(D02_Dflist, test_fisher)
Filtered_Fish = fish_filter(D02_Dflist)
trespescados = Filtered_Fish[[3]]

## This mess does PCA and histograms
set3 = D02_Dflist$Set_3
droplist = c("GLA_02-0054", "GLA_02-0032", "GLA_02-0014", "GLA_02-0057", "GLA_02-0016")
dropset = set3[(set3$Internal.LIMS.ID %in% droplist),]
dropset = cbind(dropset[1:2],dropset[row.names(set3_highlights)])
set3 = set3[!(set3$Internal.LIMS.ID %in% droplist),]
#row.names(set3)
set3_comp = comparative_statistics(set3)
set3_highlights = set3_comp[set3_comp$P_Value<0.05,]
set3_highlights = set3_highlights[order(set3_highlights$P_Value, decreasing =  TRUE),]
set3_trim = cbind(set3[1:2],set3[row.names(set3_highlights)])

hist(set3_comp$P_Val, breaks=20, labels=TRUE, xlab= "P-Value", main="Histogram of P-Values (Set 3)")
axis(side=1, at=seq(0,.2,0.1), labels=TRUE)
axis(side=1, at=seq(0,.4,0.05), labels=FALSE)

hist(set3_comp$Qval, breaks = 20, xlab = "Q Value", labels=TRUE, main="Histogram of Q Values (Set 3)")

set3_full.pca = prcomp(set3[-1:-2])
set3.pca=prcomp(set3_trim[-1:-2])
autoplot(set3.pca, data=set3_trim, colour="group", label=TRUE)
autoplot(set3_full.pca, data=set3, colour="group", label=TRUE)


dropset[-1:-2] = dropset[-1:-2][order(colMeans(dropset[-1:-2]))]
dropset_plot = melt(dropset)
repset = set3_trim[07,]
repset = melt(repset)
ggplot(data = dropset_plot, aes(x=variable, y=value, group = 1)) +
  geom_jitter(aes(color=factor(Internal.LIMS.ID))) + 
  geom_col(data = repset, aes(x= variable, y=value)) +
  labs(alpha = NULL, y="MFI Level", x="Analyte", colour = "Patient ID") +
  scale_color_manual(labels = c("0014", "0016", "0032", "0054", "0057", "Normal Healthy Representative"), values = c("Red", "Blue", "gold2", "seagreen", "Purple4")) +
  scale_y_continuous(trans='log2') +
  theme(axis.text.x=element_blank()) 


a = colnames(set3_trim)
b = colnames(trespescados)

a[a %in% b]
unique(c(a, b))


doubletrim = set3[unique(c(colnames(set3_trim), colnames(trespescados)))]
doubletrim_log = cbind(doubletrim[1:2], log(doubletrim[-1:-2], 2)) 
scatter_plot(doubletrim_log)

set = set3_trim
set = doubletrim
set$HPRA045601=NULL
set$HPRA005621[29]=mean(set$HPRA005621[set$group==1][set$HPRA005621 != set$HPRA005621[29]])

subset = as.matrix(set[-1:-2])
rownames(subset) = paste0(set$group,"-",set$Internal.LIMS.ID)
subset = subset[,order(colMeans(subset), decreasing = TRUE)]
subset = subset[,-1]
subset = log2(subset)
subset = scale(subset)
subset = t(subset)
pheatmap(subset, cluster_rows = FALSE, cluster_cols = TRUE)

pheatmap(subset, scale = "column", cluster_rows=FALSE, cluster_cols=TRUE, show_colnames = TRUE)
pheatmap(subset, scale = "row", cluster_cols = TRUE, cluster_rows=FALSE, show_colnames = TRUE)
pheatmap(subset, scale = "none", show_colnames = TRUE)


### Okay, back to the organized stuff

D03_Dflist <- define_doubletrim(D02_Dflist)
names(D03_Dflist) = names(D02_Dflist)
sumdoc = generate_summary_doc(D03_Dflist, I00_Antigens$SBA01_Antigen_list.xlsx)
names(sumdoc) = names(D03_Dflist)
export_excel(sumdoc)

### For GSEA
full_genes = unique(I00_Antigens$SBA01_Antigen_list.xlsx$Gene.name)
dput(full_genes, "fullgenes.txt")
set3_comp = set3_comp[order(set3_comp$Log2FC),]
row.names(set3_comp)[row.names(set3_comp) %in% antigens$Antigen.name]

for (rname in row.names(set3_comp)[row.names(set3_comp) %in% antigens$Antigen.name]){
  set3_comp[rname, "Gene.name"] = antigens[antigens$Antigen.name == rname, "Gene.name"]
}
ordered_genes = set3_comp$Gene.name
o2 = unlist(strsplit(ordered_genes, ","))
sink("orderedgenes.txt")
writeLines(unlist(lapply(o2, paste, collapse=" ")))
sink()

### Return To Order
mean(set3$HPRA045605[set$group==0])
mean(set3$HPRA045605[set$group==1])

## Pretty Pictures
FF01 = Filtered_Fish[[1]]
FF02 = Filtered_Fish[[2]]
FF03 = Filtered_Fish[[3]]
pp_analyte_plots(FF01)
pp_boxplot_set(set3_trim)
pp_boxplot_set(trespescados)
pp_boxplot_set(doubletrim)
pp_boxplot_set(D02_Dflist$Set_3)
bar_plot(set3_trim)
bar_plot(trespescados)

mean(trespescados$HPRA045184)
mean(trespescados$HPRA045184[trespescados$group==1])
scatter_plot(trespescados)
scatter_plot(D02_Dflist$Set_3)

pp_analyte_plots(set3_trim)
pp_analyte_plots(trespescados)
pp_analyte_plots(doubletrim)
scatter_plot(doubletrim_log)
pp_analyte_plots(D01_Dflist$Set_3)
pp_analyte_plots(set3)
pp_analyte_plots(D03_Dflist$Set_3$doubletrim)


set = D01_Dflist$`SBA02_Data Intensity.xlsx`
set

emptyset = set[set$Internal.LIMS.ID=="EMPTY-0001",]
inset = emptyset[-1:-2][colMeans(emptyset[-1:-2])<(median(colMeans(emptyset[-1:-2]))+(2*sd(colMeans(emptyset[-1:-2]))))]
if (length(inset) < 0.95*length(emptyset)){
  calset = emptyset[-1:-2]
}else{
  calset = inset
}
max(colMeans(calset))+2*sd(colMeans(calset))
mean(colMeans(calset))
i = 4
emptyset
max(colMeans(emptyset[-1:-2]))
sd(colMeans(emptyset[-1:-2]))
colnames(set[i])
antigens= I00_Antigens$SBA02_Antigen_list.xlsx
antigens[antigens$analyte == analytenum,]$Antigen.name
gsub("^.*?_","_","ATGAS_1121")
as.numeric(strsplit(colnames(set[i]), split='.', fixed = TRUE)[[1]][2])
set[,!(names(set) %in% controls)]


set = I01_Import$`SBA01_Data Intensity.xlsx`
compress_duplicates(set, layout="sample_data/03. AP0211 GLA02_SampleLayout.xlsx")
antigens = I00_Antigens$SBA01_Antigen_list.xlsx

sumdoc$Set_3$pval
pvals = sumdoc$Set_3$pval
p.adjust(pvals)
pvals
newq
length(set3[set3>0])
length(set3[set$group==0,][set3>1000])
length(set3[set$group==1,][set3>1000])
set3c = set3[set3$group==0,]
set3e = set3[set3$group==1,]
length(set3c[set3c>2500])
length(set3e[set3e>2500])

set = D02_Dflist$Set_3
settrim = set[-1:-2]
settrim[settrim > 2500] = NA
settrim = cbind(set[1:2], settrim)
setlog = cbind(set[1:2],log2(set[-1:-2]))
setlogsig = sig_test(settrim)
setlogcomp = comparative_statistics(settrim)
settrim_highlights = setlogcomp[setlogcomp$P_Value<0.05,]
settrim_highlights = settrim_highlights[order(settrim_highlights$P_Value, decreasing =  TRUE),]
set_trim = cbind(settrim[1:2],settrim[row.names(settrim_highlights)])
antigens = I00_Antigens$SBA01_Antigen_list.xlsx

rownames(settrim_highlights)
for (row in rownames(settrim_highlights)){
  print(antigens$Gene.name[antigens$Antigen.name==row])
  settrim_highlights[row,"Gene.name"] = antigens$Gene.name[antigens$Antigen.name==row]
}
settrim_highlightsq = settrim_highlights[settrim_highlights$Qval<0.05,]


set_trim_plugged = plug_holes(set_trim)
settrim.pca=prcomp(set_trim_plugged[-1:-2])
autoplot(settrim.pca, data=set_trim_plugged, colour="group", label=FALSE)

subset = as.matrix(set_trim[-1:-2])
rownames(subset) = set$group
subset = subset[,order(colMeans(subset), decreasing = TRUE)]
subset = subset[,-1]
subset = log2(subset)
subset = scale(subset)
subset = t(subset)
pheatmap(subset, cluster_rows = FALSE, cluster_cols = TRUE)

pheatmap(subset, scale = "column", cluster_rows=FALSE, cluster_cols=TRUE, show_colnames = TRUE)
pheatmap(subset, scale = "row", cluster_cols = TRUE, cluster_rows=FALSE, show_colnames = TRUE)
pheatmap(subset, scale = "none", show_colnames = TRUE)


dropset = set_trim
dropset[-1:-2] = set_trim[-1:-2][order(colMeans(set_trim[-1:-2]))]
colMeans(dropset[-1:-2])
dropset_plot = melt(dropset)
repset = set_trim[13,]
repset = melt(repset)
ggplot(data = dropset_plot, aes(x=variable, y=value, group = 1)) +
  geom_jitter(aes(color=factor(group))) + 
  geom_col(data = repset, aes(x= variable, y=value, alpha=0.2)) +
  scale_y_continuous(limits=c(0,300)) +
  labs(alpha = "Representative") +
  theme(axis.text.x=element_text(angle=270)) 


scatter_plot(set_trim)
set=set_trim
pp_boxplot_set(set_trim)
bar_plot(set_trim)
pp_analyte_plots(set_trim)

scatter_plot(set_trim)
junk = c(5,5,5,5,40,60,40,60, 1200)
summary(junk)
sd(junk)
t.test(junk)
