dataset = "GLA02"
controls = c("Anti-human IgG", "EBNA1", "Bare-bead","His6ABP")
source(paste(getwd(),"/scripts/init3.R", sep=""))
library(limma)
### Workflow
import(dataset)

# Adjust for empty values, remove some noise
D01_Dflist = lapply(I01_Import, initial_clean)
# Generate Merged Set
D02_Dflist = generate_combined_set(D01_Dflist)

Limma_set = D02_Dflist$Set_3
L2 = Limma_set

droplist = c("GLA_02-0054", "GLA_02-0032", "GLA_02-0014", "GLA_02-0057", "GLA_02-0016")
L3 = Limma_set[!(Limma_set$Internal.LIMS.ID %in% droplist),]

# design <- model.matrix(~0 + group, data=Limma_set)
#design <- Limma_set$group
#design <- matrix(as.numeric(design))
#colSums(design)

#pset <- Limma_set[-1:-2]

limma_funct <- function(set) {
  t_set=t(set[-1:-2])
  design <- model.matrix(~0 + group, data=set)
  colnames(design) <- c("case", "control")
  contrasts = makeContrasts(Diff= control - case, levels=design)
  fit<-lmFit(t_set, design, method="robust", maxit=1000)
  contrast_fit <- contrasts.fit(fit,contrasts)
  ebay_fit <- eBayes(contrast_fit)
  DE_results <- topTable(ebay_fit, n=ncol(L2), adjust.method = "fdr", confint = TRUE)
  print(summary(decideTests(ebay_fit)))
  return(DE_results)
}

full_lima = limma_funct(L2)
trim_lima = limma_funct(L3)
limma_funct(D02_Dflist$Set_1)
limma_funct(D02_Dflist$Set_2)
# colnames(design) <- c("case", "control")
# contrasts = makeContrasts(Diff= control - case, levels=design)
# print(head(contrasts))
# fit<-lmFit(L2, design, method="robust", maxit=1000)
# contrast_fit <- contrasts.fit(fit,contrasts)
# ebay_fit <- eBayes(contrast_fit)
# print(summary(decideTests(ebay_fit)))
# print(summary(decideTests(contrast_fit)))
# DE_results <- topTable(ebay_fit, n=ncol(L2), adjust.method = "fdr", confint = TRUE)
# 
# 
# lookatset = cbind(Limma_set[1:2], Limma_set$HPRA015538)
# 
# run_limma <- function(set) {
#   #unnecessary, just following the model given
#   dataset_group = set
#   print(table(dataset_group$group))
#   design <- model.matrix(~as.factor(dataset_group$group))
#   colnames(design)<-c("case", "control")
#   #what to compare
#   contrast <- makeContrasts(Diff = control-case)
# }
# run_limma(Limma_set)
# 
# set = Limma_set
# lmFit()
# design[,2]
