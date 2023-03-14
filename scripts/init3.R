library(xlsx)
library(readxl)
library(ggplot2)
library(ggfortify)
library(reshape2)
library(plyr)
library(dplyr)
library(broom)
library(qvalue)
library(pheatmap)



### Workflow
import <- function(dataset){
  sample_data=paste(getwd(),"/sample_data/", sep="")
  names = list.files(path=sample_data, pattern = ".xlsx", recursive=TRUE)
  for (file in names){
    shortindex = gregexpr(pattern=dataset,file)[[1]][1]
    shortname = substring(file,shortindex+6)
    assign(paste0(shortname), read_xlsx(paste(sample_data,file,sep="")))
  }
  assign(paste(dataset, "_A00", sep=""),filterclean("_Data"))
  assign(paste(dataset, "_Antigens", sep=""), filterclean("_Antigen"))
  
  A01_Input = lapply(get(paste(dataset, "_A00", sep="")),function(df){
    df = data.frame(df)
    names(df)[2] = "group" 
    df$group =ifelse(grepl("GC",df$group), 1, 
                     ifelse(grepl("HD",df$group), 0, "NA"))
    df
  })
  AA_Antigens = lapply(get(paste(dataset, "_Antigens", sep="")),function(df){
    df = data.frame(df)
    names(df)[1] = "analyte"
    df
  })
  assign("I01_Import", A01_Input, env=globalenv())
  assign("I00_Antigens", AA_Antigens, env=globalenv())
}
generate_combined_set <- function(dflist){
  dflist = set_colname_adapter(dflist)
  dflist = lapply(dflist, remove_controls)
  set3 = mergedown(dflist)
  dflist = c(dflist, list(set3))
  nameslist = namelist(dflist)
  names(dflist) = nameslist
  return(dflist)
}
initial_clean <- function(set, cutoff=NULL, layout = "sample_data/layout.xlsx"){
  if (is.null(cutoff)){
    cutoff = detect_background_noise(set)
  }
  print(cutoff)
  set = purgebackground(set, cutoff)
  set = compress_duplicates(set, layout)
  set = emptyadjust(set)
  return(set)
}
fish_filter <-function(dflist){
  lapply(1:length(dflist), function(setid){
    set = dflist[[setid]]
    fishreport = test_fisher(set)
    fishlist = row.names(fishreport)[fishreport$`Fisher's P Val`<0.05]
    filtered = cbind(set[1:2], set[fishlist])
    return(filtered)
  })
}

define_doubletrim <- function(dflist){
  Filtered_Fish = fish_filter(dflist)
  lapply(1:length(dflist), function(setid){
    set = dflist[[setid]]
    dropset = set[(set$Internal.LIMS.ID %in% droplist),]
    set = set[!(set$Internal.LIMS.ID %in% droplist),]
    set_comp = comparative_statistics(set)
    set_highlights = set_comp[set_comp$P_Value<0.05,]
    set_highlights = set_highlights[order(set_highlights$P_Value, decreasing =  TRUE),]
    set_trim = cbind(set[1:2],set[row.names(set_highlights)])
    doubletrim = set[unique(c(colnames(set_trim), colnames(Filtered_Fish[[setid]])))]
    # doubletrim = set3[unique(c(colnames(set_trim), colnames(trespescados)))]
    return (list(doubletrim = doubletrim, dropset = dropset, set_comp= set_comp))
  })
}

generate_summary_doc <- function(dflist, antigens){
  lapply(1:length(dflist), function(setid){
    set = dflist[[setid]][[1]]
    set_comp = dflist [[setid]][[3]]
    set_comp = subset(set_comp, P_Value < 0.05)
    fishreport = test_fisher(set)
    fishreport = subset(fishreport, `Fisher's P Val` < 0.05 )
    summarydoc = data.frame(row.names = colnames(set[-1:-2]))
    for (rname in row.names(summarydoc)[row.names(summarydoc) %in% row.names(fishreport)]){
      summarydoc[rname, "pval"] = fishreport[rname, 1]
      summarydoc[rname, "Gene"] = antigens[antigens$Antigen.name == rname, "Gene.name"]
    }
    for (rname in row.names(summarydoc)[row.names(summarydoc) %in% row.names(set_comp)]){
      summarydoc[rname, "pval"] = set_comp[rname, "P_Value"]
      summarydoc[rname, "qval"] = set_comp[rname, "Qval"]
      summarydoc[rname, "FC"] = set_comp[rname, "FC"]
      summarydoc[rname, "Log2FC"] = set_comp[rname, "Log2FC"]
      summarydoc[rname, "Gene"] = antigens[antigens$Antigen.name == rname, "Gene.name"]
    }
    return(summarydoc)
  })
}

### Corrections ###
emptyadjust <- function(set){
  emptyset = set[set$Internal.LIMS.ID=="EMPTY-0001",]
  emptyvector = colMeans(emptyset[-1:-2])
  set = set[set$Internal.LIMS.ID!="EMPTY-0001" & set$Internal.LIMS.ID!="MIX_2-0029",]
  labels = set[1:2]
  set = cbind(set[1:2],sweep(set[-1:-2],2,FUN="-",emptyvector))
  ### This sets anything with less read than the empty to 0, then adds one to everything, so that it doesn't break the log()
  set[-1:-2][set[-1:-2]<0] <- 0
  set[-1:-2] = set[-1:-2]+1
  return(set)
}
purgebackground <- function(set, cutoff){
  bgmap = set[-1:-2]> cutoff
  fset = set[-1:-2][,colSums(bgmap)>1]
  set = cbind(set[1:2],fset)
  return(set)
}

compress_duplicates <- function(set, layout){
  outlist = c()
  ### Reads the layout from excel file
  slayout <- read_excel(layout)
  ## Looks through the column named Tube Label for any items containing a hyphen, then uses any name preceding a hyphen as the for item
  for (n in strsplit(slayout$`Tube label`[grepl("-",slayout$`Tube label`)],"-")){
    # Filters the layout to only hold items either ending in or containing the for item immediately before the hyphen
    mergerows = ((slayout %>% filter(grepl(paste0(n[1],"$|",n[1],"-"),`Tube label`)))$`Sample id_LIMS`)
    # Returns the sample id for those samples, which is present in the main data set
    mergeset = set[grepl(paste(mergerows, collapse = "|"), set$Internal.LIMS.ID),]
    # Checks that both technical replicates are close enough together that their deviation is not greater than either item (like if one was 1000 and one was 10, this would print notifications)
    outlierflag = colSums(sweep(mergeset[-1:-2],2,apply(mergeset[-1:-2], 2, sd ), '-')<0)
    if (sum(outlierflag)>0){
      print(mergerows)
      print(outlierflag[outlierflag>0])
      #tupsum = c(mergerows, colnames(outlierflag[outlierflag>0]))
      #print(tupsum)
    }
    # Reassigns the first row in the set of matches so that it is equal to the mean of all matching sets, for each variable
    set[grepl(paste(mergerows, collapse = "|"), set$Internal.LIMS.ID),][1,][-1:-2] = colMeans(set[grepl(paste(mergerows, collapse = "|"), set$Internal.LIMS.ID),][-1:-2])
    # Eliminates all matching samples except the first (which is now the average of all matching) from the dataset
    for (i in 2:length(mergerows)){
      set = set[row.names(set) != row.names(set[grepl(paste(mergerows, collapse = "|"), set$Internal.LIMS.ID),][i,]),]
    }}
  return(set)
}

remove_controls <- function(set, controls=get("controls", globalenv())){
  set = set[,!(names(set) %in% controls)]
}

### Tests
test_fisher <- function(set, cutoff=cutoff_val){
  #name = deparse(substitute(set))
  #path = pathing(name,"mosaic")
  fishreport = data.frame(row.names = colnames(set[-1:-2]))
  set0 = set[set$group==0,][-1:-2]>cutoff
  set1 = set[set$group==1,][-1:-2]>cutoff
  for (i in colnames(set[-1:-2])){
    fishstack = data.frame(
      "healthy" = c(sum(set0[,i]), length(set0[,i])-sum(set0[,i])),
      "disease" = c(sum(set1[,i]), length(set1[,i])-sum(set1[,i])), 
      row.names = c("Above Cutoff", "Below Cutoff")
    )
    #mosaicplot(fishstack, main = i, color = TRUE)
    #ggsave(paste0(path, i, ".png"))
    #fishreport[i,"P_Val"] = fisher.test(fishstack)$p.value
    fishreport[i,"Fisher's P Val"] = fisher.test(fishstack)$p.value
  }
  #fishreport[fishreport$`Fisher's P Val` >= 1,] = 1
  #qobj <- qvalue(p = fishreport$`Fisher's P Val`)
  #fishreport["Q_val"] = qobj$qvalues
  return(fishreport)
}
comparative_statistics <- function(set, mode=2){
  mode=mode
  sigtest = sig_test(set)
  sigtest$bh_p.value = p.adjust(sigtest$p.value, method="BH")
  qobj <- qvalue(p = sigtest$p.value)
  #qobj <- qvalue(p = tlist$p.value, lambda=0)
  sigtest$q.value <- qobj$qvalues
  setC = set[set$group == 0,]
  setE = set[set$group == 1,]
  comp = data.frame(matrix(nrow = ncol(set)-2, ncol = 8))
  rownames(comp) <- colnames(set[-1:-2])
  colnames(comp) = c( "FC", "Log2FC", "P_Value", "Method", "BH_P-value", "Qval", "ShapP-E", "ShapP-C")
  for (i in 3:ncol(set)){
    #print(names(set[i]))
    FC = mean(setE[[i]], na.rm = TRUE)/mean(setC[[i]], na.rm=TRUE)
    #tlist$p.value[tlist$analyte==colnames(set)[i]]
    #tlist$analyte
    #colnames(set)[i]
    #tt = tlist[tlist$analyte==colnames(set)[i],]
    sigtest$p.value
    Pval = sigtest[sigtest$analyte==colnames(set)[i],]$p.value
    Method = sigtest[sigtest$analyte==colnames(set)[i],]$method
    BH_Pval = sigtest[sigtest$analyte==colnames(set)[i],]$bh_p.value
    Qval = sigtest[sigtest$analyte==colnames(set)[i],]$q.value
    ShapE = sigtest[sigtest$analyte==colnames(set)[i],]$'ShapE$p.value'
    ShapC = sigtest[sigtest$analyte==colnames(set)[i],]$'ShapC$p.value'
    vec <- list(FC, log2(FC),Pval,Method,BH_Pval, Qval, ShapE, ShapC)
    comp[names(set[i]),] <- vec
  }
  return(comp)
}

sig_test <- function(set){
  sigtestlist = data.frame()
  for (i in colnames(set[-1:-2])){
    setC = set[set$group==0,][[i]]
    setE = set[set$group==1,][[i]]
    ShapE = shapiro.test(setE)
    ShapC = shapiro.test(setC)
    if (ShapE$p.value >0.05 & ShapC$p.value >0.05){
      testrow = cbind(tidy(t.test(setC, setE))[c("statistic", "p.value", "method")],ShapC$p.value, ShapE$p.value)
    } else {
      testrow = cbind(tidy(wilcox.test(setC, setE))[c("statistic", "p.value", "method")],ShapC$p.value, ShapE$p.value)
    }
    testrow$analyte = i
    sigtestlist = rbind(sigtestlist, testrow)
  }
  return(sigtestlist)
}
detect_background_noise <- function(set){
  emptyset = set[set$Internal.LIMS.ID=="EMPTY-0001",]
  inset = emptyset[-1:-2][colMeans(emptyset[-1:-2])<(median(colMeans(emptyset[-1:-2]))+(1*sd(colMeans(emptyset[-1:-2]))))]
  if (length(inset) < 0.95*length(emptyset)){
    calset = emptyset[-1:-2]
  }else{
    calset = inset
  }
  return(max(colMeans(calset))+sd(colMeans(calset)))
}

### Utility Functions ###
filterclean <- function(codeword, e = parent.frame()) {
  dflist = Filter(function(x) is (x, "data.frame"),
                  mget(ls(e),envir= e))
  dflist = dflist[grepl(codeword,ls(dflist))]
  return (dflist)
}
set_colname_adapter <- function (dflist, antigens = I00_Antigens, mode=1){
  lapply(1:length(dflist), function(setid){
    set = dflist[[setid]]
    antigens = antigens [[setid]]
    for (i in 3:length(set)){
      analytenum = as.numeric(strsplit(colnames(set[i]), split='.', fixed = TRUE)[[1]][2])
      if (mode == 1){
        colnames(set)[i] = antigens[antigens$analyte == analytenum,]$Antigen.name
      } else if (mode == 2){
        colnames(set)[i] = antigens[antigens$analyte == analytenum,]$Gene.name
      }}
    return(set)
  })
} 
mergedown <- function (dflist){
  outputdf = dflist[[1]]
  for (setid in 2:length(dflist)){
    originlist = colnames(dflist[[1]])[-1:-2]
    mergelist = colnames(dflist[[setid]])[-1:-2]
    uniquelist = mergelist[!(mergelist %in% originlist)]
    mergelist = mergelist[mergelist %in% originlist]}
  for (name in mergelist) {
    #print(name)
    outputdf[,name] = rowMeans(data.frame(dflist[[1]][,name], dflist[[setid]][,name]), na.rm = TRUE)
  }
  if (length(uniquelist) != 0){
    for (name in uniquelist) {
      #print(name)
      outputdf[name] = dflist[[setid]][,name]
    }
  }
  return(outputdf)
}
namelist <- function (dflist){
  namelist = list()
  for (i in 1:length(dflist)){
    namelist[i] = paste0("Set_", i)
  }
  return(namelist)
}
pathing <- function (name, term){
  #name = deparse(substitute(name))
  dir.create(paste0("plots\\",dataset, "\\", name))
  dir.create(paste0("plots\\",dataset, "\\", name, "\\", term))
  path = paste0("plots\\",dataset, "\\", name, "\\", term, "\\")
  return(path)
}
export_excel <-function(dflist, name = deparse(substitute(dflist))){
  lapply(1:length(dflist), function(setid){
    print(setid)
    if (!file.exists(paste(getwd(),"\\stats\\",dataset, "\\", name,".xlsx", sep=""))){
      print("File does not exist")
      write.xlsx(x=dflist[[setid]], file=paste(getwd(),"\\stats\\",dataset, "\\", name,".xlsx", sep=""), sheetName = paste("Dataset_", setid, sep=""), append=FALSE)
    }
    else {
      print("File Exists")
      write.xlsx(x=dflist[[setid]], file=paste(getwd(),"\\stats\\",dataset, "\\", name,".xlsx", sep=""), sheetName = paste("Dataset_", setid, sep=""), append=TRUE)
    }
  })
}

### Pretty Pictures ###
pp_analyte_plots <- function(set){
  name = deparse(substitute(set))
  dir.create(paste0("plots\\",dataset, "\\", name))
  dir.create(paste0("plots\\",dataset, "\\", name, "\\analyte_by_group"))
  path = paste0(getwd(), "\\plots\\",dataset,"\\", name, "\\analyte_by_group\\")
  for (n in 3:length(set)){
    plotwork = melt(set[c(1:2,n)])
    plotwork.summary <- plotwork %>% group_by(group) %>% summarise (
      #   iqr = IQR(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      value = mean(value, na.rm = TRUE)
    )
    #print(plotwork.summary)
    #names(plotwork)[1] = "name"
    #print(names(plotwork))
    ggplot(plotwork, aes(x=group, y=value)) + 
      geom_col(data=plotwork.summary, fill = NA, color = "black") + 
      geom_jitter( position = position_jitter(0.2), color = "black") +
      #geom_text(aes(label=ifelse(value>(max(plotwork.summary$value)+3*max(plotwork.summary$iqr)),as.character(name),'')),hjust=0,vjust=0) + 
      #geom_errorbar(aes(ymin = value-iqr, ymax = value+3*iqr), data = plotwork.summary, width = 0.2) + 
      geom_errorbar(aes(ymin = max(c(0, (value-sd))), ymax = value+sd), data = plotwork.summary) + 
      scale_x_discrete(labels=c("0" = "Control", "1" = "Diseased")) +
      labs(x="Group", y="MFI Value") +
      #geom_hline(yintercept = 250, colour = "red") +
      ggtitle (colnames(set)[n])  
    ggsave(paste0(path,colnames(set)[n],"_by_group.png", sep=""))}
}
pp_box_plot <- function(set){
  melted = melt(set)
  melted$variable = as.factor(melted$variable)
  pplot = ggplot(melted, aes(x=variable, y=value)) + 
    geom_boxplot() + 
    theme(axis.text.x=element_blank()) + 
    ggtitle("Box plot of patient antibody expression levels by analyte") +
    scale_y_continuous(limits=c(1,2500)) +
    geom_hline(yintercept = 250, colour = "red") +
    xlab("Analyte") + ylab("Median Flourescent Intensity (MFI)")
  pplot
}
pp_boxplot_set <- function(set){
  name = deparse(substitute(set))
  path = paste(getwd(), "/plots/", dataset, "/", "/boxplot_", name, ".png", sep="")
  plot(pp_box_plot(set))
  ggsave(path, width = 30, height = 20, units="cm")
}
pp_boxplot_dflist <- function(dflist){
  name = deparse(substitute(dflist))
  lapply(1:length(dflist), function(setid){
    set = dflist[[setid]]
    dir.create(paste(getwd(), "/plots/", dataset, "/", name, sep=""))
    path = paste(getwd(), "/plots/", dataset, "/", name, "/boxplot_set", setid, ".png", sep="")
    pp_box_plot(set)
    ggsave(path, width = 30, height = 20, units="cm")
  })
}
bar_plot <- function(set){
  name = deparse(substitute(set))
  dir.create(paste("plots\\",dataset, "\\", name, sep=""))
  dir.create(paste("plots\\",dataset, "\\", name, "\\barplot", sep=""))
  path = paste(getwd(), "\\plots\\",dataset,"\\", name, "\\barplot\\",sep="")
  melted = melt(set)
  melted = na.omit(melted)
  means <- ddply(melted, c("group", "variable"), summarise,
                 mean=mean(value))
  colnames(melted)[1] = "ID"
  melted$ID = as.factor(melted$ID)
  ggplot(data=means,aes(x=variable, y=mean, fill=group))+
    geom_col(position=position_dodge()) + theme(axis.text.x=element_text(angle = 270, vjust = 0.5, hjust = 0)) + 
    ggtitle( "Analytes  - Set 3")
  ggsave(paste(path, name, "barplot.png", sep=""))
}
scatter_plot <- function(set){
  name = deparse(substitute(set))
  dir.create(paste("plots\\",dataset, "\\", name, sep=""))
  dir.create(paste("plots\\",dataset, "\\", name, "\\scatterplot", sep=""))
  path = paste(getwd(), "\\plots\\",dataset,"\\", name, "\\scatterplot\\",sep="")
  melted = melt(set)
  melted = na.omit(melted)
  colnames(melted)[1] = "ID"
  melted$ID = as.factor(melted$ID)
  melted$group = as.factor(melted$group)
  ggplot(data=melted, aes(x=ID, y=value))+
    #geom_point(aes(color=factor(group)))+
    geom_jitter(aes(color=factor(group))) + 
    #geom_hline(yintercept = 250) +
    scale_y_continuous(trans='log10') + 
    labs(title="Scatterplot of all MFI measurements by patient",x="Patient ID", y="MFI Level", fill="Group", colour="Group") +
    scale_color_manual(labels = c("Healthy", "Diseased"), values = c("Blue","Red")) +
    theme_bw() +
    #ggtitle( "Scatterplot of all Analytes  - Set 3") +
    theme(axis.text.x=element_text(angle=270)) 
}

