### Create matrix of counts 
#########################################################
# github: rodriguezmDNA
# j rodriguez medina
# Brady lab @ ucdavis
# Generate plots on mapping efficiency (%) from alignment scripts
# Last modified
# 2017.07.jrm


# From htseq or cDNA obtained individual library counts 
setwd("~/Desktop/junk/")

option="bwt1_cDNA" #bwt1_genome
#option="bwt1_genome" #

#### Load libraries
library(reshape)
library(ggplot2)
library(RColorBrewer)


###########

filesDir = paste0("04_Counts/",option)
countsDir <- as.list(list.files(filesDir,pattern = ".counts.txt",full.names = T))
countsList <- lapply(countsDir, read.table,row.names = 1)
names(countsList) <- lapply(countsDir, basename)



#########
source("Scripts/functions_jrm.R")

countsMatrix <- condenseListTables(countsList)
sampleNames <- gsub(".counts.txt*$","",basename(unlist(countsDir)))
colnames(countsMatrix) <- sampleNames
head(countsMatrix)
dim(countsMatrix)
totalCounts <- colSums(countsMatrix)
# --

# -- Check htseq info if used genome alignment
if (option=="bwt1_genome"){
  cat("Plotting htseq stats \n")
  ###
  logCounts <- as.list(list.files("logs/Counts/",pattern = ".counts",full.names = T))
  summaryList <- lapply(logCounts, read.table,row.names = 1)
  summaryMatrix <- do.call(cbind,summaryList)
  colnames(summaryMatrix) <- sampleNames
  rownames(summaryMatrix) <- gsub("_","",rownames(summaryMatrix))
  
  # --
  head(summaryMatrix)
  summaryMatrix <- rbind(summaryMatrix,"totalWithFeature"=totalCounts)
  #summaryMatrix <- summaryMatrix[rowSums(summaryMatrix) != 0,]
  
  ## Write to file
  titulo <- paste0("images/htseq_stats_",option,".txt")
  write.table(file = titulo,summaryMatrix,quote = F,col.names = NA,sep="\t")
  
  ## Plot
  titulo <- paste0("images/htseq_stats_",option,".pdf")
  pdf(titulo,paper = "a4r")
  mar.default <- c(15,4,8,6) + 0.1
  
  ### rplot
  #barplot(as.matrix(summaryMatrix),
  #        col=rainbow(3),legend.text = T,width = 3,las=2,cex.names = 0.6)
  
  ## Use ggplot2
  summaryMatrix$Description <- rownames(summaryMatrix)
  meltSummary <- melt(summaryMatrix)
  
  x <- ggplot(meltSummary, aes(variable, value)) +   
    geom_bar(aes(fill = Description), position = "dodge", stat="identity") +
    labs(title = "HTSeq Annotation") + theme_gray() +
    scale_fill_manual(values = rev(brewer.pal(nrow(summaryMatrix),"Accent"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~variable,nrow=2,scales = 'free_x')
  print (x)
  
  ## Close pdf
  dev.off()
  
} else {
  cat("Skipping htseq stats \n")
}


## Save everything
outDir <- "05_RawCounts/"
dir.create(outDir,showWarnings = F)
#
outCounts <- "RawCounts_"
outCounts <- paste(outDir,outCounts,option,".csv",sep="")
write.csv(countsMatrix,file = outCounts,quote = F,row.names = T)
dim(read.csv(outCounts,row.names = 1)) #Check it got written correctly.
