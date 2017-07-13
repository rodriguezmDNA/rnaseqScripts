setwd("~/Desktop/RootSeq_NitrogenPioneer/")

library(RColorBrewer)
library(ggplot2)
outFile  <- "images/00_Raw_LibrariesSize.pdf"

outPath <- "images"
dir.create(outPath, showWarnings=FALSE)

pdf(outFile, width=12, height=12)


### 
# Before Trimming: 
TotalCounts <- read.table("00b_LibrarySizes/Raw_TotalSequences.txt", sep=" ")
TotalCounts <- TotalCounts[-grep("U",TotalCounts$V1),] #Unknown libraries away
colnames(TotalCounts) <- c("Library","Counts")

TotalCounts$"Counts" <- as.numeric(sub("Total Sequences ","",TotalCounts$"Counts"))


TotalCounts$"Counts" <- TotalCounts$"Counts"/1000000 #million reads


rownames(TotalCounts) <- sub("_fastqc.zip$","",TotalCounts[,1])
TotalCounts <- TotalCounts[,2,drop=F]

#
TotalCounts$Group <- gsub("^str.[0-9]_|_set.[0-9]$","\\1",rownames(TotalCounts))
TotalCounts$Time <- gsub("_tr.[ab]$","\\1",TotalCounts$Group)
TotalCounts$Treat <- gsub("^tm.[0-9]*_","\\1",TotalCounts$Group)
TotalCounts$Library <- gsub("^str.[0-9]_|set|[.]|tr|tm|_","\\1",rownames(TotalCounts))

TotalCounts <- TotalCounts[order(TotalCounts$Group),]





###
par(oma = c(0,0,2,0), mar = c(6.1, 6, 2.1, 1.0),mgp = c(5, 0.5, 0))
## Combined size of libraries  OK
ggplot(data=TotalCounts, aes(x=Treat, y=Counts, fill=Time)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  scale_fill_brewer(palette="Set3") +
  xlab("Treatment") + ylab("Counts (millions)") +
  ggtitle("Total library sizes (Grouped by time)")

## Individual libraries OK

ggplot(data=TotalCounts, aes(x=TotalCounts$Library, y=Counts ,color=Time)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  #geom_point(stat="identity") +  
  scale_fill_brewer(palette="Set3") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Group, scales = "free", space="free") +
  xlab("Library") + ylab("Counts (millions)") +
  ggtitle("Individual library sizes (Time~Treatment)")


dev.off()
