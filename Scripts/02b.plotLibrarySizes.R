# setwd("~/Desktop/junk/")

library(RColorBrewer)
library(ggplot2)


outPath <- "images"
dir.create(outPath, showWarnings=FALSE)

outFile  <- paste0(outPath,"/Raw_and_Trimmed_LibrariesSize.pdf")

## Read data
Raw <- read.table("02_LibrarySizes/Raw_TotalSequences.txt", sep="",as.is = T,header = T,stringsAsFactors = F)
Trim <- read.table("02_LibrarySizes/Trimmed_TotalSequences.txt", sep="",as.is = T,header = T,stringsAsFactors = F)


# Process data frame
Raw$Type <- "Raw"
Trim$Type <- "Trimmed"
TotalCounts <- rbind(Raw,Trim)

## Get as million counts
TotalCounts$"Counts" <- TotalCounts$Counts/1000000 #million reads

# Use regex to get genotypes and treatments
TotalCounts$Library <- TotalCounts$File
TotalCounts$Group <- gsub("_S[0-9].*$","\\1",TotalCounts$Library)
#TotalCounts$Time <- gsub("_tr.[ab]$","\\1",TotalCounts$Group)
TotalCounts$Treat <- sapply(strsplit(gsub("Col|SALK-","",TotalCounts$Group),"-"),"[",2)
TotalCounts$Genotype <- sapply(strsplit(TotalCounts$Group,"-"),"[",1)
TotalCounts <- TotalCounts[order(TotalCounts$Group),]


### Save to PDF
pdf(outFile, width = 6,height = 8)
#pdf(outFile, paper = "a4")
### ggplot2
par(oma = c(0.2,0.2,0.2,0.2), mar=c(3.5, 3.5, 2, 1),mgp = c(5, 0.5, 0))
## Combined size of libraries  OK
ggplot(data=TotalCounts, aes(x=Library, y=Counts, fill=Type,width=0.7)) +
  geom_bar(stat = "identity", position = position_dodge(width=0.7)) + 
  scale_fill_brewer(palette="Accent") +
  coord_flip() + theme_minimal() +
  theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),axis.text=element_text(size=5)) +
  xlab("Library") + ylab("Counts (millions)") + #facet_grid(Genotype~Treat) + 
  ggtitle("Raw Library Sizes") 

dev.off()
