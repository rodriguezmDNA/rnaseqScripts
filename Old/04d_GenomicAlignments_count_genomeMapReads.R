#library(GenomicRanges)
#library(Rsamtools)

#source("https://bioconductor.org/biocLite.R")
#biocLite(c("Rsamtools","GenomicFeatures","GenomicAlignments"))
#biocLite("GenomicFeatures")

#setwd("../")

library("ggplot2")
library("GenomicFeatures")
library("Rsamtools")
library("GenomicAlignments")
library("reshape")


option="bwt1_genome" #or bwt2_genome

## Where the BAM files are
pathFile <- paste0("03_alignment/",option)
files <- list.files(pathFile,pattern = "Col-1-1_S109_L006",full.names = T) ## After bowtie alignment
#files <- list.files("03b_deduplicated/",pattern = "bam",full.names = T) ## After dedupliating with overAmp

## Read GFF object with annotation data
gffFile <- list.files("meta/",pattern = "gff",full.names = T)
TxDB <- makeTxDbFromGFF(file = gffFile,circ_seqs = character())
TxDB
seqlevels(TxDB)
## Manually fix annotation. In the TAIR genome sequence files plastid genomes are different than in the annotation.
#seqlevels(TxDB)[6] <- "chloroplast"
#seqlevels(TxDB)[7] <- "mitochondria"


## Select annotation model
exons <- exonsBy(TxDB, by="gene")
gene <- transcriptsBy(TxDB, by="gene")
cds <- cdsBy(TxDB, by="gene")


## Outdir
outDir <- "05_RawCounts/"
dir.create(outDir,showWarnings = F)

imgDir <- "images/"
dir.create(imgDir,showWarnings = F)

countDir <- paste0("04_Counts/","GenomicAlignments_",option,"/")
dir.create(countDir,showWarnings = F,recursive = T)


listSummary <- list()

for (each in files){
  name <- gsub(".bam","",basename(each))
  cat ("Getting counts for", name,"\n")
  ##
  bamfiles <- BamFileList(each)
  
  seExon <- summarizeOverlaps(features=exons, 
                              reads=bamfiles,
                              mode="Union",
                              singleEnd=T,
                              ignore.strand=F,
                              inter.feature = T,
                              fragments=F)
  head(assay(seExon),12)
  
  # Save to list
  listSummary[[name]] <- assay(seExon)
  ## Save individual files
  # Skip and put everything into a single table.
  titulo = paste0(countDir,name,".countsGA.txt")
  write.table(titulo,x = assay(seExon),quote = F,row.names = T,col.names = NA)
}

length(listSummary)
sapply(listSummary,nrow)


## Create matrix
totalSizes <- as.data.frame(sapply(listSummary,colSums))
rownames(totalSizes) <- gsub(".bam","",rownames(totalSizes))
colnames(totalSizes) <- "Total Reads"

######## Melt and Plot with ggplot2
totalSizesMelt <- melt(t(totalSizes))

titulo <- paste0("images/GenomicAlignments_libSizes_genome",".pdf")
pdf(titulo,paper = "a4r")
mar.default <- c(15,4,8,6) + 0.1

ggplot(totalSizesMelt, aes(X2, value)) +   
  geom_bar(aes(fill = X1), position = "dodge", stat="identity") +
  labs(title = "GenomicAlignments Annotation") + theme_gray() +
  #scale_fill_manual(values = rev(brewer.pal(nrow(totalSizesMelt),"Accent"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ facet_wrap(~variable,nrow=2,scales = 'free_x')
dev.off()
####################################################################################


######## Write to file
titulo <- paste0("images/GenomicAlignments_libSizes_genome",".txt")
write.table(file = titulo,totalSizes,quote = F,col.names = NA,sep="\t")

####################################################################################

source("Scripts/functions_jrm.R")
RawCounts <- condenseListTables(listSummary)

#RawCounts <- do.call("cbind",listSummary)
head(RawCounts)
dim(RawCounts)

## Genes witout counts
table(rowSums(RawCounts) != 0)


titulo = paste0(outDir,"GenomicAlignments_countsExon.csv")
write.csv(titulo,x = RawCounts)




