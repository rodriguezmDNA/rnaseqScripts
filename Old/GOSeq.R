setwd("~/Desktop/ath_GxT/06_DGE/")


shortname <- "TFmutants_wtNtreat"
load("genome_GAcounts_wtTreatment_filterd/DEListgenome_GAcounts_wtTreatment_filterd.RData")

sapply(DEList,nrow)


## Read files
##GO
GFFfile = "../meta/TAIR10_withTransposons.gff"
GOFile = "../meta/ATH_GO_GOSLIM.txt"

source("../Scripts/metaFunctions_forNetworkAnalysis.R")

## Filter logFC of significant genes                             
logfcExprs <- lapply(DEList,function(x){x[x[,grep("adj.P.Val",colnames(x))] < 0.05,grep("logFC",colnames(x)),drop=F]})
deGenesList <- lapply(logfcExprs, rownames)


### GO terms of DE genes in data
## Get GO terms enriched within each group
background <- Reduce(union,lapply(DEList,rownames))
par(mfrow=c(3,3))
GOenrich_DE <- findGOTerms(deGenesList,background,GFFfile,GOFile)
par(mfrow=c(1,1))
lapply(GOenrich_DE,dim)
## --

# Filter by pValue and only Biological Process
sapply(GOenrich_DE, nrow)
signGOList <- lapply(GOenrich_DE, function(x){  x[x$over_represented_pvalue < 0.05 & x$ontology == "BP",] })
## Use only overrepresentated pValue
for (each in names(signGOList)){
  rownames(signGOList[[each]]) <- signGOList[[each]]$term
  signGOList[[each]] <- signGOList[[each]][,"over_represented_pvalue",drop=F]
  colnames(signGOList[[each]]) <- each
}
sapply(signGOList, nrow)

## Heatmap of GO categories per cluster
library(RColorBrewer)
library(gplots)
## --
titulo <- "GOEnrichment_from_DE_genes"
colors <- colorRampPalette(c("skyblue","steelblue2","steelblue4"))
path <- paste("",titulo,"_",shortname,".pdf",sep="")
pdf(path,paper = "a4")
#
hmData <- as.matrix(condenseListTables(signGOList))
head(hmData)
hmData[hmData==0] <- NaN
head(hmData)
hmData <- -(log10(hmData))
#######
heatmap.2(hmData,col=colors(120),
          keysize = 1.5,
          symkey = F,
          #key.par=list(mar=c(3.5,0,3,0)),
          na.color = "white",
          margins = c(2,6),
          
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(3, 8, 2), 
                     c(7, 4, 9),
                     c(6, 5, 1)), 
          #lhei=c(1.5, 0,5), 
          #lwid=c(2.5, 4, 1),
          #lmat = lmat, 
          lhei=c(0.15,0.22,0.9),
          lwid=c(0.8,0.4,0.3),
          ##
          scale="none",
          density.info = "none", 
          key.xlab = "-log10(pVal)", 
          trace = "none",
          dendrogram = "none",
          cexRow = 0.064,
          cexCol = 0.15,
          Rowv = F,Colv = F,
          main="GO Enrichment",cex.main=0.5)
dev.off()


titulo <- paste("",titulo,"_",shortname,".csv",sep="")
write.csv(titulo,x = hmData)
