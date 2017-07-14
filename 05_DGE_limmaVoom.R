# RNA limma

## http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
# Chapter 15

##In the limma approach to RNA-seq, read counts are converted to log2-counts-per-million (logCPM)
#mean-variance relationship is modelled either with precision weights or with an empirical Bayes prior trend.
#The precision weights approach is called “voom” and the prior trend approach is called “limma-trend”

setwd("~/Desktop/junk/")
library(edgeR)
library(gplots)
library(calibrate)
library(RColorBrewer)
library(limma)

pValCut <- 0.05
#### Start PDF



## Created a folder with a name to identify the analysis. Empty by default
shortName <- "genome_GAcounts"
## Pre
outDir = paste0("06_DGE/",shortName)
dir.create(outDir,showWarnings = F,recursive = T)
doFilter <- F
#if (doFilter) {
#  removeSamples <- c("ARID5_10mM_rep3_str12","HAT22_10mM_rep3_str12")
#}


imgPath <- paste0(outDir,"/DE_limmaVoom_Analysis_",shortName,".pdf")
pdf(imgPath,paper = "USr")


source("Scripts/limma_DEG_Functions.R")

######## Read Data
#GeneCounts <- read.csv("05_RawCounts/RawCounts_bwt1_cDNA.csv",row.names = 1) #1 gene
#GeneCounts <- read.csv("05_RawCounts/RawCounts_bwt1_genome.csv",row.names = 1) #1 gene
GeneCounts <- read.csv("05_RawCounts/GenomicAlignments_countsExon.csv",row.names = 1) #1 gene
rownames(GeneCounts) <- gsub("\\.[0-9].*$","",rownames(GeneCounts))

###
dim(GeneCounts)

dataInfo <- data.frame(do.call("rbind",strsplit(colnames(GeneCounts),"_")),row.names = 1)
colnames(dataInfo) <- c("S","Lane","R","other1")
head(dataInfo)

## Rename columns using the metadata file.
#meta <- read.csv("meta/phenoInfoFull.csv",as.is = T,row.names = 1)
meta <- read.csv("meta/phenoInfoFull_swapped_12-1_-12-2.csv",as.is = T,row.names = 1)
## Add extra info
meta <- cbind(meta[rownames(dataInfo),],dataInfo)

#
# Rename columns of counts matrix.
tmp <- gsub("_.*$","",colnames(GeneCounts))
meta <- meta[tmp,,drop=F] #Reorder
head(meta)

head(cbind(colnames(GeneCounts),meta$SampleNameFull,meta$SampleLine,meta$LineName))
colnames(GeneCounts) <- gsub("_str.*$","",meta$SampleNameFull) #Assign sample names to counts table
#head(meta)
#head(GeneCounts)


## Remove samples that looked bad in the MDS plots
removeSamples <- c("ARID5_1mM_rep2_str3", #MDS
                   "NAC102_1mM_rep2_str3", #MDS
                   "LBD4_1mM_rep3_str7", #MDS
                   "ARID5_10mM_rep3_str12", #Prob Switched
                   "HAT22_10mM_rep3_str12" #Prob switched
)

## Remove samples
if (doFilter){
  rmIDX <- which(meta$SampleNameFull %in% removeSamples)
}


ncol(GeneCounts)
if (doFilter){
  cat("Samples to remove:\n",
      paste(colnames(GeneCounts)[rmIDX],"\n",sep=""),"\n"
  );
  GeneCounts <- GeneCounts[,-rmIDX]
  meta <- meta[-rmIDX,]
} else {
  cat("Not removing any samples \n")
}
cat(ncol(GeneCounts),"samples in the data frame \n")
#
ncol(GeneCounts) == nrow(meta)
dim(GeneCounts)
##

cat("Removing genes with 0 counts on all conditions \n")
cat("Initial number of genes:",nrow(GeneCounts),"\n")
rmIDX <- which(rowSums(GeneCounts) == 0)
if (length(rmIDX) != 0){
  cat("Removing",length(rmIDX),"genes \n")
  GeneCounts <- GeneCounts[-rmIDX,]
  cat("Remaining number of genes:",nrow(GeneCounts),"\n")
} else {
  cat("No genes with 0 counts on all  \n")
}

## Linear modeling
Genotype <- relevel(as.factor(meta$GenotypeName),ref="Col0")
Treatment <- relevel(as.factor(meta$Treatment),ref="10mM")
Group <- relevel(as.factor(meta$GroupName),ref="Col0_10mM")
Strip <- relevel(as.factor(meta$Strip),ref="str1")
Full <- as.factor(meta$SampleNameFull)
Set <- as.factor(meta$Set)
Lane <- as.factor(meta$Lane)


meta$GenotypeName <- relevel(as.factor(meta$GenotypeName),ref="Col0")
meta$Treatment <- as.factor(meta$Treatment)

## Design matrix

#design <- model.matrix(~GenotypeName*Treatment+Set+Lane,data = meta)#+Lane
design <- model.matrix(~0+GroupName+Set+Lane, data =meta)#+Lane
colnames(design) <- gsub("meta|Name|Genotype|Treatment|:|-|/|Group","",colnames(design))
head(design)


### Use cpms to uncover lowly expressed genes
dge <- DGEList(counts=GeneCounts,remove.zeros = T)
# Filter: Genes with total counts more than 
# minCount <- 50
# A <- rowSums(dge$counts)
# isexpr <- A > minCount
# table(isexpr)
# # OR
range(table(meta$GroupName))
sampleMin <- min(table(meta$GroupName))
minCPM <- 1
isexpr <- rowSums(cpm(dge) > minCPM) >= sampleMin #Make sure to use the minimum number of reps
table(isexpr)
cat("Removing lowly expressed genes ( <",minCPM,"cpm on at least",sampleMin,"samples) \n")

cat("Removing",table(isexpr)[1],"genes \n")
cat("Remaining number of genes:",table(isexpr)[2],"\n")
dge <- dge[isexpr,,keep.lib.size = FALSE]

## Do TMM normalization after filtering (otherwise lib sizes are unbalanced)
dge <- calcNormFactors(dge)

#In the limma-trend approach, the counts are converted to logCPM values using edgeR’s cpm function:
#logCPM <- cpm(dge, log=TRUE, prior.count=1) #The prior count is used here to damp down the variances of logarithms of low counts
#dim(logCPM)

cpmNormalizedExpression  <- cpm(dge)
#save(cpmNormalizedExpression,file = "06_DGE/cpmNormalizedExpression.RData")

### Use voom directly on 
#v <- voomWithQualityWeights(dge, design, plot=TRUE)
#v <- voom(dge, design, plot = TRUE,normalize.method = "quantile") #"scale","cyclicloess" or "quantile"
v <- voomWithQualityWeights(dge, design=design, normalization="quantile", plot=TRUE)
dim(v)


library(RColorBrewer)


## -- Produce graphics
Targets <- meta

col.concent <- c("darkgreen", "darkblue")
#col.Genotype <- c("#000000",RColorBrewer::brewer.pal(12,"Set3"))
col.Group <- c("#ba4343",RColorBrewer::brewer.pal(12,"Set3"))

col.Genotype <- c(
  "#ff7527",
  "#2d66a0",	
  "#ba4343",
  "#84ac38",
  "#267501",	
  "#513267",
  "#0f005b",
  "#3fbfcd",
  "#660099",
  "#1466a2",
  "#f49502",
  "#da0000",
  "#e8a3d2")

col.Genotype <- unlist(lapply(col.Genotype,function(x){rep(x,2)}))
#col.Group <- as.vector(sapply(col.Group, function (x) rep(x,2)))

## Boxplot
boxplot(v$E, range=0,col=col.Group[as.factor(meta$GroupName)], 
        ylab="log2[counts]", xlab="sample", main="voom normalized counts",
        cex.axis=0.5,las=2)


Targets$Treatment <- factor(Targets$Treatment)
Targets$GenotypeName <- factor(Targets$GenotypeName)
Targets$GroupName <- factor(Targets$GroupName)
Targets$Set <- factor(Targets$Set)
Targets$Strip <- factor(Targets$Strip)

###


#pcDimensions <- list(c(1,2),c(3,4),c(5,6))
## Treatment
# par(mfrow=c(2,2))
# lapply(pcDimensions,function(x){
#   plotMDS(v, labels=meta$SampleNameFull, 
#           col=col.concent[Targets$Treatment],
#           cex=0.7,dim.plot = x)
#     })
# plot.new()
# legend("center",legend = unique(Targets$Treatment),fill = unique(col.concent[Targets$Treatment]),bty="n")
## Do the same for Genotype, group and other factors...



par(mfrow=c(1,1))
#Or just use glimma
library(Glimma)
glMDSPlot(v, labels=meta$SampleNameFull, groups=Targets[,c("GenotypeName","Treatment","GroupName","Set","Strip","Lane")],
          folder=paste0(outDir,"/glimma"),
          launch=T)

## --


### Do contrasts

#### Fit and do Diff Expression
#cat("Using Set as blocking factor")
cor <- duplicateCorrelation(v, design, block = meta$Set)
cor$consensus

## Fit
v2 <- lmFit(v, design, block = meta$Set, correlation = cor$consensus)

### do contrasts
contrastMatrix <- makeContrasts(
  #ANAC032
  "nac032_int"=(ANAC032_1mM-ANAC032_10mM)-(Col0_1mM-Col0_10mM),
  "nac032_10mM"=ANAC032_10mM-Col0_10mM,
  "nac032_1mM"=ANAC032_1mM-Col0_1mM,
  #"nac032"=(ANAC032_1mM+ANAC032_10mM)-(Col0_10mM+Col0_1mM),
  ##ARF9
  # "arf9_int"=(ARF9_1mM-ARF9_10mM)-(Col0_1mM-Col0_10mM),
  # "arf9_10mM"=ARF9_10mM-Col0_10mM,
  # "arf9_1mM"=ARF9_1mM-Col0_1mM,
  # #"arf9"=((ARF9_1mM+ARF9_10mM))-((Col0_10mM+Col0_1mM)),
  # ##
  # ##ARF181
  # "arf181_int"=(ARF181_1mM-ARF181_10mM)-(Col0_1mM-Col0_10mM),
  # "arf181_10mM"=ARF181_10mM-Col0_10mM,
  # "arf181_1mM"=ARF181_1mM-Col0_1mM,
  # #"arf181"=((ARF181_1mM+ARF181_10mM))-((Col0_10mM+Col0_1mM)),
  # ##ARF182
  # "arf182_int"=(ARF182_1mM-ARF182_10mM)-(Col0_1mM-Col0_10mM),
  # "arf182_10mM"=ARF182_10mM-Col0_10mM,
  # "arf182_1mM"=ARF182_1mM-Col0_1mM,
  # #"arf182"=((ARF182_1mM+ARF182_10mM))-((Col0_10mM+Col0_1mM)),
  # # "ARID5"   
  # "arid5_int"=(ARID5_1mM-ARID5_10mM)-(Col0_1mM-Col0_10mM),
  # "arid5_10mM"=ARID5_10mM-Col0_10mM,
  # "arid5_1mM"=ARID5_1mM-Col0_1mM,
  # #"arid5"=((ARID5_1mM+ARID5_10mM))-((Col0_10mM+Col0_1mM)),
  # # "ERF107"
  # "erf107_int"=(ERF107_1mM-ERF107_10mM)-(Col0_1mM-Col0_10mM),
  # "erf107_10mM"=ERF107_10mM-Col0_10mM,
  # "erf107_1mM"=ERF107_1mM-Col0_1mM,
  # #"erf107"=((ERF107_1mM+ERF107_10mM))-((Col0_10mM+Col0_1mM)),
  # # "HAT22"   
  # "hat22_int"=(HAT22_1mM-HAT22_10mM)-(Col0_1mM-Col0_10mM),
  # "hat22_10mM"=HAT22_10mM-Col0_10mM,
  # "hat22_1mM"=HAT22_1mM-Col0_1mM,
  # #"hat22"=((HAT22_1mM+HAT22_10mM))-((Col0_10mM+Col0_1mM)),
  # # "HMG"
  # "hmg_int"=(HMG_1mM-HMG_10mM)-(Col0_1mM-Col0_10mM),
  # "hmg_10mM"=HMG_10mM-Col0_10mM,
  # "hmg_1mM"=HMG_1mM-Col0_1mM,
  # #"hmg"=((HMG_1mM+HMG_10mM))-((Col0_10mM+Col0_1mM)),
  # # "LBD4" 
  # "lbd4_int"=(LBD4_1mM-LBD4_10mM)-(Col0_1mM-Col0_10mM),
  # "lbd4_10mM"=LBD4_10mM-Col0_10mM,
  # "lbd4_1mM"=LBD4_1mM-Col0_1mM,
  # #"lbd4"=((LBD4_1mM+LBD4_10mM))-((Col0_10mM+Col0_1mM)),
  # # "MYB29"    
  # "myb29_int"=(MYB29_1mM-MYB29_10mM)-(Col0_1mM-Col0_10mM),
  # "myb29_10mM"=MYB29_10mM-Col0_10mM,
  # "myb29_1mM"=MYB29_1mM-Col0_1mM,
  # #"myb29"=((MYB29_1mM+MYB29_10mM))-((Col0_10mM+Col0_1mM)),
  # # "NAC102"
  # "nac102_int"=(NAC102_1mM-NAC102_10mM)-(Col0_1mM-Col0_10mM),
  # "nac102_10mM"=NAC102_10mM-Col0_10mM,
  # "nac102_1mM"=NAC102_1mM-Col0_1mM,
  # #"nac102"=((NAC102_1mM+NAC102_10mM))-((Col0_10mM+Col0_1mM)),
  # # "RAV2EDF2"
  # "rav2_int"=(RAV2EDF2_1mM-RAV2EDF2_10mM)-(Col0_1mM-Col0_10mM),
  # "rav2_10mM"=RAV2EDF2_10mM-Col0_10mM,
  # "rav2_1mM"=RAV2EDF2_1mM-Col0_1mM,
  # #"rav2"=((RAV2EDF2_1mM+RAV2EDF2_10mM))-((Col0_10mM+Col0_1mM)),
  levels = design)


fit2 <- contrasts.fit(v2, contrastMatrix)
fit2 <- eBayes(fit2)


#######################
results <- decideTests(fit2)#,lfc = logCut,p.value = pValCut)
summary(results)

if (ncol(results) <= 5){
  vennDiagram(results,include = c("up","down"))
}

DESummary <- t(summary(decideTests(fit2)))
# --

titulo <- paste0(outDir,"/DESummaryInteractions_",shortName,".csv")
write.csv(x=DESummary,titulo,quote = F,row.names = T)

###
# 
# plotData <- t(DESummary)#prop.table(DESummary,margin = 1))
# plotData <- plotData[-2,]
# yMax <- max(colSums(plotData))
# rownames(plotData) <- c("Down","Up")
# barplot(cbind(plotData,NA,NA),legend.text = rownames(plotData),col=c("tan","steelblue4"),
#         horiz = T,beside = F, las=2,#ylim = c(0,yMax*1.02),las=2,cex.names = 0.6, #border = F, bty="n",
#         main="log2 of total Up/Down regulated genes")
#dev.off()


cat("decideTests() function results: \n")
print(DESummary)

# --

voomNormalizedExpression <- v$E
titulo <- paste0(outDir,"/voomNormalizedExpression",shortName,".RData")
save(voomNormalizedExpression,file = titulo)


# --
geneAlias <- read.delim("meta/gene_aliases.txt",header = T, fill = NA,stringsAsFactors = F)
geneAlias <- split(geneAlias,geneAlias$locus_name)

AGI2Symbol <- lapply(geneAlias, function(x){
  idx <- grepl("^(?!A[Tt])",x[,"symbol"],perl = T); #Negates genes starting with At or AT
  if (all(idx == F)){
    paste(x[,"symbol"],collapse = "/")
  } else {
    symbols <- paste(x[idx,"symbol"],collapse = "/") }
})

#head(AGI2Symbol)
AGI2Symbol <- unlist(AGI2Symbol)
head(AGI2Symbol)

## --
logCut <- 1
pValCut <- 0.05

### Set list of contrasts to drop:
## For genotypes
#uniqContrasts <- unique(gsub("meta|Name|Genotype|Treatment|:|-|/|Group","",meta$GenotypeName)) #Each genotype
#uniqContrasts <- uniqContrasts[!uniqContrasts=="Col0"] #Remove controls
# For treatment
#treatIndex <- which(colnames(fit2$coefficients) %in% levels(meta$Treatment)) #Add treatment

uniqContrasts <- colnames(contrastMatrix)
# For treatment
#treatIndex <- which(colnames(fit2$coefficients) %in% levels(meta$Treatment)) #Add treatment

cat("Genotypes that will be analyzed:" ,paste0(uniqContrasts),"\n")
#cat("Along with:" ,colnames(fit2$coefficients)[treatIndex],"\n")

### Prepare lists
DEList <- list()
DESignificant <- list()

uniqContrasts <- uniqContrasts[grep("_",uniqContrasts)]
#
for (contrast in uniqContrasts){
  cat(" - - - \n")
  #
  # Get indices that correspond to treatment, genotype and interaction
  
  #genoIndex <- grep(contrast,colnames(fit2$coefficients))
  dropContrasts <- contrast#c(genoIndex)#,treatIndex)
  ##
  contrastNames <- contrast#colnames(fit2$coefficients)[dropContrasts]
  cat("Contrast:", paste0(contrastNames),"\n")
  #
  tmp <- topTable(fit2, coef=dropContrasts,number = Inf,sort.by = "none")
  tmp[,"Symbol"] <- rownames(tmp)
  #-- Add gene symbols
  Genes <- rownames(tmp)
  idx <- intersect(names(AGI2Symbol),Genes)
  tmp[idx,"Symbol"] <- AGI2Symbol[idx]
  
  
  
  ###Change names of columns
  # Get indices of columns that correspond to logFC
  #tmpIDX <- grep(paste(contrastNames,collapse="|"),colnames(tmp))
  colnames(tmp) <- paste(contrast,colnames(tmp),sep = ".")
  #colnames(tmp)[tmpIDX] <- paste(colnames(tmp)[tmpIDX],"logFC",sep = ".")
  #
  DEList[[contrast]] <- tmp
  
  ## Filter
  tmpIDX <- grep("adj.P.Val",colnames(tmp))
  tmpSign <- tmp[tmp[,tmpIDX] < pValCut,] 
  nrow(tmpSign)
  DESignificant[[contrast]] <- tmpSign
  #
  cat ("Number of DEGs on",contrast,":",nrow(tmpSign),"\n")
  #cat("\n")
  #
}  


###
sapply(DEList,function(x){table(x[,grep("adj.P.Val",colnames(x))] < 0.05)})
significantDE <- lapply(DEList,function(x){ x[x[,grep("adj.P.Val",colnames(x))] < 0.05,] })
sapply(significantDE,nrow)
##


#
### Save full list 
titulo <- paste0(outDir,"/DEList",shortName,".RData")
save(DEList,file = titulo)
### Save list of significant genes
titulo <- paste0(outDir,"/Significant_DElist",shortName,".RData")
save(significantDE,file = titulo)

### Save each contrast individually
dirContrasts <- paste0(outDir,"/contrasts")
dir.create(dirContrasts)
for (each in names(significantDE))
{
  cat ("Saving contrast:", each,"\n")
  titulo <- paste0(dirContrasts,"/",shortName,"_",each,".RData")
  write.csv(file = titulo,significantDE[[each]])
}
##################


# pdf("06_DGE/Heatmaps_Significant.pdf",paper = "US")
# lapply(list, function)contrast){
#   titulo <- contrast[x]
#   data <-listDE[[x]][,grep("logFC",colnames(listDE[[x]]))]
#   #head(data)
#   NMF::aheatmap(data,main =titulo)
# },listDE=DESignificant,contrast=names(DESignificant))
# dev.off()


dev.off()
###########################

## -- Make tables
source("Scripts/functions_jrm.R")

## Convert everything to a single table
DE_All <- condenseListTables(DEList) ## Use a custom function
dim(DE_All)
head(DE_All)

################################
## Get all symbols
Symbols <- Reduce(unique,lapply(DEList,function(x){x[,grep("Symbol",colnames(x))]}))
GeneNames <- Reduce(unique,lapply(DEList,function(x){rownames(x)}))
names(Symbols) <- GeneNames
################################

### This won't work in this case, we have more logFC columns than adjPVals columns
# ## Get logFC values and adjPvals of all genes and contrasts
All_logFC <- DE_All[,grep("logFC",colnames(DE_All))]
#All_logFC <- All_logFC[,-grep("_",colnames(All_logFC))] #Only interaction FC
All_adjpVal <- DE_All[,grep("adj.P.Val",colnames(DE_All))]
head(All_adjpVal)
head(All_logFC)

## Filter out non significant logFCs
DE_Significant <- All_logFC
DE_Significant[All_adjpVal > pValCut] <- NA

#Create a binary matrix of significance
signMatrix <- All_adjpVal
signMatrix[All_adjpVal < pValCut] <- "*"
signMatrix[All_adjpVal > pValCut] <- NA
## --

#write.csv(x = All_logFC,"06_DGE/All_DEG_Interactionss.csv",quote = F,row.names = T)
#write.csv(x = All_adjpVal,"06_DGE/All_adjPVals_Interactionss.csv",quote = F,row.names = T)
#write.csv(x = DE_Significant,"06_DGE/Significant_Unfiltered_Interactionss.csv",quote = F,row.names = T)
## -- 


### Check the expression values due to genotype of the genes where the TDNA insertion is
lines <- c(
           "AT1G77450", #nac32
           "AT4G23980", #arf9
           "AT3G61830", #arf18
           "AT1G76510", #arid5
           "AT5G61590", #erf107
           "AT4G37790", #hat22
           "AT1G04880", #hmg
           "AT1G31320", #lbd4
           "AT5G07690", #myb29,
           "AT5G63790",#nac102
           "AT1G68840") #rav2/edf2

linesSymbols <- c(
  "NAC032", #nac32
  "ARF9", #arf9
  "ARF18", #arf18
  "ARID5", #arid5
  "ERF107", #erf107
  "HAT22", #hat22
  "HMGB15", #hmg
  "LBD4", #lbd4
  "MYB29", #myb29,
  "NAC102",#nac102
  "RAV2") #rav2/edf2


TDNA <- cbind(lines,linesSymbols)

#Get only values for the genotype.
All_logFC <- DE_All[,grep("logFC",colnames(DE_All))]
#All_logFC <- All_logFC[,-grep("1mM",colnames(All_logFC))] 
# get p values
All_adjpVal <- DE_All[,grep("adj.P.Val",colnames(DE_All))]
signMatrix <- All_adjpVal 
signMatrix[All_adjpVal < pValCut] <- "*"
signMatrix[All_adjpVal > pValCut] <- NA


hmData <- as.matrix(All_logFC[lines,])
# Change AGIs to Symbols
tmp <- cbind("AGI"=lines,
             "Symbol"=linesSymbols)

idxNA <- is.na(tmp[,"Symbol"])
tmp[idxNA,"Symbol"] <- tmp[idxNA,"AGI"]
rownames(hmData) <- tmp[,"Symbol"]

hmData[idxNA,] <- 0

# colors <- colorRampPalette(
#   rev(c("indianred4",
#         "orange2",
#         "goldenrod1",
#         "white",
#         "powderblue",
#         "skyblue4",
#         "steelblue4"))
# )

colors <- colorRampPalette(
  RColorBrewer::brewer.pal(20,"RdYlBu")
)


titulo <- paste0(outDir,"/Heatmap_TDNA_",shortName,".pdf")
pdf(titulo,paper = "USr")
#####
lmat=rbind(c(9, 4, 3), 
           c(2, 7, 8),
           c(6, 5, 1))
lhei=c(0.15,0.2,0.5)
lwid=c(0.2,0.3,0.7)


colnames(hmData) <- tolower(gsub(".logFC","",colnames(hmData)))
# colnames(hmData)[3] <- "arf18-1"
# colnames(hmData)[4] <- "arf18-2"
# colnames(hmData)[8] <- "hmgb15"
# colnames(hmData)[12] <- "rav2"
heatmap.2(hmData, trace = "none", na.rm = T,  
          col=rev(colors(125)), density.info = "none",
          margins = c(8,8), keysize =1.5,
          scale = "col",
          Rowv = F, Colv = F, dendrogram = "none",
          lmat = lmat, lhei = lhei, lwid = lwid,
          cexRow = 0.6, cexCol = 0.6, na.color = "grey98",
          cellnote = signMatrix[tmp[,"AGI"],], notecol = "black",
          main = "Change in Expression of \ngenes with T-DNA Insertion")

dev.off()
############################################################


############################################################
############################## Enrichment tests
############################ Pooled
uniqGenotypes <- unique(gsub("_.*$","",colnames(contrastMatrix)))
deGenesList2 <- list()
for (singleGeno in uniqGenotypes){
  deGenesList2[[singleGeno]] <- Reduce(union,sapply(DESignificant[grep(singleGeno,names(DESignificant))], rownames))
}

############################ Individual
## Filter logFC of significant genes    
pValCutOff <- 0.05
cat ("Filtering by pValue  <", pValCutOff , " \n")
logfcExprs <- lapply(DEList,function(x){x[x[,grep("adj.P.Val",colnames(x))] < pValCutOff ,grep("logFC",colnames(x)),drop=F]})
deGenesList <- lapply(logfcExprs, rownames)


sapply(deGenesList,length)
sapply(deGenesList2,length)


#### Read network
## Some functions are stored in a different script
source("Scripts/metaFunctions_forNetworkAnalysis.R")
# Read interaction data
cat ("Reading interaction data:","\n")
network <- read.csv("meta/NUENetwork_01_May_17_AGaudinier.csv",header = T,stringsAsFactors = F)

### Get unique gene names and symbols
TF_AGI <- network$TF_Gene_Name
names(TF_AGI) <- network$TF_AGI
#
Target_AGI <- network$Promoter_Gene_Name
names(Target_AGI) <- network$Promoter_AGI  
## Combine and get unique genes
GenesInNetwork <- c(TF_AGI,Target_AGI)
GenesInNetwork <- GenesInNetwork[!duplicated(GenesInNetwork)]
#########
cat (nrow(network),"interactions \n")
cat (length(GenesInNetwork)," unique genes in network \n")
# ----


#sapply(logfcExprs,nrow)
DE_Full <- condenseListTables(DEList)

Universe <- rownames(DE_Full)
length(Universe)

contingencyTable2 <- lapply(deGenesList,function(x){
  DEgenes <- x;
  print( length(DEgenes));
  nonDE <- Universe[!Universe %in% (x)];
  #
  cbind(
    'DE_inSet'=length(DEgenes[DEgenes%in%names(GenesInNetwork)]),
    'DE_notinSet'=length(DEgenes[!DEgenes%in%names(GenesInNetwork)]),
    'nonDE_inSet'=length(DEgenes[nonDE%in%names(GenesInNetwork)]),
    'nonDE_notinSet'=length(DEgenes[!nonDE%in%names(GenesInNetwork)]))
  
})

contingencyTable2 <- as.data.frame(do.call(rbind,contingencyTable2))
rownames(contingencyTable2) <- names(logfcExprs)
head(contingencyTable2)

## Sums of the elements should equal to the number of genes analyzed (universe)
if (all(apply(contingencyTable2,1,function(x){sum(x)}) == length(Universe))){
  cat("Contingency table for Fisher's exact test: OK ")
} else{
  cat("Contingency table for Fisher's exact test: Check ")
}

## Perform fisher exact test
fisherResults2 <- apply(contingencyTable2,1,function(x){fisher.test(matrix(x,nrow = 2))})
names(fisherResults2)
fisherPVals2 <- sapply(fisherResults2,function(x){x$p.value})
fisherPVals2 < 0.05 #Which are < 0.05
signifFisher <- ifelse(fisherPVals2 < 0.05,"*","NS")
contingencyTable2 = cbind(contingencyTable2,"pVals"=fisherPVals2,"signif"=signifFisher)

titulo <- paste0(outDir,"/FisherTest_IndividualContrasts_",shortName,".csv")
write.csv(x=contingencyTable2,file = titulo,quote = F,row.names = T)


########## Pooled
contingencyTable1 <- lapply(deGenesList2,function(x){
  DEgenes <- x;
  print( length(DEgenes));
  nonDE <- Universe[!Universe %in% (x)];
  #
  cbind(
    'DE_inSet'=length(DEgenes[DEgenes%in%names(GenesInNetwork)]),
    'DE_notinSet'=length(DEgenes[!DEgenes%in%names(GenesInNetwork)]),
    'nonDE_inSet'=length(DEgenes[nonDE%in%names(GenesInNetwork)]),
    'nonDE_notinSet'=length(DEgenes[!nonDE%in%names(GenesInNetwork)]))
  
})

contingencyTable1 <- as.data.frame(do.call(rbind,contingencyTable1))
rownames(contingencyTable1) <- names(deGenesList2)
head(contingencyTable1)

## Sums of the elements should equal to the number of genes analyzed (universe)
if (all(apply(contingencyTable1,1,function(x){sum(x)}) == length(Universe))){
  cat("Contingency table for Fisher's exact test: OK ")
} else{
  cat("Contingency table for Fisher's exact test: Check ")
}

## Perform fisher exact test
fisherResults1 <- apply(contingencyTable1,1,function(x){fisher.test(matrix(x,nrow = 2))})
names(fisherResults1)
fisherPVals1 <- sapply(fisherResults1,function(x){x$p.value})
fisherPVals1 < 0.05 #Which are < 0.05
##
signifFisher1 <- ifelse(fisherPVals1 < 0.05,"*","NS")
contingencyTable1 = cbind(contingencyTable1,"pVals"=fisherPVals1,"signif"=signifFisher1)


titulo <- paste0(outDir,"/FisherTest_PooledGenes_perGenotype_",shortName,".csv")
write.csv(x=contingencyTable1,file = titulo,quote = F,row.names = T)
