

###
volcanoPlots <- function(DEList,pValCut=0.05,logCut=2,fcColumn=NA){ 
  
  for (contrast in names(DEList)) {
    print (contrast)#}
    
    # --
    DE <- DEList[[contrast]]
    head(DE)
    
    if (is.na(fcColumn)){
      cat ("Using first column as logFC values \n")
      res <- DE[,colnames(DE)[c(1,grep("adj.P.Val|Symbol",colnames(DE)))]] #Use first column as logFC data
    } else {
      cat ("Using",fcColumn,"column as logFC values \n")
      res <- DE[,c(fcColumn,colnames(DE)[grep("adj.P.Val|Symbol",colnames(DE))])]
    }
    
    
    ## Temporary subset and change names
    colnames(res) <- c('logFC','adjpVal','Symbol')
    #res$Symbol <- gsub("\\|.*$","",res$Symbol) #If symbol has multiple names, use only the first one.
    #head(res)
    
    # Make a basic volcano plot
    with(res, plot(logFC, -log10(adjpVal), pch=".", 
                   main=paste("Volcano plot of",contrast,sep=" "), 
                   ylim=c(0,max(-log10(res$adjpVal))+5),  col="lightgray", 
                   xlim=c(min(res$logFC)-4,max(res$logFC)+4)))
    
    colors <- c("darkslategray4","darksalmon","lightseagreen")
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(res, adjpVal<pValCut), points(logFC, -log10(adjpVal), pch=".", col=colors[1]))
    with(subset(res, abs(logFC)>logCut), points(logFC, -log10(adjpVal), pch=".", col=colors[2]))
    with(subset(res, adjpVal<pValCut & abs(logFC)>logCut), points(logFC, -log10(adjpVal), pch="o", col=colors[3]))
    # Label points with the textxy function from the calibrate plot
    
    with(subset(res, adjpVal<pValCut & abs(logFC)>logCut), 
         textxy(logFC, -log10(adjpVal), labs=Symbol, cex=.5))
    legend("topright", c(paste0("AdjpVal <",pValCut),
                         paste0("logFC >",logCut),
                         "pass pVal & logFC"),
           fill = colors,border = F,
           cex = 0.8,bty="n")
    abline(v = c(2,-2),col="skyblue",lty = 2)
  }
}