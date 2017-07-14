### Mapping efficiency
#########################################################
# github: rodriguezmDNA
# j rodriguez medina
# Brady lab @ ucdavis
# Generate plots on mapping efficiency (%) from alignment scripts
# Last modified
# 2017.07.jrm


#install.packages("reshape")
#install.packages("ggplot2")
library(ggplot2)
library(reshape)

setwd ("~/Desktop/junk/")

### User defined options
option <- "genomebwt1" #cDNAbwt1,genomebwt1,genomebwt2,cDNAbwt2. Choose one.


############### 
logDir <- paste0("logs/bowStats_",option)

logFiles <- list.files(logDir,pattern = ".log",full.names = T)
length(logFiles)
dir.create("images",showWarnings = F)


#Names <- gsub("^str.[0-9]_|_set.[0-9]$|set|_|tr|[.]|bowtie.logo","\\1",basename(logFiles))
Names <- basename(logFiles)
tempIndex = order(gsub("tm|[ab][0-9]$","\\1",Names))
logFiles[tempIndex]


## Save stats into a list
MapStatsList <- list()
for (each in logFiles[tempIndex]){
  
  print (each)
  #
  name <- basename(each)
  name <- gsub(".log","\\1",name)
  
  #
  temp <- read.table(each,header = F,fill=T,comment.char = "",stringsAsFactors = F,as.is = T)
  #head(temp)
  total   <- as.numeric ( temp[1,4]) # Total number of reads
  aligned <- as.numeric ( temp[2,9]) # Aligned reads
  un      <- as.numeric (temp[3,7]) # Unaligned reads
  supressed      <- as.numeric (temp[4,9]) # Unaligned reads
  MapStatsList[[name]] <- c(total, aligned ,un,supressed)
}

MapStats <- do.call("rbind",MapStatsList)
colnames(MapStats) <- c("total", "aligned" ,"unaligned","supressed")

### Barplot
toGraph <- as.data.frame(rbind(MapStats[,c(2,3,4),drop=F]))
#time = gsub('^tm|[ab].*$', '\\1',rownames(toGraph))
#group = gsub('^tm|[0-9]*$', '\\1',rownames(toGraph))
#treatment = gsub('[^ab]', '\\1',toGraph$group)


## Calculate percentage
percentGraph = as.data.frame(prop.table(as.matrix(toGraph),margin = 1))

## 
SaveTable <- cbind(toGraph,"",percentGraph)
colnames(SaveTable) <- c(colnames(toGraph),"",paste0("%",colnames(percentGraph)))

# Melt: 
#
percentGraph$Library <- row.names(percentGraph)
percentGraph = melt(percentGraph,id.vars = 'Library')
#
toGraph$Library <- row.names(toGraph)
toGraph = melt(toGraph,id.vars = 'Library')

# Set other variables to facet
#toGraph$time = gsub('^tm|[ab].*$', '\\1',toGraph$category)
#toGraph$group = gsub('^tm|[0-9]*$', '\\1',toGraph$category)
#toGraph$treatment = gsub('[^ab]', '\\1',toGraph$group)
#
#percentGraph$time = gsub('^tm|[ab].*$', '\\1',percentGraph$category)
#percentGraph$group = gsub('^tm|[0-9]*$', '\\1',percentGraph$category)
#percentGraph$treatment = gsub('[^ab]', '\\1',percentGraph$group)

titulo <- paste0("images/MappingStats_",option,".pdf")
pdf (titulo, paper = "a4r")#width=20, height=12)
mar.default <- c(15,4,8,6) + 0.1


ggplot(toGraph, aes(Library, value)) +   
  geom_bar(aes(fill = variable), position = "stack", stat="identity") +
  labs(title = "Mapping efficiency (counts)") +
  theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=6)) # +
  #facet_wrap(~group,nrow=2,scales = 'free_x')

ggplot(percentGraph, aes(Library, value)) +   
  geom_bar(aes(fill = variable), position = "stack", stat="identity") +
  labs(title = "Mapping efficiency (%)") +
  theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text=element_text(size=6)) # +
  #facet_wrap(~group,nrow=2,scales = 'free_x')

dev.off()

# Write to file 
titulo <- paste0("images/MappingStats_",option,".txt")
write.table(file = titulo,SaveTable,quote = F,col.names = NA,sep="\t")



