#############################################
### Create a larger table with all columns 
condenseListTables <- function(listDFs) {
  cat ("-- make condenseListTables function called \n")
  # First make an empty table
  uniqRows <- Reduce(union,lapply(listDFs,rownames))
  uniqCols <- unlist(sapply(listDFs, colnames))
  zeroTable <- as.data.frame(matrix(0,
                                    nrow = length(uniqRows),
                                    length(uniqCols)),
                             row.names =uniqRows)
  colnames(zeroTable) = uniqCols
  # Then Fill it
  for (each in names(listDFs)){
    cat ("Filling binary table:", each,"\n")
    zeroTable[rownames(listDFs[[each]]),
              colnames(listDFs[[each]])] <- listDFs[[each]]
  }
  #zeroTable <- apply(zeroTable,c(1,2),as.numeric)
  #zeroTable[is.na(zeroTable)] <- 0
  return(zeroTable)
  cat ("-- binary table done \n")
  cat ("\n")
}