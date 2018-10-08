install.packages("rJava")
source("https://bioconductor.org/biocLite.R")
biocLite("BridgeDbR")

library("BridgeDbR", lib.loc="~/R/win-library/3.4")
library("RCurl", lib.loc="~/R/win-library/3.4")

#download .bridge file from Figshare
mbmaps = loadDatabase("metabolites_20180201.bridge")
setwd("C:/Users/marvin.martens/Documents/R")
CASRNs <- read.table("~/R/ListofCasrns.txt", quote="\"", comment.char="")

mapperFunction = function(x) {
  mappings = paste(map(mbmaps, "Ca", x, "Ce"),collapse = ",")
  return(mappings)
}

CHEBI = cbind(sapply(as.character(CASRNs[,1]), mapperFunction))
write.table(CHEBI, "Mapping_results.txt", sep="\t",col.names = FALSE, quote = FALSE)
