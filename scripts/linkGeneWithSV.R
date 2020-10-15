options(stringsAsFactors = F)
args = commandArgs(TRUE)
geneName = args[1]
geneChr = args[2]
geneStart = as.numeric(args[3])
geneEnd = as.numeric(args[4])
fileName = args[5]

getSVs <- function(geneName,geneChr,geneStart,geneEnd,fileName) {
  gzFile = gzfile(fileName,'rt')  
  SVTable = read.table(file = gzFile, header = F, sep = '\t', quote = "")
  
  SVChrTable = SVTable[SVTable[,1]==geneChr,]
  SVgeneTable = SVChrTable[(((SVChrTable[,2]>geneEnd)|(SVChrTable[,3]<geneStart))==0),]
  if (nrow(SVgeneTable)!=0) {
    outputTSV = data.frame(variant_id = SVgeneTable[,4], gene_id = rep(geneName,nrow(SVgeneTable)))
    write.table(SVgeneTable, file=paste(fileName,geneName,'tsv',sep = '.'), quote=FALSE, sep='\t', col.names = F, row.names = F)
  }
}

getSVs(geneName,geneChr,geneStart,geneEnd,fileName)
