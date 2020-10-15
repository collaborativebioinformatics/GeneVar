options(stringsAsFactors = F)
library(stringr)
args = commandArgs(TRUE)
geneName = args[1]
geneChr = args[2]
geneStart = as.numeric(args[3])
geneEnd = as.numeric(args[4])
SVFile = args[5]
annotationFile = args[6]

getSVAndAnnotation <- function(geneName,geneChr,geneStart,geneEnd,SVFile,annotationFile) {
  cat(geneName,'\n')
  geneChr = paste('chr',geneChr,sep='')
  SVFileName = paste(SVFile,geneChr,'POS','tsv','gz',sep = '.')
  gzFile = gzfile(SVFileName,'rt')  
  SVTable = read.table(file = gzFile, header = F, sep = '\t', quote = "")
  close(gzFile)
  SVTable[,2] = as.numeric(SVTable[,2])
  SVTable[,3] = as.numeric(SVTable[,3])
  SVTable = SVTable[!is.na(SVTable[,2]+SVTable[,3]),]
  minPos = sapply(1:nrow(SVTable), function(x) min(SVTable[x,2:3]))
  maxPos = sapply(1:nrow(SVTable), function(x) max(SVTable[x,2:3]))
  SVTable[,1] = paste('chr',SVTable[,1],sep='')
  SVTable[,2] = minPos; SVTable[,3] = maxPos;
  gzFile = gzfile(annotationFile,'rt')  
  annotationTable = read.table(file = gzFile, header = F, sep = '\t', quote = "")
  close(gzFile)
  
  SVChrTable = SVTable[SVTable[,1]==geneChr,]
  SVgeneTable = SVChrTable[(((SVChrTable[,2]>=geneEnd)|(SVChrTable[,3]<=geneStart))==0),]
  
  geneAnnoTab = annotationTable[which(str_detect(annotationTable[,9], geneName)==T),]
  
  eltType = rep('',nrow(SVgeneTable))
  for (i in 1:length(eltType)) {
    startTmp = SVgeneTable[i,2];endTmp = SVgeneTable[i,3]
    geneAnnoTabForSV = geneAnnoTab[(((geneAnnoTab[,4]>=endTmp)|(geneAnnoTab[,5]<=startTmp))==0),]
    if (nrow(geneAnnoTabForSV)>0) {eltType[i] = paste(unique(geneAnnoTabForSV[,3]),collapse = ';')}
  }
  
  if (nrow(SVgeneTable)!=0) {
    outputTSV = data.frame(variant_id = SVgeneTable[,5], gene_id = rep(geneName,nrow(SVgeneTable)), elt_type = eltType)
    write.table(outputTSV, file=paste(SVFile,geneName,'tsv',sep = '.'), quote=FALSE, sep='\t', col.names = T, row.names = F)
  }
}

getSVAndAnnotation(geneName,geneChr,geneStart,geneEnd,SVFile,annotationFile)
