# load relevant libraries
library(matrixStats)
library(readr)

# gene annotation file
gene_anno = read.delim("/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/LCL.featureCounts.features.tsv",as.is=T)
testCombinations = NULL

# make chunks
nGenes = 50
startPos = 0
endOffset = 1000000000
sink("/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/chunks.txt")
for(chr in unique(gene_anno$chromosome)){
  annotationRel = gene_anno[which(gene_anno$chromosome==chr),]
  annotationRel = annotationRel[order(annotationRel$start,annotationRel$end),]
  # First go through the list to fix genes to they all touch.
  annotationRel$start[1] = startPos
  for(i in 2:nrow(annotationRel)){
    if(i == nrow(annotationRel)){
      annotationRel$end[i] = annotationRel$end[i]+endOffset
    }
    # If "overlapping" than we don't need to do anything.
    if((annotationRel$start[i]>annotationRel$end[i-1])){
      distance = (annotationRel$start[i]-annotationRel$end[i-1])/2
      annotationRel$start[i] = ceiling(annotationRel$start[i]-distance)
      annotationRel$end[i-1] = floor(annotationRel$end[i-1]+distance)
    }
  }
  chunks = seq(1, nrow(annotationRel),nGenes)
  # Need to add the last one as a separate entry.
  if(chunks[length(chunks)] < nrow(annotationRel)){
    chunks = c(chunks,(nrow(annotationRel)+1))
  }
  for(i in 1:(length(chunks)-1)){
    print(paste(chr,":",annotationRel$start[chunks[i]],"-",annotationRel$end[(chunks[i+1]-1)],sep=""))
  }
}
sink()


df = read.csv("/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/chunks.txt",sep="\t",header=F)
df$V1 = gsub("\\[1\\] ","",df$V1)
write.table(df,"/share/ScratchGeneral/anncuo/OneK1K/Sex_interactions/chunks.txt",sep="\t",col.names=F,row.names=F,quote=F)
