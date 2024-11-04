saveCnvInput <- function(mb, resDir,sId,  targColumn="") {
  saveRDS(mb, paste0(resDir,sId,"_obj.RDS" ))
  
  # for InferCNV
  
  annTable <- mb@meta.data
  if (nchar(targColumn) == 0) {
    annTable2 <- annTable[,"seurat_clusters",drop=F]
    annTable2$seurat_clusters <- paste0("cl",annTable2$seurat_clusters)
  } else {
    annTable2 <- annTable[,targColumn,drop=F]
  }
  write.table(annTable2, paste0(resDir,sId,"_cnv_ann.txt") ,col.names = F,sep="\t",quote=F)
  
  gz1 <- gzfile(paste0(resDir,sId,"_raw_counts.txt.gz"), "w")
  rawCounts <- as.matrix(mb@assays$ATAC@counts)
  
  rawCounts <- rawCounts[ grep("chr",rownames(rawCounts)),]
  write.table(rawCounts,gz1,sep="\t",quote=F)
  close(gz1)
  
  peakRegions <- GRanges(sub("-",":",rownames(rawCounts), fixed=TRUE))
  seqlevels(peakRegions) <- paste0("chr",c(1:22,"X","Y")) # fix order
  peakRegions <- sort(peakRegions)
  
  
  peakDf <- data.frame(peakRegions)[,1:3]
  rownames(peakDf) <- paste0(  peakDf[,1],"-",peakDf[,2],"-",peakDf[,3] )
  peakDf$seqnames <- gsub("chr","", peakDf$seqnames)
  summary(rownames(rawCounts) %in% rownames(peakDf))
  write.table(peakDf, paste0(resDir,sId,"_cnv_ref.txt") ,col.names = F,sep="\t",quote=F)
}
