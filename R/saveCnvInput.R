saveCnvInput <- function(mb, resDir,sId,  targColumn="") {

   # for InferCNV

  annTable <- mb@meta.data
  if (nchar(targColumn) == 0) {
    annTable2 <- annTable[,"seurat_clusters",drop=FALSE]
    annTable2$seurat_clusters <- paste0("cl",annTable2$seurat_clusters)
  } else {
    annTable2 <- annTable[,targColumn,drop=FALSE]
  }
  write.table(annTable2, paste0(resDir,sId,"_cnv_ann.txt") ,col.names = FALSE,sep="\t",quote=FALSE)

  gz1 <- gzfile(paste0(resDir,sId,"_raw_counts.txt.gz"), "w")
  rawCounts <- as.matrix(mb@assays$ATAC@counts)

  rawCounts <- rawCounts[ grep("chr",rownames(rawCounts)),]
  write.table(rawCounts,gz1,sep="\t",quote=FALSE)
  close(gz1)

  peakRegions <- GRanges(sub("-",":",rownames(rawCounts), fixed=TRUE))
  seqlevels(peakRegions) <- paste0("chr",c(seq_len(22),"X","Y")) # fix order
  peakRegions <- sort(peakRegions)


  peakDf <- data.frame(peakRegions)[,seq_len(3)]
  rownames(peakDf) <- paste0(  peakDf[,1],"-",peakDf[,2],"-",peakDf[,3] )
  peakDf$seqnames <- gsub("chr","", peakDf$seqnames)
  summary(rownames(rawCounts) %in% rownames(peakDf))
  write.table(peakDf, paste0(resDir,sId,"_cnv_ref.txt") ,col.names = FALSE,sep="\t",quote=FALSE)
}
