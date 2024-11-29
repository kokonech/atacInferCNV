


prepareAnalysis <- function(mb, resDir, sId, annData, targColumn, ctrlObj) {
  if (!is.null(annData)) {
    cIds = intersect(rownames(annData), rownames(mb@meta.data))
    if (length(cIds) == 0) {
      print("Error: no overlapping cells IDs in annotation")
    }
    print("Adjust for custom annotation")
    mb <- subset(mb, cells= cIds)
    #print(summary(rownames(mb@meta.data) == cIds))
    mb@meta.data <- cbind(mb@meta.data,  annData[rownames(mb@meta.data),targColumn])
    colnames(mb@meta.data)[ncol(mb@meta.data)] <- targColumn
  } else {
    mb <- NucleosomeSignal(mb)
    mb <- TSSEnrichment(mb)

    resName = paste0(resDir, sId,"_VlnPlot_QC.pdf")
    pdf(resName, width=8,height=6)
    VlnPlot(
      object = mb,
      features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
      pt.size = 0
    )
    dev.off()

    # filter out low quality cells
    mb <- subset(
      x = mb,
      subset = nCount_ATAC < 100000 &
        nCount_ATAC > 1000 &
        nucleosome_signal < 2 &
        TSS.enrichment > 1
    )

    pdf(paste0(resDir, sId,"_VlnPlot_QC.after_filter.pdf"), width=8, height=6)
    VlnPlot(mb, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),     pt.size=0)
    dev.off()
  }

  if (!is.null(ctrlObj)) {
    print("Merge with external control object")
    # adjusted number of control cells should no go over 33% of number of tumor cells
    expNumCtrlCells <- round( 0.33 * ncol(mb))
    if (ncol(ctrlObj) > expNumCtrlCells) {
      print(paste("Adjust external control, decrease num cells to",expNumCtrlCells))
      ctrlObj <- ctrlObj[ , 1:expNumCtrlCells]
    }
    #print(ctrlObj)
    #print(ctrlObj@assays$ATAC@ranges)
    #print(mb@assays$ATAC@ranges)
    # NOTE : not most optimal merging - regions of mb used as a reference
    mb.merged <- merge(mb, y=ctrlObj)
    mb.merged@meta.data[ , targColumn] <- c(mb@meta.data[,targColumn],
                                            rep("ExtControl",ncol(ctrlObj)) )
    mb <- mb.merged
  }


  print("Normalization...")

  mb <- RunTFIDF(mb)
  mb <- FindTopFeatures(mb, min.cutoff = 'q0')
  mb <- RunSVD(mb)

  print("Dimensional reduction...")

  ndim = 30 # default 30

  mb <- RunUMAP(object = mb, reduction = 'lsi', dims = 2:ndim)
  mb <- FindNeighbors(object = mb, reduction = 'lsi', dims = 2:ndim)
  mb <- FindClusters(object = mb, verbose = FALSE, algorithm = 3)

  pdf(paste0(resDir,sId,"_UMAP.pdf"),width = 8, height = 6)
  print(DimPlot(object = mb, pt.size=1, label=T))
  if (nchar(targColumn) > 0) {
    print(DimPlot(mb, reduction = "umap", label = TRUE,group.by = targColumn))
  }
  dev.off()

  mb
}


saveCnvInput <- function(mb,resDir, sId, targColumn) {
  saveRDS(mb, paste0(resDir,sId,"_obj.RDS" ))

  # for InferCNV

  annTable <- mb@meta.data
  annTable2 <- annTable[,targColumn,drop=F]

  if (targColumn == "seurat_clusters" ) {
      annTable2$seurat_clusters <- paste0("cl",annTable2$seurat_clusters)
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

#' Prepare analysis for the CNV calling from ATAC data
#'
#' @param dataPath Path to the input data in 10X format
#' @param annData Path to annotation of the cells. Should have column
#' @param resDir Path to the result directory
#' @param sId Result name. Default: "Sample"
#' @param targColumn Name of the target column in annotation. Default: "CellType"
#' @param ctrlGrp Name for the reference control cell type. Default: "Normal"
#' @param ctrlObj Seurat/Signac object to use as non-tumor control. Default: NULL
#' @param meta Set TRUE to use meta cells, default FALSE
#'
#' @return NULL
#' @export
#'
prepareAtacInferCnvInput <- function(dataPath,
                                     annPath,
                                     resDir, sId = "sample",
                                     targColumn = "CellType",
                                     ctrlGrp = "Normal",
                                     ctrlObj = NULL,
                                     metaCells = FALSE) {

  print("Loading input...")
  if (!(file.exists(dataPath))) {
    stop("Input data folder is not found:",dataPath)

  }
  countsPath = paste0(dataPath,"/filtered_feature_bc_matrix.h5")
  if (!(file.exists(countsPath))) {
    stop("Input feature counts matrix is not found:",countsPath)
  }
  counts <- Read10X_h5(countsPath)

  fragpath = paste0(dataPath, "/atac_fragments.tsv.gz")
  if (!(file.exists(fragpath))) {
    stop("Input fragments loci is not found:",countsPath)
  }
  #annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  #seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

  # create ATAC assay and add it to the object
  chrom_assay  <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath
  )

  mb <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    project = sId
  )

  resDir = paste0(resDir,"/") # make sure subfolder usage
  if (!(dir.exists(resDir))) {
    print(paste("Creating result directory:", resDir))
    dir.create(resDir)
  }

   if (!(file.exists(annPath))) {
    stop(paste("Annotation file is not found:",annPath))

  }
  annData <- read.delim(annPath)

  if (! (targColumn %in% colnames(annData) ) ) {
    stop(paste("Required annotation column is not available:", targColumn))
  } else {
    print(paste("Using target annotation column:",targColumn))
    #print(class(annData[, targColumn]))
    annInfo = summary(as.factor(annData[, targColumn]))
    print(annInfo)
    print(ctrlObj)
    if (! (ctrlGrp %in% names(annInfo)) ) {
        stop(paste("Non-tumor control group is not found in annotation:", ctrlGrp))
    }
    if (!is.null(ctrlObj)) {
      if (!inherits(ctrlObj, "Seurat")) {
        stop(paste0("Non-tumor external control input is not Seurat object!"))
      }
      print("Using external control:")
      print(ctrlObj)
      ctrlGrp = "ExtControl"
    }

  }

  print("Prepare input data...")
  mb <- prepareAnalysis(mb, resDir, sId, annData, targColumn, ctrlObj)


  print("Save signal...")
  saveCnvInput(mb, resDir, sId, targColumn)

  if (metaCells) {
    print("Forming meta-cells...")
    mb$seurat_clusters <- as.factor(mb@meta.data[, targColumn])
    print(summary(mb$seurat_clusters))
    extractMetacells(resDir, sId, mb)
  }

  print("Write configuration...")
  writeConfig(resDir, sId, ctrlGrp, metaCells)
  print("Prepared input")

}
