
mergeWithControl <- function(mb, ctrlObj, targColumn) {
  # adjusted number of control cells should no go over 33% of number of tumor cells
  expNumCtrlCells <- round( 0.33 * ncol(mb))
  if (ncol(ctrlObj) > expNumCtrlCells) {
    print(paste("Adjust external control, decrease num cells to",expNumCtrlCells))
    ctrlObj <- ctrlObj[ , 1:expNumCtrlCells]
  }
  #print(ctrlObj)
  #print(ctrlObj@assays$ATAC@ranges)
  #print(mb@assays$ATAC@ranges)
  # NOTE : not most optimal merging - regions of x are used as main

  mb.merged <- merge(mb, y=ctrlObj)

  mb.merged@meta.data[ , targColumn] <- c(mb@meta.data[,targColumn],
                                          rep("ExtControl",ncol(ctrlObj)) )
  mb.merged

}


prepareAnalysis <- function(mb, resDir, sId, annData, targColumn, ctrlObj) {
  if (!is.null(annData)) {
    cIds = intersect(rownames(annData), rownames(mb@meta.data))
    if (length(cIds) == 0) {
      print("Error: no overlapping cells IDs in annotation")
    }
    print("Adjust for custom annotation")
    mb <- subset(mb, cells= cIds)
    #mb <- mb[ ,cIds]
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
    mb <- mergeWithControl(mb, ctrlObj, targColumn )
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

  saveRDS(mb, paste0(resDir,sId,"_obj.RDS" ))

  mb
}


saveCnvInput <- function(mb,resDir, sId, targColumn) {

  # for InferCNV

  annTable <- mb@meta.data
  annTable2 <- annTable[,targColumn,drop=F]

  #if (targColumn == "seurat_clusters" ) {
  #    annTable2$seurat_clusters <- paste0("cl",annTable2$seurat_clusters)
  #}
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
#' @param annPath Path to annotation of the cells in the target column
#' @param resDir Path to the result directory
#' @param inObj Precomputed Seurat/Signac object with required input data
#' @param sId Result name. Default: "Sample"
#' @param targColumn Name of the target column in annotation. Default: "CellType"
#' @param ctrlGrp Name for the reference control cell type. Could be several names, separated by comma. Default: "Normal"
#' @param ctrlObj Seurat/Signac object to use as non-tumor control. Default: NULL
#' @param binSize Apply custom bin size to combine signals in windows for CNV calling
#' e.g. 500000 for 500 KBp. Default: NULL (not use this option)
#' @param chromLength Numeric vector of chromosome sizes, specific for genome. Default: NULL
#' @param metaCells Set TRUE to use meta cells (n=5 cells will be used to merge) or assign a number of cells. Default: FALSE
#'
#' @return NULL
#' @export
#'
prepareAtacInferCnvInput <- function(dataPath = "",
                                     annPath = "",
                                     resDir = "", inObj = NULL, sId = "sample",
                                     targColumn = "CellType",
                                     ctrlGrp = "Normal", ctrlObj = NULL,
                                     binSize = NULL, chromLength = NULL,
                                     metaCells = FALSE) {

  print("Loading input...")
  if (nchar(resDir) == 0) {
    stop("Path to result is not provided!")
  }
  if (is.null(inObj)) {
    if (!(file.exists(dataPath))) {
      stop("Input data folder is not found:",dataPath)

    }
    if (dir.exists(dataPath)) {
      print("Input is directory, assuming 10X data folder...")
      # scMulti-omics
      countsPath = paste0(dataPath,"/filtered_feature_bc_matrix.h5")
      fragpath = paste0(dataPath, "/atac_fragments.tsv.gz")
      # scAtac
      countsPath2 = paste0(dataPath,"/filtered_peak_bc_matrix.h5")
      if (!(file.exists(countsPath))) {
        if (file.exists(countsPath2)) {
          print("10X scATAC format identified")
          countsVals = Read10X_h5(countsPath2)
          fragpath = paste0(dataPath, "/fragments.tsv.gz")
        } else {
          stop(paste("Input feature counts matrix is not found in path:",dataPath,
             "Expected formats: filtered_feature_bc_matrix.h5 or filtered_peak_bc_matrix.h5"))
        }
      } else {
        print("10X scMulti-omics format identified")
        countsVals <- Read10X_h5(countsPath)$Peaks
      }


      if (!(file.exists(fragpath))) {
        stop("Input fragments loci is not found:",fragpath)
      }

      # create ATAC assay and add it to the object
      chrom_assay  <- CreateChromatinAssay(
        counts = countsVals,
        sep = c(":", "-"),
        fragments = fragpath
      )
    } else {
      print("Input is a file, assuming peaks matrix...")
      countsVals <- read.delim(dataPath,check.names = F)
      chrom_assay  <- CreateChromatinAssay(
        counts = countsVals,
        sep = c(":", "-")
      )
    }

    mb <- CreateSeuratObject(
      counts = chrom_assay,
      assay = "ATAC",
      project = sId
    )

    if (!(file.exists(annPath))) {
      stop(paste("Annotation file is not found:",annPath))

    }
    annData <- read.delim(annPath)

    if (! (targColumn %in% colnames(annData) ) ) {
      stop(paste("Required annotation column is not available:", targColumn))
    } else {
      print(paste("Using target annotation column:",targColumn))
      #print(class(annData[, targColumn]))
      if (is.null(ctrlObj)) {
        annInfo = summary(as.factor(annData[, targColumn]))
        print(annInfo)
        ctrlStatus = as.numeric(str_split_1(ctrlGrp,pattern = ",") %in% names(annInfo))
        if (any(ctrlStatus == 0) ) {
          stop(paste("Non-tumor control group is not found in annotation:", ctrlGrp))
        }
      } else {
        if (!inherits(ctrlObj, "Seurat")) {
          stop(paste0("Non-tumor external control input is not Seurat object!"))
        }
        print("Using external control (assigned as ExtControl):")
        print(ctrlObj)
        ctrlGrp = "ExtControl"
      }

    }
  } else {
      print("Using existing Signac/Seurat object")
      if (!inherits(inObj, "Seurat")) {
        stop(paste0("Pre-computed input object is not Seurat object!"))
      }
      #if (!is.null(ctrlObj)) {
      #  stop(paste0("Usage of external reference is not supported if pre-computed object is provided."))
      #}
      mb <- inObj
      annData <- inObj@meta.data
      if (! (targColumn %in% colnames(annData) ) ) {
        stop(paste("Required annotation column is not available in pre-computed input object:", targColumn))
      }
      if (is.null(ctrlObj)) {
        annInfo = summary(as.factor(annData[, targColumn]))
        print(annInfo)
        ctrlStatus = as.numeric(str_split_1(ctrlGrp,pattern = ",") %in% names(annInfo))
        if (any(ctrlStatus == 0) ) {
          stop(paste("Non-tumor control group is not found in annotation:", ctrlGrp))
        }
      } else {
         print("Merge with external control object")
         mb <- mergeWithControl(mb,ctrlObj,targColumn)
         ctrlGrp = "ExtControl"
         #print(head(mb@meta.data))
      }
  }

  resDir = paste0(resDir,"/") # make sure subfolder usage
  if (!(dir.exists(resDir))) {
    print(paste("Creating result directory:", resDir))
    dir.create(resDir)
  }


  if (!is.null(binSize)) {
    if (is.null(chromLength)) {
      stop("The chromosomes length vector is not provided for the bin size adjustment!")
    }
  }

  if (is.null(inObj)) {
    print("Prepare input data...")
    mb <- prepareAnalysis(mb, resDir, sId, annData, targColumn, ctrlObj)
  }

  print("Save signal...")
  saveCnvInput(mb, resDir, sId, targColumn)
  if (!is.null(binSize)) {
    print(paste("Re-format input signal matrix for bins of size", binSize))
    mb <- aggregateBins(mb, resDir, sId, binSize, chromLength)
  }
  #print(head(mb@meta.data))
  if (metaCells) {
    print("Forming meta-cells...")
    print(targColumn)
    if (is.numeric(metaCells)) {
      print("Using custom meta-cell count...")
      metaCount = metaCells
    } else {
      print("Using default meta-cell count...")
      metaCount = 5
    }
    extractMetacells(resDir, sId, mb, targColumn, metacell_content = metaCount )
  }

  print("Write configuration...")
  writeConfig(resDir, sId, ctrlGrp, binSize, metaCells)
  print("Prepared input")

}
