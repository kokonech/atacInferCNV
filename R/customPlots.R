
#' Function to plot CNV blocks
#'
#' This function creates a plot for CNV assigned/identified subclones
#' @param resDir Path to the result directory with input
#' @return Invisibly returns NULL.
#' @examples
#' resPath = tempfile()
#' inPath = system.file("extdata", "MB183_ATAC_subset.tsv.gz", package = "atacInferCnv")
#' sAnn = system.file("extdata", "MB183_ATAC_subset.CNV_blocks_ann.txt", package = "atacInferCnv" )
#' prepareAtacInferCnvInput(inPath,sAnn,resPath, targColumn = "cnvBlock", ctrlGrp = "Normal")
#' runAtacInferCnv(resPath)
#' plotCnvBlocks(resPath)
#' @export

plotCnvBlocks <- function( resDir) {

  cnvDir <- paste0(resDir,"/sample_infercnv")
  if (!(dir.exists(cnvDir))) {
    stop("The result directory with InferCNV input does not exist:",resDir)
  }

  infercnv_obj <- readRDS(paste0(cnvDir, "/run.final.infercnv_obj"))

  message("Load InferCNV result...")

  obsFile <- paste0(cnvDir,"/infercnv.observations.txt")
  if (file.exists(obsFile)) {
    cnvMtx <- read.table(obsFile,check.names = FALSE)
  } else {
    warning("Output file infercnv.observations.txt does not exist in result folder, using expr matrix.\n
         Update InferCnv version (>= 1.3.3 ) to support normalized output.")
    cIds <- unlist(infercnv_obj@observation_grouped_cell_indices)
    cnvMtx <- infercnv_obj@expr.data[ ,cIds]
  }
  gnMtx <- infercnv_obj@gene_order

  cMax <- max(cnvMtx)
  cMin <- min(cnvMtx)
  cMean <- mean(rowMeans(cnvMtx))

  # draw lines
  curChr <- "1"
  breaks <- c()
  for (i in seq_len(nrow(gnMtx))) {
    chrId <- as.character(gnMtx[i,1] )
    if (chrId != curChr) {
      #print(chrId)
      breaks <- c(breaks,i)
      curChr <- chrId
    }
  }

  chrBorders <- c(1,breaks,nrow(gnMtx))
  names(chrBorders) <- c(seq_len(22),"")

  rbPal <- colorRampPalette(c('red','green'))
  chrBreaks <- breaks
  if (infercnv_obj@options$k_obs_groups  == 1) {
    blocks <- c("Full", names(infercnv_obj@tumor_subclusters$subclusters))
  } else {
    hcBlocks <- cutree(infercnv_obj@tumor_subclusters$hc$all_observations, infercnv_obj@options$k_obs_groups )
    hcBlocksAdj <- paste0("C",hcBlocks)
    names(hcBlocksAdj) <- names(hcBlocks)
    blocks <- c("Full", unique(hcBlocksAdj))
  }
  resName <- paste0(cnvDir,"/subclone_CNV_plot.pdf")
  pdf(resName, width = 14, height= 6)

  for (targ in blocks) {
    message(targ)
    if (targ == "Full") {
        vals<-rowMeans(cnvMtx)
    } else {
      if (infercnv_obj@options$k_obs_groups  == 1) {
        cellIds <-names(infercnv_obj@tumor_subclusters$subclusters[[targ]][[1]])
        if (sum(cellIds %in% colnames(cnvMtx)) == 0) {
          message("Skipping...")
          next
        }
      } else {
        cellIds <- names(hcBlocksAdj)[hcBlocksAdj == targ]
      }
      vals <- rowMeans(cnvMtx[,cellIds])
    }
    adjCols <- rbPal(10)[as.numeric(cut(vals,breaks = 10))]

    p<- plot(vals,pch=20,cex=0.3,
       main=paste0("inferCNV signal:",targ), col=adjCols,ylim=c(cMin,cMax),
       xlab="Chromosome", ylab="Modified expr",  axes=FALSE)
    p <- p + abline(h=1, col="grey",lwd=1, lty=2)
    if (length(chrBreaks) > 0) {
      p <- p + abline(v=c(1,breaks,nrow(gnMtx)),col="black",lwd=1,lty=2)
      p <- p + axis(side=1,at=chrBorders,labels=names(chrBorders),cex.axis=0.45 )
    }

    # print REQUIRED for ggplot2 output
    #abline(v=markerLoci, col="red", lwd=0.2, lty=2)
    p <- p + axis(side=2,seq(round(cMin,2),round(cMax,2),0.01) )
    print(p)
  }

  dev.off()

  invisible(NULL)
}
