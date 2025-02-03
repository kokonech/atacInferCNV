
#' Function to plot CNV blocks
#'
#' This function creates a plot for CNV assigned/identified subclones
#' @param resDir Path to the result directory with input
#'
#' @export
plotCnvBlocks <- function( resDir) {

  cnvDir = paste0(resDir,"/sample_infercnv")
  if (!(dir.exists(cnvDir))) {
    stop("The result directory with InferCNV input does not exist:",resDir)
  }

  infercnv_obj <- readRDS(paste0(cnvDir, "/run.final.infercnv_obj"))

  print("Load InferCNV result...")

  cnvMtx <- read.table(paste0(cnvDir,"/infercnv.observations.txt"),check.names = F)
  gnMtx <- infercnv_obj@gene_order

  cMax = max(cnvMtx)
  cMin = min(cnvMtx)
  cMean = mean(rowMeans(cnvMtx))

  # draw lines
  curChr <- "1"
  breaks <- c()
  for (i in 1:nrow(gnMtx)) {
    chrId <- as.character(gnMtx[i,1] )
    if (chrId != curChr) {
      #print(chrId)
      breaks <- c(breaks,i)
      curChr <- chrId
    }
  }

  chrBorders <- c(1,breaks,nrow(gnMtx))
  names(chrBorders) <- c(1:22,"")

  rbPal <- colorRampPalette(c('red','green'))
  chrBreaks = breaks
  if (infercnv_obj@options$k_obs_groups  == 1) {
    blocks <- c("Full", names(infercnv_obj@tumor_subclusters$subclusters))
  } else {
    hcBlocks <- cutree(infercnv_obj@tumor_subclusters$hc$all_observations, infercnv_obj@options$k_obs_groups )
    hcBlocksAdj <- paste0("C",hcBlocks)
    names(hcBlocksAdj) = names(hcBlocks)
    blocks <- c("Full", unique(hcBlocksAdj))
  }
  resName = paste0(cnvDir,"/subclone_CNV_plot.pdf")
  print(targ)
  pdf(resName, width = 14, height= 6)

  for (targ in blocks) {
    print(targ)
    if (targ == "Full") {
        vals = rowMeans(cnvMtx)
    } else {
      if (infercnv_obj@options$k_obs_groups  == 1) {
        cellIds <-names(infercnv_obj@tumor_subclusters$subclusters[[targ]][[1]])
        if (sum(cellIds %in% colnames(cnvMtx)) == 0) {
          print("Skipping...")
          next
        }
      } else {
        cellIds = names(hcBlocksAdj)[hcBlocksAdj == targ]
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

    #abline(v=markerLoci, col="red", lwd=0.2, lty=2)
    p <- p + axis(side=2,seq(round(cMin,2),round(cMax,2),0.01) )
    print(p)
  }

  dev.off()

}
