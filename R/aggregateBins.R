#' Combine ATAC signals across bins of specific size
#'
#' @param sample Input Seurat object
#' @param resDir Path to the result directory
#' @param sId Result name
#' @param bin_size Size of the bin e.g. 100000 for 100Kbp
#' @param chrom_lengths Numeric vector of chromosome sizes, specific for genome
#'
#' @return NULL
#' @export
#'
aggregateBins <- function(sample, resDir, sId, bin_size, chrom_lengths) {

  #print("Get regions")
  regions <- sample@assays$ATAC@ranges
  # adjust for chromosome size and extract signals matrix
  valid_chromosomes <- intersect(seqlevels(regions), names(chrom_lengths))
  filtered_regions <- keepSeqlevels(regions, valid_chromosomes, pruning.mode = "coarse")
  seqlengths(filtered_regions) <- chrom_lengths[seqlevels(filtered_regions)]
  #print(filtered_regions)
  signal_matrix <- sample@assays$ATAC@counts[ gsub(":","-",as.character(filtered_regions)), ]

  # tile into bins
  bins <- tileGenome(seqlengths(filtered_regions), tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
  # overlap with bins
  region_bins <- findOverlaps(filtered_regions, bins)

  query_hits <- queryHits(region_bins)
  subject_hits <- subjectHits(region_bins)
  targ_bins <- unique(subject_hits)
  num_bins <- length(bins)
  # use cpp-compiled function
  binned_signal <- aggregate_bins_cpp_sparse(targ_bins, query_hits, subject_hits,
                                              signal_matrix, num_bins)

  # save result
  rownames(binned_signal) <- as.character(bins)
  colnames(binned_signal) <- colnames(signal_matrix)

  # re-write raw counts matrix
  gz2 <- gzfile(sprintf("%s/%s_raw_counts.txt.gz", resDir,sId), "w")
  write.table(binned_signal,gz2,sep="\t",quote=F)
  close(gz2)

  peakDf <- data.frame(bins)[,1:3]
  rownames(peakDf) <- paste0(  peakDf[,1],":",peakDf[,2],"-",peakDf[,3] )
  peakDf$seqnames <- gsub("chr","", peakDf$seqnames)
  summary(rownames(binned_signal) == rownames(peakDf))
  write.table(peakDf, sprintf("%s/%s_cnv_ref.binsize_%d.txt", resDir, sId,bin_size),
              col.names = F,sep="\t",quote=F)

  # assign and return novel matrix for possible further analysis: meta-cells
  sample@assays$ATAC@counts <- binned_signal
  sample

}
