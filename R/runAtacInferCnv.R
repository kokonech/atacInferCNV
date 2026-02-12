#' Wrapper function to run InferCNV calling
#'
#' This function calls InferCNV from generated input. It has support for all
#' original inferCnv options.
#' @param resDir Path to the result directory with input
#' @param configFile Name of configuration file with InferCnv input data
#' @param numClusters Number of clusters for hier. clustering. If equals one (by default)
#' then no clustering is performed and provided annoitation used for the formation of CNV groups.
#' @param chrToExclude Chromosomes to exclude. Default: Y,MT
#' @param addDenoise Activate denoise (InferCNV param). Deafult: TRUE
#' @param clusterRefs Cluster also reference (InferCNV param). Default: FALSE
#' @param smoothMethod Method for smoothing (InferCNV param). Default: runmeans
#' @param ... Other parameters to provide for infercnv::run, more details in documentation of this function
#' @return Invisibly returns NULL.
#' @examples
#' resPath = tempfile()
#' inPath = system.file("extdata", "MB183_ATAC_subset.tsv.gz", package = "atacInferCnv")
#' sAnn = system.file("extdata", "MB183_ATAC_subset.CNV_blocks_ann.txt", package = "atacInferCnv" )
#' prepareAtacInferCnvInput(inPath,sAnn,resPath, targColumn = "cnvBlock", ctrlGrp = "Normal")
#' runAtacInferCnv(resPath)
#'
#' @export
#'
runAtacInferCnv <- function(resDir, configFile = "infercnv_config.yml",
                            numClusters = 1, chrToExclude = c("Y","MT"),
                            addDenoise = TRUE, clusterRefs = FALSE,
                            smoothMethod = "runmeans",
                            ...) {

  if (!(dir.exists(resDir))) {
    stop("The result directory with InferCNV input does not exist:",resDir)
  }

  # Ensure the working directory is restored when the function exits
  originalDir <- getwd()
  on.exit(setwd(originalDir))
  setwd(resDir)

  message("Loading InferCNV configuration from: ",configFile)
  cfg <- config::get(file=configFile)

  message("Processing ",cfg$resName)
  message("Input: ",cfg$countsFile)
  message("Annotation: ",cfg$annFile)
  message("Normal clusters: ",cfg$refGroup)
  message("Cut off: ",cfg$cutOff)

  groupUsage <- ifelse(numClusters > 1,FALSE,TRUE)
  message("Num tumor clusters: ",numClusters)

  refGroups <- str_split(cfg$refGroups,",")[[1]]

  message("Assign custom reference: ",cfg$customRef)
  geneOrderRef <- cfg$customRef


  # specific:  selected cluster as reference
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=cfg$countsFile,
                                      annotations_file=cfg$annFile,
                                      delim="\t",
                                      gene_order_file=geneOrderRef,
                                      ref_group_names=refGroups,
                                      chr_exclude=chrToExclude
  )

  # this is required to allow plots, otherwise there is a fail in plotting
  assign("infercnv_obj", infercnv_obj, envir = .GlobalEnv)

  infercnv_obj <- infercnv::run(infercnv_obj ,
                               # cutoff: 1 for SmartSeq, 0.1 for 10x Genomics, mean to meta
                               cutoff=cfg$cutOff,
                               out_dir=cfg$resName,
                               cluster_by_groups=groupUsage,
                               k_obs_groups =  numClusters,
                               output_format = "pdf", # issue with atac
                               # further already custom params to play with
                               denoise=addDenoise,
                               cluster_references = clusterRefs,
                               smooth_method=smoothMethod,
                               ...

  )

  invisible(NULL)
}
