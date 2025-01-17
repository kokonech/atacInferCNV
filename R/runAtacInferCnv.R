#' Wrapper function to run InferCNV calling using generated input data
#'
#' @param resDir Path to the result directory with input
#' @param configFile Name of configuration file with InferCnv input data
#' @param numClusters Number of clusters for hier. clustering. If equals one then no clustering performed.
#' @param chrToExclude Chromosomes to exclude. Default: Y,MT
#' @param addDenoise Activate denoise (InferCNV param). Deafult: TRUE
#' @param clusterRefs Cluster also referecne (InferCNV param). Default: FALSE
#' @param smoothMethod Method for smoothing (InferCNV param). Default: runmeans
#' @param ... Other parameters to provide for infercnv::run, more details in documentation of this function
#' @return NULL
#' @export
#'
runAtacInferCnv <- function(resDir, configFile = "infercnv_config.yml",
                            numClusters = 3, chrToExclude = c("Y","MT"),
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

  print(paste("Loading InferCNV configuration from:",configFile))
  cfg <- config::get(file=configFile)

  print(paste("Processing",cfg$resName))
  print(paste("Input:",cfg$countsFile))
  print(paste("Annotation:",cfg$annFile))
  print(paste("Normal clusters:",cfg$refGroup))
  print(paste("Cut off:",cfg$cutOff))

  groupUsage <- ifelse(numClusters > 1,FALSE,TRUE)
  print(paste("Num tumor clusters:",numClusters))

  refGroups <- str_split(cfg$refGroups,",")[[1]]

  print(paste("Assign custom reference:",cfg$customRef))
  geneOrderRef = cfg$customRef


  # specific:  selected cluster as reference
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=cfg$countsFile,
                                      annotations_file=cfg$annFile,
                                      delim="\t",
                                      gene_order_file=geneOrderRef,
                                      ref_group_names=refGroups,
                                      chr_exclude=chrToExclude
  )

  # this is required to allow plots, otherwise there is a fail in plotting
  assign("infercnv_obj", infercnv_obj, envir = .GlobalEnv)

  infercnv_obj = infercnv::run(infercnv_obj ,
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

}
