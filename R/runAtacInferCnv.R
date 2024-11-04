#' Run InferCNV calling using generated input config
#'
#' @param resDir Path to the result directory with input
#' @param configFile Name ofconfiguration file with InferCnv params
#' @return NULL
#' @export
#'
runAtacInferCnv <- function(resDir, configFile = "infercnv_config.yml") {

  setwd(resDir)

  print(paste("Loading InferCNV configuration from:",configFile))
  cfg <- config::get(file=configFile)

  print(paste("Processing",cfg$resName))
  print(paste("Input:",cfg$countsFile))
  print(paste("Annotation:",cfg$annFile))
  print(paste("Normal clusters:",cfg$refGroup))
  print(paste("Cut off:",cfg$cutOff))

  numClusters <- 1
  groupUsage <- TRUE
  if (!is.null(cfg$numClusters)) {
    numClusters <- cfg$numClusters
    groupUsage <- FALSE
  }
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
                                      chr_exclude=c("Y","MT") # default
                                      #chr_exclude=c("MT")
  )

  # this is required to allow plots, otherwise there is a fail in plotting
  assign("infercnv_obj", infercnv_obj, envir = .GlobalEnv)

  infercnv_obj = infercnv::run(infercnv_obj ,
                               # cutoff: 1 for SmartSeq, 0.1 for 10x Genomics, mean to meta
                               cutoff=cfg$cutOff,
                               out_dir=cfg$resName,
                               cluster_by_groups=groupUsage,
                               cluster_references = FALSE,
                               k_obs_groups =  numClusters,
                               smooth_method="runmeans",
                               # analysis_mode="subclusters", # verification
                               output_format = "pdf", # issue with atac
                               denoise=T,
                               plot_steps=F,
                               HMM=F,
                               no_plot = F
  )

}
