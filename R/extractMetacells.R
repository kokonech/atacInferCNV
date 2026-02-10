#' Compute meta-cells out of the main
#'
#' @param resDir Path to the result directory
#' @param sId Result name
#' @param sample Input Seurat object to split
#' @param targColumn Name of annotation column to split
#' @param metacell_content Amount of cells for adjustment, default n=5
#' @return Invisibly returns NULL.
#'
extractMetacells <- function(resDir, sId, sample, targColumn, metacell_content = 5) {

  # Checking list:
  #Seurat::Idents(sample) <- sample$cluster_names
  message("Number of cells for meta-cell formation:",metacell_content)
  #print(targColumn)
  targAnn <- factor(sample[[targColumn,drop=TRUE]])
  #print(targAnn)
  rawCounts <- sample@assays$ATAC@counts
  #print(dim(rawCounts))

  # Will store all the metacells. The test column will be removed at the end.
  whole_metacells <- data.frame(test = rownames(rawCounts), row.names = rownames(rawCounts))


  # Will store the complete annotation for the metacells.
  whole_annotation <- data.frame(cluster_names = "test", row.names = "test")
  meta_counter <- 0 # To keep a count of the metacells that are created.

  # Generate a new metadata column storing the mapping cell-metacell.
  sample[["metacell_mapping"]] <- "not_mapped"


  # for restored object
  for (cluster_id in levels(targAnn)){
    # custom version
    #for (cluster_id in levels(sample$)){
    message("Computing metacells for cluster", cluster_id)
    # Will store the metacells per cluster.
    metacells <- data.frame(test = rownames(rawCounts), row.names = rownames(rawCounts))
    #chunksample <- sample[, sample$cluster_names == cluster_id]
    # Subset the sample by each cluster ID.
    #chunksample <- sample[, targAnn  == cluster_id]
    # Get the count data as a data frame and transpose it so columns are GENES and rows are CELLS.
    #countdata <- t(as.data.frame(Seurat::GetAssayData(chunksample, slot = "counts")))
    countdata <- t(as.data.frame(rawCounts[ ,targAnn  == cluster_id]))
    #print(head(countdata))
    # Get the possible amount of metacells.
    times <- trunc(dim(countdata)[1] / metacell_content)
    for (i in seq(1,times)){
      meta_counter <- meta_counter + 1
      # Generate slice points for each metacell. i.e: 1-5, 6-10, 11-15...
      start <- ((i -1) * metacell_content + 1)
      end <- i * metacell_content
      # Compute the slice as a data frame containing the sum of the subsetted cells. dims = 1 row (metacell), X columns (genes)
      slice <- as.data.frame(colSums(countdata[start:end, ]))
      # Get the name of the cells merged.
      cell_names <- rownames(countdata[start:end, ])
      # Add the metacell.
      col_name <- sprintf("metacell_%s", meta_counter)
      metacells[[col_name]] <- slice[,1]
      # Add the metacell mapping to the cells in the sample.
      sample$metacell_mapping[colnames(sample) %in% cell_names] <- col_name
    }
    # Delete the test column as we already have more than 1 column in our data frame.
    metacells[["test"]] <- NULL
    # Will contain the annotation of the generated metacells. Columns: cluster identities. Rows: each metacell.
    annotation <- data.frame(cluster_names = colnames(metacells), row.names = colnames(metacells))
    # Replace the dummy cluster_names column's values for the actual label for the cluster.
    annotation$cluster_names <- cluster_id
    # Add the annotation data and the metacell data to the "whole" dataframe. \
    # In the end: Number of Columns for metacell object = Number of rows for annotation object.
    whole_metacells <- cbind(whole_metacells, metacells)
    whole_annotation <- rbind(whole_annotation, annotation)
  }

  # Delete the test row from the global annotation data.
  whole_annotation <- whole_annotation[!rownames(whole_annotation) %in% c("test"), , drop = FALSE]

  # Delete the test column from the global metacell data.
  whole_metacells$test <- NULL

  # Path and name of the annotation file that will be used in the inferCNV call.
  #cnv_analysis_folder <- "" # Path to the folder that will store the annotation file.
  #dir.create(cnv_analysis_folder, recursive = TRUE)

  annotation_file <- sprintf("%s/%s_annotation_metacells.tsv", resDir, sId)

  # Save the annotation object.
  utils::write.table(whole_annotation,
                     file = annotation_file,
                     sep = "\t",
                     row.names = TRUE,
                     col.names = FALSE,
                     quote = FALSE)

  # Return the metacell object as a matrix (required for running inferCNV).
  whole_metacells <- as.matrix(whole_metacells)

  # It would be wise to save the Seurat object with the metacell mapping.

  ann_file2 <- sprintf("%s/%s_metacells_mapping.tsv", resDir , sId)
  write.table(sample@meta.data, ann_file2, sep="\t", quote=FALSE)

  mtxFile <- paste0(resDir,"/", sId, "_raw_matrix_metacells.txt.gz")
  gz1 <- gzfile(mtxFile, "w")
  write.table(whole_metacells, gz1, quote=FALSE,sep="\t")
  close(gz1)

  invisible(NULL)

}
