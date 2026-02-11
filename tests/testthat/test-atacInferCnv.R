test_that("prepare input function parses data correctly", {

  ctrlId = "Normal"
  sId = "MB183_ATAC_test"
  inPath = system.file("extdata", "MB183_ATAC_subset.tsv.gz", package = "atacInferCnv")
  sAnn = system.file("extdata", "MB183_ATAC_subset.CNV_blocks_ann.txt", package = "atacInferCnv" )
  resPath = "incorrect_path"
  # no result path
  expect_error(
    prepareAtacInferCnvInput("bad_input", "meta"),
    regexp = "Path to result is not provided!"
  )

  # input mtx incorrect
  expect_error(
    prepareAtacInferCnvInput("bad_input", "meta", resPath),
    regexp = "Input data folder is not found"
  )

  # ann mtx incorrect
  expect_error(
    prepareAtacInferCnvInput(inPath, "meta", resPath),
    regexp = "Annotation file is not found"
  )

  # ann column
  expect_error(
    prepareAtacInferCnvInput(inPath, sAnn, resPath),
    regexp = "Required annotation column is not available"
  )

  # normal cell types
  expect_error(
    prepareAtacInferCnvInput(inPath, sAnn, resPath,
                             targColumn = "cnvBlock", ctrlGrp = "boo"),
    regexp = "Non-tumor control group is not found in annotation"
  )


})

test_that("atacInferCnv works on toy data", {

  ctrlId = "Normal"
  sId = "MB183_ATAC_test"
  inPath = system.file("extdata", "MB183_ATAC_subset.tsv.gz", package = "atacInferCnv")
  sAnn = system.file("extdata", "MB183_ATAC_subset.CNV_blocks_ann.txt", package = "atacInferCnv" )
  resPath = tempfile()

  expect_message(prepareAtacInferCnvInput(inPath, sAnn, resPath,
                                          targColumn = "cnvBlock",ctrlGrp = ctrlId),
                 regexp = "Prepared input")

  expect_true(file.exists(paste0(resPath,"/sample_raw_counts.txt.gz")))
  expect_true(file.exists(paste0(resPath,"/sample_cnv_ref.txt")))
  expect_true(file.exists(paste0(resPath,"/sample_cnv_ann.txt")))

  expect_message( runAtacInferCnv(resPath),
                  regexp = "Making the final infercnv heatmap" )
  expect_true(file.exists(paste0(resPath,"/sample_infercnv/run.final.infercnv_obj")))


})

