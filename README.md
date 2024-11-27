The R package focuses on the copy number variance profiling of the scATAC-seq data from 10X protocols ATAC and snMultiomics. 
It is based on the usage of [InferCNV tool](https://github.com/broadinstitute/infercnv) after initial custom preparation of input data. 
This approach was applied also in our [study](https://www.biorxiv.org/content/10.1101/2024.02.09.579690v1.full). 

For questions/comments: please create a github issue or contact the author directly.  


## Installation ##

Installation can be performed via devtools. For this purpose download the repo to a folder. 
Afterwards run command: 

```
library(devtools)
install("path/to/the/source/code")
```

## How to run the analysis ##

The example input data is located [here](https://drive.google.com/drive/folders/1okTZxc4yeuv1U2BsSEEccn1tng-sLqu7?usp=drive_link). 

To launch it, download the dataset and enter the corresponding folder path as mainDir variable. 
By default, the required input are 
1) full signals matrix 
2) annotation of cells to distinguish tumor and normal.   

The path to input (raw signals matrix and cells annotation) should be prepared via the following code: 

```
mainDir = "/path/to/test/data"
inDir = paste0(mainDir, "MB183_ATAC_data")
sId="MB183_ATAC_sub"
sAnn = paste0(mainDir, "MB183_ATAC.CNV_blocks_ann.subsample_filtered.txt" )
```

Afterwards, the name of the result should be stated and the input data formatted for CNV calling. 
Corresponding function processes input data, outputs it formatted to infercnv and also creates custom configuration.
It has a set of various options including for example meta-cells activation or usage of external reference control. 

```
resPath=paste0(mainDir,sId,"_result")
prepareAtacInferCnvInput(inDir,sAnn,resPath,targColumn = "cnvBlock")
```

Finally, the wrapper function to launch InferCnv is applied. It uses the generated configuration and outputs full result as additional subfolder:
```
# run InferCNV
runAtacInferCnv(resPath)

```






