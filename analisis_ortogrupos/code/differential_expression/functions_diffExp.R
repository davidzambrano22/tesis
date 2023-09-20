require(DESeq2)

# differential expression -------------------------------------------------
get_DEGs <- function(dataset, cond.1.vector, cond.2.vector, treatments){ #p.value, covariable=FALSE){
  #Check that the first column in dataset has the HOGs names
  tryCatch(
    {check.first.col(dataset=dataset)},
    warning=function(w){
      print("The get_DEGs function will return row numbers instead of HOGS names")
    },
    finally = {
      HOGs.names <- dataset[[1]]
      expression.set <- dataset[2:dim(dataset)[2]]
      #Create DESeq object
      coldata <- data.frame(antocianina=treatments)
      dds0 = DESeq2::DESeqDataSetFromMatrix(
        countData=dataset[c(cond.1.vector, cond.2.vector)],
        colData = coldata,
        design = antocianina
      )
      dds0=DESeq2::estimateSizeFactors(dds0)
    }
  )
}
