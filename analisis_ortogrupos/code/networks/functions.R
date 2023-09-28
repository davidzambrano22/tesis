# ========================= Creating dirs ===============================
create.dir <- function(dir){
  if (!dir.exists(dir)){
    dir.splitted <- strsplit(dir, "/")[[1]]
    dir <- ""
    for (subDir in dir.splitted){
      dir <- paste0(dir, paste0(subDir, "/"))
      if (dir.exists(dir)){next}
      else {dir.create(dir)}
      }
    }
}


# ======================== Filtering =====================================
# Obtain filtered dataframe -----------------------------------------------
get.filtered.matrix <- function(x, dataframe){
  if (is.vector(x) && length(x) > 0 && is.character(x)){
    x <- names(table(x))
    print(cat("Vector es character --> ", is.character(x)))
    filtered.frame <- dataframe[dataframe$`Row.Labels` %in% x,]
    return(filtered.frame)
  }else {
    print("Input is not a vector or is not a character vector")
  }
}

# Filtering by correlation -------------------------------------------------
get.highCorr.genes <- function(cor.matrix, threshold){
  high.corr.genes <- c()
  pos <- 1
  for (i in rownames(cor.matrix)){
    print(pos)
    print(i)
    corr.names <- rownames(cor.matrix[rownames(cor.matrix) != i,])
    correlated.hogs <- colnames(cor.matrix[,cor.matrix[i, corr.names] >= threshold])
    if (length(correlated.hogs) != 0){
      print(correlated.hogs)
      high.corr.genes <- c(high.corr.genes, c(i, correlated.hogs))
    }
    pos <- pos +1 
  }
  return(names(table(high.corr.genes)))
}

#==================================================================
