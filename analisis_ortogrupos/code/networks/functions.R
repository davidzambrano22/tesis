library(dplyr)

# ====================================== CREATE DIRS ===============================
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


# ====================================== FILTERING =====================================
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
    corr.names <- rownames(cor.matrix[rownames(cor.matrix) != i,])
    correlated.hogs <- colnames(cor.matrix[,cor.matrix[i, corr.names] >= threshold])
    if (length(correlated.hogs) != 0){
      high.corr.genes <- c(high.corr.genes, c(i, correlated.hogs))
    }
    pos <- pos +1
  }
  return(names(table(high.corr.genes)))
}


getHighcorGenes <- function(expression.matrix, threshold){
  high.corr.genes <- c()
  
  nrows <- 1:dim(expression.matrix)[1]
  for (i in nrows){
    for (j in nrows){
      print(cat(i, j))
      if (i == j) next
      corr <- abs(cor(unlist(expression.matrix[i,]), unlist(expression.matrix[j,])))
      print(cat("Corr is = ", corr))
      if (is.na(corr)) next
      else if (corr >= threshold) {high.corr.genes <- c(high.corr.genes, rownames(expression.matrix[i,]), rownames(expression.matrix[i+1,]))} 
    }
  }
  return(high.corr.genes)
}


# ====================================== NETWORK ANALYSIS =================================

# ====================== DEGREE ======================

view.hogs.degree <- function(degree.frame, n = F){
  # (This function return a list with info about nodes degree) #
  if (is.numeric(n)){
    highest.degrees <- degree.frame %>% arrange(desc(Degree))
    View(highest.degrees[1:n,])
    return(list(Names = unlist(highest.degrees[1:n,]$`X`), Degrees = unlist(highest.degrees[1:n,]$`Degree`)))
  } else {
    highest.degrees <- degree.frame %>% arrange(desc(Degree))
    View(highest.degrees)
  }
}

get.central.HOGS.forNode <- function(modules.file, degree.list){
  modules.frame <<- read.csv(modules.file)
  message("Modules...")
  modules <<- table(modules.frame[,2])
  print(modules)
  n.modules <- length(modules)
  message(sprintf("Network has %s modules; most connected HOGS are in modules: ", n.modules))
  
  for (h in degree.list$Names){
    print(cat(h, "= Node ", modules.frame[modules.frame[,1] == h, 2]))
  }
}

connected.Hogs.forNode <- function(degree.frame, module, n){
  hogs <- modules.frame[modules.frame[,2] == module, 1]
  degree.frame.filtered <- degree.frame[degree.frame$X %in% hogs,] %>% arrange(desc(Degree))
  degree.frame.filtered <-  degree.frame.filtered[1:n,]
  print(" HOG       DEGREE")
  for (x in unlist(degree.frame.filtered$`X`)){
    print(cat(x, degree.frame.filtered[degree.frame.filtered$`X` == x, "Degree"]))
  }
  return(degree.frame.filtered)
}


# ====================== MODULES ANALYSIS ===========
draw.flavonoidInModules <- function(){
  modules.vector <- c()
  modules.HOGs.vector <- c()
  flav.HOGS.in.vector <- c()
  
  for (m in as.numeric(names(modules))){
    modules.vector <- c(modules.vector, m)
    
    module.HOGs <- modules.frame[modules.frame[,2] == m, 1]
    n.module.HOGs <- length(module.HOGs)
    modules.HOGs.vector <- c(modules.HOGs.vector, n.module.HOGs)
    
    flav.HOGS.in <- module.HOGs[module.HOGs %in% floavonoid.path.HOGS]
    n.flav.HOGS.in <- length(flav.HOGS.in)
    flav.HOGS.in.vector <- c(flav.HOGS.in.vector, n.flav.HOGS.in)
  }
  print(length(modules.vector))
  print(length(modules.HOGs.vector))
  print(length(flav.HOGS.in.vector))
  
  bar.data <- matrix(c(modules.HOGs.vector, flav.HOGS.in.vector), nrow = 2, byrow = TRUE)
  barplot(bar.data, beside = TRUE, col = c("blue", "red"), names.arg = modules.vector,
          legend.text = c("Module", "Flavonoid Hogs"), 
          main = "Grouped Barplot of Modules and Flavonoid HOGs")
}
