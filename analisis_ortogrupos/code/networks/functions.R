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
    cat("Vector es character --> ", is.character(x), "\n")
    cat("No. de.novo dataset HOGs -->", length(dataframe$`Row.Labels` %in% x), "\n")
    cat("No. Network HOGs in de.novo dataset-->", sum(dataframe$`Row.Labels` %in% x), "\n")
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

view.hogs.degree <- function(degree.frame, n = F, modules.t){
  modules.f <- data.frame(modules.t) %>% mutate(ID = rownames(.))
  #colnames(modules.frame) <- c("", "")
  # View(modules.frame)
  
  colnames(modules.frame) <- c("ID", "Module")
  # (This function return a list with info about nodes degree) #
  if (is.numeric(n)){
    print("n is numeric")
    highest.degrees <- degree.frame %>% arrange(desc(Degree))
    View(highest.degrees[1:n,])
    return(list(Names = unlist(highest.degrees[1:n,]$`X`), Degrees = unlist(highest.degrees[1:n,]$`Degree`)))
    
  } else {
    highest.degrees <- degree.frame %>% arrange(desc(Degree)) %>%
      left_join(annot.dataframe, by = c("X" = "HOG_ID"), keep = T) %>% dplyr::select(X, Degree, Keywords, Annot_pathway) %>%
      left_join(modules.f, by = c("X" = "ID")) %>% dplyr::select(X, Degree, modules.t, Keywords, Annot_pathway)
    print("n is False")
    View(highest.degrees)
    print(getwd())
    if (partial.corr) {
      write.csv(highest.degrees, sprintf("/home/andres/Documents/Tesis/codigos_tesis/analisis_ortogrupos/results/networks/big/raw/partial.corr/%s/most_connected_hogs_annot.csv",TAO))
    }
    else {
      write.csv(highest.degrees, sprintf("/home/andres/Documents/Tesis/codigos_tesis/analisis_ortogrupos/results/networks/big/raw/pearson_corr/%s/most_connected_hogs_annot.csv",TAO))
    }
    
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
    cat(h, "= Node ", modules.frame[modules.frame[,1] == h, 2], "\n")
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
  diff.HOGS.in.vector <- c()
  
  total.flav.HOGS.in.names <- c()
  total.DiffExp.HOGS.in.names <- c()
  
  for (m in as.numeric(names(modules))){
    modules.vector <- c(modules.vector, m)
    
    module.HOGs <- modules.frame[modules.frame[,2] == m, 1]
    n.module.HOGs <- length(module.HOGs)
    modules.HOGs.vector <- c(modules.HOGs.vector, n.module.HOGs)
    
    # Flavonoid HOGS
    flav.HOGS.in <- module.HOGs[module.HOGs %in% floavonoid.path.HOGS]
    total.flav.HOGS.in.names <- c(total.flav.HOGS.in.names, flav.HOGS.in)
    n.flav.HOGS.in <- length(flav.HOGS.in)
    flav.HOGS.in.vector <- c(flav.HOGS.in.vector, n.flav.HOGS.in)
    
    
    # Differentially expressed HOGs
    diff.HOGS.in <-  module.HOGs[module.HOGs %in% D.E.HOGS]
    total.DiffExp.HOGS.in.names <- c(total.DiffExp.HOGS.in.names, diff.HOGS.in)
    n.diff.HOGS.in <- length(diff.HOGS.in)
    diff.HOGS.in.vector <- c(diff.HOGS.in.vector, n.diff.HOGS.in)
    
  }
  cat("Number of total most variable HOGs = ", length(most.var.HOGS.names), "\n")
  cat("Total number of Flavonoid path HOGs in Network = ", sum(most.var.HOGS.names %in% modules.frame$X), "\n")
  cat("Number of removed Flavonoid path HOGs = ", length(most.var.HOGS.names) - sum(most.var.HOGS.names %in% modules.frame$X), "\n\n")
  
  
  cat("Number of total HOGs from flavonoids path = ", length(floavonoid.path.HOGS), "\n")
  cat("Total number of Flavonoid path HOGs in Network = ", sum(floavonoid.path.HOGS %in% modules.frame$X), "\n")
  cat("Number of removed Flavonoid path HOGs = ", length(floavonoid.path.HOGS) - sum(floavonoid.path.HOGS %in% modules.frame$X), "\n\n")
  
  
  cat("Number of total HOGs from flavonoids path = ", length(D.E.HOGS), "\n")
  cat("Total number of FDIfferentially expressed  HOGs in Network = ", sum(D.E.HOGS %in% modules.frame$X), "\n")
  cat("Number of removed DIfferentially expressed  HOGs  = ", length(D.E.HOGS) - sum(D.E.HOGS %in% modules.frame$X))
  
  
  # cat("Total number of Flavonoid path HOGs in Network is ", length(total.flav.HOGS.in.names), "\n")
  # cat("Number of removed Flavonoid path HOGs is = ", length(floavonoid.path.HOGS) - length(total.flav.HOGS.in.names), "\n")
  # 
  # cat("Number of flavonoid HOGs in DE. HOGs =", sum(total.DiffExp.HOGS.in.names %in% total.flav.HOGS.in.names ), "\n")
  # 
  # cat("Total number of FDIfferentially expressed  HOGs in Network is ", length(total.DiffExp.HOGS.in.names), "\n")
  # cat("Number of removed DIfferentially expressed  HOGs is = ", length(D.E.HOGS) - length(total.DiffExp.HOGS.in.names))
  
  bar.data <- matrix(c(modules.HOGs.vector, flav.HOGS.in.vector, diff.HOGS.in.vector), nrow = 3, byrow = TRUE)
  barplot(bar.data, beside = TRUE, col = c("blue", "red", "green"), names.arg = modules.vector,
          legend.text = c("Module", "Flavonoid Hogs", "Differentially Expressed"), 
          main = "Grouped Barplot of Modules and Flavonoid HOGs")
  
  bar.data <- matrix(c(flav.HOGS.in.vector, diff.HOGS.in.vector), nrow = 2, byrow = TRUE)
  barplot(bar.data, beside = TRUE, col = c("red", "green"), names.arg = modules.vector,
          legend.text = c("Flavonoid Hogs", "Differentially Expressed"), 
          main = "Grouped Barplot of Flavonoid and DIfferentially expressed HOGs ")
}
