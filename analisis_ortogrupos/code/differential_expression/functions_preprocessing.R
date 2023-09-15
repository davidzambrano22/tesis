

# ================================== get.accessions.by.quantile() ======================================
"Returns a list of vectors containing levels of anthocyanin for each accession and each quartile"

get.accessions.by.quantile <- function(anthocyanins.vector,
                                   accessions.vector, 
                                   accesions.dataset, 
                                   quantiles.vector){
  #In the following loop all anthocyanin values for each accession are extracted into vectors and saved in anthocyanin.levels List.
  anthocyanin.levels <- list()
  for (a in anthocyanins.vector){
    anthocyanin <- a
    level.values <-  accesions.dataset[accesions.dataset[,1] %in% accessions.vector, anthocyanin]
    names(level.values) <- accesions.dataset[accesions.dataset$Taxa %in% accessions.vector, "Taxa"]
    anthocyanin.levels[[anthocyanin]] <- level.values
  }
  
  #For every anthocyanin values vector obtain the value corresponding to each quantile in quantiles.vector function argument
  quantiles.values <- list()
  for (a in names(anthocyanin.levels)){
    anth.values <- anthocyanin.levels[[a]]
    quantiles. <- list()
    for (q in quantiles.vector){
      ql.value <- quantile(anth.values, probs=q/100)
      quantiles.[[q]] <- ql.value
    }
    quantiles.values[[a]] <- quantiles.
  }
  
  #Obtain accession names per quantile filter
  accesions.per.antho <- list()
  for (a in names(quantiles.values)){
    accesions.per.filter <- list()
    for (q in quantiles.vector){
     low.ammount <- c()
     high.ammount <- c()
     ql.value <- quantiles.values[[a]][[q]]
     for (n in accessions.vector){
       anth.value <- accesions.dataset[accesions.dataset$Taxa %in% n, a]
       if (anth.value <= ql.value) low.ammount <- c(low.ammount, n)
       else high.ammount <- c(high.ammount, n)
     }
     accesions.per.filter[[q]] <- list(low.ammount=low.ammount, high.ammount=high.ammount)
    }
    accesions.per.antho[[a]] <- accesions.per.filter
  }
  return(accesions.per.antho)
}

# ================================== get.contrast.vector() ======================================
get.contrast.vector <- function(accessions.contrast, below_quantile, above_quantile){
  accessions.contrast[accessions.contrast %in% below_quantile] <- 0
  accessions.contrast[accessions.contrast %in% above_quantile] <- 1
  return(accessions.contrast)
}
# ================================== get.clade.vector() ======================================
"Returns a vector containing information about clade for each accession"

get.clade.vector <- function(dataset, contrast, mode){
  if (!sum(colnames(dataset) %in% "Taxa") == 1 && (!sum(colnames(dataset) %in% "Clado_filogenia") == 1 || (!sum(colnames(dataset) %in% "PC1") == 1))){
    stop("There is no Taxa or CLado_filogenia column in accessions info dataset")
  }
  clade.vector <- c()
  if (mode == "clade"){
    for (acc in 1:length(contrast)){
      clade.vector[acc] <- dataset[as.vector(dataset$Taxa) %in% contrast[acc],]$`Clado_filogenia`
    } 
    print(cat("Clade vector is ==> ", clade.vector))
    print(cat("Length of clade vector is ==> ", length(clade.vector)))
    
  } else if (mode == "PCA"){
    for (acc in 1:length(contrast)){
      clade.vector[acc] <- dataset[as.vector(dataset$Taxa) %in% contrast[acc],]$`PC1`
    } 
    print(cat("Clade vector is ==> ", clade.vector))
    print(cat("Length of clade vector is ==> ", length(clade.vector)))
  }
  
  return(clade.vector)
}




#Control de calidad ------------------------------------------------------

# sub-functions -----------------------------------------------------------
check.first.col <- function(dataset){
  if (class(dataset[[1]]) != "character"){
    warning("The first column of the dataset doesn't have HOGs names")
  }
  return(TRUE)
}

