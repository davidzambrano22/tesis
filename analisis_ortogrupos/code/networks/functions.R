
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
  for (i in rownames(cor.matrix)){
    print(i)
    corr.names <- rownames(cor.matrix)[colnames(cor.matrix) != i]
    print(length(corr.names))
    corr.mean <- mean(cor.matrix[i, corr.names])
    print(sprintf("Corr mean for gene %1s is %2s==> ", i, corr.mean))
    if (corr.mean <= threshold){
      high.corr.genes <- c(high.corr.genes,i)
    }
  }
  return(high.corr.genes)
}

# ======================== Network making ===============================

# ---- For HARD THRESHOLDING ----

#Get.A.matrix() function gets the adjacency graph for a adjacency matrix
Get.A.matrix <- function(n.genes, t,corrMatrix){
  A <- matrix(ncol = n.genes, nrow = n.genes)
  #Fill Adjacency matrix
  for (g in n.genes){
    A[which(corrMatrix[,g] >= t),g] <- 1
    A[which(corrMatrix[,g] < t),g] <- 0
  }
  #Rename adjacency matrix names with HOGS names
  colnames(A) <- rownames(A) <- rownames(cor.matrix)
  
  #Convert diagonal to zero (no dirigida):
  #diag(A)<-0
  
  # Transform A matrix into a igraph object
  
  return(A)
}


# Functions get.network.coefficients() gets expected and observed
# clustering coefficients to get best threshold
get.network.coefficients <- function(corrMatrix){
  #starting threshold vector, degree and transitivity matrix for all genes  
  n.genes <- dim(corrMatrix)[1]
  C <- matrix(nrow = n.genes, ncol = 100)
  K <- matrix(nrow = n.genes, ncol = 100)
  threshold <- seq(0.01,0.99,0.01)
  
  Cr=Co=rep(0,100)
  #Compute Adjacency matrix for all genes for each threshold value
  for (t in threshold){
    A <- Get.A.matrix(n.genes = n.genes, t = t,corrMatrix = corrMatrix) 
    A <- graph.adjacency(A, mode = "undirected", diag = F)

    #Extract graph features to populate K and C matrix
    Cv=transitivity(A,type="local")
    Kv=degree(A,loops=FALSE)
    # K[,round(t*100,0)]<-Kv
    # C[,round(t*100,0)]<-Cv
    
    #------------------------Section to comment --------------------------------
    print(sprintf("-- %s ----------------------------", t))
    if (t == 0.5){
      gn <- which(Kv>=1) #Get positions of nodes in the degree vector which value is greater than 1
      print(cat("Number of nodes with degree >= 1", Kv[gn]))
      kn <- length(gn) #Get the number of nodes in the network
      print(cat("Kn for degree >= 1", kn))
      k1=1/kn*sum(Kv[gn])
      k2=1/kn*sum(Kv[gn]^2)
      Co[round(t*100,0)]=((k2-k1)^2)/(kn*k1^3) #Expected clustering coefficient for a random generated network
      if(kn==0){Co[round(t*100,0)]=0}
      
      gn<-which(Kv>1)
      print(cat("gn for degree > 1", Kv[gn]))
      kn=length(gn)
      print(cat("Kn for degree > 1", kn))
      Cr[round(t*100,0)]=1/kn*sum(Cv[gn]) #Observed clustering coefficient of network
      if(kn==0){Cr[round(t*100,0)]=0}
      next
    }
    gn <- which(Kv>=1) #Get positions of nodes in the degree vector which value is greater than 1
    #print(cat("gn for degree >= 1", Kv[gn]))
    kn <- length(gn) #Get the number of nodes in the network
    #print(cat("Kn for degree >= 1", kn))
    k1=1/kn*sum(Kv[gn])
    k2=1/kn*sum(Kv[gn]^2)
    Co[round(t*100,0)]=((k2-k1)^2)/(kn*k1^3) #Expected clustering coefficient for a random generated network
    if(kn==0){Co[round(t*100,0)]=0}
    
    gn<-which(Kv>1)
    #print(cat("gn for degree > 1", Kv[gn]))
    kn=length(gn)
    #print(cat("Kn for degree > 1", kn))
    Cr[round(t*100,0)]=1/kn*sum(Cv[gn]) #Observed clustering coefficient of network
    if(kn==0){Cr[round(t*100,0)]=0}
    #---------------------------------------------------------------------------
  }
 # ----------------------------------------------------------------------------
  # for (t in round(threshold * 100, 0)){
  #   print(sprintf("-- %s ----------------------------", t))
  #   gn <- which(K[,t]>=1)
  #   kn <- length(gn)
  #   if (t == 3) {
  #     View(K)
  #     print(cat("gn => ", gn))
  #     print(cat("Kn => ", kn))
  #     print(cat("sum(K[gn, t] = ",sum(K[gn, t])))
  #     k1=1/kn*sum(K[gn, t])
  #     print(cat("K1 =>", k1))
  #     k2=1/kn*sum(K[gn, t]^2)
  #     print(cat("K2 =>", k2))
  #     Co[t]=((k2-k1)^2)/(kn*k1^3) #Expected clustering coefficient for a random generated network
  #     if(kn==0){Co[t]=0}
  #     gn<-which(Kv>1)
  #     kn=length(gn)
  #     Cr[t]=1/kn*sum(C[gn,t]) #Observed clustering coefficient of network
  #     if(kn==0){Cr[t]=0}
  #   }
  #   print(cat("Kn for degree >= 1", kn))
  #   k1=1/kn*sum(K[gn, t])
  #   print(cat("K1 for degree >= 1", k1))
  #   k2=1/kn*sum(K[gn, t]^2)
  #   print(cat("K2 for degree >= 1", k2))
  #   Co[t]=((k2-k1)^2)/(kn*k1^3) #Expected clustering coefficient for a random generated network
  #   if(kn==0){Co[t]=0}
  #   gn<-which(Kv>1)
  #   print(cat("Kn for degree > 1", kn))
  #   kn=length(gn)
  #   Cr[t]=1/kn*sum(C[gn,t]) #Observed clustering coefficient of network
  #   if(kn==0){Cr[t]=0}
  # }
  coefficients. <- list(Cr,Co)
  return(coefficients.)
}

# --Function get.network.graph() gets the ipgraph object for the elements of the Adjacency matrix
get.network.graph <- function(corrMatrix, thresholds){
  n.genes <- dim(corrMatrix)[1]
  thresholds <- c(thresholds) # Definte the thresholds based on the previous plot
  for (t in thresholds) {
    #Define Adjacency Matrix:
    A <- Get.A.matrix(n.genes = n.genes, t = t, corrMatrix = corrMatrix)
    #Eliminate unconnected nodes
    #A <- A[Kv > 0, Kv > 0]
    
    #Get adjacency graph
    A <- graph.adjacency(A, mode = "undirected", diag = F)
    Kv <- degree(A,loops=FALSE)
    print(cat("Length of Degree vector ==> ", length(Kv)))
    
    
    
    print(cat("Clase del objeto matriz de adjacencia es => ",class(A)))
    print(cat("Dimensiones del objeto Matriz de adjacencia son ==> ", dim(A)))
    A=graph.adjacency(A,mode="undirected",add.colnames=NULL,diag=FALSE)
  }
  return(A)
}
