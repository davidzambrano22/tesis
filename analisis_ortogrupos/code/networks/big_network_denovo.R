
library(dplyr)
library(corrplot)
library(igraph)
library(infotheo)
source("./functions.R")

SHRINKAGE <- F
log_transform <- T
mutual_info <- F

for (dir. in list.dirs("../../data/")){
  if (grepl("denovo$", dir.) || grepl("REFs$", dir.)){
    ### Entering potato_denovo expression files (Universidad Nacional)
    if (grepl("denovo$", dir.)){
      print(dir.)
      for (file. in list.files(dir.)){
        if (grepl("^A_TPM", file.)){
          print(file.)
          de.novo.table <- read.delim(paste0(dir.,"/", file.), sep="\t", header=TRUE, na.strings="") %>% select("Row.Labels", starts_with("A")) 
          
          # ====================== Loading gene names of interest  ================================================
          # --Load genes from Falvonoids path
          message("Loading Anthocyanin annotated  HOGS...")
          floavonoid.path.HOGS <- read.csv("../../results/filtered_genes/flav.HOGS.csv")$`annot.info` %>% as.vector()
          
          # --Load most variable HOGS
          message("Loading most variable HOGS...")
          most.var.HOGS.names <- read.csv("../../results/filtered_genes/6000_most_variable_genes.csv6000.csv")$`Var1` %>% as.vector()
          
          # --Load DEHS (differentially expressed HOGS)
          if (!SHRINKAGE){
            message("Loading differentially expressed HOGS...")
            D.E.HOGS <- c()
            DEHs.Cyanidin <-  read.csv("../../results/DEGS/DEGS_denovo/model/Cyanidin.csv")$`total.DEH` %>% as.vector() 
            DEHs.Delphinidin <-  read.csv("../../results/DEGS/DEGS_denovo/model/Delphinidin.csv")$`total.DEH` %>% as.vector()
            DEHs.Peonidin <-  read.csv("../../results/DEGS/DEGS_denovo/model/Peonidin.csv")$`total.DEH` %>% as.vector()
            
            D.E.HOGS <- c(D.E.HOGS, DEHs.Cyanidin, DEHs.Delphinidin, DEHs.Peonidin) 
            D.E.HOGS <- names(table(D.E.HOGS)) #Define unique HOGS
            message(cat("Number of unique differentially expressed HOGS is ==>", length(D.E.HOGS)))
            
          } else if (SHRINKAGE){
            message("Loading differentially expressed HOGS from SHRINKAGE workflow...")
            D.E.HOGS <- c()
            DEHs.Cyanidin <-  read.csv("../../results/DEGS/DEGS_denovo/shrinkage/Cyanidin.csv")$`total.DEH` %>% as.vector() 
            DEHs.Delphinidin <-  read.csv("../../results/DEGS/DEGS_denovo/shrinkage/Delphinidin.csv")$`total.DEH` %>% as.vector()
            DEHs.Peonidin <-  read.csv("../../results/DEGS/DEGS_denovo/shrinkage/Peonidin.csv")$`total.DEH` %>% as.vector()
            
            D.E.HOGS <- c(D.E.HOGS, DEHs.Cyanidin, DEHs.Delphinidin, DEHs.Peonidin) 
            D.E.HOGS <- names(table(D.E.HOGS)) #Define unique HOGS
            message(cat("Number of unique differentially expressed HOGS is ==>", length(D.E.HOGS)))
          }
          
          # ====================== Obtaining expression matrix for all genes of interest ==========================
          if (log_transform){
            message("Obtaining log 2 expression matrix for all genes of interest...")
            network.genes <- c(most.var.HOGS.names, D.E.HOGS, floavonoid.path.HOGS) # --> Enter here as a vector all genes loaded from previous step
            expression.data <- get.filtered.matrix(names(table(c(most.var.HOGS.names, D.E.HOGS, floavonoid.path.HOGS))), de.novo.table)
            #stop("Determining expression data matrix properties")
            #View(expression.data)
            rownames(expression.data) <- expression.data$Row.Labels #Rename expression data with gene names
            expression.data <- expression.data %>% select(-Row.Labels) %>% log(base = 2)
            View(expression.data)
            #stop("Determining log matrix properties")
            print(cat("Filtered dataframe dimenstions are --> ", dim(expression.data)))
          } 
          
          else {
            message("Obtaining expression matrix for all genes of interest...")
            network.genes <- c(most.var.HOGS.names, D.E.HOGS, floavonoid.path.HOGS) # --> Enter here as a vector all genes loaded from previous step
            expression.data <- get.filtered.matrix(names(table(c(most.var.HOGS.names, D.E.HOGS, floavonoid.path.HOGS))), de.novo.table)
            rownames(expression.data) <- expression.data$Row.Labels #Rename expression data with gene names
            expression.data <- expression.data %>% select(-Row.Labels)
            print(cat("Filtered dataframe dimenstions are --> ", dim(expression.data)))
          }
          
          # ====================== Obtaining similarity matrix for all genes of interest ==========================
          if (mutual_info){
            # --Obtaining mutual information
            message("Obtaining mututal information for all genes of interest...")
            discretized.matrix <-discretize(t(expression.data), disc = "equalfreq", nbins = ceiling(1+log(ncol(expression.data))))
            sim.matrix <- mutinformation(discretized.matrix, method = "shrink")
            
            
            S1 <- sqrt(1/exp(-2*sim.matrix))
            S1[which(is.na(S1))] <- 0
            hist(S1, xlim=c(0,1))
            stop()
            
          } else {
            # --Obtaining correlation matrix for all genes
            message("Obtaining correlation matrix for all genes of interest...")
            sim.matrix <- abs(cor(t(expression.data), use = "pairwise.complete.obs"))
            #Rename correlation matrix with HOGS names
            colnames(cor.matrix) <- rownames(expression.data)
            rownames(cor.matrix) <- rownames(expression.data)
            stop("Checking corr matrix")
          }
          
          # ====================== Check point for correlated genes ==========================
          # --Filtering out low correlated genes
          # High.corr.genes <- get.highCorr.genes(cor.matrix = cor.matrix, threshold = 0.5)
          # stop("Stopping code for mantenaince")
          
          # ===================== Determining the best threshold to work with on getting the network =============
          message("Getting clustering coefficients...")
          # net.coefficients <- get.network.coefficients(cor.matrix) #Get clustering coefficients
          # stop("Stop for mantainence")
          
          n <- nrow(cor.matrix)
          C <- matrix(nrow = n, ncol = 100)
          K <- matrix(nrow=n, ncol=100)
          ltaos <- seq(0.01,0.99,by=0.01)
          
          for(tao in ltaos){
            print(tao)
            ##Matriz de adyacencia:
            A=matrix(0,nrow=n,ncol=n)    
            ##Completa la matriz de adyacencia usando la funci?n de adyacencia:
            for(i in 1:n){  
              A[which(cor.matrix[,i]>=tao),i]<-1
              A[which(cor.matrix[,i]<tao),i]<-0
            }
            ##Transforma la matriz A en un objeto igraph:
            A=graph.adjacency(A,mode="undirected",diag=FALSE)
            ##Calcula el Ci de los nodos:
            Cv=transitivity(A,type="local")
            ##Calcula el ki de los nodos:
            Kv=degree(A,loops=FALSE)
            ##Guarda Ci y ki en los vectores C y K respectivamente:
            K[,round(tao*100,0)]<-Kv
            C[,round(tao*100,0)]<-Cv
          }
          
          Cr=Co=rep(0,100)
          ##Para cada valor de tao:
          for(i in round(ltaos*100,0))
          {  
            gn<-which(K[,i]>=1)#Posici?n de los genes conectados en la red
            kn=length(gn)#N?mero de nodos en la red
            k1=1/kn*sum(K[gn,i])#Variable en ecuaci?n 3 (Ver Elo et.al. 2007)
            k2=1/kn*sum(K[gn,i]^2)#Variable en ecuaci?n 3 (Ver Elo et.al. 2007)
            Co[i]=((k2-k1)^2)/(kn*k1^3) #Coeficiente de agrupamiento esperado para una red aleatoria
            if(kn==0){Co[i]=0}#Si no hay nodos conectados: Co=0
            gn<-which(K[,i]>1)#Posici?n de los genes con k1>1.
            kn=length(gn)#N?mero de genes con m?s de una arista en la red
            Cr[i]=1/kn*sum(C[gn,i])#Coeficiente de agrupamiento observado en la red.
            if(kn==0){Cr[i]=0}#Si no hay nodos conectados: Cr=0
          }
          
          hist(K[,99],xlab = "Degree", main = "Threshold 0.99") ##### Plot de grado
          plot(ltaos,abs(Cr-Co)[ltaos*100],t="l",xlab="Threshold",ylab="|C-Co|")
          
          dif=runmed(abs(Cr-Co),k=3,endrule="constant")[1:100]
          plot(ltaos,dif[ltaos*100],t="l",xlab="Threshold",ylab="|C-Co|")
          identify(ltaos,dif[ltaos*100],n=1) / 100
          
          ## ========================================CREATE NETWORK
          TAO <- 0.92
          if (!dir.exists(sprintf("../../results/networks/big/%s/", TAO))){
            dir.create(sprintf("../../results/networks/big/%s/", TAO))
          }
          if (mutual_info){
            if (!dir.exists(sprintf("../../results/networks/big/%s/mutual_info", TAO))){
              dir.create(sprintf("../../results/networks/big/%s/mutual_info", TAO))
            }
            if (!mutual_info){
              if (!dir.exists(sprintf("../../results/networks/big/%s/corr", TAO))){
                dir.create(sprintf("../../results/networks/big/%s/corr", TAO))
              }
              
              red_info <- list()
              for (tao in c(TAO)){
                red_n <- 1
                #Define matriz de adyacencia:
                A=matrix(0,nrow=n,ncol=n)   
                
                #Completa la matriz de adyacencia usando la función de adyacencia:
                for(i in 1:n)
                {  
                  A[which(cor.matrix[,i]>=tao),i]<-1
                  A[which(cor.matrix[,i]<tao),i]<-0
                }
                
                #Trabajando con los nombres de los genes
                colnames(A) <- rownames(A) <- rownames(cor.matrix)
                
                #Convierte la diagonal en ceros (red no dirigida):
                diag(A)<-0
                
                #Elimina nodos no conectados:
                A=A[which(K[,round(tao*100,0)]>0),which(K[,round(tao*100,0)]>0)]
                
                class(A)
                dim(A)
                stop()
                
                A=graph.adjacency(A,mode="undirected",add.colnames=NULL,diag=FALSE)
                #plot(A)
                
                vertex_colors <- rep("lightblue", 1181)
                vertex_colors[1:924] <- "red"
                
                plot.igraph(A,layout=layout_with_kk,vertex.color = vertex_colors,vertex.size=5,edge.color="darkgreen",vertex.label.font=1,vertex.label.cex=0.6, height=800, width=800, rescale = T)
              }
              
              #Análisis de red
              
              #save graph object in the list info
              red_info[[as.character(tao)]][["A"]] <- A
              
              #Número de nodos
              red_info[[as.character(tao)]][["Length"]] <-  length(V(A))
              
              #Extrae clusters
              red_info[[as.character(tao)]][["Clusters"]] <- clusters(A, mode=c("weak"))
              write.csv(data.frame(red_info[[as.character(tao)]][["Clusters"]]$membership), sprintf("../../results/Networks/big/%s/node_genes/node_info.csv", tao))
              
              #Tamaño de clusters
              red_info[[as.character(tao)]][["Clusters size"]] <- clusters(A, mode=c("weak"))$csize
              
              #Transitividad
              red_info[[as.character(tao)]][["Transitividad"]] <- transitivity(A, type="global") 
              
              #Diámetro
              red_info[[as.character(tao)]][["Diametro"]] <- diameter(A)
              
              ## ==================== SAVE NETWORK AND INFORMATION ABOUT DEGREE DEPENDING ON THE NATURE OF DATA (LOG OR NOT) =======================
              if (log_transform){
                if (!dir.exists("../../results/networks/big/log_transformed")){
                  dir.create("../../results/networks/big/log_transformed")
                }
                #Guarda la red
                if (mutual_info){
                  if (!dir.exists(sprintf("../../results/networks/big/%s/mutual_info/graph/", TAO))){
                    dir.create(sprintf("../../results/networks/big/%s/mutual_info/graph/", TAO))
                    write.graph(A,sprintf("../../results/networks/big/%s/mutual_info/graph/edges_%s.txt", tao, tao),format="ncol")
                  }
                } else {
                  if (!dir.exists(sprintf("../../results/networks/big/%s/corr/graph/", TAO))){
                    dir.create(sprintf("../../results/networks/big/%s/corr/graph/", TAO))
                    write.graph(A,sprintf("../../results/networks/big/%s/corr/graph/edges_%s.txt", tao, tao),format="ncol")
                  }
                }
                
                #degree
                deg <- degree(A, mode="all")
                red_info[[as.character(tao)]][["Grado"]] <- deg
                #Guarda tabla con información de grado
                degree.table <- data.frame(Degree = red_info[["0.92"]][["Grado"]])
                
                if (mutual_info){
                  if (!dir.exists(sprintf("../../results/networks/big/%s/mutual_info/degree/", TAO))){
                    dir.create(sprintf("../../results/networks/big/%s/mutual_info/degree/", TAO))
                    write.graph(A,sprintf("../../results/networks/big/%s/mutual_info/degree/edges_%s.txt", tao, tao),format="ncol")
                  }
                } else {
                  if (!dir.exists(sprintf("../../results/networks/big/%s/corr/degree/", TAO))){
                    dir.create(sprintf("../../results/networks/big/%s/corr/degree/", TAO))
                    write.graph(A,sprintf("../../results/networks/big/%s/corr/degree/edges_%s.txt", tao, tao),format="ncol")
                  }
                }
                
                
              } else{
                #Guarda la red
                if (mutual_info){
                  if (!dir.exists(sprintf("../../results/networks/big/%s/mutual_info/graph/", TAO))){
                    dir.create(sprintf("../../results/networks/big/%s/mutual_info/graph/", TAO))
                    write.graph(A,sprintf("../../results/networks/big/%s/mutual_info/graph/edges_%s.txt", tao, tao),format="ncol")
                  }
                } else {
                  if (!dir.exists(sprintf("../../results/networks/big/%s/corr/graph/", TAO))){
                    dir.create(sprintf("../../results/networks/big/%s/corr/graph/", TAO))
                    write.graph(A,sprintf("../../results/networks/big/%s/corr/graph/edges_%s.txt", tao, tao),format="ncol")
                  }
                }
                
                #degree
                deg <- degree(A, mode="all")
                red_info[[as.character(tao)]][["Grado"]] <- deg
                #Guarda tabla con información de grado
                degree.table <- data.frame(Degree = red_info[["0.92"]][["Grado"]])
                
                if (mutual_info){
                  if (!dir.exists(sprintf("../../results/networks/big/%s/mutual_info/degree/", TAO))){
                    dir.create(sprintf("../../results/networks/big/%s/mutual_info/degree/", TAO))
                    write.graph(A,sprintf("../../results/networks/big/%s/mutual_info/degree/edges_%s.txt", tao, tao),format="ncol")
                  }
                } else {
                  if (!dir.exists(sprintf("../../results/networks/big/%s/corr/degree/", TAO))){
                    dir.create(sprintf("../../results/networks/big/%s/corr/degree/", TAO))
                    write.graph(A,sprintf("../../results/networks/big/%s/corr/degree/edges_%s.txt", tao, tao),format="ncol")
                  }
                }
                
              }
              #============================================================================================================================
              
              barplot(table(deg),las=2)
              plot(A, vertex.size=deg)
              plot.igraph(A,layout=layout_with_kk,vertex.color = 'lightblue',vertex.size=deg,edge.color="darkgreen",vertex.label.font=1,vertex.label.cex=0.6 )
              
              hist(deg, breaks=1:vcount(A)-1, main="Histogram of node degree")
              
              
              deg.dist <- degree_distribution(A, cumulative=T, mode="all")
              plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
                    xlab="Degree", ylab="Cumulative Frequency")
              
              #Guarda información de centralización
              red_info[[as.character(tao)]][["centr_degree"]] <- centr_degree(A, mode="in", normalized=T)
              # cd<-centr_degree(A, mode="in", normalized=T)
              # cd$centralization
              # cd$res
              
              #Eigen centrality
              red_info[[as.character(tao)]][["eigen_centrality"]] <- eigen_centrality(A, directed=F, weights=NA)
              # eigc<-eigen_centrality(A, directed=F, weights=NA)
              
              #Hub score
              red_info[[as.character(tao)]][["hub_score"]] <- hub_score(A, weights=NA)$vector
              # hs <- hub_score(A, weights=NA)$vector
              # hs
              # sort(hs,decreasing = TRUE)
              # sort(deg,decreasing = TRUE)
              # 
              
              }
            }
          }
        }
      }
    }
  }
}