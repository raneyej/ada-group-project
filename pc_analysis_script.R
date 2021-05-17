# Principal Component Analysis Package
#
# Here are the functions for the pcanalysis package
#
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

pcanalysis.pca <- function(x, center = TRUE, scale = FALSE, svd = TRUE) {
  x <- as.matrix(x) #input as matrix
  k <- ncol(x) #number of samples/variables
  x <- scale(x, center = center, scale = scale)
  if(svd) {        #SVD method (default)
    sv <- svd(x)
    pcs <- sv$v
    d <- sv$d
    pcs.t <- t(pcs)
    x.t <- t(x)
    pcx <- pcs.t %*% x.t
    pcx <- t(pcx)
    colnames(pcs) <- paste('PC', c(1:ncol(pcs)), sep = '')
    rownames(pcs) <- colnames(x)
    pca <- list(x=head(x), pcx=pcx, k=k, PCs=pcs,    #final list for output
                sdev= d / sqrt(max(1, nrow(x) - 1)))
  } else {
    R <- cov(x)   #covarance matrix
    eig <- eigen(R)
    pcs <- eig$vectors[,1:k] #get PCs
    pcs.t <- t(pcs)
    x.t <- t(x)
    pcx <- pcs.t %*% x.t
    pcx <- t(pcx)
    colnames(pcs) <- paste('PC', c(1:ncol(pcs)), sep = '')
    rownames(pcs) <- colnames(x)
    eigenvals <- eig$values[1:k]
    pca <- list(x=head(x), pcx=pcx, k=k, PCs=pcs,     #final list for output
                Eigenvalues=eigenvals)
  }
  return(pca)
}

#Plotting function
pcanalysis.plot <- function(pca){
  plot(pca$PCs[,1], pca$PCs[,2], #plots PCs 1 and 2
       main = paste("Principal Components 1 and 2"),
       xlab = paste("PC 1"),
       ylab= paste("PC 2"),
       pch=19)
  plot(pca$PCs[,1], pca$PCs[,3], #Plots PCs 1 and 3
       main = paste("Principal Components 1 and 3"),
       xlab = paste("PC 1"),
       ylab= paste("PC 3"),
       pch=19)
  p <- 100*pca$sdev^2 / sum(pca$sdev^2)  #calculates variance for each PC and take the percentage
  percent <- data.frame(p = p, PC = as.factor(1:length(p))) #creates a data frame for easy plotting
  barplot(percent$p, ylim=c(0,100), #Plots the percent of variance explained by each PC
          main="Percent of variance explained",
          xlab="PC",
          ylab="Percent",
          names.arg=percent$PC)
}

#Creates a scree plot for cross-validation to find optimal k
pcanalysis.scree <- function(data, kmax){
  WSS <- as.vector(kmax)
  for(i in 1:kmax){
    WSS[i] <- kmeans(data, centers = i)$tot.withinss
  }
  scree_p <- plot(1:kmax, WSS, type="b",
                  xlab = "Number of clusters",
                  ylab = "Within-group SS",
                  main = "K-means Scree Plot",
                  cex.names=.5)
  
  return(scree_p)
}

pcanalysis.cluster <- function(x, pca, k) {
  #assign samples to clusters based on PCS
  #input should be the raw data matrix, the output from pcanalysis.pca, and the optimal k
  data <- t(x)
  pca <- pca$PCs
  z <- kmeans(data, k)
  princomps <- as.data.frame(pca)
  princomps$Sample <- rownames(princomps)
  clusters <- data.frame("Sample" = z[1])
  clusters$Sample <- rownames(clusters)
  rownames(clusters) <- NULL
  nd <- merge(princomps, clusters, by="Sample")
  plot(nd$PC1, nd$PC2, col=factor(nd$cluster), main = "Clustered PCA", xlab ="PC1", ylab="PC2", pch=19)
  return(clusters)
}
