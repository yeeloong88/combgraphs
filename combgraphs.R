#############
## GENERIC ##
#############
combgraphs <-  function (x, ...) UseMethod("combgraphs")



#############
## DEFAULT ##
#############
combgraphs.default <- function(x, ngrp=NULL, plotAll=FALSE, outputAll=FALSE, plot=TRUE, layout=layout.circle, method=c("union","intersect"), col.pal=funky, truenames=TRUE, ...){
  stop(paste("No method for objects of class",class(x)))
} # end gengraph.default



###############################
## LIST ('dist' OR 'matrix') ##
###############################
##
## this is the basic method
##
combgraphs.list <- function(x, ngrp=NULL, plotAll=FALSE, outputAll=FALSE, plot=TRUE, layout=layout.circle, method=c("union","intersect"), col.pal=funky ,truenames=TRUE, ...){
  ## CHECK1: IGRAPH ##
  if(!require("igraph")) stop("igraph is required")
  if(!require("adegenet")) stop("adegenet is required")
  
  ## CHECK2: LIST OF MATRICES ##
  if(!is.list(x)) stop("List of distances or distance matrices required")
  
  ## CHECK3: LIST CONTAINS ONLY MATRICES ##
  K <- length(x)  # K is the number of distance matrices
  for(k in 1:K){if(class(x[[k]])=="dist") x[[k]] <- as.matrix(x[[k]])}
  for(k in 1:K){if(!is.matrix(x[[k]])) stop(paste("No method for objects of class",class(x[[k]])))}
  
  ## CHECK4: DISTANCE MATRICES ARE INFORMATIVE ##
  del.index <- c()
  for(k in 1:K){
    if(sum(x[[k]])==0){
      del.index <- c(del.index,k) 
      warning(paste("Data from distance matrix",k,"is uninformative and is discarded from analysis"),immediate.=T)
    } 
  }
  if(length(del.index)!=0){
    x <- x[-del.index]  #discard uninformative distance matrices
    K <- length(x)
  }
  
  
  ## GET GENGRAPH OUTPUTS IN A LIST ##
  g <- list()
  for(k in 1:K){
    g[[k]] <- gengraph.matrix(x[[k]],ngrp,truenames=truenames)
    #g[[k]]$graph$name <- paste("Data Type",k)
  }
  
  l <- layout(g[[1]]$graph)
  igraph.options(sparsematrices=TRUE)  #Use sparse matrices for plotting igraph
  
  if(plotAll){
    chooseAgain <- TRUE
    while(chooseAgain){
      cat("\nPlease enter the distance matrix would you like to plot:  ")
      ans <- NA
      while(is.null(ans) || is.na(ans)) suppressWarnings(ans <- as.numeric(readLines(n = 1)))
      plot.igraph(g[[ans]]$graph, main=paste("Distance matrix",ans), layout=l, edge.label.cex=0.8)
      if(truenames){
        V(g[[ans]]$graph)$label <- rownames(x[[k]])
      }
      ans <- ""  
      while(!ans %in% c("c","q")){
       cat("\nPlease enter 'c' to plot another distance matrix or 'q' to quit: ")
       ans <- tolower(readLines(n = 1))
      }
      if(ans=="q") chooseAgain <- FALSE
    }
  }
  
  if(outputAll){
    return(g)
  }
 
  ## GET INDIVIDUAL ADJACENT MATRICES ##
  x.adj <- x
  for (k in 1:K){
    x.adj[[k]][x.adj[[k]] <= g[[k]]$cutoff] <- 1
    x.adj[[k]][x.adj[[k]] > g[[k]]$cutoff] <- 0
  } 
  
  ## GET WEIGHTS FOR COMBINED GRAPHS ##
  x.wts <- x
  for (k in 1:K){
    if(max(x.adj[[k]]!=0)){
      x.wts[[k]][x.wts[[k]]==0] <- 0.0000001  # Assign small non-zero value to 'perfectly' connected individuals
      x.wts[[k]][x.wts[[k]] > g[[k]]$cutoff] <- 0 # Assign 0 weight to non-connected individuals
      x.wts[[k]] <- x.wts[[k]]/max(x.wts[[k]]) # Standardised weights (i.e. w^k_ij)
    }else warning("The maximum distance in the distance matrix is 0")
  } 
  comb.wts <- round((Reduce('+',x.wts))/K,2)    # (i.e. omega_ij) 
  
  method <- match.arg(method)
  switch(method,
         ## COMBINING BY UNION ##
         union = {
           comb.adj <- Reduce('+',x.adj)   # union: adding up adjacency distance matrices 
           comb.adj[comb.adj!=0] <-1
           comb.adj <- comb.adj*comb.wts   # getting weighted adjacency matrix with omega_ij
           comb.g <- graph.adjacency(comb.adj, mode="undirected", weighted=TRUE, diag=FALSE)
           comb.g$name <- "UNION"
           clust <- clusters(comb.g)   #get number of clusters
           V(comb.g)$color <- col.pal(clust$no)[clust$membership]
           col <- col.pal(clust$no)[1:clust$no]
           V(comb.g)$label <- rownames(comb.adj)
           if(length(E(comb.g))>0){
             E(comb.g)$label <- E(comb.g)$weight
           } #assign labels to edges
           if(plot){
             plot.igraph(comb.g, main="Union of graphs",layout=l, edge.label.cex=0.8)   #plot graph
           }
           res <- list(graph=comb.g, clust=clusters(comb.g), col=col)    #store results
           return(res)
         },
         ## COMBINING BY INTERSECT ##
         intersect= {   
           comb.adj <- Reduce('*',x.adj)   # intersect: multiplying adjacency distance matrices
           comb.adj[comb.adj!=0] <-1
           comb.adj <- comb.adj*comb.wts   # getting weighted adjacency matrix with omega_ij
           comb.g <- graph.adjacency(comb.adj, mode="undirected", weighted=TRUE, diag=FALSE)
           comb.g$name <- "INTERSECT"
           clust <- clusters(comb.g)   #get number of clusters
           V(comb.g)$color <- col.pal(clust$no)[clust$membership]
           col <- col.pal(clust$no)[1:clust$no]
           V(comb.g)$label <- rownames(comb.adj)
           if(length(E(comb.g))>0){
             E(comb.g)$label <- E(comb.g)$weight
           } #assign labels to edges
           if(plot){
             plot.igraph(comb.g, main="Intersect of graphs",layout=l, edge.label.cex=0.8)   #plot graph 
           }
           res <- list(graph=comb.g, clust=clusters(comb.g), col=col)    #store results
           return(res)
         }
    )
}#end of method

 