# This script contains helper functions to generate the figures from 
# Linnenbrink, Milà, Ludwig and Meyer (2023)


#' Function to plot integral between Gj and Gij
#' @details
#' Plotting of Gj and Gij modified from https://github.com/HannaMeyer/CAST/blob/master/R/plot.R,
#' calculation of the integral between Gj and Gij taken from https://github.com/cdowd/twosamples
#' @param x kNNDM-object.
plot.knndm.integral <- function(x, ...){
  
  
  # Prepare data for plotting: Gij function
  Gij_df <- data.frame(r=x$Gij[order(x$Gij)])
  Gij_df$Function <- "1_Gij(r)"
  
  # Prepare data for plotting: G function
  Gj_df <- data.frame(r=x$Gj[order(x$Gj)])
  Gj_df$Function <- "2_Gj(r)"
  
  # Merge data for plotting
  Gplot <- rbind(Gij_df, Gj_df)
  
  n = length(Gij_df$r)
  dfx = data.frame(sort(Gij_df$r),cumsum(rep(1/n,n)),"x")
  n1=length(Gj_df$r)
  dfy = data.frame(sort(Gj_df$r),cumsum(rep(1/n1,n1)),"y")
  colnames(dfx) = colnames(dfy) = c("point","ecdf","sample")
  df = rbind(dfy,dfx)
  
  combined = sort(c(Gij_df$r,Gj_df$r))
  
  midpoints = combined[-(2*n)]+diff(combined)/2
  fhat = ehat = dhat = width = numeric(2*n-1)
  for (i in 1:(2*n-1)) {
    point = midpoints[i]
    width[i] = (point-combined[i])*2
    fhat[i] = mean(Gj_df$r < point)
    ehat[i] = mean(Gij_df$r < point)
    dhat[i] = mean(combined < point)
  }
  gap = fhat-ehat
  
  
  x1 = midpoints-0.5*width
  x2 = midpoints+0.5*width
  y1 = ehat
  y2 = fhat
  block = 1:(2*n-1)
  alpha = 1/(dhat*(1-dhat))
  
  df3.1 = cbind(x1,y1,block,alpha)
  df3.2 = cbind(x1,y2,block,alpha)
  df3.3 = cbind(x2,y2,block,alpha)
  df3.4 = cbind(x2,y1,block,alpha)
  
  df3 = rbind(df3.1,df3.2,df3.3,df3.4)
  colnames(df3) = c("point","ecdf","block","weights")
  df3 = as_tibble(df3) |> arrange(block)
  
  # Plot
  ggplot2::ggplot(data=Gplot) +
    ggplot2::geom_vline(xintercept=0, lwd = 0.1, ggplot2::aes_string(x="r", group="Function", col="Function")) +
    ggplot2::geom_hline(yintercept=0, lwd = 0.1, ggplot2::aes_string(x="r", group="Function", col="Function")) +
    ggplot2::geom_hline(yintercept=1, lwd = 0.1, ggplot2::aes_string(x="r", group="Function", col="Function")) +
    geom_polygon(data=df3,aes(x=point,y=ecdf,group=block),alpha=0.7, fill="orange") +
    ggplot2::stat_ecdf(geom = "step", lwd = 1, ggplot2::aes_string(x="r", group="Function", col="Function")) +
    ggplot2::scale_colour_manual(values=c("#000000", "#56B4E9"),
                                 labels=c(expression(hat(G)[ij](r)),
                                          expression(hat(G)[j]^"*"*"(r,L)"),
                                          expression(hat(G)[j](r)))) +
    ggplot2::ylab(expression(paste(hat(G)[ij](r), ", ",
                                   hat(G)[j]^"*"*"(r,L)", ", ",
                                   hat(G)[j](r))))
}

#' Function that keeps every split of kNNDM to reproduce fig 2
#' @details
#' Manipulated kNNDM-function to keep every split.
#' @param k integer. Number of folds.
#' @param maxp numeric. Maximum fold size allowed.
#' @param tpoints sf or sfc point object. training points.
#' @param ppoints sf or sfc point object. prediction points.
k_plot <- function(k, maxp, tpoints, ppoints) {
  
  distclust <- function(tr_coords, folds){
    alldist <- rep(NA, length(folds))
    for(f in unique(folds)){
      alldist[f == folds] <- c(FNN::knnx.dist(query = tr_coords[f == folds,,drop=FALSE],
                                              data = tr_coords[f != folds,,drop=FALSE], k = 1))
    }
    alldist
  }
  
  # Gj: NN distance function for a cluster per point, i.e. LOO CV
  tcoords <- sf::st_coordinates(tpoints)[,1:2]
  Gj <- c(FNN::knn.dist(tcoords, k = 1))
  
  # Gij: prediction to training NN distances
  Gij <- c(FNN::knnx.dist(query = sf::st_coordinates(ppoints)[,1:2],
                          data = tcoords, k = 1))
  
  # Build grid of number of clusters to try - we sample low numbers more intensively
  clustgrid <- data.frame(nk = as.integer(round(exp(seq(log(k), log(nrow(tpoints)-1),
                                                        length.out = 6)))))
  clustgrid$W <- NA
  clustgrid <- clustgrid[!duplicated(clustgrid$nk),]
  
  clustgroups=Gjstar_i <- list()
  
  # Compute 1st PC for ordering clusters
  pcacoords <- stats::prcomp(tcoords, center = TRUE, scale. = FALSE, rank = 1)
  
  j = 0
  for(nk in clustgrid$nk){
    
    # Cut nk clusters
    clust_nk <- stats::kmeans(tcoords, nk)$cluster
    tabclust <- as.data.frame(table(clust_nk))
    tabclust <- tabclust[order(tabclust$Freq, decreasing=T),]
    tabclust$clust_k <- NA
    
    # compute cluster centroids and apply PC loadings to shuffle along the 1st dimension
    centr_tpoints <- sapply(tabclust$clust_nk, function(x){
      centrpca <- matrix(apply(tcoords[clust_nk %in% x, , drop=FALSE], 2, mean), nrow = 1)
      colnames(centrpca) <- colnames(tcoords)
      return(predict(pcacoords, centrpca))
    })
    tabclust$centrpca <- centr_tpoints
    tabclust <- tabclust[order(tabclust$centrpca),]
    
    # We don't merge big clusters
    clust_i <- 1
    for(i in 1:nrow(tabclust)){
      if(tabclust$Freq[i] >= nrow(tpoints)/k){
        tabclust$clust_k[i] <- clust_i
        clust_i <- clust_i + 1
      }
    }
    rm("clust_i")
    
    # And we merge the remaining into k groups
    clust_i <- setdiff(1:k, unique(tabclust$clust_k))
    tabclust$clust_k[is.na(tabclust$clust_k)] <- rep(clust_i, ceiling(nk/length(clust_i)))[1:sum(is.na(tabclust$clust_k))]
    tabclust2 <- data.frame(ID = 1:length(clust_nk), clust_nk = clust_nk)
    tabclust2 <- merge(tabclust2, tabclust, by = "clust_nk")
    tabclust2 <- tabclust2[order(tabclust2$ID),]
    clust_k <- tabclust2$clust_k
    
    # Compute W statistic if not exceeding maxp
    j <- j+1
    if(!any(table(clust_k)/length(clust_k)>maxp)){
      Gjstar_i[[j]] <- distclust(tcoords, clust_k)
      clustgrid$W[clustgrid$nk==nk] <- twosamples::wass_stat(Gjstar_i[[j]], Gij)
      clustgroups[[paste0("nk", nk)]] <- clust_k
    }
  }
  list(clustgroups,clustgrid, Gjstar_i, Gij, Gj, clustgrid$nk)
}

#' Function to plot ECDFs of Gi, Gij and Gij*
#' @details
#' see CAST (https://github.com/HannaMeyer/CAST/blob/master/R/plot.R)
#' @param x kNNDM-object.
#' @param overlay are two lines overlaying?
plot.knndm <- function(x, overlay=FALSE, ...){
  
  # Prepare data for plotting: Gij function
  Gij_df <- data.frame(r=x$Gij[order(x$Gij)])
  Gij_df$Function <- "1_Gij(r)"
  
  # Prepare data for plotting: Gjstar function
  Gjstar_df <- data.frame(r=x$Gjstar[order(x$Gjstar)])
  Gjstar_df$Function <- "2_Gjstar(r)"
  
  # Prepare data for plotting: G function
  Gj_df <- data.frame(r=x$Gj[order(x$Gj)])
  Gj_df$Function <- "3_Gj(r)"
  
  # Merge data for plotting
  Gplot <- rbind(Gij_df, Gjstar_df, Gj_df)
  
  # Plot
  if (isFALSE(overlay)) {
    ggplot2::ggplot(data=Gplot, ggplot2::aes_string(x="r", group="Function", col="Function")) +
      ggplot2::geom_vline(xintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=1, lwd = 0.1) +
      ggplot2::stat_ecdf(geom = "step", lwd = 1) +
      ggplot2::scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9"),
                                   labels=c(expression(hat(G)[ij](r)),
                                            expression(hat(G)[j]^"*"*"(r,L)"),
                                            expression(hat(G)[j](r)))) +
      ggplot2::ylab(expression(paste(hat(G)[ij](r), ", ",
                                     hat(G)[j]^"*"*"(r,L)", ", ",
                                     hat(G)[j](r))))
  } else if (isTRUE(overlay)) {
    ggplot2::ggplot(data=Gplot, ggplot2::aes_string(x="r", group="Function", col="Function", size="Function")) +
      ggplot2::geom_vline(xintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=0, lwd = 0.1) +
      ggplot2::geom_hline(yintercept=1, lwd = 0.1) +
      ggplot2::stat_ecdf(geom = "step") +
      ggplot2::scale_size_manual(values = c(1,1,0.5),
                                 labels=c(expression(hat(G)[ij](r)),
                                          expression(hat(G)[j]^"*"*"(r,L)"),
                                          expression(hat(G)[j](r)))) +
      ggplot2::scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9"),
                                   labels=c(expression(hat(G)[ij](r)),
                                            expression(hat(G)[j]^"*"*"(r,L)"),
                                            expression(hat(G)[j](r)))) +
      ggplot2::ylab(expression(paste(hat(G)[ij](r), ", ",
                                     hat(G)[j]^"*"*"(r,L)", ", ",
                                     hat(G)[j](r))))
  }
  
}
