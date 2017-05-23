set.seed(1)
### three-dimensional case
### p = 3 and n = 100
x<-rbind(matrix(rnorm(90,sd=0.1),ncol=3),
         matrix(rnorm(90,mean=1,sd=0.2),ncol=3),
         matrix(rnorm(90,mean=5,sd=0.1),ncol=3),
         matrix(rnorm(90,mean=7,sd=0.2),ncol=3))

x <- x[,-3]

set.seed(8)
x.miss <- x
x.miss[sample(1:(1*120),20),1] <- NA
x <- x.miss

plot((hclust(dist(x))))

hcluster <- hclust(dist(x,method="euclidean"),method = "complete")

# defining some clusters
mycl <- cutree(hcluster, k=2:10)

res <- NbClust(x, distance = "euclidean", min.nc=2, max.nc=8, 
             method = "complete", index = "ch")
res$All.index
res$Best.nc
res$Best.partition

ch.index <- function(x, cl) {
  n <- length(cl)
  k <- max(cl)
  centers <- matrix(nrow = k, ncol = ncol(x))
  for (i in 1:k) {
    if (ncol(x) == 1) {
      centers[i, ] <- mean(x[cl == i, ])
    }
    if (is.null(dim(x[cl == i, ]))) {
      bb <- matrix(x[cl == i, ], byrow = FALSE, nrow = 1, 
                   ncol = ncol(x))
      centers[i, ] <- apply(bb, 2, mean)
    }
    else {
      centers[i, ] <- apply(x[cl == i, ], 2, mean,na.rm=TRUE) ### mean for cluster i
    }
  }
  allmean <- apply(x, 2, mean,na.rm=TRUE) ### overall mean
  dmean <- sweep(x, 2, allmean, "-")
  ## dmean2 <- t(t(x) - allmean)
  ## dmean2 == dmean
  allmeandist <- sum(dmean^2,na.rm=TRUE) ## total sum of squares
  withins <- rep(0, k)
  x.2 <- (x - centers[cl, ])^2
  for (i in 1:k) {
    withins[i] <- sum(x.2[cl == i, ],na.rm=TRUE)
  }
  wgss <- sum(withins)
  bgss <- allmeandist - wgss
  
  ### calculating bgss ## between cluster variation
  bgss_star <- 0
  for(i in 1:k){
    bgss_star <- bgss_star + (sum(cl==i)*((sum((centers[i,] - allmean)^2,na.rm=TRUE))))
  }
  
  ### bgss_star is the same as bgss if there's no missing data!
  
  ch <- (bgss/(k-1))/(wgss/(n-k))
  
  results <- ch
  return(results)
}


indices <- apply(mycl,2,ch.index,x=x)
plot(2:10,indices)


