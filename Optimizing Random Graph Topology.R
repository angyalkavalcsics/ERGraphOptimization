################################################################################
# Problem #1
################################################################################
# Erdos-Renyi Random Graph
ER <- function(p, beta){
  # p choose 2 possible edges e without self-loops
  e <- (p*(p-1))/2
  # make adjacency matrix X
  X <- matrix(0, p, p)
  # create e possible random edges, X is a Bernoulli(beta) random variable
  edges <- rbinom(e, 1, beta)
  # Fill X
  X[upper.tri(X, diag = FALSE)] <- edges
    # make X symmetric
  X[lower.tri(X)] = t(X)[lower.tri(X)]
  return(list("graph" = X, "p" = p, "beta" = beta))
}

# Tests
A <- ER(p = 3, beta = 0)
A$graph
B <- ER(p = 3, beta = 1)
B$graph
C <- ER(p = 5, beta = 0.2)
C$graph

################################################################################
# Problem #2
################################################################################
# Probability distribution associated with the 
# degree sequence of an ER random graph.
degSeq <- function(X){
  # d is the row sums of the adjacency matrix X
  p <- X$p
  beta <- X$beta
  X <- X$graph
  d <- rowSums(X)
  hist(d, main = paste("p = ", p, ", beta = ", beta))
}

# Tests:
par(mfrow = c(2,2))

X1 <- ER(p = 1000, beta = 0.05)
degSeq(X1)

X2 <- ER(p = 50, beta = 0.95)
degSeq(X2)

X3 <- ER(p = 1000, beta = 0.95)
degSeq(X3)

X4 <- ER(p = 50, beta = 0.05)
degSeq(X4)

################################################################################
# Problem #3
################################################################################
# Dijkstra returns shortest path lengths 
# from source vertex s to all other vertices

# X is an adjacency matrix
dijkstra <- function(X, s){
  A <- X$graph
  # number of vertices p
  p <- X$p
  #set of vertex indices
  Q <- c(1:p)
  # initialize distances to source
  dist <- rep(Inf, p)
  # shortest path from a vertex to itself is 0
  dist[s] <- 0
  while (length(Q) > 0) {
    # Find index of minimum dist
    u <- which.min(dist[Q[1]:Q[length(Q)]])
    # Remove that element from the set of indices
    Q <- Q[-u]
    # Find the neighborhood of that element i.e,
    # the indices for which the adjacency matrix has a 1
    v <- which(A[,u] == 1)
    # For each element in the neighborhood, we calculate the length
    # of the path from the source to that element and return
    # only the shortest path
    for (neighbor in v) {
      temp <- dist[u] + A[u, neighbor]
      if(temp < dist[neighbor]){
        dist[neighbor] <- temp
      }
    }
  }
  dist[dist == Inf] <- 0
  return(dist)
}
A <- ER(p = 10, beta = 0.7)

for(i in 1:A$p){
  print(dijkstra(A, i))
}
FW(A)
################################################################################
# Problem #4
################################################################################
# Floyd-Warshall algorithm for computing the lengths 
# of the shortest paths between all pairs of vertices.
FW <- function(X){
  A <- X$graph
  p <- X$p
  D <- A
  # set all elements in D == 0 to Inf
  for (i in 1:p) {
    for (j in 1:p) {
      if(A[i, j] == 0){
        D[i, j] = Inf
        D[j, i] = Inf
      }
    }
  }
  # shortest path from a vertex to itself is 0
  diag(D) = 0
  # compare all possible paths through the graph 
  # between each pair of vertices
  for (k in 1:p) {
    for (i in 1:p) {
      for (j in 1:p) {
        if(D[i,j] > (D[i, k] + D[k, j])){
          D[i,j] <- D[i,k] + D[k,j]
          D[j,i] <- D[i,j]
        }
      }
    }
  }
  return(D)
}

################################################################################
# Problem #5
################################################################################
createCirclePlot <- function(X){
  A <- X$graph
  p <- X$p# number of points you want on the unit circle
  pts.circle <- t(sapply(1:p,function(r)c(cos(2*r*pi/p),sin(2*r*pi/p))))
  plot(pts.circle, col='red', pch=19, xlab='x', ylab='y', main = "Circle Plot")
  vertices.label <- c(seq(1:p))
  for (vertex in vertices.label) {
    vertices.label[vertex] <- paste0("V", vertex)
  }
  
  text(pts.circle, labels = vertices.label, cex=0.9, pos=4)
  
  for (v in 1:p) {
    adj <- which(A[,v] == 1)
    segments(pts.circle[v, 1], pts.circle[v, 2], pts.circle[adj, 1], pts.circle[adj, 2])
  }
  return(list(init.x = pts.circle[,1], init.y = pts.circle[,2]))
}

################################################################################
# Problem #6
################################################################################
# A function for the computing the value of the loss function
loss.value <- function(x, y, L){
  l <- 0
  for (i in 1:p) {
    for (j in 1:p) {
      l <- l + (sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2) - L[i, j])
    }
  }
  return(l)
}
# A function for the partial derivative of L with 
# respect to xk
partial.xk <- function(x, y, L, k){
  partial <- 0
  for (j in 1:p) {
    if(j != k){
      dkj <- sqrt((x[k]-x[j])^2 + (y[k]-y[j])^2)
      partial <- partial + (2*(x[k]-x[j])*(dkj - L[k, j]))/(dkj)
    }
  }
  return(partial)
}
# A function for the partial derivative of L with 
# respect to yk
partial.yk <- function(x, y, L, k){
  partial <- 0
  for (j in 1:p) {
    if(j != k){
      dkj <- sqrt((x[k]-x[j])^2 + (y[k]-y[j])^2)
      partial <- partial + (2*(y[k]-y[j])*(dkj - L[k, j]))/(dkj)
    }
  }
  return(partial)
}

# A generalized plotting function for plotting the graph
plot.GD <- function(pts, X){
  A <- X$graph
  p <- X$p
  plot(pts, col='red', pch=19, xlab='x', ylab='y', 
       main = "Optimize Loss by Gradient Descent")
  vertices.label <- c(seq(1:p))
  for (vertex in vertices.label) {
    vertices.label[vertex] <- paste0("V", vertex)
  }
  
  text(pts, labels = vertices.label, cex=0.9, pos=4)
  
  for (v in 1:p) {
    adj <- which(A[,v] == 1)
    segments(pts[v, 1], pts[v, 2], pts[adj, 1], pts[adj, 2])
  }
}
# Gradient Descent
minLoss.GD <- function(A, eps = 1e-4, alpha = 0.01){
  # Get Adjacency graph and p = # of vertices
  X <- A$graph
  p <- A$p
  
  # Matrix of shortest path lengths
  L <- FW(A)

  # Get initial x and y vectors from circle plot
  init.pts <- t(sapply(1:p,function(r)c(cos(2*r*pi/p),sin(2*r*pi/p))))
  x1 <- init.pts[,1]
  y1 <- init.pts[,2]
  # plot initial x and y coordinates
  pts <- cbind(x1, y1)
  
  plot.GD(pts, A)
  
  x0 <- x1 + 1
  y0 <- y1 + 1
  
  ct <- 0
  while (mean(abs(x0-x1)) > eps || mean(abs(y0-y1)) > eps) {
    ct <- ct + 1
    #print(ct)
    Dx <- rep(0,p)
    Dy <- rep(0,p)
    # Find partial derivatives
    for (k in 1:p) {
      Dx[k] <- partial.xk(x1, y1, L, k)
      Dy[k] <- partial.yk(x1, y1, L, k)
    }
    # Update x and y
    x0 <- x1
    x1 <- x0 - alpha*Dx
    
    y0 <- y1
    y1 <- y0 - alpha*Dy
    
    # Plot each update
    pts <- cbind(x1, y1)
    plot.GD(pts, A)
    Sys.sleep(0.05)
  }
  loss <- loss.value(x1, y1, L)
  return(list("opt.x" = x1, "opt.y" = y1, "iterations" = ct, "loss" = loss))
}

plot.SGD <- function(pts, X){
  A <- X$graph
  p <- X$p
  plot(pts, col='red', pch=19, xlab='x', ylab='y', 
       main = "Optimize Loss by Stochastic Gradient Descent")
  vertices.label <- c(seq(1:p))
  for (vertex in vertices.label) {
    vertices.label[vertex] <- paste0("V", vertex)
  }
  
  text(pts, labels = vertices.label, cex=0.9, pos=4)
  
  for (v in 1:p) {
    adj <- which(A[,v] == 1)
    segments(pts[v, 1], pts[v, 2], pts[adj, 1], pts[adj, 2])
  }
}
# Stochastic Gradient Descent
minLoss.SGD <- function(A, eps = 1e-4, alpha = 0.01){
  # Get Adjacency graph and p = # of vertices
  X <- A$graph
  p <- A$p
  
  # Matrix of shortest path lengths
  L <- FW(A)
  
  # Get initial x and y vectors from circle plot
  init.pts <- t(sapply(1:p,function(r)c(cos(2*r*pi/p),sin(2*r*pi/p))))
  x1 <- init.pts[,1]
  y1 <- init.pts[,2]
  # plot initial x and y coordinates
  pts <- cbind(x1, y1)
  
  plot.SGD(pts, A)
  
  x0 <- x1 + 1
  y0 <- y1 + 1
  
  ct <- 0
  while (mean(abs(x0-x1)) > eps || mean(abs(y0-y1)) > eps) {
    ct <- ct + 1
    print(ct)
    k <- sample.int(p, 1)
    # Find partial derivatives
    dx <- partial.xk(x1, y1, L, k)
    dy <- partial.yk(x1, y1, L, k)
    
    # Update x and y
    x0[k] <- x1[k]
    x1[k] <- x0[k] - alpha*dx
    
    y0[k] <- y1[k]
    y1[k] <- y0[k] - alpha*dy
    
    # Plot each update
    pts <- cbind(x1, y1)
    plot.SGD(pts, A)
    Sys.sleep(0.05)
  }
  loss <- loss.value(x1, y1, L)
  return(list("opt.x" = x1, "opt.y" = y1, "iterations" = ct, "loss" = loss))
}

B <- ER(p = 10, beta = 0.8)
out2 <- minLoss.SGD(B, eps = 1e-4, alpha = 0.01)
out2

# Coordinate Descent
# A generalized plotting function for plotting the graph
plot.CD <- function(pts, X){
  A <- X$graph
  p <- X$p
  plot(pts, col='red', pch=19, xlab='x', ylab='y', 
       main = "Optimize Loss by Coordinate Descent")
  vertices.label <- c(seq(1:p))
  for (vertex in vertices.label) {
    vertices.label[vertex] <- paste0("V", vertex)
  }
  
  text(pts, labels = vertices.label, cex=0.9, pos=4)
  
  for (v in 1:p) {
    adj <- which(A[,v] == 1)
    segments(pts[v, 1], pts[v, 2], pts[adj, 1], pts[adj, 2])
  }
}
minLoss.CD <- function(A, eps = 1e-4, alpha = 0.01){
  # Get Adjacency graph and p = # of vertices
  X <- A$graph
  p <- A$p
  
  # Matrix of shortest path lengths
  L <- FW(A)
  
  # Get initial x and y vectors from circle plot
  init.pts <- t(sapply(1:p,function(r)c(cos(2*r*pi/p),sin(2*r*pi/p))))
  x1 <- init.pts[,1]
  y1 <- init.pts[,2]
  # plot initial x and y coordinates
  pts <- cbind(x1, y1)
  plot.CD(pts, A)
  
  x0 <- x1 + 1
  y0 <- y1 + 1
  
  ct <- 0
  for(k in 1:p){
    while(mean(abs(x0[k]-x1[k])) > eps) {
      dx <- partial.xk(x1, y1, L, k)
      x0[k] <- x1[k]
      x1[k] <- x0[k] - alpha*dx
      pts <- cbind(x1, y1)
      ct <- ct + 1
      #print(ct)
      plot.CD(pts, A)
      Sys.sleep(0.05)
    }
    while(mean(abs(y0[k]-y1[k])) > eps){
      dy <- partial.yk(x1, y1, L, k)
      y0[k] <- y1[k]
      y1[k] <- y0[k] - alpha*dy
      pts <- cbind(x1, y1)
      plot.CD(pts, A)
      Sys.sleep(0.05)
    }
  }
  loss <- loss.value(x1, y1, L)
  return(list("opt.x" = x1, "opt.y" = y1, "iterations" = ct, "loss" = loss))
}

A <- ER(p = 10, beta = 0.8)
minLoss.GD(A)
minLoss.SGD(A)
minLoss.CD(A)

createCirclePlot(A)


