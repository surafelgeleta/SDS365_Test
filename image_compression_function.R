library(readr)

image1 <- read.csv("~/SDS365/image_compression/HW2imagefiles/image1.csv")
image2 <- read.csv("~/SDS365/image_compression/HW2imagefiles/image2.csv")
image3 <- read.csv("~/SDS365/image_compression/HW2imagefiles/image3.csv")
image4 <- read.csv("~/SDS365/image_compression/HW2imagefiles/image4.csv")

image1 <- data.matrix(image1)
image2 <- data.matrix(image2)
image3 <- data.matrix(image3)
image4 <- data.matrix(image4)

image_error <- function(M, M_k){
  # Frobenius norm of difference X - X_k
  error_k <- norm(M - M_k, type = "F")
  return (error_k)
}

compress_image <- function(k, M){
  # Function takes n x p matrix M and integer k
  # Returns: rank k approximation of image matrix using PCA
  # and, error of rank k approximation w/Frobenius norm
  
  # Cross product of M
  V <- crossprod(M, M)
  
  # Eigendecomposition of V
  eigendecomp_V <- eigen(V)
  
  # Eigenvector matrix
  U <- eigendecomp_V$vectors
  
  # PC Score matrix
  Z <- M %*% U
  
  # Initialize rank approx
  M_k <- matrix(0, nrow = nrow(M), ncol = ncol(M))
  
  # Rank approximation
  for (i in 1:k){
    M_k <- M_k + (Z[, i] %o% U[, i])
  }
  
  # Frobenius norm of difference M - M_k
  error_k <- image_error(M, M_k)
  
  return (list(M_k = M_k, error_k = error_k, U = U, Z = Z))
}

error_k_plot <- function(M, k){
  
  k_seq <- seq(1, k, by = 1)
  
  Xk_grid <- sapply(k_seq, compress_image, M = M)
  
  return (Xk_grid)
}

image3_10 <- compress_image(image3, 345) 
image(image3,)
image(image3_10$X_k, axes = FALSE)

image1_error <- error_k_plot(image3)
