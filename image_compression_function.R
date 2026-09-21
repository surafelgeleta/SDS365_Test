library(readr)
library(tidyverse)

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

error_k <- function(M, k){
  
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
  
  errors <- numeric(k)
  
  for (i in 1:k){
    M_k <- M_k + (Z[, i] %o% U[, i])
    errors[i] <- image_error(M, M_k)
  }
  
  error_data <- data.frame(k = 1:k,
                           error = errors)
  
  return (error_data)
}

# Plot first 10 k images
## Image 1
for (i in 1:10){
  image1_k <- compress_image(i, image1)
  png(paste0("~/SDS365/image_compression/image_approximations/image1_", i, ".png"))
  image(image1_k$M_k, axes = FALSE)
  dev.off()
}

## Image 2
for (i in 1:10){
  image2_k <- compress_image(i, image2)
  png(paste0("~/SDS365/image_compression/image_approximations/image2_", i, ".png"))
  image(image2_k$M_k, axes = FALSE)
  dev.off()
}

## Image 3
for (i in 1:10){
  image3_k <- compress_image(i, image3)
  png(paste0("~/SDS365/image_compression/image_approximations/image3_", i, ".png"))
  image(image3_k$M_k, axes = FALSE)
  dev.off()
}

## Image 4
for (i in 1:10){
  image4_k <- compress_image(i, image4)
  png(paste0("~/SDS365/image_compression/image_approximations/image4_", i, ".png"))
  image(image4_k$M_k, axes = FALSE)
  dev.off()
}

# Plotting eigenvector and pcscore
## Image 1
image_1_eigenvector <- compress_image(1, image1)$U[,1]
image_1_pc <- compress_image(1, image1)$Z[,1]



# Plotting errors against k
image1_errors <- error_k(image1, ncol(image1))
image2_errors <- error_k(image2, ncol(image2))
image3_errors <- error_k(image3, ncol(image3))
image4_errors <- error_k(image4, ncol(image4))

# Image 1
png("~/SDS365/image_compression/plots/image1_kplot.png")
ggplot(image1_errors, aes(x = k, y = error)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, colour = "red") +
  labs(x = "K",
       y = "Frobenius Error",
       title = "Error vs. K for Image 1")
dev.off()

# Image 2
png("~/SDS365/image_compression/plots/image2_kplot.png")
ggplot(image2_errors, aes(x = k, y = error)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, colour = "red") +
  labs(x = "K",
       y = "Frobenius Error",
       title = "Error vs. K for Image 2")
dev.off()

# Image 3
png("~/SDS365/image_compression/plots/image3_kplot.png")
ggplot(image3_errors, aes(x = k, y = error)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, colour = "red") +
  labs(x = "K",
       y = "Frobenius Error",
       title = "Error vs. K for Image 3")
dev.off()

# Image 4
png("~/SDS365/image_compression/plots/image4_kplot.png")
ggplot(image4_errors, aes(x = k, y = error)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, colour = "red") +
  labs(x = "K",
       y = "Frobenius Error",
       title = "Error vs. K for Image 4")
dev.off()
