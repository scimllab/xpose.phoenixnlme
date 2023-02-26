#' Covariance parameterization
#' @param n number of samples to draw
#' @param mu p-vector mean
#' @param Sigma covariance matrix (p x p)
#' @return matrix of dimension p x n of samples
rMVNormC <- function(n, mu, Sigma){
  p <- length(mu)
  Z <- matrix(rnorm(p*n), p, n)
  L <- t(chol(Sigma)) # By default R's chol fxn returns upper cholesky factor
  X <- L%*%Z
  X <- sweep(X, 1, mu, FUN=`+`)
  return(X)
}

#' Precision parameterization
#' @param n number of samples to draw
#' @param mu p-vector mean
#' @param Omega precision matrix (p x p)
#' @return matrix of dimension p x n of samples
rMVNormP <- function(n, mu, Omega){
  p <- length(mu)
  Z <- matrix(rnorm(p*n), p, n)
  U <- chol(Omega) # By default R's chol fxn returns upper cholesky factor # A= R* R'
  X <- backsolve(U, Z) # more efficient and stable than acctually inverting
  X <- sweep(X, 1, mu, FUN=`+`)
  return(X)
}

# http://www.statsathome.com/2018/10/19/sampling-from-multivariate-normal-precision-and-covariance-parameterizations/

# Cholesky decomposition means factorizes symmetric positive definite matrix A into an upper triangular R that satisfies A = R'*R. 
# If A is nonsymmetric , then chol treats the matrix as symmetric and uses only the diagonal and upper triangle of A.
R_S=chol(Sigma)
Equal_Omega<-chol2inv(chol(Sigma))

Equal_Sigma<-chol2inv(chol(Omega))


set.seed(153)
rWishart(2, 4, diag(4))

n <- 150
mu <- 1:4
Sigma <- rWishart(1, 10, diag(4))[,,1] # random covariance matrix 
# https://en.wikipedia.org/wiki/Wishart_distribution
# it is probability distribution of pxp or var-cov matrix (XX') S=GG' S ~ Wp(df,V); V is variance
Omega <- solve(Sigma)
x1 <- rMVNormC(n, mu, Sigma)
x2 <- rMVNormP(n, mu, Omega)

# My example
n <- 15
mu <- log(c(2.6,42,22,57))
Sigma <- rWishart(1, 10, diag(4))[,,1]
Omega <- solve(Sigma)
x2 <- rMVNormP(n, mu, Omega)

