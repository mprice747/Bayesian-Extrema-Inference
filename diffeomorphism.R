library(pracma)
library(microbenchmark)

SQRT_2 <- sqrt(2)


check_matrix <- function(X) {
  
  if (!is.matrix(X)) {
    
    X <- matrix(X, ncol = 1)
  }
  
  return(X)
}

turn_1d_into_matrix <- function(x) {
  
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1)
    return(x)
  }
  return(x)
}


cosine_basis <- function(X, Beta) {
  
  # X - nx1 Matrix of observations 
  # Beta - mxj Matrix of different beta vectors
  
  # Output- nxm Matrix of representing the n observations evaluated at m 
  # different cosine bases
  
  Cos_Mat <- cos(pi * (X %*% matrix(seq(1, ncol(Beta)), nrow = 1)))
  
  result <- tcrossprod(Cos_Mat, (SQRT_2 * Beta))
  
  return(result)
}

int_cosine_basis <- function(X, Beta) {
  
  j_vec <- seq(1, ncol(Beta))
  
  Sin_Mat <- sin(pi * (X %*% matrix(j_vec, nrow = 1)))
  
  result <- Sin_Mat %*% (t((SQRT_2 * Beta))/(j_vec * pi))
  
  return(result)
}


cos_2_integral <- function(x, j) {
  
  # Integral of cos^2(x * pi * j)
  
  two_pi_j <- 2 * pi * j 
  result <- sin(two_pi_j * x)/(2 * two_pi_j) + x/2
  
  return(result)
  
}

cos_cos_integral <- function(x, j1, j2){
  
  # Integral of cos(x * pi * j1)cos(x * pi * j2) where j1 != j2
  
  pi_j1_j2 <- pi * (j1 + j2)
  pi_j2_sub_j1 <- pi * (j2 - j1)
  
  result <- sin(pi_j1_j2 * x)/(2 * pi_j1_j2) +
    sin(pi_j2_sub_j1 * x)/(2 * pi_j2_sub_j1)
  
  return(result)
}

different_beta_product <- function(beta) {
  
  beta_comb <- combn(beta, 2)
  
  return(beta_comb[1, ] * beta_comb[2, ])
}


int_cosine_basis_2 <- function(X, Beta) {
  
  
  j_seq <- seq(1, ncol(Beta))
  j_comb <- combn(j_seq, 2)
  
  
  cos_cos_ints <- turn_1d_into_matrix(apply(X, 1, cos_cos_integral, 
                         j1 = j_comb[1, ], j2 = j_comb[2, ]))
  
  diff_betas_mult <- turn_1d_into_matrix(apply(Beta, 1, different_beta_product))
  
  cos_cos_all <- 4 * crossprod(cos_cos_ints, diff_betas_mult)
  
  cos_2_ints <- turn_1d_into_matrix(apply(X, 1, cos_2_integral, j = j_seq))
  
  cos_2_all <- 2 * t(Beta^2 %*% cos_2_ints)
  
  result <- cos_cos_all + cos_2_all
  
  return(result)
}


gamma_function <- function(X, Beta) {
  
  Beta_norms <- sqrt(rowSums(Beta^2))
  
  int_1 <- t(int_cosine_basis(X, Beta))
  
  int_2 <- t(int_cosine_basis_2(X, Beta))
  
  cos_beta <- cos(Beta_norms)
  sin_beta_div <- sin(Beta_norms)/Beta_norms
  
  cos_beta_X <- outer(cos(Beta_norms)^2, X[, 1])
  
  exp_map <- cos_beta_X + ((2 * cos_beta * 
                   sin_beta_div) * int_1) + sin_beta_div^2 * int_2
  
  return(t(exp_map))
 
}

