
source('diffeomorphism.R')


linear_interp_area <- function(X, Y){
  
  n <- length(X)
  
  base1 <- Y[1:(n - 1)]
  base2 <- Y[2:n]
  height <- X[2:n] - X[1:(n - 1)]
  
  areas <- (base1 + base2)/2 * height
  int_answer <- sum(areas)
  
  return(int_answer)
}

