# Estatística 2 - Lista 7
# Luciano Fabio Busatto Venturim

# First, we will generate the data. w1 and w2 are the estimates for the simulated
#Wald statistics
set.seed(1)

B <- 1000
n <- 100
beta0 <- c(1,1,1)

w1_simulated <- numeric(B)
w2_simulated <- numeric(B)

for (i in 1:B) {
  x <- matrix(runif(3*n), nrow = n, ncol = 3)
  colnames(x) <- c('x1','x2','x3')
  v <- rnorm(n)
  y <- x %*% beta0 + v

  
  #Linear projection coefficients
  beta <- solve(t(x)%*%x,t(x)%*%y)
  e <- y - x%*%beta
  
  #Asymptotic variance estimator
  xx_inv <- solve(t(x)%*%x)
  xe <- x*(e%*%matrix(1,1,3))
  sigma_hat <- t(xe)%*%xe
  var_beta <- xx_inv%*%sigma_hat%*%xx_inv
  
  #Wald test statistics
  r1 <- matrix(c(1,-1,0,0,1,-1), nrow = 3, ncol = 2);
  r2 <- matrix(c(exp(beta[1]),-exp(beta[2]),0,0,exp(beta[2]),-exp(beta[3])),
               nrow = 3, ncol = 2)
  r_beta <- matrix(c(exp(beta[1])-exp(beta[2]),exp(beta[2])-exp(beta[3])),
                   nrow = 2, ncol = 1)
  w1 <- t(beta)%*%r1%*%solve(t(r1)%*%var_beta%*%r1)%*%t(r1)%*%beta
  w2 <- t(r_beta)%*%solve(t(r2)%*%var_beta%*%r2)%*%r_beta
  
  w1_simulated[i] <- w1
  w2_simulated[i] <- w2
}

# Item a). For both cases, the test statistic converge in distribution to a
#chi-squared distribution with degree 2, the number of restrictions. The critical
#value c and the estimates p1 and p2 are given below.

c <- qchisq(p = 0.1,df = 2, lower.tail = FALSE)

p1 <- sum(w1_simulated>c)/B
p2 <- sum(w2_simulated>c)/B

