#Estatística 2 - Lista 8 - Luciano Fabio Venturim

# First, we read the data.
data <- read.table("Estatística 2 - Lista 8 - Data.txt", header = TRUE, sep = ",")
colnames(data) <- c('z', 'choice')

n <- nrow(data)


# Lets define Qn, the log likelihood criterion function.

qn <- function(lambdas) {
  q <- numeric(n)
  for (i in 1:n) {
    qi <- (1-data[i,'choice'])*(lambdas[1]+lambdas[2]*data[i,'z'])-
      log(1+exp(lambdas[1]+lambdas[2]*data[i,'z']))
    q[i] <- qi
  }
  return(-mean(q))
}

# To maximize Qn, we will use the function optim.

optmizedqn <- optim(par = c(1,1), fn = qn, method = "Nelder-Mead")

lambda_hat <- optmizedqn$par


# The MM estimate of -H^{-1} is as below.

h <- matrix(c(0,0,0,0),nrow = 2,ncol = 2)

h11 <- numeric(n)
h12 <- numeric(n)
h21 <- numeric(n)
h22 <- numeric(n)
  
for (i in 1:n) {
  h11[i] <-  -exp(lambda_hat[1]+lambda_hat[2]*data[i,'z'])/(1+exp(lambda_hat[1]+lambda_hat[2]*data[i,'z']))^2
  h12[i] <-  -data[i,'z']*exp(lambda_hat[1]+lambda_hat[2]*data[i,'z'])/(1+exp(lambda_hat[1]+lambda_hat[2]*data[i,'z']))^2
  h21[i] <-  -data[i,'z']*exp(lambda_hat[1]+lambda_hat[2]*data[i,'z'])/(1+exp(lambda_hat[1]+lambda_hat[2]*data[i,'z']))^2
  h22[i] <-  -(data[i,'z'])^2*exp(lambda_hat[1]+lambda_hat[2]*data[i,'z'])/(1+exp(lambda_hat[1]+lambda_hat[2]*data[i,'z']))^2
}

h[1,1] <- mean(h11)
h[1,2] <- mean(h12)
h[2,1] <- mean(h21)
h[2,2] <- mean(h22)

var_asy <- -solve(h)

# The variances of lambda_hat are approximately (-h)^{-1}/n

var_lambdas <- var_asy/n

#and the standard errors are

std_errors_lambdas <- sqrt(diag(var_lambdas))

# Lets construct the asymptotic variance of P.

grad_g <- c(-exp(lambda_hat[1]+11/3*lambda_hat[2])/(1+exp(lambda_hat[1]+11/3*lambda_hat[2]))^2,
            -11/3*exp(lambda_hat[1]+11/3*lambda_hat[2])/(1+exp(lambda_hat[1]+11/3*lambda_hat[2]))^2)

std_errors_p <- sqrt(t(grad_g)%*%var_lambdas%*%grad_g)

p <- 1/(1+exp(lambda_hat[1]+11/3*lambda_hat[2]))

ci1 <- p - 1.96*std_errors_p
ci2 <- p + 1.96*std_errors_p

# So the CI is [0.258,0.427]
