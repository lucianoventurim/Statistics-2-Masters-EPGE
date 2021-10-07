# Estatística 2 - Luciano Fabio Busatto Venturim

# Loading the data.
data <- read.csv(file = "Estatística 2 - LIsta 4 - Data.csv",
                 header = TRUE, row.names = 1)

# Since there are 5 groups, we will create 4 dummy variables.
data$group1 <- ifelse(data$Groups == 1, 1, 0)
data$group2 <- ifelse(data$Groups == 2, 1, 0)
data$group3 <- ifelse(data$Groups == 3, 1, 0)
data$group4 <- ifelse(data$Groups == 4, 1, 0)

# We define the matrix X and y:
X <- cbind(rep(1,length(data$x)), 
           as.matrix(data[,c('x','group1','group2','group3','group4')]))
colnames(X) <- c('cte','x','group1','group2','group3','group4')
y <- as.matrix(data$y)

# The OLS estimator is:
Beta <- solve((t(X)%*%X),t(X)%*%y)

# The residuals are:
e <- y - X%*%Beta

# Assuming homoskedasticity, the unbiased estimator of the variance of the error
#is sig2:

n <- length(y)
k <- length(Beta)

sig2 <- as.numeric((t(e)%*%e)/(n-k))
xx <- solve(t(X)%*%X)

# The estimator of the variance of the Beta is:
v <- xx*sig2

# The t-statistic for the homoskedastic normal model to test the null hypothesis
#Beta1 = 0 is the following:
t1 <- (Beta[2,1]-0)/sqrt(v[2,2])
print(t1)


# Lets plot the residuals to see if they seem to follow a standard normal distribution:
library(ggplot2)

ggplot(data = data.frame(e)) +
  geom_density(mapping=aes(x = e))

# We can also compare the distribution of y for each group 
#to see if the homoskedastic assumption is plausible.
data[,"Groups"] <- factor(data[,"Groups"], labels=c("1","2","3","4","5"))

library(ggridges)

ggplot(data = data, mapping = aes(x = y, y = Groups, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3) +
  scale_fill_viridis_c(name = "y", option = "C")

# The densities seem to be indeed equally spread around their means.

# The t-statistic is -28.18537 and the critical value for a 5% test is 1.96,
#since n-k = 4994 is big enough for the distribution t with n-k degrees of
#freedom be approximate by a standard normal distribution. Thus
#|t1|>1.96 and we reject the null hypothesis.

# Inverting the pivotal t-statistic, we can create the 95% confidence interval:
l <- Beta[2,1] - 1.96*sqrt(v[2,2])
u <- Beta[2,1] + 1.96*sqrt(v[2,2])

c(l,u)

# Thus, the 95% confidence interval is [-1.0702240, -0,9310558].