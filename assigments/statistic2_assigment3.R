# Loading ggplot2 and ggridges.
library(ggplot2)
#install.packages("ggridges")
library(ggridges)

# Loading the data.
data <- read.csv(file = "https://raw.githubusercontent.com/lucianoventurim/Statistics-2-Masters-EPGE/main/data/statistic2_assigment3_data.csv",
                 header = TRUE, row.names = 1)

# Question 5 item a. Lets plot the data for each group and a regression line 
#for each:
data[,"Groups"] <- factor(data[,"Groups"], labels=c("1","2","3","4","5"))

ggplot(data = data, mapping = aes(x = x, y = y, fill = Groups)) +
  geom_point(shape = 21, size = 1.5) +
  geom_smooth(color = "black", method = lm, formula = 'y ~ x', se = TRUE)

# The dispersion of the dependent variable y on group 5 seems to be higher, 
#which might indicate that the data is not homoskedastic.

# Question 5 item b. Now, we create a ridge density plot:
ggplot(data = data, mapping = aes(x = y, y = Groups, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3) +
  scale_fill_viridis_c(name = "y", option = "C")

# We can see that the mean and variance within the groups increase from group 1
#to group 5.

# Question 5 item c. Now we want to create dummy variables for the groups. Since
# there are 5 groups, we will create 4 dummy variables.
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
b1 <- Beta[2]

# The residuals are:
e <- y - X%*%Beta

# Now, we calculate the estimators of the variance of Beta:
n <- length(y)
k <- length(Beta)
a <- n/(n-k)
sig2 <- as.numeric((t(e)%*%e)/(n-k))
xx <- solve(t(X)%*%X)
M <- diag(n)-X%*%(xx)%*%t(X)
h <- diag(M)
hsqrt <- sqrt(h)

# The estimators are of the form Vl = (X'X)^(-1)(X'DlX)(X'X)^(-1), where Dl is a
#diagonal matrix where the diagonal entries are the estimators of ei^(2). Each
#matrix Dl, l=0,2,3, considers different estimators. In order to simplify the 
#calculation, we make X'DlX = Xl'Xl where Xl is construct to satisfy this
#equation. Since Dl is a diagonal matrix, Xl is construct by multiplying the ith
#row of X by êi.
x0 <- X*(e%*%matrix(1,1,k))
x2 <- X*((e/hsqrt)%*%matrix(1,1,k))
x3 <- X*((e/h)%*%matrix(1,1,k))

# The variance matrices are:
hom <- xx*sig2
hc0 <- xx%*%(t(x0)%*%x0)%*%xx
hc1 <- a*(xx%*%(t(x0)%*%x0)%*%xx)
hc2 <- xx%*%(t(x2)%*%x2)%*%xx
hc3 <- xx%*%(t(x3)%*%x3)%*%xx

# The standard errors are:
se_hom <- sqrt(diag(hom))
se_hc0 <- sqrt(diag(hc0))
se_hc1 <- sqrt(diag(hc1))
se_hc2 <- sqrt(diag(hc2))
se_hc3 <- sqrt(diag(hc3))

table <- matrix(c(Beta,se_hom,se_hc0,se_hc1,se_hc2,se_hc3), 6,6)
colnames(table) <- c('estimates', 'hom', 'hco', 'hc1', 'hc2', 'hc3')
rownames(table) <- c('Beta0', 'Beta1', 'Delta1', 'Delta2', 'Delta3', 'Delta4')

table <- as.table(table)
table

# Question 5 item d. As expected, hc0<hc2<hc3. Moreover, the standard error is 
#larger under the Homoskedasticity assumption, which is expected since the 
#errors seem to be heteroskedastic.