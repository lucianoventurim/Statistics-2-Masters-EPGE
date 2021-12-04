# Loading ggplot2.
library(ggplot2)

# Loading the data.
data <- read.csv(file = "https://raw.githubusercontent.com/lucianoventurim/Statistics-2-Masters-EPGE/main/data/statistic2_assignment2_data.csv",
                 header = TRUE, row.names = 1)

# Now we want to create dummy variables for the groups. Since there are 5 groups,
#we will create 4 dummy variables.
data$group1 <- ifelse(data$Groups == 1, 1, 0)
data$group2 <- ifelse(data$Groups == 2, 1, 0)
data$group3 <- ifelse(data$Groups == 3, 1, 0)
data$group4 <- ifelse(data$Groups == 4, 1, 0)

# Question 3, item a. We will first do a regression of y on x and a constant. First,
#we will create our matrices y and X = [1 x]'
y <- as.matrix(data$y)
X_a <- as.matrix(data$x)
X_a <- cbind(rep(1,length(data$x)),X_a)

colnames(y) <- c('y')
colnames(X_a) <- c('cte','x')

# Thus, we have that the OLS estimator is B such that (X'X)B = X'y and the residuals
#are e = y - XB. We then solve for this system of linear equations for B (Beta_a).
A_a <- t(X_a) %*% X_a
c_a <- t(X_a) %*% y
Beta_a <- solve(A_a,c_a)
print(Beta_a)

# The estimator is B = [7.29 1.55]. The residuals e = y - XB are:
e_a = y - X_a%*%Beta_a

# The SSE is e'e and equals 1244385 and the R2 is 1-SSE/TSS equals 0.6197:
SSE_a = t(e_a) %*% e_a
TSS_a = t(y-rep(1,length(y))*mean(y)) %*% (y-rep(1,length(y))*mean(y))
R2_a = 1-(SSE_a/TSS_a)
print(R2_a)

# Now, we use the function lm from R to compare the results, and a simple inspection
# show that the values are the same in both cases:
model_a <- lm(data$y ~ data$x)

model_coeff_a <- as.matrix(model_a$coefficients)
print(model_coeff_a)
model_residuals_a <- as.matrix(model_a$residuals)
model_R2_a <- as.matrix(summary(model_a)$r.squared)
print(model_R2_a)

# We plot the data and the regression line:
x_seq = seq(min(data$x), max(data$x), length.out = 5000)
fitted_values_a = Beta_a[1,1] + Beta_a[2,1]*x_seq
manual_reg_line_a <- data.frame('x' = x_seq, 'fitted' = fitted_values_a)

ggplot(data = data) +
  geom_point(mapping = aes(x = x, y = y), color = 'gray') +
  geom_line(data = manual_reg_line_a, mapping = aes(x = x, y = fitted)) +
  geom_smooth(mapping = aes(x = x, y = y), method = lm, formula = 'y ~ x',
              color = 'red')


# Question 3, item b. To regress y on x and the dummy variables, we will
#redefine X and Beta to include the other variables and do the regression just
#like before.
X_b <- cbind(X_a,as.matrix(data[,c('group1','group2','group3','group4')]))
A_b <- t(X_b) %*% X_b
c_b <- t(X_b) %*% y
Beta_b <- solve(A_b,c_b)
print(Beta_b)

# See that Beta_b = [124.85 -1.00 -99.92 -74.98 -49.69 - 25.35], thus the
#coefficient of xi changed. Again, the new residuals, R2, and SSE are:
e_b = y - X_b%*%Beta_b
SSE_b = t(e_b) %*% e_b
TSS_b = t(y-rep(1,length(y))*mean(y)) %*% (y-rep(1,length(y))*mean(y))
R2_b = 1-(SSE_b/TSS_b)
print(R2_b)

#i.e, R2_b = 0.847 and SSE_b = 500360.3. The results from using the lm function
#are:
model_b <- lm(data$y ~ data$x + data$group1 + data$group2 + data$group3 +
                data$group4)

model_coeff_b <- as.matrix(model_b$coefficients)
print(model_coeff_b)
model_residuals_b <- as.matrix(model_b$residuals)
model_R2_b <- as.matrix(summary(model_b)$r.squared)
print(model_R2_b)

# As we can see, the values are the same.

# Question 3, item c. In the FWL approach we partition X = [X1 X2], where X1 
# (X1_c) is x and X2 (X2_c) is the matrix with the dummy variables for the groups and constant.
# We first regress y on X2:
X2_c = X_b[, c('cte', 'group1', 'group2', 'group3', 'group4')]
X1_c = as.matrix(X_b[, 'x'])
colnames(X1_c) <- 'x'

A1_c <- t(X2_c) %*% X2_c
c1_c <- t(X2_c) %*% y
Beta1_c <- solve(A1_c,c1_c)
print(Beta1_c)

# The residualised y, y_tilda is:
y_tilda = y - X2_c%*%Beta1_c

# Then, we regress X on X2:
c2_c <- t(X2_c) %*% X1_c
Beta2_c <- solve(A1_c,c2_c)
print(Beta2_c)

# The residualised x, x_tilda is:
x_tilda = X1_c - X2_c%*%Beta2_c

# Finally, we do the regression of y_tilda on x_tilda:
A_c <- t(x_tilda) %*% x_tilda
c_c <- t(x_tilda) %*% y_tilda
Beta_c <- solve(A_c,c_c)
print(Beta_c)

# As we would expect, it yields the same estimator of the coefficient of x as
#in the OLS model.

# Now, we plot y_tilda and x_tilda and the regression line:
residualized_variables = data.frame('x_tilda' = x_tilda, 'y_tilda' = y_tilda)
x_seq_c = seq(min(x_tilda), max(x_tilda), length.out = 5000)
fitted_values_c = Beta_c[1]*x_seq_c
manual_reg_line_c <- data.frame('x' = x_seq_c, 'fitted' = fitted_values_c)


ggplot(data = residualized_variables) +
  geom_point(mapping = aes(x = x_tilda, y = y_tilda), color = 'gray') +
  geom_line(data = manual_reg_line_c, mapping = aes(x = x, y = fitted))


# Question 3, item d. Lets plot the data again highlighting the groups and
#the regression lines for each one:
data[,"Groups"] <- factor(data[,"Groups"], labels=c("1","2","3","4","5"))

ggplot(data = data, mapping = aes(x = x, y = y, fill = Groups)) +
  geom_point(shape = 21, size = 1.5) +
  geom_smooth(color = "black", method = lm, formula = 'y ~ x', se = FALSE)
