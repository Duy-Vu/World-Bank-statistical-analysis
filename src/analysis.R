#### Data processing ####
library(foreign)
data_dir <- "../data/vaestotieto.sav"
path <- file.path(data_dir)
dataset <- read.spss(path, to.data.frame=TRUE)

# Remove data of regions, or groups of country (Valtioiden muodostama alue)  
dataset.analysis <- subset(dataset, V3 != "Valtioiden muodostama alue")
summary(dataset.analysis)

# Explanatory variables: 
# 1. V20-GDP per capita, 
# 2. V23 healthcare expense, 
gdp.per.capita <- dataset.analysis$V20
healthcare.expense <- dataset.analysis$V23

# Dependent variables: 
# 1. V27-life expectancy
# 2. V14 infant mortality
# 3. V13-mortality
life.expectancy <- dataset.analysis$V27
infant.mortality <- dataset.analysis$V14
mortality <- dataset.analysis$V13

# Create new dataframe to analyse and remove rows with NA value
countries <- dataset.analysis$V1
dataset.analysis <-na.omit(data.frame(countries,
                                      gdp.per.capita, 
                                      healthcare.expense,
                                      life.expectancy,
                                      infant.mortality,
                                      mortality,
                                      stringsAsFactors=FALSE))
row.names(dataset.analysis) <- NULL
write.csv(dataset.analysis, file="../output/data/clean_data.csv", row.names=FALSE)
summary(dataset.analysis)



#### Exploratory data analysis ####
# Calculate correlation coefficients between all variables
# Get correlation matrix between 5 columns 
cor_matrix <- cor(dataset.analysis[,2:6])

library(corrplot)
corrplot(cor_matrix, method="number", type="lower", order="hclust")

plot(log(dataset.analysis$gdp.per.capita/1000),
     dataset.analysis$life.expectancy,
     main="gdp/capita vs life expectancy: 0.65")
# Strong relationship between life expectancy and log of gdp/capita
cor(dataset.analysis$life.expectancy, 
    dataset.analysis$gdp.per.capita, 
    method="spearman")

plot(dataset.analysis$gdp.per.capita/1000,
     dataset.analysis$healthcare.expense,
     main="gdp/capita vs healthcare expense: 0.35")
# Weak relationship

# Strongest negative linear relationship
plot(dataset.analysis$infant.mortality,
     dataset.analysis$life.expectancy,
     main="infant mortality vs life expectancy: -0.93")



#### Sampling ####
n <- 35
set.seed(n)
random_id <- sample(1:174, n)
dataset.sample <- dataset.analysis[random_id,]
write.csv(dataset.sample, file="../output/data/sample_data.csv", row.names=FALSE)
summary(dataset.sample)

X <- dataset.sample$infant.mortality
Y <- dataset.sample$life.expectancy
plot(X, Y,
     xlab='Infant mortality (per 1,000 live births)',
     ylab='Life expectancy', 
     pch = 16, col = "blue")


#### Fit linear line to the data ####
# Without R
X_1 <- cbind(1,X)
BM<-solve(t(X_1)%*%X_1)%*%t(X_1)%*%Y
BM

# With R
lmTemp = lm(Y~X)
summary(lmTemp)
abline(lmTemp)  # Draw the fit line


#### Cross-validation OLS ####
library(bootstrap)
# The fit-function for crossval(), n is the degree of the model
theta.fitn <- function (x,y,n){
  a <- cbind(x)
  if (n>1) {
    for (i in 2:n) {
      a <- cbind(a, I(x^i))
    }
  }
  lsfit(a, y)
}
# The predict-function for crossval()
theta.predictn<- function (fit,x){
  a <- cbind(1, x)
  n <- length(fit$coef) - 1
  if (n>1) {
    for (i in 2:n) {
      a <- cbind(a, I(x^i))
    }
  }
  a%*%fit$coef
}

# Go through models up to the 3rd degree
deg <- 3
qs <- vector(length=deg)
for (i in 1:deg) {
  results<-crossval(dataset.sample$infant.mortality,
                    dataset.sample$life.expectancy, 
                    theta.fitn, theta.predictn, n=i)
  # Mean of the loss-function
  Q <- sum((dataset.sample$life.expectancy-results$cv.fit)^2)/52
  qs[i] <- Q
}
# Choose the degree with the smallest mean of the loss-function
# Degree
which.min(qs)

# Smallest Q
qs
qs[which.min(qs)]



#### Hypothesis testing ####
set.seed(n)
r.obt <- cor(X, Y)
cat("The obtained correlation is ", r.obt,'\n')
nreps <- 1e4
r.random <- numeric(nreps)
for (i in 1:nreps) {
  X.permute <- sample(X, length(X), replace=FALSE)
  r.random[i] <- cor(X.permute, Y)
}
plot(density(r.random),
     main = "Distribution of correlation coefficient r", 
     xlab = "r from randomized samples",
     xlim = c(-1,1))
r.obt <- round(r.obt, digits = 2)
abline(v=r.obt)
abline(v=-r.obt)
legend(-1, .2, r.obt, bty = "n")
legend(0.66, .2, -r.obt, bty = "n")
r.obt

prob <- sum(r.random[-r.obt <= r.random] + r.random[r.random >= r.obt])/nreps
prob



#### Estimate correlation coefficients using Bootstrap and Jackknife method ####
set.seed(n)
library(bootstrap)
library(boot)

data <- cbind(X,Y)
# Jackknife
theta1 <- function(x,m){
  cor(m[x,1], m[x,2])
}
l1 <- jackknife(1:n, theta1, data)
# Pseudo-values and the Jackknife-estimator
bjack <- n*cor(X,Y) - (n-1)*l1$jack.values 
theta.est <- sum(bjack) / n
theta.est

# Confidence interval
t.value <- qt(0.975, n-1)
dev <- t.value*sqrt(var(bjack) / n)
lower <- theta.est - dev
upper <- theta.est + dev
CI <- c(lower, upper)
CI


# Bootstrap
theta2 <- function(m,x) {
  cor(m[x,1], m[x,2])
}
l2 <- boot(data, theta2, R=1000)
l2$t0

boot.ci(l2, type = "perc")
