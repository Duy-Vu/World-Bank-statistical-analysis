#### Data processing ####
library(foreign)
data_dir <- "vaestotieto.sav"
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

# Summary
summary(dataset.analysis)
#  countries         gdp.per.capita   healthcare.expense life.expectancy infant.mortality    mortality     
# Length:174         Min.   :   293   Min.   : 2.270     Min.   :52.24   Min.   : 1.70   Min.   : 1.169  
# Class :character   1st Qu.:  1862   1st Qu.: 4.466     1st Qu.:67.02   1st Qu.: 6.20   1st Qu.: 5.849  
# Mode  :character   Median :  5413   Median : 6.421     Median :73.75   Median :14.40   Median : 7.236  
#                    Mean   : 13570   Mean   : 6.510     Mean   :72.15   Mean   :22.04   Mean   : 7.700  
#                    3rd Qu.: 15826   3rd Qu.: 8.135     3rd Qu.:77.57   3rd Qu.:34.27   3rd Qu.: 9.200  
#                    Max.   :107627   Max.   :17.004     Max.   :84.10   Max.   :85.90   Max.   :15.500 

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
# [1] 0.8631761

plot(dataset.analysis$gdp.per.capita/1000,
     dataset.analysis$healthcare.expense,
     main="gdp/capita vs healthcare expense: 0.35")
# Weak relationship

plot(dataset.analysis$infant.mortality,
     dataset.analysis$life.expectancy,
     xlab='Infant mortality (per 1,000 live births)',
     ylab='Life expectancy')
# Strongest negative linear relationship
cor(dataset.analysis$infant.mortality,
    dataset.analysis$life.expectancy)

##### Sampling ####
n <- 35
set.seed(n)
random_id <- sample(1:174, n)
dataset.sample <- dataset.analysis[random_id,]
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
n=3
qs <- vector(length=n)
for (i in 1:n) {
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
# [1] 3

# Smallest Q
qs[which.min(qs)]
qs



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
legend(-1, .3, r.obt, bty = "n")
legend(0.7, .3, -r.obt, bty = "n")
r.obt
prob <- sum(r.random[-r.obt <= r.random] + r.random[r.random >= r.obt])/nreps
prob



#### Estimate correlation coefficients using Bootstrap and Jackknife method ####  
set.seed(1000)
data <- cbind(X,y)
# Jackknife
theta1 <- function(x,m){
  cor(m[x,1],m[x,2])
}
l1 <- jackknife(1:10, theta1, data)
# Pseudo-values and the Jackknife-estimator
bjack<-10*cor(x,y)-9*l1$jack.values
theta.est <- sum(bjack) / 10
theta.est
## [1] 0.9272874
# Confidence interval
t.value <- qt(0.975, 9)
dev <- t.value*sqrt(var(bjack) / 10)
lower <- theta.est - dev
upper <- theta.est + dev
CI <- c(lower, upper)
CI
## [1] 0.8704885 0.9840864
4# Bootstrap
theta2 <- function(m,x) {
  cor(m[x,1], m[x,2])
}
l2 <- boot(data, theta2, R=1000)
l2$t0
## [1] 0.9328462
boot.ci(l2, type = "perc")



#### 95% Confidence interval ####























