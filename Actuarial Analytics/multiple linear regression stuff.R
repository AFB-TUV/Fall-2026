#Multiple Linear Regression R Demo

library(ISLR)
library(broom)
data("Carseats")
summary(Carseats)

mlr = lm(Sales ~ ., data = Carseats)
summary(mlr)
confint(mlr, level = 0.95)

mlr2 = lm(Sales ~ CompPrice + Income + Advertising + Price + ShelveLoc + Age, data = Carseats)
summary(mlr2)

glance(mlr2) #model level
tidy(mlr2) # summary coefficients
augment(mlr2) # observation-level data

#what hurtw Sales? bad shelf location, price, age (higher is worse for both)
?Carseats

# F-test
ftest = anova(mlr, mlr2)
ftest
#fail to reject H0; more complex not significantly better

