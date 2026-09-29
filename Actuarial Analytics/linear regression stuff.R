install.packages("tidyverse")
install.packages("broom")
install.packages("ggplot2")

library(tidyverse)
library(broom)
library(ggplot2)

data("anscombe")
help(anscombe)

summary(anscombe)

#run function through every column in dataset
sapply(anscombe, FUN = mean)

#create linear models
lm1 = lm(y1 ~ x1, data = anscombe)
summary(lm1)
lm1$coefficients

lm2 = lm(y2 ~ x2, data = anscombe)
summary(lm1)
lm2$coefficients

lm3 = lm(y3 ~ x3, data = anscombe)
summary(lm1)
lm3$coefficients

lm4 = lm(y4 ~ x4, data = anscombe)
summary(lm1)
lm4$coefficients

#individual model statistics
glance(lm1)
glance(lm1)$r.squared
glance(lm1)$sigma #residual standard error (s)

glance(lm2)$r.squared
glance(lm2)$sigma #residual standard error (s)

glance(lm3)$r.squared
glance(lm3)$sigma #residual standard error (s)

glance(lm4)$r.squared
glance(lm4)$sigma #residual standard error (s)

#correlation function
cor(anscombe)
cor(anscombe$x1, anscombe$y1) # positive correlation, expect b1 > 0
cor(anscombe$x2, anscombe$y2)
cor(anscombe$x3, anscombe$y3)
cor(anscombe$x4, anscombe$y4)

# Predictions

predict(lm1, newdata = data.frame(x1 = 10))
predict(lm2, newdata = data.frame(x2 = 10))
predict(lm3, newdata = data.frame(x3 = 10))
predict(lm4, newdata = data.frame(x4 = 10))

#Prediction Interval
predict(lm1, newdata = data.frame(x1 = 10), interval = "prediction")


# Confidence Intervals

confint(lm1, level = 0.95)
confint(lm2, level = 0.95)
confint(lm3, level = 0.95)
confint(lm4, level = 0.95)
predict(lm1, newdata = data.frame(x1 = 10), interval = "confidence")


# Visualize
ggplot(data = anscombe, aes(x = x1, y = y1)) + geom_point() #ideal for linear regress9ion
ggplot(data = anscombe, aes(x = x2, y = y2)) + geom_point() #non-linear! (bad)
ggplot(data = anscombe, aes(x = x3, y = y3)) + geom_point() # maybe; outlier
ggplot(data = anscombe, aes(x = x4, y = y4)) + geom_point() #maybe? "high leverage point" to investigate (drives whole model)
