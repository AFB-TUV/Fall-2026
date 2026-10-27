library(broom)
library(car)
library(ggplot2)
library(ISLR)
library(lmtest)
library(VGAM)
library(ggfortify)

help(College)
summary(College)
head(College)

fullmod = lm(Apps ~ ., data = College)
summary(fullmod)

# VIF
vif(fullmod)
#Enroll & F. undergrad severely collinear

diagnosis = lm(Apps ~ Accept + Enroll + Top10perc + Outstate + Expend, data = College)
summary(diagnosis)
vif(diagnosis)

#Standardized residuals
stdresid = rstandard(diagnosis)
which(abs(stdresid) > 2)

# Leverage
leverage = hatvalues(diagnosis)
p = 5
threshold = 3*(p+1)/nrow(College)
threshold
which(leverage > threshold)

#Cook's Distance
cooks = cooks.distance(diagnosis)
head(sort(cooks, decreasing = TRUE))

#Breusch-Pagan test for homoscedasticity
bptest(diagnosis) #rejects H0 => conclude heteroscedasticity

#visual diagnostic plots
autoplot(diagnosis) #only gives 4
autoplot(diagnosis, which = 1) #residuals vs fitted - heteroscedasticity
autoplot(diagnosis, which = 2) #normal QQ plot - very non-normal
autoplot(diagnosis, which = 3) #scale location (fitted vs standardized residuals) - similar to 1st but easier to spot outliers
autoplot(diagnosis, which = 4) #Cook's Distance
autoplot(diagnosis, which = 5) #Residuals vs Leverage - high leverage, outliers (maybe combine into influence), heteroscedasticity
autoplot(diagnosis, which = 6) #Cook's vs leverage - high leverage, high influence
