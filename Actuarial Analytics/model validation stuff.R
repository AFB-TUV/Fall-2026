library(ISLR)
library(MASS)
library(leaps)

data(Credit)
head(Credit)
help(Credit)

cred = Credit
cred$ID = NULL

#Boundary Models

nullmod = lm(Balance ~ 1, data = cred)
summary(nullmod)

allmod = lm(Balance ~ ., data = cred)
summary(allmod)

#Backward Selection

aic.back = stepAIC(allmod, direction = "backward")
summary(aic.back)
#Balance ~ Income + Limit + Rating + Cards + Age + Student
AIC(aic.back)

aic.fwd = stepAIC(nullmod, direction = "forward", scope = list(upper = allmod, lower = nullmod))
summary(aic.fwd)
#Balance ~ Rating + Income + Student + Limit + Cards + Age

#Mixed stepwise
aic.mixed = stepAIC(nullmod, direction = "both", scope = list(upper = allmod, lower = nullmod))
summary(aic.mixed)

# Balance ~ Rating + Income + Student + Limit + Cards + Age

#use BIC as criteria
bic.mixed = stepAIC(nullmod, direction = "both", k = log(nrow(cred)), scope = list(upper = allmod, lower = nullmod))
summary(bic.mixed)
BIC(bic.mixed)

# Balance ~ Income + Student + Limit + Cards


# Best subset

bestlin = regsubsets(Balance ~ ., data = cred, method = "exhaustive", nvmax = 11)
#nvmax: max # of predictors (defaults to 9)
bestlinsumm = summary(bestlin)
bestlinsumm

#does not bother calculating metric for null model, so 1st entry = 1 predictor model
bestlinsumm$adjr2
which.max(bestlinsumm$adjr2)
coef(bestlin, 7)

which.min(bestlinsumm$cp)
coef(bestlin, 6)

which.min(bestlinsumm$bic)
coef(bestlin, 4)

plot(bestlin, scale = "adjr2")
plot(bestlin, scale = "Cp")
plot(bestlin, scale = "bic")

