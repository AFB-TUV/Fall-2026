install.packages("tidyverse")
install.packages("ggplot2")

library(tidyverse)
library(ggplot2)

data(mtcars)
?mtcars

#overall syntax
#ggplot(data = , aes(x=, y=, col=, fill=, shape=)) + geom_point() + ...

ggplot(data = mtcars, aes(x = cyl)) + geom_histogram()
ggplot(data = mtcars, aes(x = cyl)) + geom_bar()
ggplot(data = mtcars, aes(x = cyl)) + geom_boxplot()

ggplot(data = mtcars, aes(x = cyl, y = mpg)) + geom_point()
ggplot(data = mtcars, aes(x = cyl, y = mpg)) + geom_line()

ggplot(data = mtcars, aes(x = cyl, y = mpg)) + geom_point() + geom_smooth()
ggplot(data = mtcars, aes(x = cyl, y = mpg, color = disp)) + geom_point() + geom_smooth(method = "lm")

ggplot(data = mtcars, aes(x = hp, y = mpg, color = factor(cyl))) + geom_point() + geom_smooth(method = "lm")

ggplot(data = mtcars, aes(x = hp, y = mpg, color = factor(am))) + geom_point() + geom_smooth(method = "lm")
# want to label 0,1 to something descriptive

mtcars$am.fact = factor(mtcars$am, levels = c(0,1), labels = c("Automatic", "Manual"))
ggplot(data = mtcars, aes(x = hp, y = mpg, color = factor(am.fact))) + geom_point() + geom_smooth(method = "lm")

ggplot(data = mtcars, aes(x = mpg, fill = am.fact)) + geom_histogram()
#filled each bar w/ different color by category

ggplot(data = mtcars, aes(x = mpg, color = am.fact)) + geom_histogram()
#outside coloring instead

ggplot(data = mtcars, aes(x = mpg, fill = am.fact)) + geom_bar()
#bar splits auto vs. manual and doesn't stack

savep = ggplot(data = mtcars, aes(x = hp, y = mpg, color = factor(am.fact))) + geom_point() + geom_smooth(method = "lm")
savep
ggsave(filename = "scatter.png", plot = savep, width = 6, height = 4, dpi = 300)
wd()
