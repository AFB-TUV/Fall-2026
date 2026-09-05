#loading packages

library(ggplot2)

myclass = "Actuarial Analytics"
# <- is also the assignment operator
myclass

enrollment = 20
enrollment + 5
enrollment 

bool = TRUE
bool

?length
library(dslabs)
data("murders")

numbers = 1:10
numbers

length(numbers)
select = numbers > 5
select
#functions use ()
#accessing elements uses []

numbers[select]
murders[1,2] #row-column notation

murders$population # datasetname#columnname
ifelse(numbers > 5, "Big", "Small")
