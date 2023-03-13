
install.packages("fitdistrplus")

library(fitdistrplus)

# set working directory to where the CSV file is located
setwd("/path/to/your/directory")

# read the CSV file into R as a data frame
mydata <- read.csv("filename.csv")

# view the first few rows of the data frame
head(mydata)
