#####################################
### Implied volatility

#####################################
## Read-in data
getwd()
setwd("/Users/thijsbenschop/Dropbox/Implied volatility CO2/Data")
data15 <- read.csv("co2_futures15min.csv")
dim(data15)
names(data15)

## Create roll-over data of 
# Roll over on November 30