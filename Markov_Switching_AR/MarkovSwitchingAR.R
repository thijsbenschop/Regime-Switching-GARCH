## Clear history
rm(list = ls(all = TRUE))
graphics.off()

## Install and load packages
libraries = c("DEoptim")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Set working directory
setwd("")

# Estimation of Markov Switching AR(4) model

# Function to evaluate the expectation stage of the EM algorithm
condLogLikPar <- function(theta, series) {
  # calculates the conditional loglikelihood for a AR Markov Switching model with 2
  # states arguments - p11 transition probability to stay in state 1 - p22
  # transition probability to stay in state 2 - c1, c2, constant term for AR
  # process in resp. state 1 and 2 - phi - sigma, standard deviation of Gaussian
  # disturbances returns : value of conditional log-likelihood
  p11    <- theta[1]
  p22    <- theta[2]
  c1     <- theta[3]
  c2     <- theta[4]
  phi11  <- theta[5]
  phi12  <- theta[6]
  phi21  <- theta[7]
  phi22  <- theta[8]
  phi31  <- theta[9]
  phi32  <- theta[10]
  phi41  <- theta[11]
  phi42  <- theta[12]
  sigma1 <- theta[13]
  sigma2 <- theta[14]
  n     <- length(series)  # number of obs.
  p     <- cbind(c(p11, 1 - p22), c(1 - p11, p22))  # transition matrix
  dzeta <- cbind(rep(0, n), rep(0, n))
  f     <- cbind(rep(0, n))
  eta   <- cbind(rep(0, n), rep(0, n))
  
  dzetaInit <- c((1 - p[2, 2])/(2 - p[1, 1] - p[2, 2]), (1 - p[1, 1])/(2 - p[2, 
    2] - p[1, 1]))
  # startvalue for iterations, assuming the Markov chain is ergodic alternative
  # dzetaInit <- c(0.5, 0.5)
  
  for (i in 5:n) {
    # Evaluate the densities under the two regimes create if else for different
    # functional forms
    mean1 <- c1 + phi11 * series[i - 1] + phi21 * series[i - 2] + phi31 * series[i - 
      3] + phi41 * series[i - 4]
    mean2 <- c2 + phi12 * series[i - 1] + phi22 * series[i - 2] + phi32 * series[i - 
      3] + phi42 * series[i - 4]
    eta[i, 1] <- dnorm(x = series[i], mean = mean1, sd = sigma1)
    eta[i, 2] <- dnorm(x = series[i], mean = mean2, sd = sigma2)
    
    # Evaluate the conditional density of the ith observation
    if (i == 5) {
      f[i] <- t(p %*% c(eta[i, 1], eta[i, 2])) %*% dzetaInit
    } else {
      f[i] <- t(p %*% c(eta[i, 1], eta[i, 2])) %*% c(dzeta[i - 1, 1], dzeta[i - 
        1, 2])
    }
    # Evaluate the state probabilities
    if (i == 5) {
      dzeta[i, 1] <- dzetaInit[1]
      dzeta[i, 2] <- dzetaInit[2]
    } else {
      dzeta[i, 1] <- (eta[i, 1] * (p[, 1] %*% c(dzeta[i - 1, 1], dzeta[i - 
        1, 2])))/f[i]
      dzeta[i, 2] <- (eta[i, 2] * (p[, 2] %*% c(dzeta[i - 1, 1], dzeta[i - 
        1, 2])))/f[i]
    }
  }
  logf <- sum(log(f[5:n]))
  # cat(logf, '\n') flush.console()
  output <- cbind(eta, f, dzeta)
  return(list(dzeta))
}

# Read-in data
data <- read.table("dataCO2.txt", header = TRUE)

names(data)
dim(data)

# Starting values for llh function
startVal <- c(0.5, 0.5, -6e-04, -6e-04, 0.0988, 0.0988, -0.1391, -0.1391, 0.0795, 
  0.0795, 0.0609, 0.0609, 0.02397082, 0.02397082)

# Optimize the conditional log likelihood function iteratively
# test1 <- optim(f = condLogLik, p <- startVal, series = data[2:725, 5])
# test2 <- optim(f = condLogLik, p <- test1$par, series = data[2:725, 5])
# test3 <- optim(f = condLogLik, p <- test2$par, series = data[2:725, 5])
# test4 <- optim(f = condLogLik, p <- test3$par, series = data[2:725, 5])
# test5 <- optim(f = condLogLik, p <- test4$par, series = data[2:725, 5])
# test6 <- optim(f = condLogLik, p <- test5$par, series = data[2:725, 5])
# test7 <- optim(f = condLogLik, p <- test6$par, series = data[2:725, 5])
# test8 <- optim(f = condLogLik, p <- test7$par, series = data[2:725, 5])
# test9 <- optim(f = condLogLik, p <- test8$par, series = data[2:725, 5])
# test10 <- optim(f = condLogLik, p <- test9$par, series = data[2:725, 5])
# test11 <- optim(f = condLogLik, p <- test10$par, series = data[2:725, 5])
# test12 <- optim(f = condLogLik, p <- test11$par, series = data[2:725, 5])
# test13 <- optim(f = condLogLik, p <- test12$par, series = data[2:725, 5])

# test2 <- nlm(condLogLik, p <- c(0.7, 0.8, 1.2, 0.8, 0.6, 1), series <-
# data[,1], iterlim =500) test2

# parValues <- test13$par
# use parValues to reduce computation time
parValues <- c(0.981773697, 0.969827092, 0.001729956, -0.003293143, -0.059681647, 
  0.164715677, -0.066217239, -0.194683826, 0.0086037, 0.11156463, -0.086799545, 
  0.107836808, 0.015979023, 0.032438789)

# unconditonal mean state 1
parValues[3]/(1 - parValues[5] - parValues[7] - parValues[9] - parValues[11])
# unconditonal mean state 2
parValues[4]/(1 - parValues[6] - parValues[8] - parValues[10] - parValues[12])

# Plot dzeta, prob in high regime and plot of series in-sample period
dzeta     <- condLogLikPar(theta <- parValues, series <- data[2:725, 5])  #calculate dzeta
dzeta1    <- dzeta[[1]][, 1]
dzeta2    <- dzeta[[1]][, 2]
par(mfrow = c(2, 1))
jan08     <- as.Date("01/01/08", "%d/%m/%y")
dec10     <- as.Date("01/01/11", "%d/%m/%y")
DatePlot  <-   as.Date(data$Date[6:725], ,"%d/%m/%y")

par(mar = c(4, 4, 2, 2))
plot(dzeta1[5:724] ~ DatePlot, type = "l", ylim = c(0, 1), lwd = 1.2, xlab = "Date", 
  ylab = "Regime probabilities", xaxs = "i", yaxs = "i", xlim = c(jan08, dec10), 
  yaxt = "n")
axis(side = 2, at = c(0, 0.5, 1), las = 1)

par(mar = c(4, 4, 0.5, 2))
plot(data[6:725, 5] ~ DatePlot, type = "l", lwd = 1.2, xaxs = "i", yaxs = "i", xlim = c(jan08, 
  dec10), xlab = "Date", ylab = "Log returns", ylim = c(-0.15, 0.15), yaxt = "n")
axis(side = 2, at = c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2), las = 1)
