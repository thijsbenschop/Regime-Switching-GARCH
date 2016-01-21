condLogLikPar <- function(theta, series){
	# calculates the conditional loglikelihood for a Markov Switching model 
	# with 2 states
	# Arguments
	# p11 transition probability to stay in state 1
	# p22 transition probability to stay in state 2
	# mu1, mu2, mean of distribution in resp. state 1 and 2
	# sigma1, sigma2, variance of distribution in resp. state 1 and 2
	# returns : value of conditional log-likelihood
	theta[1] -> p11
	theta[2] -> p22
	theta[3] -> c1
	theta[4] -> c2 
	theta[5] -> alpha01
	theta[6] -> alpha02
	theta[7] -> alpha11
	theta[8] -> alpha12
	theta[9] -> beta11
	theta[10] -> beta12
	theta[11] -> phi11
	theta[12] -> phi12
	theta[13] -> phi21
	theta[14] -> phi22
	theta[15] -> phi31
	theta[16] -> phi32
	theta[17] -> phi41
	theta[18] -> phi42
	
	if ((alpha11 + beta11) >= 1 || (alpha12 + beta12) >= 1){logf <- -999999999}
	# add restrictions on alpha 1 and beta 1
	else{
	n <- length(series)	# number of obs.
	p <- cbind(c(p11, 1-p22), c(1-p11, p22)) # transition matrix
	dzeta <- cbind(rep(0, n), rep(0, n))
	f <- cbind(rep(0, n))
	eta <- cbind(rep(0,n), rep(0,n))
	
	dzetaInit <- c((1-p[2,2])/(2-p[1,1]-p[2,2]), (1-p[1,1])/(2-p[2,2]-p[1,1]))
	# startvalue for iterations, assuming the Markov chain is ergodic
	# alternative  dzetaInit <- c(0.5, 0.5)	
	
	# creating sigma and epsilon vectors
	sigma1 <- cbind(rep(0,n))
	sigma2 <- cbind(rep(0,n))
	sigma1[5] <- sd(series)
	sigma2[5] <- sd(series)
	epsilon21 <- (series-c1)^2
	epsilon22 <- (series-c2)^2
		#for (i in 2:n){
		#sigma1[i] <- sqrt(alpha01 + (alpha11 * epsilon21[i-1]) + (beta11 * (sigma1[i-1])^2))
		#sigma2[i] <- sqrt(alpha02 + (alpha12 * epsilon22[i-1]) + (beta12 * (sigma2[i-1])^2))
		#}
	#sigma1 <- sqrt(sigma21)
	#sigma2 <- sqrt(sigma22)
	
	for (i in 5:n){
		
		# Evaluate the densities under the two regimes
		############### create if else for different functional forms ##############
		mean1 <- c1+phi11*series[i-1]+phi21*series[i-2] +phi31*series[i-3] +phi41*series[i-4] 
		mean2 <- c2+phi12*series[i-1]+phi22*series[i-2] +phi32*series[i-3] +phi42*series[i-4]
		
		eta[i, 1] <- dnorm(x=series[i], mean=mean1, sd=sigma1[i])
		eta[i, 2] <- dnorm(x=series[i], mean=mean2, sd=sigma2[i])
		
		# Evaluate the conditional density of the ith observation
		if (i == 5){
		f[i] <- t(p %*% c(eta[i,1], eta[i,2])) %*% dzetaInit
		}
		else{
		f[i] <- t(p %*% c(eta[i,1], eta[i,2])) %*% c(dzeta[i-1, 1], dzeta[i-1, 2])
		}
		# Evaluate the state probabilities
		if(i==5){
		dzeta[i, 1] <- dzetaInit[1]
		dzeta[i, 2]	<- dzetaInit[2]
		}
		else{
		dzeta[i, 1] <- (eta[i,1] * (p[,1] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		dzeta[i, 2] <- (eta[i,2] * (p[,2] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		}

		# Calculating sigma2
		if(i == 5){
		sigma1[i+1] <- sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		sigma2[i+1] <- sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		}	
		else{
		sigma1[i+1] <- sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		sigma2[i+1] <- sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		}	
	}
	logf <- sum(log(f[5:n]))
	}
	if(is.nan(logf)==TRUE){
		cat("Error : Returned not a number ", "\n")
		flush.console()
		logf <- -999999999
		}
	#output <- cbind(eta, f, dzeta)
	cat(logf, "\n")
	#cat(p11, " ", p22, " ", c1, " ", c2, "\n ", alpha01, " ", alpha02, " ", alpha11, " ", alpha12, "\n ", beta11, " ", beta12, "\n ", logf, "\n")
	flush.console()
	return(list(dzeta, sigma1, sigma2, epsilon21, epsilon22))	
}


# Function to evaluate the loglikelihood, 2 states, both AR(4)-GARCH(1,1)
condLogLik <- function(theta, series){
	# calculates the conditional loglikelihood for a Markov Switching model 
	# with 2 states
	# Arguments
	# p11 transition probability to stay in state 1
	# p22 transition probability to stay in state 2
	# mu1, mu2, mean of distribution in resp. state 1 and 2
	# sigma1, sigma2, variance of distribution in resp. state 1 and 2
	# returns : value of conditional log-likelihood
	theta[1] -> p11
	theta[2] -> p22
	theta[3] -> c1
	theta[4] -> c2 
	theta[5] -> alpha01
	theta[6] -> alpha02
	theta[7] -> alpha11
	theta[8] -> alpha12
	theta[9] -> beta11
	theta[10] -> beta12
	theta[11] -> phi11
	theta[12] -> phi12
	theta[13] -> phi21
	theta[14] -> phi22
	theta[15] -> phi31
	theta[16] -> phi32
	theta[17] -> phi41
	theta[18] -> phi42
	
	if ((alpha11 + beta11) >= 1 || (alpha12 + beta12) >= 1){logf <- -999999999}
	# add restrictions on alpha 1 and beta 1
	else{
	n <- length(series)	# number of obs.
	p <- cbind(c(p11, 1-p22), c(1-p11, p22)) # transition matrix
	dzeta <- cbind(rep(0, n), rep(0, n))
	f <- cbind(rep(0, n))
	eta <- cbind(rep(0,n), rep(0,n))
	
	dzetaInit <- c((1-p[2,2])/(2-p[1,1]-p[2,2]), (1-p[1,1])/(2-p[2,2]-p[1,1]))
	# startvalue for iterations, assuming the Markov chain is ergodic
	# alternative  dzetaInit <- c(0.5, 0.5)	
	
	# creating sigma and epsilon vectors
	sigma1 <- cbind(rep(0,n))
	sigma2 <- cbind(rep(0,n))
	sigma1[5] <- sd(series)
	sigma2[5] <- sd(series)
	epsilon21 <- (series-c1)^2
	epsilon22 <- (series-c2)^2
		#for (i in 2:n){
		#sigma1[i] <- sqrt(alpha01 + (alpha11 * epsilon21[i-1]) + (beta11 * (sigma1[i-1])^2))
		#sigma2[i] <- sqrt(alpha02 + (alpha12 * epsilon22[i-1]) + (beta12 * (sigma2[i-1])^2))
		#}
	#sigma1 <- sqrt(sigma21)
	#sigma2 <- sqrt(sigma22)
	
	for (i in 5:n){
		
		# Evaluate the densities under the two regimes
		############### create if else for different functional forms ##############
		mean1 <- c1+phi11*series[i-1]+phi21*series[i-2] +phi31*series[i-3] +phi41*series[i-4] 
		mean2 <- c2+phi12*series[i-1]+phi22*series[i-2] +phi32*series[i-3] +phi42*series[i-4]
		
		eta[i, 1] <- dnorm(x=series[i], mean=mean1, sd=sigma1[i])
		eta[i, 2] <- dnorm(x=series[i], mean=mean2, sd=sigma2[i])
		
		# Evaluate the conditional density of the ith observation
		if (i == 5){
		f[i] <- t(p %*% c(eta[i,1], eta[i,2])) %*% dzetaInit
		}
		else{
		f[i] <- t(p %*% c(eta[i,1], eta[i,2])) %*% c(dzeta[i-1, 1], dzeta[i-1, 2])
		}
		# Evaluate the state probabilities
		if(i==5){
		dzeta[i, 1] <- dzetaInit[1]
		dzeta[i, 2]	<- dzetaInit[2]
		}
		else{
		dzeta[i, 1] <- (eta[i,1] * (p[,1] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		dzeta[i, 2] <- (eta[i,2] * (p[,2] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		}

		# Calculating sigma2
		if(i == 5){
		sigma1[i+1] <- sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		sigma2[i+1] <- sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		}	
		else{
		sigma1[i+1] <- sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		sigma2[i+1] <- sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		}	
	}
	logf <- sum(log(f[5:n]))
	}
	if(is.nan(logf)==TRUE){
		cat("Error : Returned not a number ", "\n")
		flush.console()
		logf <- -999999999
		}
	#output <- cbind(eta, f, dzeta)
	#cat(logf, "\n")
	#cat(p11, " ", p22, " ", c1, " ", c2, "\n ", alpha01, " ", alpha02, " ", alpha11, " ", alpha12, "\n ", beta11, " ", beta12, "\n ", logf, "\n")
	#flush.console()
	return(-logf)	
}


#########################################################
getwd()
setwd("/Users/thijsbenschop/Desktop/MasterThesis/Data/MasterThesis3")
setwd("E:/MasterThesis/Data")
setwd("H:/MTBackup110113")
getwd()
data <- read.table("dataCO2.txt", header=TRUE)
names(data)
dim(data)
########################################################

#Coefficient(s) of GARCH(1,1) without regime switching:
# constant, alpha0, alpha1, beta1
#         mu        omega       alpha1        beta1  
#-2.6460e-04   4.6171e-06   7.2568e-02   9.1988e-01  

# Test function with values of GARCH(1,1)
test <- condLogLik(theta <- c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988, 0.0988, 0.0988, -0.1391, -0.1391, 0.0795, 0.0795, 0.0609 , 0.0609), series <- data[2:725, 5])
test
test <- condLogLik(theta <- parDE, series <- data[2:725, 5])
test

#Plot state probs
dzetaTest <- condLogLikDzeta(theta <- c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988, 0.0988, 0.0988, -0.1391, -0.1391, 0.0795, 0.0795, 0.0609 , 0.0609), series <- data[2:725, 5])
plot(dzetaTest[[1]][,1], type="l")

## Defining constaints
# Matrix with lienar combinations of parameters (p11. p22, c1, c2, alpha01, alpha02, alpa11, alpha12, beta11, beta12, phi11, phi12, phi21, phi22, phi31, phi32, phi41, phi42)
# (k x p) matrix, k constraints, p parameters 
#					
constraintMat <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0,
						  1, 0, 0, 0, 0, 0, 0, 0,
						  0, 1, 0, 0, 0, 0, 0, 0,
						  0, 0, 1, 0, 0, 0, -1, 0,
						  0, 0, 0, 1, 0, 0, 0, -1,
						  0, 0, 0, 0, 1, 0, -1, 0,
						  0, 0, 0, 0, 0, 1, 0, -1,
						  0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0, 						  0, 0, 0, 0, 0, 0, 0, 0,
						  0, 0, 0, 0, 0, 0, 0, 0
						  ), 8, 18)
constraintMat
constraintVec <-c(0, 0, 0, 0, 0, 0, -1, -1)

logRetFirstPeriod <- data[2:725, 5]

## Other optimization alghorithms
# deoptim
library(DEoptim)
#		p11, p22, c1,   c2,   alf01, alf02, alf11, alf12, beta12, beta22, phi11,  phi12,  phi21,  phi22,  phi31,  phi32,  phi41,  phi42
lowParDE <- c(0,   0,   -0.3, -0.3, 0,     0,     0,     0,     0,    0, 	-1,	  -1,     -1,     -1, 	  -1,     -1,     -1,     -1)
upParDE <-  c(1,   1,    0.3,  0.3, 0.2,   0.2,   1,     1,     1,    1, 	 1,	   1,      1,      1, 	   1,      1,      1,     1)
controlDE <- list(NP = 240, itermax = 2500)
testDE <- DEoptim(fn=condLogLik, lower=lowParDE, upper=upParDE, series=logRetFirstPeriod, control=controlDE)
# Result
Iteration: 1908 bestvalit: -1750.935434 bestmemit:    0.974023    0.881826    0.001113   -0.009027    0.000013    0.000223    0.007800    0.195207    0.864445    0.750989   -0.033932    0.301394   -0.063675   -0.210758    0.026095    0.196542   -0.031466    0.251212
Iteration: 1638 bestvalit: -1750.834437 bestmemit:    0.974023    0.881826    0.001101   -0.009656    0.000013    0.000223    0.007800    0.204368    0.864445    0.750989   -0.022887    0.301394   -0.063675   -0.210758    0.028724    0.214878   -0.041642    0.251212
parDE <- c(0.974023, 0.881826, 0.001101, -0.009656, 0.000013 , 0.000223, 0.007800, 0.204368, 0.864445, 0.750989, -0.022887, 0.301394, -0.063675, -0.210758, 0.028724, 0.214878, -0.041642, 0.251212)

# unconditional mean and sd
# state 1
parDE[3] / (1-parDE[11]-parDE[13]-parDE[15]-parDE[17])
sqrt(parDE[5] / (1-parDE[7]-parDE[9]))
# state 2
parDE[4] / (1-parDE[12]-parDE[14]-parDE[16]-parDE[18])
sqrt(parDE[6] / (1-parDE[8]-parDE[10]))

AIC <- 2*(-1750.935434) + 2*18

# Constrained optimization
test7 <- constrOptim(c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988, 0.0988, 0.0988, -0.1391, -0.1391, 0.0795, 0.0795, 0.0609 , 0.0609), ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test8 <- constrOptim(test7$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test9 <- constrOptim(test8$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test10 <- constrOptim(test9$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test11 <- constrOptim(test10$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test12 <- constrOptim(test11$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test13 <- constrOptim(test12$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test14 <- constrOptim(test13$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test15 <- constrOptim(test14$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test16 <- constrOptim(test15$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test17 <- constrOptim(test16$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test18 <- constrOptim(test17$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)

cbind(test7$par, test8$par, test9$par, test10$par, test11$par, test12$par, test13$par, test14$par, test15$par, test16$par, test17$par, test18$par) 

####
results of test 18	
parValues <- test18$par
parValues <- c

#$value
#[1] -2492.734

#$par
# [1]  1.604017e+00  3.607635e-01 -2.251704e-05 -1.613421e-01  4.475618e-10
# [6]  1.838259e+00  7.254282e-02  7.256048e-02  9.198838e-01  9.198772e-01
#[11]  1.058183e-02  4.204175e-02  7.197604e-03  1.524742e-02  6.525468e-02
#[16]  1.096865e-01 -1.724443e-03  1.949872e-01

final <- condLogLik(theta <- parValues, series <- data[2:725, 5])
dzetaFinal <- condLogLikDzeta(theta <- parValues, series <- data[2:725, 5])
plot(dzetaFinal[,1], type="l")

#different starting values
#test9 <- constrOptim(c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.01, 0.01, 0.98, 0.98), ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)

#########################################################
# Plot dzeta, prob in high regime and plot of series
# in-sample period
dzeta <- condLogLikPar(theta <- parDE, series <- data[2:725, 5]) #calculate dzeta
dzeta1 <- dzeta[[1]][,1]
dzeta2 <- dzeta[[1]][,2]
par(mfrow=c(2,1))
jan08 <- as.Date("01/01/08", "%d/%m/%y")
dec10 <- as.Date("01/01/11", "%d/%m/%y")
DatePlot <- as.Date(data$Date[6:725], , "%d/%m/%y")
par(mar=c(4,4,2,2))
plot(dzeta1[5:724]~DatePlot, type="l", ylim=c(0.0,1.0), lwd=1.2, xlab="Date", ylab="Regime probabilities", xaxs="i", yaxs="i", xlim=c(jan08, dec10), yaxt="n")
axis(side=2, at=c(0, 0.5, 1), las=1)
par(mar=c(4,4,0.5,2))
plot(data[6:725,5]~DatePlot, type="l", lwd=1.2, xaxs="i", yaxs="i", xlim=c(jan08, dec10), xlab="Date", 
ylab="Log returns", ylim=c(-0.15, 0.15), yaxt="n")
axis(side=2, at=c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15), las=1)


par(mfrow=c(2,1))
jan08 <- as.Date("01/01/08", "%d/%m/%y")
dec10 <- as.Date("01/01/11", "%d/%m/%y")
DatePlot <- as.Date(data$Date[6:725], , "%d/%m/%y")
par(mar=c(4,4,2,2))
plot(dzeta1[5:724]~DatePlot, type="l", ylim=c(0.0,1.0), lwd=1.2, xlab="Date", ylab="Regime probabilities", xaxs="i", yaxs="i", xlim=c(jan08, dec10), yaxt="n")
axis(side=2, at=c(0, 0.5, 1), las=1)
par(mar=c(4,4,0.5,2))
plot(data[6:725,5]~DatePlot, type="l", lwd=1.2, xaxs="i", yaxs="i", xlim=c(jan08, dec10), xlab="Date", ylab="Log returns", ylim=c(-0.15, 0.15), yaxt="n")
axis(side=2, at=c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15), las=1)



prob <- cbind(dzeta1, dzeta2)
write.table(prob, "probMSAR4GARCH.txt")

#compute the fitted values
series <- data[2:725, 5])
fitted.y <- rep(0,724)
fitted.y[1:4]<-NA
for(i in 5:724){
		mean1 <- parDE[3]+parDE[11]*series[i-1]+parDE[13]*series[i-2] +parDE[15]*series[i-3] +parDE[17]*series[i-4] 
		mean2 <- parDE[4]+parDE[12]*series[i-1]+parDE[14]*series[i-2] +parDE[16]*series[i-3] +parDE[18]*series[i-4]		
		fitted.y[i] <- dzeta1[i] * mean1 + dzeta2[i] * mean2
	} 
par(new=TRUE)
plot(fitted.y, , type="l", ylim=c(-0.12, 0.12), col="red")


#########################################################
####################### Forecasting #####################
#########################################################
### point forecast for period 2011-2012 (726-1183)
#1# static approach
# create forecasts and calculate errors
#matrix with estimate, error, estimated variance
static <- cbind(rep(0, 458), rep(0, 458))
# names(static) <- c("estimate", "error")

estimate <- constrOptim(c(0.8, 0.8, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988), ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)

for (i in 726:1183){
	static[i-725, 1] <- estimate					# estimate
	static[i-725, 2] <- data$logRet[i] - estimate   	# estimation error
	cat("iteration: ", i)
	flush.console()
	}
staticMSE <- (1/length(static))*sum(static[,2]^2)
staticMAE <- (1/length(static))*sum(abs(static[,2]))
staticMSE
staticMAE

#2# reestimation, recursive
#matrix with estimate, error, estimated variance
recursivePar <- parValues
prevPar <- parValues
# reestimate the parameters for the longer time series
for (i in 726:1183){
	reestMSGARCH11 <- constrOptim(prevPar, ui=constraintMat, ci=constraintVec, condLogLik, series = data[2:i, 5], grad=NULL)
	recursivePar <- cbind(recursivePar, reestMSGARCH11$par) 
	prevPar <- reestMSGARCH11$par
	cat("iteration: ", i, "/n")
	}

# write reestimation to file
write.table(recursivePar, "reestRecurMSAR4GARCH11.txt")

# load recursivePar
recursivePar <- read.table("reestRecurMSAR4GARCH11.txt", header=TRUE)

# calculate point and density estimates, confidence intervals
pointForecastRec <- rep(0, 458)
sigma2Forecast <- rep(0, 458)
p1Forecast <- rep(0, 458)
p2Forecast <- rep(0, 458)
c1Forecast <- rep(0, 458)
c2Forecast <- rep(0, 458)
dzeta1 <- rep(0,458)
dzeta2 <- rep(0,458)

for (i in 726:1183){	
	# Forecasting values for i
	# calculate dzeta[i], sigma[i]
	info <- condLogLikPar(theta <- recursivePar[,(i-725)], series <- data[2:(i-1), 5]) #calculate dzeta and sigma
	curPar <- recursivePar[,(i-725)]	# estimated parameters
	dzeta1Cur <- info[[1]][(i-2),1]
	dzeta2Cur <- info[[1]][(i-2),2]
	p1 <- dzeta1Cur * curPar[1] + dzeta2Cur * (1-curPar[2]) # prob of being in state 1 in forecast
	p2 <- dzeta2Cur * curPar[2] + dzeta1Cur * (1-curPar[1]) # prob of being in state 2 in forecast
	p1Forecast[(i-725)] <- p1
	meanCur1 <- curPar[3] + curPar[11] * data[(i-1),5] + curPar[13] * data[(i-2),5]+ curPar[15] * data[(i-3),5]+ curPar[17] * data[(i-4),5]
	meanCur2 <- curPar[4] + curPar[12] * data[(i-1),5] + curPar[14] * data[(i-2),5]+ curPar[16] * data[(i-3),5]+ curPar[18] * data[(i-4),5]
	c1Forecast[(i-725)] <-meanCur1
	c2Forecast[(i-725)] <-meanCur2
	#point estimate
	pointForecastRec[i-725] <- p1 * curPar[3] + p2 * curPar[4]
	sigma2Cur <- dzeta1Cur * (info[[2]][(i-2)])^2 + dzeta2Cur * (info[[3]][(i-2)])^2
	sigma2Forecast[i-725] <- p1 * (curPar[5] + (curPar[7] * info[[4]][(i-2)]) + curPar[9]*sigma2Cur) + p2 * (curPar[6] + (curPar[8] * info[[5]][(i-2)]) + curPar[10]*sigma2Cur)
	}

estErrorRecur <- pointForecastRec - data$logRet[726:1183]
resPlot <- cbind(estErrorRecur, sigma2Forecast)
write.table(resPlot, "resMSAR4GARCH11.txt")

recursiveMSE <- (1/length(estErrorRecur)) * sum(estErrorRecur^2)
recursiveMAE <- (1/length(estErrorRecur)) * sum(abs(estErrorRecur))
recursiveMSE
recursiveMAE

# logreturns and predicted 95 percent confidence intervals
# calculate 2,5 and 97,5 quantiles
confInterval <- cbind(rep(0, 458), rep(0, 458))
for (i in 1:458){
	confInterval[i,1] <- qnorm(0.025, mean=pointForecastRec[i], sd=sqrt(sigma2Forecast[i]))
	confInterval[i,2] <- qnorm(0.975, mean=pointForecastRec[i], sd=sqrt(sigma2Forecast[i]))
	}

conf <- cbind(confInterval, pointForecastRec)
write.table(conf, "confMSAR4GARCH.txt")

plot(data[726:1183,5], type="l", ylim=c(-0.5, 0.5))
par(new=TRUE)
plot(confInterval[2:458,1], type="l", ylim=c(-0.5, 0.5))
par(new=TRUE)
plot(confInterval[2:458,2], type="l", ylim=c(-0.5, 0.5))
par(new=TRUE)
plot(pointForecastRec, type="l", ylim=c(-0.3, 0.3), col="red")

## Checking univformity of density forecasts
# Density transformation
u <- c(rep(0, 458))
for (i in 1:length(u)){
	u[i] <- pnorm(data$logRet[i+725], mean=pointForecastRec[i], sd=sqrt(sigma2Forecast[i]))
	}
hist(u, breaks=20)

write.table(u, "uMSAR4GARCH11.txt")

# test on uniformity
#Kolmogorov Smirnov
ks.test(u, "punif")




#3# reestimation, rolling window
rolwindow12 <- cbind(rep(0, 458), rep(0, 458))	# 12 month

for (i in 726:1183){
	reestARfit12  <- arima(data$logRet[(i-241):(i-1)], order=c(order,0,0))
	estimate12  <- reestARfit3$coef[5] + reestARfit3$coef[1]*data$logRet[i-1] + reestARfit3$coef[2]*data$logRet[i-2] + reestARfit3$coef[3]*data$logRet[i-3] + reestARfit3$coef[4]*data$logRet[i-4] 
	
	rolwindow12[i-725, 1] <- estimate12		# estimate
	rolwindow12[i-725, 2]  <- data$logRet[i] - estimate12   # estimation error
	cat("iteration: ", i, "\n")
	flush.console()
	}

rolwindowMSE12  <- (1/length(rolwindow12))*sum(rolwindow12[,2]^2)
rolwindowMAE12  <- (1/length(rolwindow12))*sum(abs(rolwindow12[,2]))

########
MSEMAE <- cbind(c(staticMSE, recursiveMSE, rolwindowMSE), c(staticMAE, recursiveMAE, rolwindowMAE))
colnames(MSEMAE) <- c("MSE", "MAE")
rownames(MSEMAE) <- c("static", "recursive", "rolwin6", "rolwin12","rolwin18","rolwin24")


##############################################################
# Density forecasts
# One-day ahead forecast AR(4)-GARCH(1,1)
order <- 4
densityForecast <- cbind(rep(0, 458), rep(0, 458)) # mean and var
for (i in 726:1183){
	# mean
	reestimate <- garchFit(formula = ~arma(4,0) + garch(1,1), data=data$logRet[2:(i-1)])
	densityForecast[i-725, 1] <- estimate <- reestimate$coef[1] + reestimate$coef[2]*data$logRet[i-1] + reestimate$coef[3]*data$logRet[i-2] + reestimate$coef[4]*data$logRet[i-3] + reestimate$coef[5]*data$logRet[i-4] 
	# standard deviation TO BE CHECKED ABOUT FORECASTING
	densityForecast[i-725, 2] <- sqrt(reestimate$coef[6] + reestimate$coef[7] * (data$logRet[] - densityForecast[i-726, 1]-)^2+ reestimate$coef[8] * )
	}
reestimate@fit$coef[5]
?garchFit

# Density transformation
u <- c(rep(0, 458))
for (i in 1:length(u)){
	u[i] <- pnorm(data$logRet[i+725], mean=densityForecast[i, 1], sd=densityForecast[i,2] )
	}

hist(u, breaks=20)
?hist
# test on uniformity!!!!

# logreturns and predicted 95 percent confidence intervals
# calculate 2,5 and 97,5 quantiles
confInterval <- cbind(rep(0, 458), rep(0, 458))
for (i in 1:458){
	confInterval[i,1] <- qnorm(0.025, mean=densityForecast[i, 1], sd=densityForecast[i,2])
	confInterval[i,2] <- qnorm(0.975, mean=densityForecast[i, 1], sd=densityForecast[i,2])
	}

plot(data$logRet[726:1183], type="l", ylim=c(-0.15, 0.2))
par(new=TRUE)
plot(confInterval[,1], type="l", ylim=c(-0.15, 0.2))
par(new=TRUE)
plot(confInterval[,2], type="l", ylim=c(-0.15, 0.2))
confInterval

#################################################################################
#################################################################################
#################################################################################
test


test3 <- condLogLik(theta <- c(0.5, 0.5, mean(logRetFirstPeriod), mean(logRetFirstPeriod), var(logRetFirstPeriod), var(logRetFirstPeriod), 0.5, 0.5, 0.5, 0.5), series <- data[2:725, 5])

logRetFirstPeriod <- data[2:725,5]
class(logRetFirstPeriod)

# Starting values max of GARCH(1,1) --> max = 1750, gives error, non-finite value 
test4 <- nlm(condLogLik, p <- c(0.7, 0.8, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988), data <-logRetFirstPeriod, iterlim=10000)

# Starting value 
test2 <- nlm(condLogLik, p <- c(0.8, 0.8, mean(logRetFirstPeriod), mean(logRetFirstPeriod), var(logRetFirstPeriod), var(logRetFirstPeriod), 0.072568, 0.072568, 0.91988, 0.91988), data <-logRetFirstPeriod, iterlim=10000)

test3 <- nlm(condLogLik, p <- test2, data <-logRetFirstPeriod, iterlim=10000)

k <- 6
aic <- 2*(test2$minimum) + 2*k

# Contrained optimization
test5  <- nlminb(c(0.7, 0.8, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988), condLogLik, series = logRetFirstPeriod)

?nlminb
#####
logRetFirstPeriod <- data[2:725,5]
test11  <- optim(c(0.6, 0.6, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988), condLogLik, series = logRetFirstPeriod)
test12 <- optim(f=condLogLik, p <- test11$par , series = data[2:725,5])
# further optimizing with parameters from test3
test13 <- optim(f=condLogLik, p <- test12$par , series = data[2:725,5])
# further optimizing with parameters from test4
test14 <- optim(f=condLogLik, p <- test13$par , series = data[2:725,5])
# further optimizing with parameters from test5
test15 <- optim(f=condLogLik, p <- test14$par , series = data[2:725,5])
# further optimizing with parameters from test6
test16 <- optim(f=condLogLik, p <- test15$par , series = data[2:725,5])
# further optimizing with parameters from test7
test17 <- optim(f=condLogLik, p <- test16$par , series = data[2:725,5])
# further optimizing with parameters from test8
test18 <- optim(f=condLogLik, p <- test17$par , series = data[2:725,5])
# further optimizing with parameters from test3
test19 <- optim(f=condLogLik, p <- test18$par , series = data[2:725,5])
# further optimizing with parameters from test4
test20 <- optim(f=condLogLik, p <- test19$par , series = data[2:725,5])
# further optimizing with parameters from test5
test21 <- optim(f=condLogLik, p <- test20$par , series = data[2:725,5])
# further optimizing with parameters from test6
test22 <- optim(f=condLogLik, p <- test21$par , series = data[2:725,5])
# further optimizing with parameters from test7
test23 <- optim(f=condLogLik, p <- test22$par , series = data[2:725,5])
# further optimizing with parameters from test8
test24 <- optim(f=condLogLik, p <- test23$par , series = data[2:725,5])
cbind(test11$par, test12$par, test13$par, test14$par, test15$par, test16$par, test17$par, test18$par, test19$par, test20$par, test21$par, test22$par, test23$par, test24$par)
c(test11$value, test12$value, test13$value, test14$value, test15$value, test16$value, test17$value, test18$value, test19$value, test20$value, test21$value, test22$value, test23$value, test24$value)


# Check constraints
(test6$par[7]+test6$par[9] < 1)
(test6$par[8]+test6$par[10] < 1)
 

test7  <- constrOptim(c(0.6, 0.6, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988), ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)



?constrOptim

p11	Init <- c(0.2, 0.3, 0.4, 0.6, 0.8)
p22Init <- c(0.2, 0.3, 0.4, 0.6, 0.8)
likMat <- matrix(rep(0,25), 5, 5)
p11Mat <- matrix(rep(0,25), 5, 5)
p22Mat <- matrix(rep(0,25), 5, 5)

for (i in 1:5){
	for (j in 1:5)
	reestimate <- optim(c(p11Init[i], p22Init[j], -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988), condLogLik, series = data[2:725,5])
	likMat[i,j] <- reestimate$value
	p11Mat[i,j] <- reestimate$par[1]
	p22Mat[i,j] <- reestimate$par[2]
	}
reestimate <- optim(c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988), condLogLik, series = data[2:725,5])

####################################################################################
####################### Forecasting #####################
### point forecast for period 2011-2012 (726-1183)
#1# static approach
# create forecasts and calculate errors
static <- cbind(rep(0, 458), rep(0, 458))
# names(static) <- c("estimate", "error")
reestimate <- optim(c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988), condLogLik, series = data[2:725,5])
# estimated transition matrix
p11 <- reestimate$par[1]
p22 <- reestimate$par[2]
p <- cbind(c(p11, 1-p22), c(1-p11, p22)) # transition matrix

for (i in 726:1183){
	estimate <- ARfit$coef[5] + ARfit$coef[1]*data$logRet[i-1] + ARfit$coef[2]*data$logRet[i-2] + ARfit$coef[3]*data$logRet[i-3] + ARfit$coef[4]*data$logRet[i-4] 
	static[i-725, 1] <- estimate					# estimate
	static[i-725, 2] <- data$logRet[i] - estimate   # estimation error
	cat("iteration: ", i)
	flush.console()
	}
staticMSE <- (1/length(static))*sum(static[,2]^2)
staticMAE <- (1/length(static))*sum(abs(static[,2]))
staticMSE
staticMAE

#2# reestimation, recursive
recursive <- cbind(rep(0, 458), rep(0, 458))
for (i in 726:1183){
	reestARfit <- arima(data$logRet[2:(i-1)], order=c(order,0,0))
	reest <- ARfit
	estimate <- reestARfit$coef[5] + reestARfit$coef[1]*data$logRet[i-1] + reestARfit$coef[2]*data$logRet[i-2] + reestARfit$coef[3]*data$logRet[i-3] + reestARfit$coef[4]*data$logRet[i-4] 
	recursive[i-725, 1] <- estimate					# estimate
	recursive[i-725, 2] <- data$logRet[i] - estimate   # estimation error
	cat("iteration: ", i)
	flush.console()
	}
recursiveMSE <- (1/length(recursive)) * sum(recursive[,2]^2)
recursiveMAE <- (1/length(recursive)) * sum(abs(recursive[,2]))
recursiveMSE
recursiveMAE

#3# reestimation, rolling window
rolwindow6 <- cbind(rep(0, 458), rep(0, 458))	# 6 month
rolwindow12 <- cbind(rep(0, 458), rep(0, 458))	# 12 month
rolwindow18 <- cbind(rep(0, 458), rep(0, 458))	# 18 month
rolwindow24 <- cbind(rep(0, 458), rep(0, 458))	# 24 months

for (i in 726:1183){
	reestARfit6  <- arima(data$logRet[(i-121):(i-1)], order=c(order,0,0))
	reestARfit12  <- arima(data$logRet[(i-241):(i-1)], order=c(order,0,0))
	reestARfit18  <- arima(data$logRet[(i-361):(i-1)], order=c(order,0,0))
	reestARfit24 <- arima(data$logRet[(i-481):(i-1)], order=c(order,0,0))
	estimate6  <- reestARfit1$coef[5] + reestARfit1$coef[1]*data$logRet[i-1] + reestARfit1$coef[2]*data$logRet[i-2] + reestARfit1$coef[3]*data$logRet[i-3] + reestARfit1$coef[4]*data$logRet[i-4] 
	estimate12  <- reestARfit3$coef[5] + reestARfit3$coef[1]*data$logRet[i-1] + reestARfit3$coef[2]*data$logRet[i-2] + reestARfit3$coef[3]*data$logRet[i-3] + reestARfit3$coef[4]*data$logRet[i-4] 
	estimate18  <- reestARfit6$coef[5] + reestARfit6$coef[1]*data$logRet[i-1] + reestARfit6$coef[2]*data$logRet[i-2] + reestARfit6$coef[3]*data$logRet[i-3] + reestARfit6$coef[4]*data$logRet[i-4] 
	estimate24 <- reestARfit12$coef[5] + reestARfit12$coef[1]*data$logRet[i-1] + reestARfit12$coef[2]*data$logRet[i-2] + reestARfit12$coef[3]*data$logRet[i-3] + reestARfit12$coef[4]*data$logRet[i-4] 
	
	rolwindow6[i-725, 1]  <- estimate6			# estimate
	rolwindow12[i-725, 1] <- estimate12		# estimate
	rolwindow18[i-725, 1] <- estimate18		# estimate
	rolwindow12[i-725, 1] <- estimate24			# estimate
	rolwindow6[i-725, 2]  <- data$logRet[i] - estimate6   # estimation error
	rolwindow12[i-725, 2]  <- data$logRet[i] - estimate12   # estimation error
	rolwindow18[i-725, 2]  <- data$logRet[i] - estimate18   # estimation error
	rolwindow24[i-725, 2] <- data$logRet[i] - estimate24  # estimation error
	cat("iteration: ", i, "\n")
	flush.console()
	}
rolwindowMSE6  <- (1/length(rolwindow6))*sum(rolwindow6[,2]^2)
rolwindowMAE6  <- (1/length(rolwindow6))*sum(abs(rolwindow6[,2]))
rolwindowMSE12  <- (1/length(rolwindow12))*sum(rolwindow12[,2]^2)
rolwindowMAE12  <- (1/length(rolwindow12))*sum(abs(rolwindow12[,2]))
rolwindowMSE18  <- (1/length(rolwindow18))*sum(rolwindow18[,2]^2)
rolwindowMAE18  <- (1/length(rolwindow18))*sum(abs(rolwindow18[,2]))
rolwindowMSE24 <- (1/length(rolwindow24))*sum(rolwindow24[,2]^2)
rolwindowMAE24 <- (1/length(rolwindow24))*sum(abs(rolwindow24[,2]))
rolwindowMSE   <- c(rolwindowMSE6, rolwindowMSE12, rolwindowMSE18, rolwindowMSE24)
rolwindowMAE   <- c(rolwindowMAE6, rolwindowMAE12, rolwindowMAE18, rolwindowMAE24)

########
MSEMAE <- cbind(c(staticMSE, recursiveMSE, rolwindowMSE), c(staticMAE, recursiveMAE, rolwindowMAE))
colnames(MSEMAE) <- c("MSE", "MAE")
rownames(MSEMAE) <- c("static", "recursive", "rolwin6", "rolwin12","rolwin18","rolwin24")

##############################################################
# Density forecasts
# One-day ahead forecast AR(4)-GARCH(1,1)
order <- 4
densityForecast <- cbind(rep(0, 458), rep(0, 458)) # mean and var
for (i in 726:1183){
	# mean
	reestimate <- garchFit(formula = ~arma(4,0) + garch(1,1), data=data$logRet[2:(i-1)])
	densityForecast[i-725, 1] <- estimate <- reestimate$coef[1] + reestimate$coef[2]*data$logRet[i-1] + reestimate$coef[3]*data$logRet[i-2] + reestimate$coef[4]*data$logRet[i-3] + reestimate$coef[5]*data$logRet[i-4] 
	# standard deviation TO BE CHECKED ABOUT FORECASTING
	densityForecast[i-725, 2] <- sqrt(reestimate$coef[6] + reestimate$coef[7] * (data$logRet[] - densityForecast[i-726, 1]-)^2+ reestimate$coef[8] * )
	}
reestimate@fit$coef[5]
?garchFit

# Density transformation
u <- c(rep(0, 458))
for (i in 1:length(u)){
	u[i] <- pnorm(data$logRet[i+725], mean=densityForecast[i, 1], sd=densityForecast[i,2] )
	}

hist(u, breaks=20)
?hist
# test on uniformity!!!!

# logreturns and predicted 95 percent confidence intervals
# calculate 2,5 and 97,5 quantiles
confInterval <- cbind(rep(0, 458), rep(0, 458))
for (i in 1:458){
	confInterval[i,1] <- qnorm(0.025, mean=densityForecast[i, 1], sd=densityForecast[i,2])
	confInterval[i,2] <- qnorm(0.975, mean=densityForecast[i, 1], sd=densityForecast[i,2])
	}

plot(data$logRet[726:1183], type="l", ylim=c(-0.15, 0.2))
par(new=TRUE)
plot(confInterval[,1], type="l", ylim=c(-0.15, 0.2))
par(new=TRUE)
plot(confInterval[,2], type="l", ylim=c(-0.15, 0.2))
confInterval
