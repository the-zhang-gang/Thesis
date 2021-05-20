### SGT plot fit
xvals = seq(-3,3,by=0.01)
SGT_mu <- SGT_result$estimate[1]
SGT_sigma <- SGT_result$estimate[2]
SGT_lambda <- SGT_result$estimate[3]
SGT_p <- SGT_result$estimate[4]
SGT_q <- SGT_result$estimate[5]
plotsgt <- plot(xvals, dsgt(xvals, mu = SGT_mu, sigma = SGT_sigma, lambda = SGT_lambda, p = SGT_p, q = SGT_q), col="red", type ="l",main = "SGT", ylab="Density")
plotsgt <- lines(density(coredata(series)))

### SGED Plot fit 
SGED_mu <- SGED_result$estimate[1]
SGED_sigma <- SGED_result$estimate[2]
SGED_lambda <- SGED_result$estimate[4]
SGED_p <- SGED_result$estimate[3]
SGED_q <- Inf 

plotsged <- plot(xvals, dsged(xvals, mean = SGED_mu, sd = SGED_sigma, xi = SGED_lambda, nu = SGED_q), col="red", type ="l", main = "SGED", ylab="Density")
plotsged <- lines(density(coredata(series)))

### GED Plot fit 
GED_mu <- GED_result$estimate[1]
GED_sigma <- GED_result$estimate[2]
GED_lambda <- 0
GED_p <- GED_result$estimate[3]
GED_q <- Inf
plotged <- plot(xvals, dged(xvals, mean = GED_mu, sd = GED_sigma, nu = GED_q), col="red", type ="l", main = "GED", ylab="Density")
plotged <- lines(density(coredata(series)))

### ST Plot fit (fitdist)
ST_mean <- ST_result$estimate[1]
ST_sigma <- ST_result$estimate[2]
ST_nu <- ST_result$estimate[3] # q
ST_xi <- ST_result$estimate[4] #lamda

plotst <- plot(xvals, dsstd(xvals, mean = ST_mean, sd = ST_sigma, nu = ST_nu, xi=ST_xi), col="red", type ="l", main = "ST", ylab = "Density")
plotst <- lines(density(coredata(series)))

### T Plot fit (fitdist)
T_mean <- T_result$estimate[1]
T_sd <- T_result$estimate[2]
T_nu <- T_result$estimate[3] #q

plott <- plot(xvals, dstd(xvals, mean = T_mean, sd = T_sd, nu = T_nu), col="red", type ="l", main = "T", ylab="Density")
plott <- lines(density(coredata(series)))

### Normal  Plot fit

Normal_mu <- Normal_result$estimate[1]
Normal_sigma <- Normal_result$estimate[2]

plotnorm <- plot(xvals, dnorm(xvals, mean = Normal_mu, sd = Normal_sigma), col="red", type ="l", main = "Normal", ylab="Density")
plotnorm <- lines(density(coredata(series)))