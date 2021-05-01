# packages
library(readxl)
library(xts)
library(PerformanceAnalytics)
library(kableExtra)

data <- read_excel("data/datastream.xlsx",col_types = c("date", rep("numeric", 6)),
                                       skip = 2)
colnames(data) <- c("Date",gsub(pattern = " ", replacement = "_",gsub(pattern = " - PRICE INDEX", replacement='' , colnames(data)[2:7])))
Price_indices <- as.xts(data[,-1], order.by = data$Date)
R <- diff(Price_indices, log = TRUE, na.pad = FALSE) * 100 
x <-  table.Stats(R[,-c(3:5)])
x <- x[c(3,5:6,9,14:16),] #selecting relevant columns only
x %>%
  kbl(format='latex',
      caption = "Summary statistics of the returns",
      label = 'dsTable',
      booktabs = T,
      position = "h!",
      digits = 3) %>%
  kable_classic(full_width = F) %>%
  footnote(general = sprintf('This table shows the descriptive statistics of the returns of a selection of indices (%s and %s) for the period from %s until %s', 
                             gsub(pattern = "_", replacement = " ", paste(colnames(x)[-length(colnames(x))],collapse=", ")), 
                             gsub(pattern = "_", replacement = " ",colnames(x)[length(colnames(x))]), 
                             data$Date[2],data$Date[nrow(data)]),
           footnote_as_chunk = T,
           general_title = 'Note: ',
           threeparttable = TRUE)

#### Eurostoxx 50 ====
plot(Price_indices$EURO_STOXX_50, major.ticks = "years", grid.ticks.on = "years", col = "blue", 
     main = "EuroStoxx 50 Price Index")
plot(R$EURO_STOXX_50, major.ticks = "years", grid.ticks.on = "years", col = "blue", 
     main = "EuroStoxx 50 Price log returns")
acf(coredata(R$EURO_STOXX_50), lag = 22, main = "Autocorrelation function EuroStoxx 50 Price Log Returns")
#In the Figure XXX we observe that autocorrelations are rather small so that it is difficult to predict future outcomes using, e.g., an AR model.
#However from Figure XXXX, there can be assumed that there is volatility clustering because of the conditional heteroscedasticity.

acf(coredata(abs(R$EURO_STOXX_50-mean(R$EURO_STOXX_50))), lag = 22, main = "Autocorrelation function EuroStoxx 50 Price Log Returns")
#positive autocorrelation reflects the presence of volatility clusters


chart.RollingPerformance(R$EURO_STOXX_50, FUN = "sd", width = 252, colorset = "red", main = "EuroStoxx 50 rolling 252-day volatility")
# unconditional value is the sd in table.Stats
e <- scale(R$EURO_STOXX_50, center = TRUE, scale = FALSE) # demean returns
u1 <- abs(e)
sd_window <- rollapply(data = R$EURO_STOXX_50, width = 252, FUN = "sd")
chart.TimeSeries(cbind(sd_window, u1), main = "Absolute returns and rolling vol")
# the rolling standard deviation jums up whenever there is a large absolute residual. Large absolute residuals cluster in time

# test for conditional heteroscedasticity: LB
m <- 22
Box.test(u1, lag = m, type = "Ljung-Box")
# or LM
regout <- lm(u1 ~ lag(u1, 1:m)) # run regression
r2 <- summary(regout)$r.squared # get the R-squared coefficient
lm <- r2 * length(resid(regout)) # multiply by the number of observations
plm <- pchisq(lm, m, lower.tail = FALSE) # compute p-value
print(c(m = m, LM = lm, p.value = plm)) 

# simpler using FinTS
library(FinTS)
ArchTest(R$EURO_STOXX_50, lags = m, demean = TRUE)
ArchTest(e, lags = m) # this is the same as previous
## Asymmetric volatility tested 
library(dynlm)
neg <- (e < 0) # is the news negative?
avol <- dynlm(u1 ~ lag(neg))
summary(avol) # estimate for lag(neg) is significantly positive: after a negative shock a higher absolute residual follows

# EWMA ? do we want this as well

# GARCH
library(rugarch)
arch1_spec <- ugarchspec(variance.model = list(garchOrder = c(1,0)), 
                         mean.model = list(armaOrder = c(0,0)),
                         distribution.model = "norm")
arch1_fit <- ugarchfit(spec = arch1_spec, data = R$EURO_STOXX_50)
print(arch1_fit)

plot(sigma(arch1_fit))
plot(arch1_fit, which = 3)
plot(arch1_fit, which = 10)
plot(arch1_fit, which = 8)
zhat <- residuals(arch1_fit, standardize = TRUE)
plot(zhat)
GMMTest(zhat)

garch11_spec <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                           mean.model = list(armaOrder = c(0, 0)))
garch11_fit <- ugarchfit(garch11_spec, R$EURO_STOXX_50)
plot(garch11_fit, which = 3)
plot(garch11_fit, which = 11)
plot(garch11_fit, which = 8)
tgarch11_spec <- ugarchspec(variance.model = list(garchOrder = c(1,1)), 
                            mean.model = list(armaOrder = c(0,0)), 
                            distribution.model = "std")
tgarch11_fit <- ugarchfit(tgarch11_spec, R$EURO_STOXX_50)
plot(tgarch11_fit, which = 3)
plot(tgarch11_fit, which = 11)
plot(tgarch11_fit, which = 8)
mean(sigma(arch1_fit))
(garch11_for <- ugarchforecast(garch11_fit, n.ahead = 10))
plot(garch11_for, which = 1)
plot(garch11_for, which = 3)

igarch11_spec <- ugarchspec(variance.model = list(garchOrder = c(1,1), model = "iGARCH"), 
                        mean.model = list(armaOrder = c(0,0)))
igarch11_fit <- ugarchfit(igarch11_spec, data = R$EURO_STOXX_50)
coef(igarch11_fit)

signbias(garch11_fit) #absolute values are printed is a missed opportunity, as it does not allow to check the direction of the effect
gjr.garch_spec <- ugarchspec(variance.model = list(model = "gjrGARCH",
                                                   garchOrder = c(1,1)), 
                             mean.model = list(armaOrder = c(0,0)), 
                             distribution.model = "norm")
gjr.garch_fit <- ugarchfit(gjr.garch_spec, R$EURO_STOXX_50)
show(gjr.garch_fit)
plot(garch11_fit, which = 12)
plot(gjr.garch_fit, which = 12)

egarch_spec <- ugarchspec(variance.model = list(model = "eGARCH",
                                                garchOrder = c(1,1)), 
                          mean.model = list(armaOrder = c(0,0)), 
                          distribution.model = "norm")
egarch_fit <- ugarchfit(egarch_spec, R$EURO_STOXX_50)
plot(egarch_fit, which = 12)

#### STOXX 600 ====

#### FTSE 100 ====



variance.models <- c("egarch","sGARCH","gjrGARCH", "iGARCH","fGARCH")
distribution.models <- c("norm", "std", "sstd","sged", "ged")