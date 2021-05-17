---
output:
  bookdown::pdf_document2:
    template: templates/brief_template.tex
    citation_package: biblatex
  #bookdown::word_document2: default
  #bookdown::html_document2: default
documentclass: book
bibliography: references.bib
---

# Empirical Findings {#analysis}

\minitoc <!-- this will include a mini table of contents-->

## Density of the returns



### MLE distribution parameters













In table \@ref(tab:disttable) we can see...
<!-- Stephane, I leave this to you :), maybe you can use in-text notation like in R programming or like filippo did in Data methodology for table1-->

```
## Warning in sqrt(diag(varcov)): NaNs produced

## Warning in sqrt(diag(varcov)): NaNs produced

## Warning in sqrt(diag(varcov)): NaNs produced
```

```
## Warning in matrix(value, n, p): data length [5] is not a sub-multiple or
## multiple of the number of columns [2]
```

![](04-Findings_files/figure-latex/MLE tables for diferent series-1.pdf)<!-- --> \begin{table}[h!]

\caption{(\#tab:dsTable)Maximum likelihood estimates of unconditional distribution functions}
\centering
\begin{threeparttable}
\begin{tabular}[t]{lllllllrr}
\toprule
 & \$\textbackslash{}mu\$ & \$\textbackslash{}sigma\$ & \$\textbackslash{}lambda\$ & \$p\$ & \$q\$ & \$\textbackslash{}nu\$ & \$L\$ & AIC\\
\midrule
SGT & 0.02 & 1.321 & -0.04 & 1.381 & 3.317 &  & -13973.01 & 27956.01\\
 & (0.013) & (0.026) & (0.012) & (0.071) & (0.534) &  &  & \\
SGED & 0.02 & 1.274 & -0.018 & 0.918 & Inf &  & -14008.18 & 27956.01\\
 & (0.01) & (0.016) & (0.008) & (0.016) &  &  &  & \\
GED & 0.032 & 1.276 & 0 & 0.913 & Inf &  & -14009.09 & 28028.17\\
\addlinespace
 & (0.005) & (0.016) &  & (0.016) &  &  &  & \\
ST & 0.019 & 1.487 & 0.949 &  &  & 2.785 & -13997.35 & 28002.70\\
 & (0.014) & (0.056) & (0.013) &  &  & (0.1) &  & \\
T & 0.056 & 1.494 &  &  &  & 2.765 & -14005.14 & 28016.29\\
 & (0.01) & (0.056) &  &  &  & (0.097) &  & \\
\addlinespace
Normal & 0.017 & 1.304 & 0 & 2 & Inf &  & -15093.32 & 30196.64\\
 & (0.014) & (0.015) &  &  &  &  &  & \\
\bottomrule
\end{tabular}
\begin{tablenotes}
\item Notes
\end{tablenotes}
\end{threeparttable}
\end{table}


```

## Results of GARCH with constant higher moments

<!--# Here comes our main part [FILIPPO] -> to do!  -->


```r
distributions <- c("norm", "std", "sstd", "ged", "sged")
#garchspec <- garchfit <- garchforecast <- stdret <- vector(mode = "list", length = length(distributions))
#names(garchspec) <- names(garchfit) <- names(garchforecast) <- names(stdret) <- distributions
Models.garch <- c("sGARCH", "eGARCH","fGARCH.AVGARCH","fGARCH.NAGARCH", "gjrGARCH", "fGARCH.TGARCH", "iGARCH", "EWMA")

for(i in 1:length(Models.garch)){
assign(paste0("garchspec.",Models.garch[i]),vector(mode = "list", length = length(distributions)))
assign(paste0("garchfit.",Models.garch[i]),vector(mode = "list", length = length(distributions)))
assign(paste0("stdret.",Models.garch[i]),vector(mode = "list", length = length(distributions)))
} 

# ls(pattern = "garchspec.")
# sapply(ls(pattern = "garchspec."), FUN = setNames, distributions)

#.sGARCH--------------------------
for(i in 1:length(distributions)){
# Specify a GARCH model with constant mean
garchspec.sGARCH[[i]] <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "sGARCH", garchOrder = c(1,1), variance.targeting = F), 
                     distribution.model = distributions[i])
# Estimate the model
garchfit.sGARCH[[i]] <- ugarchfit(data = R, spec = garchspec.sGARCH[[i]])
# Compute stdret using residuals()
stdret.sGARCH[[i]] <- residuals(garchfit.sGARCH[[i]], standardize = TRUE)
}

#.eGARCH-------------------
for(i in 1:length(distributions)){
# Specify a GARCH model with constant mean
garchspec.eGARCH[[i]] <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "eGARCH", variance.targeting = F), 
                     distribution.model = distributions[i])
# Estimate the model
garchfit.eGARCH[[i]] <- ugarchfit(data = R, spec = garchspec.eGARCH[[i]])
# Compute stdret using residuals()
stdret.eGARCH[[i]] <- residuals(garchfit.eGARCH[[i]], standardize = TRUE)
}

#.fGARCH.NAGARCH------------------------
for(i in 1:length(distributions)){
# Specify a GARCH model with constant mean
garchspec.fGARCH.NAGARCH[[i]] <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "fGARCH", submodel = "NAGARCH", variance.targeting = F),
                     distribution.model = distributions[i])
# Estimate the model
garchfit.fGARCH.NAGARCH[[i]] <- ugarchfit(data = R, spec = garchspec.fGARCH.NAGARCH[[i]])
# Compute stdret using residuals()
stdret.fGARCH.NAGARCH[[i]] <- residuals(garchfit.fGARCH.NAGARCH[[i]], standardize = TRUE)
}

#.fGARCH.AVGARCH------------------------
for(i in 1:length(distributions)){
# Specify a GARCH model with constant mean
garchspec.fGARCH.AVGARCH[[i]] <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "fGARCH", submodel = "AVGARCH", variance.targeting = F),
                     distribution.model = distributions[i])
# Estimate the model
garchfit.fGARCH.AVGARCH[[i]] <- ugarchfit(data = R, spec = garchspec.fGARCH.AVGARCH[[i]])
# Compute stdret using residuals()
stdret.fGARCH.AVGARCH[[i]] <- residuals(garchfit.fGARCH.AVGARCH[[i]], standardize = TRUE)
}

#.gjrGARCH------------------
for(i in 1:length(distributions)){
# Specify a GARCH model with constant mean
garchspec.gjrGARCH[[i]] <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "gjrGARCH", variance.targeting = F), 
                     distribution.model = distributions[i])
# Estimate the model
garchfit.gjrGARCH[[i]] <- ugarchfit(data = R, spec = garchspec.gjrGARCH[[i]])
# Compute stdret using residuals()
stdret.gjrGARCH[[i]] <- residuals(garchfit.gjrGARCH[[i]], standardize = TRUE)
}

#fGARCH.TGARCH-------------------
for(i in 1:length(distributions)){
# Specify a GARCH model with constant mean
garchspec.fGARCH.TGARCH[[i]] <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "fGARCH", submodel = "TGARCH", variance.targeting = F), 
                     distribution.model = distributions[i])
# Estimate the model
garchfit.fGARCH.TGARCH[[i]] <- ugarchfit(data = R, spec = garchspec.fGARCH.TGARCH[[i]])
# Compute stdret using residuals()
stdret.fGARCH.TGARCH[[i]] <- residuals(garchfit.fGARCH.TGARCH[[i]], standardize = TRUE)
}

#.iGARCH--------------------
for(i in 1:length(distributions)){
# Specify a GARCH model with constant mean
garchspec.iGARCH[[i]] <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "iGARCH", variance.targeting = F), 
                     distribution.model = distributions[i])
# Estimate the model
garchfit.iGARCH[[i]] <- ugarchfit(data = R, spec = garchspec.iGARCH[[i]])
# Compute stdret using residuals()
stdret.iGARCH[[i]] <- residuals(garchfit.iGARCH[[i]], standardize = TRUE)
}

#.csGARCH-----------------
# for(i in 1:length(distributions)){
# # Specify a GARCH model with constant mean
# garchspec.csGARCH[[i]] <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
#                      variance.model = list(model = "csGARCH", variance.targeting = F),
#                      distribution.model = distributions[i])
# # Estimate the model
# garchfit.csGARCH[[i]] <- ugarchfit(data = R, spec = garchspec.csGARCH[[i]])
# # Compute stdret using residuals()
# stdret.csGARCH[[i]] <- residuals(garchfit.csGARCH[[i]], standardize = TRUE)
# }


# we need EWMA
for(i in 1:length(distributions)){
# Specify a GARCH model with constant mean
garchspec.EWMA[[i]] <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "iGARCH", variance.targeting = F),
                     distribution.model = distributions[i], fixed.pars = list(omega=0))
# Estimate the model
garchfit.EWMA[[i]] <- ugarchfit(data = R, spec = garchspec.EWMA[[i]])
# Compute stdret using residuals()
stdret.EWMA[[i]] <- residuals(garchfit.EWMA[[i]], standardize = TRUE)
}


#  make the histogram
# 
# chart.Histogram(stdret.iGARCH[[1]], methods = c("add.normal","add.density" ),
#                 colorset = c("gray","red","blue"))
```


<!-- What do we do, do we do 6 tables for all the distributions?  -->

```r
table3 <- matrix(nrow = 12, ncol = 5)
colnames(table3) <- distributions

#trying a loop, maybe you can solve that @filippo?
## column loop i = normal distribution, std, sstd, ged, sged
table3[1,1] <- garchfit.sGARCH[[1]]@fit$coef[1] #first parameter estimate
table3[2,1] <- garchfit.sGARCH[[1]]@fit$se.coef[1] #first standard error
table3[3,1] <- garchfit.sGARCH[[1]]@fit$coef[2] #second parameter estimate
table3[4,1] <- garchfit.sGARCH[[1]]@fit$se.coef[2]



#...
table3 <- round(table3, 3)

# for (i in length(distributions)) {
#   for (j in nrow(table3)) {
#     table3[j,i] <- garchfit.sGARCH[[i]]@fit$coef
#     table3[j+1,i] <-garchfit.sGARCH[[i]]@fit$se.coef
#     }
# }

print("sGARCH")
garchfit.sGARCH[[1]]@fit$coef
garchfit.sGARCH[[1]]@fit$se.coef
```


```r
print("iGARCH")
garchfit.iGARCH[[1]]@fit$coef
garchfit.iGARCH[[1]]@fit$se.coef
```


```r
print("EWMA")
garchfit.EWMA[[1]]@fit$coef
c(garchfit.EWMA[[1]]@fit$se.coef[1:2],NA,garchfit.EWMA[[1]]@fit$se.coef[3], NA)
```


```r
print("eGARCH")
garchfit.eGARCH[[1]]@fit$coef
garchfit.eGARCH[[1]]@fit$se.coef
```


```r
print("gjrGARCH")
garchfit.gjrGARCH[[1]]@fit$coef
garchfit.gjrGARCH[[1]]@fit$se.coef
```


```r
print("NAGARCH")
garchfit.fGARCH.NAGARCH[[1]]@fit$coef
garchfit.fGARCH.NAGARCH[[1]]@fit$se.coef
```


```r
print("TGARCH")
garchfit.fGARCH.TGARCH[[1]]@fit$coef
garchfit.fGARCH.TGARCH[[1]]@fit$se.coef
```


```r
print("TSGARCH/AVGARCH")
garchfit.fGARCH.AVGARCH[[1]]@fit$coef
garchfit.fGARCH.AVGARCH[[1]]@fit$se.coef
```
<!-- Here comes the description -->

## Results of GARCH with time-varying higher moments


```r
require(racd)
require(rugarch)
require(parallel)
require(xts)
# ACD specification
sGARCH_ACDspec = acdspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(variance.targeting = TRUE),
distribution.model = list(model = 'jsu', skewOrder = c(1, 1, 1), shapeOrder = c(1,1,1), skewmodel = 'quad', shapemodel = 'pwl'))

# sGARCH
cl = makePSOCKcluster(10)
fit = acdfit(sGARCH_ACDspec, as.data.frame(R), solver = 'msoptim', solver.control = list(restarts = 10),cluster = cl) #very long process: starts from different starting values to find an optimum
```


```r
# plotxts comes from implementing https://stackoverflow.com/a/50051183/271616
# par(mfrow = c(2, 2), mai = c(0.75, 0.75, 0.3, 0.3))
# cm <- plot.zoo(xts(fit@model$modeldata$data, fit@model$modeldata$index), auto.grid = FALSE,minor.ticks = FALSE, main = 'Conditional Mean',yaxis.right = F, col = 'steelblue')
# cm <- lines(fitted(fit), col = 2)
# cm
# cs <- plot(xts(abs(fit@model$modeldata$data),fit@model$modeldata$index), auto.grid = FALSE,
# minor.ticks = FALSE, main = 'Conditional Sigma', yaxis.right = F,col = 'grey')
# cs <- lines(sigma(fit), col = 'steelblue')
# cs
# plot(racd::skewness(fit), col = 'steelblue',yaxis.right = F, main = 'Conditional Skewness')
# plot(racd::kurtosis(fit), col = 'steelblue', yaxis.right = F,main = 'Conditional Excess Kurtosis')

# pnl <- function(fitted(fit),xts(fit@model$modeldata$data, fit@model$modeldata$index), ...) {
#   panel.number <- parent.frame()$panel.number
# 	if (panel.number == 1) lines(fitted(fit), xts(fit@model$modeldata$data, fit@model$modeldata$index),col = "red")
# 	lines(fitted(fit),xts(fit@model$modeldata$data, fit@model$modeldata$index), col = "red")
# }
# plot(xts(fit@model$modeldata$data, fit@model$modeldata$index), auto.grid = T,minor.ticks = FALSE,major.ticks=T, yaxis.right = F, main = 'Conditional Mean', col = 'steelblue', xlab = "", screens = 1, ylab="") #panel = pnl
# # lines(fitted(fit), col = 2) + grid()
# 
# plot(xts(fit@model$modeldata$data, fit@model$modeldata$index), auto.grid = T,minor.ticks = FALSE,major.ticks=T, yaxis.right = F, main = 'Conditional Mean', col = 'steelblue', xlab = "", screens = 1, ylab="", )
```
