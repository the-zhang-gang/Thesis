# SGARCH ====
fit=garchfit.sGARCH;n.list = 1
z = rugarch::residuals(fit[[n.list]],standardize=T)
tmp.1 <- rugarch::GMMTest(z)
tmp.1

n.list = 2
z = rugarch::residuals(fit[[n.list]],standardize=T)
tmp.2 <- rugarch::GMMTest(z)
tmp.2

n.list = 3
z = rugarch::residuals(fit[[n.list]],standardize=T)
tmp.3<- rugarch::GMMTest(z)
tmp.3

n.list = 4
z = rugarch::residuals(fit[[n.list]],standardize=T)
tmp.4<- rugarch::GMMTest(z)
tmp.4

n.list = 5
z = rugarch::residuals(fit[[n.list]],standardize=T)
tmp.5<- rugarch::GMMTest(z)
tmp.5

par(mfrow=c(5,1), mar = c(1,1,1,1))
plot(fit[[1]],which=11)
plot(fit[[2]],which=11)
plot(fit[[3]],which=11)
plot(fit[[4]],which=11)
plot(fit[[5]],which=11)
# eGARCH ====
fit=garchfit.eGARCH;n.list = 1
z = rugarch::residuals(fit[[n.list]],standardize=T)
tmp.1 <- rugarch::GMMTest(z, lags = 22)
tmp.1

n.list = 2
z = rugarch::residuals(fit[[n.list]],standardize=T)
tmp.2 <- rugarch::GMMTest(z)
tmp.2

n.list = 3
z = rugarch::residuals(fit[[n.list]],standardize=T)
tmp.3<- rugarch::GMMTest(z)
tmp.3

n.list = 4
z = rugarch::residuals(fit[[n.list]],standardize=T)
tmp.4<- rugarch::GMMTest(z)
tmp.4

n.list = 5
z = rugarch::residuals(fit[[n.list]],standardize=T)
tmp.5<- rugarch::GMMTest(z)
tmp.5
par(mfrow=c(5,1), mar = c(1,1,1,1))
plot(fit[[1]],which=11)
plot(fit[[2]],which=11)
plot(fit[[3]],which=11)
plot(fit[[4]],which=11)
plot(fit[[5]],which=11)

#NAGARCH ====
fit=garchfit.fGARCH.NAGARCH
par(mfrow=c(5,1), mar = c(1,1,1,1))
plot(fit[[1]],which=11)
plot(fit[[2]],which=11)
plot(fit[[3]],which=11)
plot(fit[[4]],which=11)
plot(fit[[5]],which=11)

# AVGARCH ====
fit=garchfit.fGARCH.AVGARCH
par(mfrow=c(5,1), mar = c(1,1,1,1))
plot(fit[[1]],which=11)
plot(fit[[2]],which=11)
plot(fit[[3]],which=11)
plot(fit[[4]],which=11)
plot(fit[[5]],which=11)

# TGARCH ====
fit=garchfit.fGARCH.TGARCH
z = rugarch::residuals(fit[[2]],standardize=T)
tmp.4<- rugarch::GMMTest(z)
tmp.4
par(mfrow=c(5,1), mar = c(1,1,1,1))
plot(fit[[1]],which=11)
plot(fit[[2]],which=11)
plot(fit[[3]],which=11)
plot(fit[[4]],which=11)
plot(fit[[5]],which=11)

