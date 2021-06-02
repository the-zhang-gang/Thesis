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

# eGARCH ====
fit=garchfit.eGARCH;n.list = 1
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
