#Table.GARCH.function-------------------------------------
Table.GARCH.function <- function(GARCHfit.object = garchfit.eGARCH){
  #Making objects to fill   
  list.table.3.tmp <- vector(mode = "list", length = length(distributions)*2)
  list.table.3 <- vector(mode = "list", length = length(distributions))
  names(list.table.3) <- distributions
  ref.distr <- seq(from = 1, to = length(distributions)*2, by = 2)
  #Retrieving all the data from the original lists
  for(i in 1:length(distributions)){
    list.table.3.tmp[[ref.distr[i]]] <- GARCHfit.object[[i]]@fit$coef
    list.table.3.tmp[[ref.distr[i]+1]] <- GARCHfit.object[[i]]@fit$robust.matcoef[,2]
  }
  list.pvalues <- vector(mode = "list", length = length(distributions)*2)
  for(i in 1:length(distributions)){
    list.pvalues[[ref.distr[i]]] <- GARCHfit.object[[i]]@fit$coef
    list.pvalues[[ref.distr[i]+1]] <- GARCHfit.object[[i]]@fit$robust.matcoef[,4]
  }
  #From list of vectors of different lengths to list of matrixes with empty spaces filled with NAs
  for(i in 1:length(distributions)){
    list.table.3[[i]] <- cbind(list.table.3.tmp[[ref.distr[i]]], list.table.3.tmp[[ref.distr[i]+1]][match(names(list.table.3.tmp[[ref.distr[i]]]), names(list.table.3.tmp[[ref.distr[i]+1]]))])
  }
  list.table.3.p <- vector(mode = "list", length = length(distributions))
  for(i in 1:length(distributions)){
    list.table.3.p[[i]] <- cbind(list.pvalues[[ref.distr[i]]], list.pvalues[[ref.distr[i]+1]][match(names(list.pvalues[[ref.distr[i]]]), names(list.pvalues[[ref.distr[i]+1]]))])
  }
  list.table.3 <- lapply(list.table.3, round,3)
  for(i in 1:length(list.table.3)) {
    list.table.3[[i]][which(list.table.3.p[[i]][,2]>=0.10),2] <- paste0("(",list.table.3[[i]][which(list.table.3.p[[i]][,2]>=0.10),2],")")
    list.table.3[[i]][which(list.table.3.p[[i]][,2]<0.10&list.table.3.p[[i]][,2]>0.05),2] <- paste0("(",list.table.3[[i]][which(list.table.3.p[[i]][,2]<0.10&list.table.3.p[[i]][,2]>0.05),2],")*")
    list.table.3[[i]][which(list.table.3.p[[i]][,2]<=0.05&list.table.3.p[[i]][,2]>0.01),2] <- paste0("(",list.table.3[[i]][which(list.table.3.p[[i]][,2]<=0.05&list.table.3.p[[i]][,2]>0.01),2],")**")
    list.table.3[[i]][which(list.table.3.p[[i]][,2]<=0.01),2] <- paste0("(",list.table.3[[i]][which(list.table.3.p[[i]][,2]<=0.01),2],")***")
  }
  
  #Function to rearrange the list from a list of matrices to list of vectors
  list.restructure <- function(object = list.table.3){
    len.table <- length(object)
    len.inside.list <- rep(NA, len.table)
    for(i in 1:len.table){len.inside.list[i] <- nrow(object[[i]])}
    adj.list <- vector(mode = "list", length = len.table)
    ref.list <- vector(mode = "list", length = len.table)
    names.list <- vector(mode = "list", length = len.table)
    for(i in 1:len.table){ref.list[[i]] <- seq(from = 1, to = len.inside.list[i]*2, by = 2)}
    for(i in 1:len.table){adj.list[[i]] <- names.list[[i]] <- rep(NA, len.inside.list[i]*2)}
    for(i in 1:len.table){
      adj.list[[i]][ref.list[[i]]] <- object[[i]][,1]
      adj.list[[i]][ref.list[[i]]+1] <- object[[i]][,2]
    }
    for(i in 1:len.table){
      names.list[[i]][ref.list[[i]]] <- rownames(object[[i]])
      names.list[[i]][ref.list[[i]]+1] <- paste0("p-val ",rownames(object[[i]]))
      #names.list[[i]][ref.list[[i]]+1] <- ""
    }
    names(adj.list) <- distributions
    for(i in 1:len.table){names(adj.list[[i]]) <- names.list[[i]]}
    return(adj.list)
  }
  #Unlisting and removing NAs
  adj.list <- list.restructure(object = list.table.3)
  table.3.matrix <- matrix(unlist(lapply(adj.list, `length<-`, max(lengths(adj.list)))),ncol = length(distributions),nrow = max(lengths(adj.list)), byrow = F)
  colnames(table.3.matrix) <- names(adj.list)
  table.3.matrix[is.na(table.3.matrix)] <- ""
  table.3.matrix[table.3.matrix=="(NA)"] <- ""
  #Adjustments for std & ged distributions
  table.3.matrix[c(length(table.3.matrix[,2])-1,length(table.3.matrix[,2])),2] <- tail(table.3.matrix[table.3.matrix[,2]!="",2],2)
  table.3.matrix[c(length(table.3.matrix[,2])-3,length(table.3.matrix[,2])-2),2] <- ""
  table.3.matrix[c(length(table.3.matrix[,4])-1,length(table.3.matrix[,4])),4] <- tail(table.3.matrix[table.3.matrix[,4]!="",4],2)
  table.3.matrix[c(length(table.3.matrix[,4])-3,length(table.3.matrix[,4])-2),4] <- ""
  #Log-Likelyhoods
  LLH <- rep(NA, length(distributions))
  for(i in 1:length(distributions)){LLH[i] <- GARCHfit.object[[i]]@fit$LLH}
  names.table.3 <- revalue(names(adj.list[[3]]), c("mu"="$\\alpha_0$", "ar1"="$\\alpha_1$", "omega"= "$\\beta_0$", "alpha1"="$\\beta_1$", "beta1"="$\\beta_2$", "gamma1"="$\\gamma$", "skew"="$\\xi$", "shape"="$\\kappa$", "eta11"="$rot$", "eta21"="$shift$"), warn_missing = F)
  kappa.object <- cbind(matrix(data = "", nrow = 2, ncol = 3) ,table.3.matrix[(nrow(table.3.matrix)-1):(nrow(table.3.matrix)),4:5])
  eta <- cbind(table.3.matrix[(nrow(table.3.matrix)-1):(nrow(table.3.matrix)),1:3], matrix(data = "", nrow = 2, ncol = 2) )
  table.3.matrix[(nrow(table.3.matrix)-1):(nrow(table.3.matrix)),] <- kappa.object
  table.3.matrix <- rbind(table.3.matrix, eta, round(LLH,3))
  table.3.matrix <- cbind(c(names.table.3,"$\\eta$","p-val eta","$LLH$"), table.3.matrix)
  table.3.matrix <- as.data.frame(table.3.matrix, row.names = c(names.table.3,"$\\eta$","","LLH"), stringsAsFactors = FALSE)
  return(table.3.matrix)
}
#DistMLE------------------------------------------------
DistMLE <- function(R=R) {
  ### SGT
  X.data <- X ~ coredata(R)
  SGT_start <- list(mu=0,sigma=1, lambda = 0.5, p=2, q=8) # p = kappa, q = nu
  SGT_result <- sgt.mle(X.f = X.data, start = SGT_start)
  names(SGT_result$estimate) <- names(SGT_result$std.error) <-  c("alpha", "beta", "xi", "kappa", "eta")
  SGT_sumResult <- summary(SGT_result)
  SGT_AIC <- 2*length(SGT_result$estimate) - 2*SGT_sumResult$maximum
  SGT_pv <- 2* stats::pnorm(-abs(SGT_result$estimate/(SGT_result$std.error)))
  SGT_significance <- matrix(ncol = 5,nrow = 1)
  colnames(SGT_significance) <-  c("alpha", "beta", "xi", "kappa", "eta")
  for(i in 1:5){
    if(SGT_pv[i]>=0.10){
      SGT_significance[i] <-  ""
    }
    if(SGT_pv[i]<0.10){
      SGT_significance[i] <-  "*"  
    }
    if(SGT_pv[i]<0.05){
      SGT_significance[i] <-  "**"  
    }
    if(SGT_pv[i]<0.01){
      SGT_significance[i] <-  "***"  
    }
  }
  SGT_significance[is.na(SGT_significance)] <- ""
  
  ### SGED
  SGED_start <- list(mean=0,sd=1, nu = 2, xi = 1.5)
  SGED_result <- fitdistrplus::fitdist(data = as.vector(coredata(R)), distr = "sged", method = "mle", SGED_start)
  names(SGED_result$estimate) <- names(SGED_result$sd) <-  c("alpha", "beta", "kappa", "xi")
  SGED_sumResult <- summary(SGED_result)
  SGED_AIC <- 2*length(SGED_result$estimate-1) - 2*SGED_sumResult$loglik
  SGED_pv <- 2* stats::pnorm(-abs(SGED_result$estimate/(SGED_result$sd)))
  SGED_pv <- c(SGED_pv[1:2], SGED_pv[4], SGED_pv[3],NA)
  SGED_significance <- matrix(ncol =5,nrow = 1)
  colnames(SGED_significance) <-  c("alpha", "beta", "xi", "kappa", "eta")
  for(i in c(1:4)){
    if(SGED_pv[i]>=0.10){
      SGED_significance[i] <-  ""
    }
    if(SGED_pv[i]<0.10){
      SGED_significance[i] <-  "*"  
    }
    if(SGED_pv[i]<0.05){
      SGED_significance[i] <-  "**"  
    }
    if(SGED_pv[i]<0.01){
      SGED_significance[i] <-  "***"  
    }
  }
  SGED_significance[is.na(SGED_significance)] <- ""
  
  ### GED
  GED_start <- list(mean = 0, sd = 1, nu = 2)
  GED_result <- fitdistrplus::fitdist(data = as.vector(R),distr="ged", start = GED_start)
  names(GED_result$estimate) <- names(GED_result$sd) <-  c("alpha", "beta", "kappa")
  GED_sumResult <- summary(GED_result)
  GED_AIC <- 2*length(GED_result$estimate-2) - 2*GED_sumResult$loglik
  GED_pv <- 2* stats::pnorm(-abs(GED_sumResult$estimate/(GED_sumResult$sd)))
  GED_pv <- c(GED_pv[1:2], NA, GED_pv[3], NA)
  GED_significance <- matrix(ncol = 5,nrow = 1)
  colnames(GED_significance) <-  c("alpha", "beta", "xi", "kappa", "eta")
  for(i in c(1,2,4)){
    if(GED_pv[i]>=0.10){
      GED_significance[i] <-  ""
    }
    if(GED_pv[i]<0.10){
      GED_significance[i] <-  "*"  
    }
    if(GED_pv[i]<0.05){
      GED_significance[i] <-  "**"  
    }
    if(GED_pv[i]<0.01){
      GED_significance[i] <-  "***"  
    }
  }
  GED_significance[is.na(GED_significance)] <- ""
  
  ### ST (fitdist)
  ST_start <- list(mean=0,sd=1, nu = 5, xi=1.5)
  ST_result <- fitdistrplus::fitdist(data = as.vector(coredata(R)), distr = "sstd", method = "mle", ST_start)
  names(ST_result$estimate) <- names(ST_result$sd) <-  c("alpha", "beta", "eta","xi")
  ST_sumResult <- summary(ST_result)
  ST_sumResult$aic
  ST_pvalue <- 2*stats::pnorm(-abs(ST_sumResult$estimate/(ST_sumResult$sd)))
  ST_pvalue <-  c(ST_pvalue[1:2], ST_pvalue[4], NA, ST_pvalue[3])
  ST_significance <- matrix(ncol = 5,nrow = 1)
  colnames(ST_significance) <-  c("alpha", "beta", "xi", "kappa", "eta")
  for(i in c(1:3,5)){
    if(ST_pvalue[i]>=0.10){
      ST_significance[i] <-  ""
    }
    if(ST_pvalue[i]<0.10){
      ST_significance[i] <-  "*"  
    }
    if(ST_pvalue[i]<0.05){
      ST_significance[i] <-  "**"  
    }
    if(ST_pvalue[i]<0.01){
      ST_significance[i] <-  "***"  
    }
  }
  ST_significance[is.na(ST_significance)] <- ""
  
  ### T (fitdist)
  T_start <- list(mean = 0, sd = 1, nu = 5)
  T_result <- fitdistrplus::fitdist(data = as.vector(coredata(R)), distr = "std", method = "mle", T_start)
  names(T_result$estimate) <- names(T_result$sd) <-  c("alpha", "beta", "eta")
  T_sumResult <- summary(T_result)
  T_pvalues <-  2*stats::pnorm(-abs(T_result$estimate/(T_result$sd)))
  T_pvalues <-  c(T_pvalues[1:2], NA, NA, T_pvalues[3])
  T_significance <- matrix(ncol = 5,nrow = 1)
  colnames(T_significance) <-  c("alpha", "beta", "xi", "kappa", "eta")
  for(i in c(1,2,5)){
    if(T_pvalues[i]>=0.10){
      T_significance[i] <-  ""
    }
    if(T_pvalues[i]<0.10){
      T_significance[i] <-  "*"  
    }
    if(T_pvalues[i]<0.05){
      T_significance[i] <-  "**"  
    }
    if(T_pvalues[i]<0.01){
      T_significance[i] <-  "***"  
    }
  }
  T_significance[is.na(T_significance)] <- ""
  
  ### Normal (fitdist)
  Normal_start <- list(mean = 0, sd=1)
  Normal_result <- fitdistrplus::fitdist(as.vector(coredata(R)), "norm","mle", Normal_start)
  names(Normal_result$estimate) <- names(Normal_result$sd) <-  c("alpha", "beta")
  Normal_sumResult <- summary(Normal_result)
  Normal_AIC <- 2*length(Normal_result$estimate-3) - 2*Normal_sumResult$loglik
  Normal_pv <- 2*stats::pnorm(-abs(Normal_result$estimate/(Normal_result$sd)))
  Normal_significance <- matrix(ncol = 5,nrow = 1)
  colnames(Normal_significance) <-  c("alpha", "beta", "xi", "kappa", "eta")
  for(i in c(1:2)){
    if(Normal_pv[i]>=0.10){
      Normal_significance[i] <-  ""
    }
    if(Normal_pv[i]<0.10){
      Normal_significance[i] <-  "*"  
    }
    if(Normal_pv[i]<0.05){
      Normal_significance[i] <-  "**"  
    }
    if(Normal_pv[i]<0.01){
      Normal_significance[i] <-  "***"  
    }
  }
  Normal_significance[is.na(Normal_significance)] <- ""
  
  #maximum likelihood estimates of unconditional distribution functions
  Table2 <- matrix(nrow = 6, ncol = 7)
  colnames(Table2) <- c("alpha","beta","xi","kappa","eta","LLH","AIC")
  rownames(Table2) <- c("SGT","SGED","GED","ST","T","Normal")
  # "alpha","beta","xi","kappa","nu","L","AIC"
  Table2[1,1:5] <- SGT_result$estimate
  Table2[1,6] <- SGT_result$maximum
  Table2[1,7] <- SGT_AIC
  # "alpha","beta","xi","kappa","nu","L","AIC"
  Table2[2,1:2] <- SGED_result$estimate[1:2]
  Table2[2,3] <- SGED_result$estimate[4] - 1
  Table2[2,4] <- SGED_result$estimate[3]
  Table2[2,5] <- Inf
  Table2[2,6] <- SGED_result$loglik
  Table2[2,7] <- SGT_AIC
  # "alpha","beta","xi","kappa","nu","L","AIC"
  Table2[3,1:2] <- GED_result$estimate[1:2]
  Table2[3,3] <- 0
  Table2[3,4] <- GED_result$estimate[3]
  Table2[3,5] <- Inf
  Table2[3,6] <- GED_result$loglik
  Table2[3,7] <- GED_AIC
  # "alpha","beta","xi","kappa","nu","L","AIC"
  Table2[4,1:2] <- ST_result$estimate[1:2]
  Table2[4,3] <- ST_result$estimate[4] - 1
  Table2[4,4] <- 2
  Table2[4,5] <- ST_result$estimate[3]/2 # eta parameter is not exactly true here... See formula SGT package
  Table2[4,6] <- ST_result$loglik
  Table2[4,7] <- ST_result$aic
  # "alpha","beta","xi","kappa","nu","L","AIC"
  Table2[5,1:2] <- T_result$estimate[1:2]
  Table2[5,3] <- 0
  Table2[5,4] <- 2
  Table2[5,5] <- T_result$estimate[3]/2 # eta parameters
  Table2[5,6] <- T_result$loglik
  Table2[5,7] <- T_result$aic
  
  Table2[6,1:2] <- Normal_result$estimate[1:2]
  Table2[6,3] <- 0
  Table2[6,4] <- 2
  Table2[6,5] <- Inf
  Table2[6,6] <- Normal_result$loglik
  Table2[6,7] <- Normal_result$aic
  
  #adding SE
  Table2_SE <- matrix(nrow = 12, ncol = 7)
  
  Table2_SE <- as.data.frame(Table2_SE)
  #coefficients placement
  Table2_SE[1,] <- Table2[1,]
  Table2_SE[3,] <- Table2[2,]
  Table2_SE[5,] <- Table2[3,]
  Table2_SE[7,] <- Table2[4,]
  Table2_SE[9,] <- Table2[5,]
  Table2_SE[11,] <- Table2[6,]
  #fixing rownames, in a new column
  colnames(Table2_SE) <- c("alpha","beta","xi","kappa","eta","L","AIC")
  rownames(Table2_SE) <- 1:12
  tablenames <- c("SGT","","SGED","","GED","","ST","","T","","Normal","")
  Table2_SE <- cbind(tablenames,Table2_SE)
  Table2_SE[c(1,3,5,7,9,11),-1] <- round(Table2_SE[c(1,3,5,7,9,11),-1], 3) #round
  
  #SEs (basic)
  # "alpha","beta","xi","kappa","eta","L","AIC"
  Table2_SE[2,2:6] <- paste0("(",round(SGT_result$std.error, 3),")",SGT_significance)
  # "alpha","beta","xi","kappa","eta","L","AIC"
  #sd is correct, but the significance has to be fixed, same for GED, ST and T
  Table2_SE[4,2:6] <- paste0("(",round(c(SGED_result$sd[1:2],SGED_result$sd[4], SGED_result$sd[3], NA),3),")", SGED_significance)
  # "alpha","beta","xi","kappa","eta","L","AIC"
  Table2_SE[6,2:6] <- paste0("(",round(c(GED_result$sd[1:2],NA,GED_result$sd[3],NA), 3),")",GED_significance)
  # "alpha","beta","xi","kappa","eta","L","AIC"
  Table2_SE[8,2:6] <- paste0("(",round(c(ST_result$sd[1:2],ST_result$sd[4],NA,ST_result$sd[3]),3),")",ST_significance)
  # "alpha","beta","xi","kappa","eta","L","AIC"
  Table2_SE[10,2:6] <- paste0("(",round(c(T_result$sd[-3],NA,NA,T_result$sd[3]),3),")",T_significance)
  # "alpha","beta","xi","kappa","eta","L","AIC"
  Table2_SE[12,2:6] <- paste0("(",round(c(Normal_result$sd, NA,NA,NA), 3),")",Normal_significance)
  
  
  #Remove SE for limiting cases
  # Table2_SE[6,4] <- NA
  # Table2_SE[c(4,6),6] <- NA
  
  totalresults <- list(table = Table2_SE, results = list(SGT = SGT_sumResult,SGED = SGED_sumResult, GED = GED_sumResult, ST = ST_sumResult, T = T_sumResult, Norm =Normal_sumResult))
  return(totalresults)
}
#Table.3.function-------------------------------------------------------
Table.3.function <- function(distribution = "sstd"){
  ref <- which(names(Table.3.eGARCH)==distribution)#Reference column
  #Taking all the relevant values
  AVGARCH <- Table.3.fGARCH.AVGARCH[,c(1,ref)]
  iGARCH <- Table.3.iGARCH[,c(1,ref)]
  eGARCH <- Table.3.eGARCH[,c(1,ref)]
  sGARCH <- Table.3.sGARCH[,c(1,ref)]
  gjrGARCH <- Table.3.gjrGARCH[,c(1,ref)]
  EWMA <- Table.3.EWMA[,c(1,ref)]
  NAGARCH <- Table.3.fGARCH.NAGARCH[,c(1,ref)]
  TGARCH <- Table.3.fGARCH.TGARCH[,c(1,ref)]
  #Assigning the right names to columns
  colnames(AVGARCH)[-1] <- rep("AVGARCH", length(colnames(AVGARCH)[-1]))
  colnames(iGARCH)[-1] <- rep("iGARCH", length(colnames(iGARCH)[-1]))
  colnames(eGARCH)[-1] <- rep("eGARCH", length(colnames(eGARCH)[-1]))
  colnames(sGARCH)[-1] <- rep("sGARCH", length(colnames(sGARCH)[-1]))
  colnames(gjrGARCH)[-1] <- rep("gjrGARCH", length(colnames(gjrGARCH)[-1]))
  colnames(EWMA)[-1] <- rep("EWMA", length(colnames(EWMA)[-1]))
  colnames(NAGARCH)[-1] <- rep("NAGARCH", length(colnames(NAGARCH)[-1]))
  colnames(TGARCH)[-1] <- rep("TGARCH", length(colnames(TGARCH)[-1]))
  #Binding all the columns & cleaning & ordering data
  Table3 <- full_join(full_join(full_join(full_join(full_join(full_join(full_join(sGARCH, iGARCH),eGARCH),gjrGARCH),EWMA),NAGARCH),TGARCH),AVGARCH)
  Table3[is.na(Table3)] <- ""
  Table3 <- rbind(Table3[Table3[,1]!="$LLH$",], Table3[Table3[,1]=="$LLH$",])
  Table3[substr(Table3[,1],1,5)=="p-val",1] <- ""
  colnames(Table3) <- c("", colnames(Table3)[-1])
  return(Table3)
}
#roll---------------------------------------------------
roll <- function(x) { 
  cl = makePSOCKcluster(10)
  result = ugarchroll(x,data = R,n.ahead = 1,n.start = 2500,forecast.length=1,refit.every = 22,solver = "hybrid",refit.window="moving",calculate.VaR=TRUE,VaR.alpha=0.01, cluster = cl, keep.coef=TRUE)
  stopCluster(cl)
  return(result)
}
#getVaR----------------------------------------------
getVaR <- function(x) {
  VaRresult <- as.data.frame(x,which="VaR")$`alpha(1%)`
}
#gettests-----------------------------------------
gettests <- function(x) { 
  # modify this according to your needs
  uc <- BacktestVaR(R[2501:nrow(R)],x,0.01, 5)$LRuc
  cc <- BacktestVaR(R[2501:nrow(R)],x,0.01, 5)$LRcc
  ae <- BacktestVaR(R[2501:nrow(R)],x,0.01, 5)$AE
  dq <- BacktestVaR(R[2501:nrow(R)],x,0.01, 5)$DQ
  return(list = list(ae = ae, uc = uc, cc = cc, dq = dq))
}
#ESTesting------------------------------------------
ESTesting <- function(garch.back, var, dist){
  
  roll <- as.data.frame(garch.back, which = "density",)
  f <- function(x, skew, shape) qdist(dist, p = x, mu = 0, sigma = 1, skew = skew, shape = shape)
  #VaR in same fashion as ES below, and gives the same VaR as previously found so it is redundant.
  # VaR <- roll[,'Mu'] + qdist(dist, 0.01, mu = 0, sigma =1, skew = roll[,'Skew'], shape = roll[,'Shape']) * roll[,'Sigma']
  ES <- roll['Mu'] + roll['Sigma']*apply(roll, 1, function(x)
    integrate(f,0,0.01, skew = x['Skew'], shape = x['Shape'])$value/0.01)
  names(ES) <- "ES"
  ESTest_object <- ESTest(0.01, R[-c(1:2500),], ES[,1], var, conf.level = 0.99, boot = T)
  # length(R)
  # nrow(ES)
  # length(VaR.sstd.egarch)
  AE_ES <- as.numeric(ESTest_object[2])/as.numeric(ESTest_object[1])
  AE_ES.p <- ESTest_object[4] #bootstrapped p-value
  return(list(ES,AE_ES,AE_ES.p))
}
#Significances-----------------------------------------
uc.signific <- function(x) {
  if(x$uc[2]<0.1)  y = "*"
  if(x$uc[2]<0.05) y = "**"
  if(x$uc[2]<0.01) y = "***"
  else if(x$uc[2]>0.1) y = ""
  return(y)
}
cc.signific <- function(x) {
  if(x$cc[2]<0.1)  y = "*"
  if(x$cc[2]<0.05) y = "**"
  if(x$cc[2]<0.01) y = "***"
  else if(x$cc[2]>0.1) y = ""
  return(y)
}

dq.signific <- function(x) {
  if(x$dq$pvalue<0.1)  y = "*"
  if(x$dq$pvalue<0.05) y = "**"
  if(x$dq$pvalue<0.01) y = "***"
  else if(x$dq$pvalue>0.1) y = ""
  return(y)
}

ES.signific <- function(x) {
  if(x<0.1)  y = "*"
  if(x<0.05) y = "**"
  if(x<0.01) y = "***"
  else if(x>0.1) y = ""
  return(y)
}

signific <- function(x) {
  if(x<0.1)  y = "*"
  if(x<0.05) y = "**"
  if(x<0.01) y = "***"
  else if(x>0.1) y = ""
  return(y)
}
#ESTestingACD----------------------------------------
ESTestingACD <- function(acd.back, var, dist){
  roll <- rollacd@forecast$density
  f <- function(x, skew, shape) qdist(dist, p = x, mu = 0, sigma = 1, skew = skew, shape = shape)
  #VaR in same fashion as ES below, and gives the same VaR as previously found so it is redundant.
  # VaR <- roll[,'Mu'] + qdist(dist, 0.01, mu = 0, sigma =1, skew = roll[,'Skew'], shape = roll[,'Shape']) * roll[,'Sigma']
  ES <- roll['Mu'] + roll['Sigma']*apply(roll, 1, function(x)
    integrate(f,0,0.01, skew = x['Skew'], shape = x['Shape'])$value/0.01)
  names(ES) <- "ES"
  ESTest_object <- ESTest(0.01, R[-c(1:2500),], ES[,1], var, conf.level = 0.99, boot = T)
  # length(R)
  # nrow(ES)
  # length(VaR.sstd.egarch)
  AE_ES <- as.numeric(ESTest_object[2])/as.numeric(ESTest_object[1])
  AE_ES.p <- ESTest_object[4] #bootstrapped p-value
  return(list(ES,AE_ES,AE_ES.p))
}
#Table.tests.function--------------------------------------------
Table.tests.function <- function(function.to.apply, return.matrix = F){
  LJ.Test <- vector(mode = "list", length = length(Models.garch))
  names(LJ.Test) <- Models.garch
  for(i in 1:length(LJ.Test)){
    LJ.Test[[i]] <- matrix(NA, nrow = 2, ncol = length(distributions), dimnames = list(c("Statistic", "P Values"), names(Table.3)))
  }
  
  for(i in 1:length(distributions)){
    LJ.Test[[1]][,i] <- function.to.apply(garchmodel = garchfit.sGARCH, n.list = i)
    LJ.Test[[2]][,i] <- function.to.apply(garchmodel = garchfit.eGARCH, n.list = i)
    LJ.Test[[3]][,i] <- function.to.apply(garchmodel = garchfit.fGARCH.AVGARCH, n.list = i)
    LJ.Test[[4]][,i] <- function.to.apply(garchmodel = garchfit.fGARCH.NAGARCH, n.list = i)
    LJ.Test[[5]][,i] <- function.to.apply(garchmodel = garchfit.gjrGARCH, n.list = i)
    LJ.Test[[6]][,i] <- function.to.apply(garchmodel = garchfit.fGARCH.TGARCH, n.list = i)
    LJ.Test[[7]][,i] <- function.to.apply(garchmodel = garchfit.iGARCH, n.list = i)
    LJ.Test[[8]][,i] <- function.to.apply(garchmodel = garchfit.EWMA, n.list = i)
  }
  
  LJ.Test <- lapply(LJ.Test, round, 3)
  
  for(i in 1:length(LJ.Test)){
    LJ.Test[[i]][1,which(LJ.Test[[i]][2,]>=0.10)] <- LJ.Test[[i]][1,which(LJ.Test[[i]][2,]>=0.10)]
    LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<0.10&LJ.Test[[i]][2,]>0.05)] <- paste0(LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<0.10&LJ.Test[[i]][2,]>0.05)],"*")
    LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<=0.05&LJ.Test[[i]][2,]>0.01)] <- paste0(LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<=0.05&LJ.Test[[i]][2,]>0.01)],"**")
    LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<=0.01)] <- paste0(LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<=0.01)],"***")
    LJ.Test[[i]] <- LJ.Test[[i]][-2,]
  }
  
  LJ.test.matrix <- cbind(LJ.Test[[1]],LJ.Test[[2]],LJ.Test[[3]],LJ.Test[[4]],LJ.Test[[5]],LJ.Test[[6]],LJ.Test[[7]],LJ.Test[[8]])
  colnames(LJ.test.matrix) <- Models.garch.clean
  LJ.test.df <- as.data.frame(LJ.test.matrix)
  if(return.matrix == T){return(LJ.test.matrix)}else{return(LJ.test.df)}
}
#Lj.box.test.function.2------------------------------
Lj.box.test.function.2 <- function(garchmodel = garchfit.sGARCH, n.list = 1, n.lag=22){
  tmp.1 <- Box.test(rugarch::residuals(garchmodel[[n.list]], standardize = T)^2,n.lag, "Ljung-Box") 
  test <- c(tmp.1$statistic, tmp.1$p.value)
  return(test)
}
#archlmtest-------------------------------
archlmtest = function (garchmodel, lags = 22, demean = FALSE, n.list = 1){
  x <- rugarch::residuals(garchmodel[[n.list]], standardize = T) 
  if(any(!is.finite(x))) x[!is.finite(x)] = 0
  x = as.vector(x)
  if(demean) x = scale(x, center = TRUE, scale = FALSE)
  lags = lags + 1
  mat = embed(x^2, lags)
  arch.lm = summary(lm(mat[, 1] ~ mat[, -1]))
  STATISTIC = arch.lm$r.squared * length(resid(arch.lm))
  names(STATISTIC) = "Chi-squared"
  PARAMETER = lags - 1
  names(PARAMETER) = "df"
  PVAL = 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD = "ARCH LM-test"
  result = list(statistic = STATISTIC, parameter = PARAMETER,
                p.value = PVAL, method = METHOD)
  class(result) = "htest"
  final <- c("Statistic"=result$statistic,"P Value" = result$p.value)
  return(final)
}
#GMMTest.2-------------------------------------------
GMMTest.2 = function(garchmodel, lags = 22, skew=0, kurt=3, conf.level = 0.95, n.list = 1, col.numb = 1){
  require(rugarch)
  if(length(skew)>1) sk = skew[-c(1:lags)] else sk = skew
  if(length(kurt)>1) ku = kurt[-c(1:lags)] else ku = kurt
  z = matrix(rugarch::residuals(garchmodel[[n.list]],standardize=TRUE), ncol = 1)
  N = dim(z)[1] - lags #DEGREES OF FREEDOM
  zlag = z[-c(1:lags), , drop = FALSE]
  orthmat = matrix(NA, ncol = 4, nrow = 2)
  colnames(orthmat) = c("mean", "var", "skew", "ex.kurt")
  rownames(orthmat) = c("statistic", "p.value")
  f1 = zlag[,1]
  orthmat[1:2,1] = c(mean(f1), dt(mean(f1)/sqrt(mean(f1^2)/N),df = N)     )
  f2 = (zlag[,1]^2)-1
  orthmat[1:2,2] = c(mean(f2), dt(mean(f1)/sqrt(mean(f1^2)/N), df = N))
  f3 = zlag[,1]^3-sk
  orthmat[1:2,3] = c(mean(f3), dt(mean(f1)/sqrt(mean(f1^2)/N), df = N))
  f4 = zlag[,1]^4-ku
  orthmat[1:2,4] = c(mean(f4), dt(mean(f1)/sqrt(mean(f1^2)/N), df = N))
  M = rbind(t(f1), t(f2), t(f3), t(f4))
  # all moments
  moment.mat = orthmat[1:2,1:4]
  return(moment.mat[,col.numb])
}
#GMM.table.function----------------------------------------
GMM.table.function <- function(function.to.apply, return.matrix = F, col.numb = 1){
  LJ.Test <- vector(mode = "list", length = length(Models.garch))
  names(LJ.Test) <- Models.garch
  for(i in 1:length(LJ.Test)){
    LJ.Test[[i]] <- matrix(NA, nrow = 2, ncol = length(distributions), dimnames = list(c("Statistic", "P Values"), names(Table.3)))
  }
  
  for(i in 1:length(distributions)){
    LJ.Test[[1]][,i] <- function.to.apply(garchmodel = garchfit.sGARCH, n.list = i, col.numb = col.numb)
    LJ.Test[[2]][,i] <- function.to.apply(garchmodel = garchfit.eGARCH, n.list = i, col.numb = col.numb)
    LJ.Test[[3]][,i] <- function.to.apply(garchmodel = garchfit.fGARCH.AVGARCH, n.list = i, col.numb = col.numb)
    LJ.Test[[4]][,i] <- function.to.apply(garchmodel = garchfit.fGARCH.NAGARCH, n.list = i, col.numb = col.numb)
    LJ.Test[[5]][,i] <- function.to.apply(garchmodel = garchfit.gjrGARCH, n.list = i, col.numb = col.numb)
    LJ.Test[[6]][,i] <- function.to.apply(garchmodel = garchfit.fGARCH.TGARCH, n.list = i, col.numb = col.numb)
    LJ.Test[[7]][,i] <- function.to.apply(garchmodel = garchfit.iGARCH, n.list = i, col.numb = col.numb)
    LJ.Test[[8]][,i] <- function.to.apply(garchmodel = garchfit.EWMA, n.list = i, col.numb = col.numb)
  }
  
  LJ.Test <- lapply(LJ.Test, round, 3)
  
  for(i in 1:length(LJ.Test)){
    LJ.Test[[i]][1,which(LJ.Test[[i]][2,]>=0.10)] <- LJ.Test[[i]][1,which(LJ.Test[[i]][2,]>=0.10)]
    LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<0.10&LJ.Test[[i]][2,]>0.05)] <- paste0(LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<0.10&LJ.Test[[i]][2,]>0.05)],"*")
    LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<=0.05&LJ.Test[[i]][2,]>0.01)] <- paste0(LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<=0.05&LJ.Test[[i]][2,]>0.01)],"**")
    LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<=0.01)] <- paste0(LJ.Test[[i]][1,which(LJ.Test[[i]][2,]<=0.01)],"***")
    LJ.Test[[i]] <- LJ.Test[[i]][-2,]
  }
  
  LJ.test.matrix <- cbind(LJ.Test[[1]],LJ.Test[[2]],LJ.Test[[3]],LJ.Test[[4]],LJ.Test[[5]],LJ.Test[[6]],LJ.Test[[7]],LJ.Test[[8]])
  colnames(LJ.test.matrix) <- Models.garch.clean
  LJ.test.df <- as.data.frame(LJ.test.matrix)
  if(return.matrix == T){return(LJ.test.matrix)}else{return(LJ.test.df)}
}
#Jarque.bera.test.function-------------------------------
Jarque.bera.test.function <- function(garchmodel = garchfit.sGARCH, n.list = 1){
  tmp.1 <- jarqueberaTest(rugarch::residuals(garchmodel[[n.list]], standardize = T))
  test <- c(tmp.1@test$statistic, tmp.1@test$p.value)
  return(test)
}
#--------------------------
