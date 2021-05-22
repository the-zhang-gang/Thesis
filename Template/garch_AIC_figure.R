autoarfima(R, ar.max = 2, ma.max = 2, criterion = c("AIC","BIC","SIC","HQIC"), method = c("partial", "full"))
?autoarfima
do_single_garch <- function(x, 
                            type_model, 
                            type_dist, 
                            lag_ar, 
                            lag_ma, 
                            lag_arch, 
                            lag_garch) {
  require(rugarch)
  
  
  spec = ugarchspec(variance.model = list(model =  type_model, 
                                          garchOrder = c(lag_arch, lag_garch)),
                    mean.model = list(armaOrder = c(lag_ar, lag_ma)),
                    distribution = type_dist)
  
  message('Estimating ARMA(',lag_ar, ',', lag_ma,')-',
          type_model, '(', lag_arch, ',', lag_garch, ')', 
          ' dist = ', type_dist,
          appendLF = FALSE)
  
  try({
    my_rugarch <- list()
    my_rugarch <- ugarchfit(spec = spec, data = x)
  })
  
  if (!is.null(coef(my_rugarch))) {
    message('\tDone')
    
    AIC <- rugarch::infocriteria(my_rugarch)[1]
    BIC <- rugarch::infocriteria(my_rugarch)[2]
  } else {
    message('\tEstimation failed..')
    
    AIC <- NA
    BIC <- NA
  }
  
  est_tab <- tibble(lag_ar, 
                    lag_ma,
                    lag_arch,
                    lag_garch,
                    AIC =  AIC,
                    BIC = BIC,
                    type_model = type_model,
                    type_dist,
                    model_name = paste0('ARMA(', lag_ar, ',', lag_ma, ')+',
                                        type_model, '(', lag_arch, ',', lag_garch, ') ',
                                        type_dist) ) 
  
  return(est_tab)
}

do_single_garch2 <- function(x, 
                            type_model, 
                            sub_model,
                            type_dist, 
                            lag_ar, 
                            lag_ma, 
                            lag_arch, 
                            lag_garch) {
  require(rugarch)
  
  
  spec = ugarchspec(variance.model = list(model =  type_model, 
                                          submodel = sub_model,
                                          garchOrder = c(lag_arch, lag_garch)),
                    mean.model = list(armaOrder = c(lag_ar, lag_ma)),
                    distribution = type_dist)
  
  message('Estimating ARMA(',lag_ar, ',', lag_ma,')-',
          type_model, sub_model,'(', lag_arch, ',', lag_garch, ')', 
          ' dist = ', type_dist,
          appendLF = FALSE)
  
  try({
    my_rugarch <- list()
    my_rugarch <- ugarchfit(spec = spec, data = x)
  })
  
  if (!is.null(coef(my_rugarch))) {
    message('\tDone')
    
    AIC <- rugarch::infocriteria(my_rugarch)[1]
    BIC <- rugarch::infocriteria(my_rugarch)[2]
  } else {
    message('\tEstimation failed..')
    
    AIC <- NA
    BIC <- NA
  }
  
  est_tab <- tibble(lag_ar, 
                    lag_ma,
                    lag_arch,
                    lag_garch,
                    AIC =  AIC,
                    BIC = BIC,
                    type_model = type_model,
                    sub_model = sub_model,
                    type_dist,
                    model_name = paste0('ARMA(', lag_ar, ',', lag_ma, ')+',
                                        type_model,sub_model, '(', lag_arch, ',', lag_garch, ') ',
                                        type_dist) ) 
  
  return(est_tab)
}

#' Finds best ARMA-GARCH model 
#'
#' @param x A (numeric) vector of returns
#' @param type_models Type of models (see rugarch::rugarchspec)
#' @param dist_to_use Type of distributions to use (see rugarch::rugarchspec)
#' @param max_lag_AR Maximum lag for AR parameter
#' @param max_lag_MA Maximum lag for MA parameter
#' @param max_lag_ARCH Maximum lag for ARCH parameter
#' @param max_lag_GARCH Maximum lag for GARCH parameter
#'
#' @return A list with results
#' @export
#'
#' @examples
find_best_arch_model <- function(x, 
                                 type_models, 
                                 dist_to_use,
                                 max_lag_AR,
                                 max_lag_MA,
                                 max_lag_ARCH,
                                 max_lag_GARCH) {
  
  require(tidyr)
  
  df_grid <- expand_grid(type_models = type_models,
                         dist_to_use = dist_to_use,
                         arma_lag = 0:max_lag_AR,
                         ma_lag = 0:max_lag_MA,
                         arch_lag = 1:max_lag_ARCH,
                         garch_lag = 1:max_lag_GARCH)
  
  
  l_out <- pmap(.l = list(x = rep(list(x), nrow(df_grid)), 
                          type_model = df_grid$type_models,
                          type_dist = df_grid$dist_to_use,
                          lag_ar = df_grid$arma_lag,
                          lag_ma = df_grid$ma_lag,
                          lag_arch = df_grid$arch_lag,
                          lag_garch  = df_grid$garch_lag),
                do_single_garch)
  
  tab_out <- bind_rows(l_out)
  
  # find by AIC
  idx <- which.min(tab_out$AIC)
  best_aic <- tab_out[idx, ]
  
  # find by BIC
  idx <- which.min(tab_out$BIC)
  best_bic <- tab_out[idx, ]
  
  l_out <- list(best_aic = best_aic,
                best_bic = best_bic,
                tab_out = tab_out)
  
  return(l_out)
}

find_best_arch_model2 <- function(x, 
                                  type_models,
                                  sub_models,
                                  dist_to_use,
                                  max_lag_AR,
                                  max_lag_MA,
                                  max_lag_ARCH,
                                  max_lag_GARCH) {
  
  require(tidyr)
  
  df_grid <- expand_grid(type_models = type_models,
                         sub_models = sub_models,
                         dist_to_use = dist_to_use,
                         arma_lag = 0:max_lag_AR,
                         ma_lag = 0:max_lag_MA,
                         arch_lag = 1:max_lag_ARCH,
                         garch_lag = 1:max_lag_GARCH)
  
  
  l_out <- pmap(.l = list(x = rep(list(x), nrow(df_grid)), 
                          type_model = df_grid$type_models,
                          sub_model = df_grid$sub_models,
                          type_dist = df_grid$dist_to_use,
                          lag_ar = df_grid$arma_lag,
                          lag_ma = df_grid$ma_lag,
                          lag_arch = df_grid$arch_lag,
                          lag_garch  = df_grid$garch_lag),
                do_single_garch2)
  
  tab_out <- bind_rows(l_out)
  
  # find by AIC
  idx <- which.min(tab_out$AIC)
  best_aic <- tab_out[idx, ]
  
  # find by BIC
  idx <- which.min(tab_out$BIC)
  best_bic <- tab_out[idx, ]
  
  l_out <- list(best_aic = best_aic,
                best_bic = best_bic,
                tab_out = tab_out)
  
  return(l_out)
}
max_lag_AR <- 1
max_lag_MA <- 1
max_lag_ARCH <- 1
max_lag_GARCH <- 1
dist_to_use <- c('norm', 'std', 'ged', 'sstd', 'sged') # see rugarch::ugarchspec help for more
models_to_estimate <- c('sGARCH', 'iGARCH', 'eGARCH', 'gjrGARCH')# "fGARCH","fGARCH", "fGARCH") 
models_to_estimate2 <- "fGARCH"
sub_models_to_estimate2 <- c("NAGARCH", "AVGARCH","TGARCH")

library(tidyverse)
library(purrr)

out1 <- find_best_arch_model(x = R, 
                            type_models = models_to_estimate,
                            dist_to_use = dist_to_use,
                            max_lag_AR = max_lag_AR,
                            max_lag_MA = max_lag_MA,
                            max_lag_ARCH = max_lag_ARCH,
                            max_lag_GARCH = max_lag_GARCH)
out2 <- find_best_arch_model2(x = R,  
                             type_models = models_to_estimate2,
                             sub_models = sub_models_to_estimate2,
                             dist_to_use = dist_to_use,
                             max_lag_AR = max_lag_AR,
                             max_lag_MA = max_lag_MA,
                             max_lag_ARCH = max_lag_ARCH,
                             max_lag_GARCH = max_lag_GARCH)
# get table with estimation results
tab_out <- out1$tab_out
tab_out2 <- out2$tab_out
# clean tab_out
tab_out <- full_join(tab_out,tab_out2)

# pivot table to long format (better for plotting)
df_long <- tidyr::pivot_longer(data = tab_out %>%
                                 select(model_name,
                                        type_model,
                                        sub_model,
                                        type_dist,
                                        AIC, BIC),  cols = c('AIC', 'BIC'))

models_names <- unique(df_long$model_name)
best_models <- c(tab_out$model_name[which.min(tab_out$AIC)],
                 tab_out$model_name[which.min(tab_out$BIC)])

# figure out where is the best model
df_long <- df_long %>%
  mutate(order_model = if_else(model_name %in% best_models, 'Best Model', 'Not Best Model') ) 
df_long <- df_long %>%
  mutate(ARMA_order = substr(df_long$model_name, 1,9))
meanaic <- df_long %>% group_by(ARMA_order,type_dist, name) %>% summarise(mean_aic = mean(value)) %>% arrange(mean_aic)  
meanaicdf <- as.data.frame(meanaic)
ggplot(meanaicdf,aes(x = ARMA_order, y = mean_aic)) + geom_bar(stat="identity") + facet_wrap(~ type_dist)

# symmetric garches and symmetric dists
df_long1 <-  filter(df_long, type_model %in% c("sGARCH", "iGARCH") & type_dist %in% c("norm", "std", "ged"))
# symmetric garches and asymmetric dists with best symmetric
df_long2 <- filter(df_long, type_model %in% c("sGARCH", "iGARCH") & type_dist %in% c("std","sstd", "sged"))
# asymmetric garches and symmetric dists
df_long3 <-  filter(df_long, type_model %in% c("eGARCH", "gjrGARCH") & type_dist %in% c("norm", "std", "ged"))
# asymmetric garches and asymmetric dists with best symmetric
df_long4 <-  filter(df_long, type_model %in% c("eGARCH", "gjrGARCH") & type_dist %in% c("std", "sstd", "sged"))
# treshhold garches and symmetric dists
df_long5 <-  filter(df_long, type_model == "fGARCH" & type_dist %in% c("norm", "std", "ged") & !(sub_model %in% c("TGARCH", "AVGARCH") & type_dist == "ged"))
# treshhold garches and asymmetric dists with best symmetric
df_long6 <-  filter(df_long, type_model == "fGARCH" & type_dist %in% c("std", "sstd", "sged") & !(sub_model %in% c("TGARCH", "AVGARCH") & type_dist == "sged") )

# make table with best models
df_best_models <- df_long %>%
  group_by(name) %>%
  summarise(model_name = model_name[which.min(value)],
            value = value[which.min(value)],
            type_model = type_model[which.min(value)])

# plot results
p1 <- ggplot(df_long1 %>%
               arrange(type_model), 
             aes(x = reorder(model_name, 
                             order(type_model)),
                 y = value, 
                 shape = type_dist,
                 color = type_model)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  facet_wrap(~name, scales = 'free_x') + 
  geom_point(data = df_best_models, mapping = aes(x = reorder(model_name, 
                                                              order(type_model)),
                                                  y = value), 
             color = 'blue', size = 5, shape = 8) +
  labs(title = 'Selecting Garch Models by Fitness Criteria', 
       subtitle = 'The best model is the one with lowest AIC or BIC (with star)',
       x = '',
       y = 'Value of Fitness Criteria',
       shape = 'Type of Dist.',
       color = 'Type of Model') + 
  theme(legend.position = "right")
p1

p2 <- ggplot(df_long2 %>%
               arrange(type_model), 
             aes(x = reorder(model_name, 
                             order(type_model)),
                 y = value, 
                 shape = type_dist,
                 color = type_model)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  facet_wrap(~name, scales = 'free_x') + 
  geom_point(data = df_best_models, mapping = aes(x = reorder(model_name, 
                                                              order(type_model)),
                                                  y = value), 
             color = 'blue', size = 5, shape = 8) +
  labs(title = 'Selecting Garch Models by Fitness Criteria', 
       subtitle = 'The best model is the one with lowest AIC or BIC (with star)',
       x = '',
       y = 'Value of Fitness Criteria',
       shape = 'Type of Dist.',
       color = 'Type of Model') + 
  theme(legend.position = "right")
p2

p3 <- ggplot(df_long3 %>%
               arrange(type_model), 
             aes(x = reorder(model_name, 
                             order(type_model)),
                 y = value, 
                 shape = type_dist,
                 color = type_model)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  facet_wrap(~name, scales = 'free_x') + 
  geom_point(data = df_best_models, mapping = aes(x = reorder(model_name, 
                                                              order(type_model)),
                                                  y = value), 
             color = 'blue', size = 5, shape = 8) +
  labs(title = 'Selecting Garch Models by Fitness Criteria', 
       subtitle = 'The best model is the one with lowest AIC or BIC (with star)',
       x = '',
       y = 'Value of Fitness Criteria',
       shape = 'Type of Dist.',
       color = 'Type of Model') + 
  theme(legend.position = "right")
p3

p4 <- ggplot(df_long4 %>%
               arrange(type_model), 
             aes(x = reorder(model_name, 
                             order(type_model)),
                 y = value, 
                 shape = type_dist,
                 color = type_model)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  facet_wrap(~name, scales = 'free_x') + 
  geom_point(data = df_best_models, mapping = aes(x = reorder(model_name, 
                                                              order(type_model)),
                                                  y = value), 
             color = 'blue', size = 5, shape = 8) +
  labs(title = 'Selecting Garch Models by Fitness Criteria', 
       subtitle = 'The best model is the one with lowest AIC or BIC (with star)',
       x = '',
       y = 'Value of Fitness Criteria',
       shape = 'Type of Dist.',
       color = 'Type of Model') + 
  theme(legend.position = "right")
p4


p5 <- ggplot(df_long5 %>%
               arrange(type_model), 
             aes(x = reorder(model_name, 
                             order(type_model)),
                 y = value, 
                 shape = type_dist,
                 color = type_model)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  facet_wrap(~name, scales = 'free_x') + 
  geom_point(data = df_best_models, mapping = aes(x = reorder(model_name, 
                                                              order(type_model)),
                                                  y = value), 
             color = 'blue', size = 5, shape = 8) +
  labs(title = 'Selecting Garch Models by Fitness Criteria', 
       subtitle = 'The best model is the one with lowest AIC or BIC (with star)',
       x = '',
       y = 'Value of Fitness Criteria',
       shape = 'Type of Dist.',
       color = 'Type of Model') + 
  theme(legend.position = "right")
p5

p6 <- ggplot(df_long6 %>%
               arrange(type_model), 
             aes(x = reorder(model_name, 
                             order(type_model)),
                 y = value, 
                 shape = type_dist,
                 color = type_model)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  facet_wrap(~name, scales = 'free_x') + 
  geom_point(data = df_best_models, mapping = aes(x = reorder(model_name, 
                                                              order(type_model)),
                                                  y = value), 
             color = 'blue', size = 5, shape = 8) +
  labs(title = 'Selecting Garch Models by Fitness Criteria', 
       subtitle = 'The best model is the one with lowest AIC or BIC (with star)',
       x = '',
       y = 'Value of Fitness Criteria',
       shape = 'Type of Dist.',
       color = 'Type of Model') + 
  theme(legend.position = "right")
p6

