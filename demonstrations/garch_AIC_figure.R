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
df_long <- df_long %>% filter(value < 4)
df_long <- df_long %>% mutate(model_name = gsub(pattern = "fGARCH", "",model_name))
df_long <- df_long %>% mutate(model_name = gsub(pattern = "sstd", "ST",model_name))
df_long <- df_long %>% mutate(model_name = gsub(pattern = "std", "T",model_name))
df_long <- df_long %>% mutate(model_name = gsub(pattern = "sged", "SGED",model_name))
df_long <- df_long %>% mutate(model_name = gsub(pattern = "ged", "GED",model_name))
df_long <- df_long %>% mutate(model_name = gsub(pattern = "norm", "NORM",model_name))
df_long <- df_long %>% mutate(type_dist = gsub(pattern = "sstd", "ST",type_dist))
df_long <- df_long %>% mutate(type_dist = gsub(pattern = "std", "T",type_dist))
df_long <- df_long %>% mutate(type_dist = gsub(pattern = "sged", "SGED",type_dist))
df_long <- df_long %>% mutate(type_dist = gsub(pattern = "ged", "GED",type_dist))
df_long <- df_long %>% mutate(type_dist = gsub(pattern = "norm", "NORM",type_dist))

# meanaic <- df_long %>% group_by(ARMA_order,type_dist, name) %>% summarise(mean_aic = mean(value)) %>% arrange(mean_aic)  
# meanaicdf <- as.data.frame(meanaic)
# ggplot(meanaicdf,aes(x = ARMA_order, y = mean_aic)) + geom_bar(stat="identity") + facet_wrap(~ type_dist)
df_long <- df_long %>% filter(name != "BIC")

nrow(df_long)

# symmetric garches and symmetric dists
df_long1 <-  filter(df_long, type_model %in% c("sGARCH", "iGARCH") & type_dist %in% c("NORM", "T", "GED")) 
# symmetric garches and asymmetric dists with best symmetric
df_long2 <- filter(df_long, type_model %in% c("sGARCH", "iGARCH") & type_dist %in% c("ST", "SGED","T"))
# asymmetric garches and symmetric dists
df_long3 <-  filter(df_long, type_model %in% c("eGARCH", "gjrGARCH") & type_dist %in% c("NORM", "T", "GED"))
# asymmetric garches and asymmetric dists with best symmetric
df_long4 <-  filter(df_long, type_model %in% c("eGARCH", "gjrGARCH") & type_dist %in% c("ST", "SGED","T"))
# treshhold garches and symmetric dists
df_long5 <-  filter(df_long, type_model == "fGARCH" & type_dist %in% c("NORM", "T", "GED")) # & !(sub_model %in% c("TGARCH", "AVGARCH") & type_dist == "GED"))
# treshhold garches and asymmetric dists with best symmetric
df_long6 <-  filter(df_long, type_model == "fGARCH" & type_dist %in% c("ST", "SGED","T"))# & !(sub_model %in% c("TGARCH", "AVGARCH") & type_dist == "sged") )

# make table with best models
df_best_models <- df_long %>%
  group_by(name) %>%
  summarise(model_name = model_name[which.min(value)],
            value = value[which.min(value)],
            type_model = type_model[which.min(value)])


# plot results
p1 <- ggplot(df_long1 %>% arrange(type_model), 
             aes(x = reorder(model_name, order(type_model)), 
                 y = value, 
                 shape = type_dist)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  scale_shape_manual(values = c(15, 16, 17)) +
  theme_bw() +
  facet_wrap(~name, scales = 'free_x')+
  geom_point(mapping = aes(x = reorder(model_name, order(type_model)),y = value)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))+
  labs(title = 'Symmetric Garch Models and distributions', 
       x = '',
       y = '',
       shape = 'Type of Dist.') + 
  theme(legend.position = "right") 
pdf(file = "figures/aicfigures/symmetric aics.pdf", width =9.615, height= 9.23)
p1
dev.off()

p2 <- ggplot(df_long2 %>% arrange(type_model), 
             aes(x = reorder(model_name, order(type_model)), 
                 y = value, 
                 shape = type_dist)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  theme_bw() +
  scale_shape_manual(values = c(7, 2, 17)) +
  facet_wrap(~name, scales = 'free_x')+
  geom_point(mapping = aes(x = reorder(model_name, order(type_model)),y = value)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))+
  labs(title = 'Symmetric Garch Models and other distributions', 
       x = '',
       y = '',
       shape = 'Type of Dist.') + 
  theme(legend.position = "right") 
pdf(file = "figures/aicfigures/symmetric aics2.pdf", width =9.615, height= 9.23)
p2
dev.off()

p3 <- ggplot(df_long3 %>% arrange(type_model), 
             aes(x = reorder(model_name, order(type_model)), 
                 y = value, 
                 shape = type_dist)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  scale_shape_manual(values = c(15, 16, 17)) +
  theme_bw() +
  facet_wrap(~name, scales = 'free_x')+
  geom_point(mapping = aes(x = reorder(model_name, order(type_model)),y = value)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))+
  labs(title = 'Asymmetric Garch Models and symmetric distributions', 
       x = '',
       y = '',
       shape = 'Type of Dist.') + 
  theme(legend.position = "right")
pdf(file = "figures/aicfigures/asymmetric aics.pdf", width =9.615, height= 9.23)
p3
dev.off()


p4 <- ggplot(df_long4 %>% arrange(type_model), 
             aes(x = reorder(model_name, order(type_model)), 
                 y = value, 
                 shape = type_dist)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  theme_bw() +
  scale_shape_manual(values = c(7, 2, 17)) +
  facet_wrap(~name, scales = 'free_x')+
  geom_point(mapping = aes(x = reorder(model_name, order(type_model)),y = value)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))+
  labs(title = 'Asymmetric Garch Models and other distributions', 
       x = '',
       y = '',
       shape = 'Type of Dist.') + 
  theme(legend.position = "right")
pdf(file = "figures/aicfigures/asymmetric aics2.pdf", width =9.615, height= 9.23)
p4
dev.off()


p5 <- ggplot(df_long5 %>% arrange(type_model), 
             aes(x = reorder(model_name, order(type_model)), 
                 y = value, 
                 shape = type_dist)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  theme_bw() +
  scale_shape_manual(values = c(15, 16, 17)) +
  facet_wrap(~name, scales = 'free_x')+
  geom_point(mapping = aes(x = reorder(model_name, order(type_model)),y = value)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))+
  labs(title = 'Family Garch Models and symmetric distributions', 
       x = '',
       y = '',
       shape = 'Type of Dist.')+ 
  theme(legend.position = "right")
pdf(file = "figures/aicfigures/family aics.pdf", width =9.615, height= 9.23)
p5
dev.off()
# for the family garch model, the GED was very poorly for the TGARCH 

p6 <- ggplot(df_long6 %>% arrange(type_model), 
             aes(x = reorder(model_name, order(type_model)), 
                 y = value, 
                 shape = type_dist)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip()  + 
  scale_shape_manual(values = c(7, 2, 17)) +
  theme_bw() +
  facet_wrap(~name, scales = 'free_x')+
  geom_point(mapping = aes(x = reorder(model_name, order(type_model)),y = value)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))+
  labs(title = 'Family Garch Models and other distributions', 
       x = '',
       y = '',
       shape = 'Type of Dist.') + 
  theme(legend.position = "right")
pdf(file = "figures/aicfigures/family aics2.pdf", width =9.615, height= 9.23)
p6
dev.off()
# for the family garch models, the SGED was very poorly for the TGARCH 

