# Initial set up ----

rm(list = ls())

library(tidyverse)
library(sbt)

theme_set(theme_bw())

# Create data list ----

attach(data_csv1)
lr <- data_csv1

Data <- list(last_yr = 2022, age_increase_M = 25,
             length_m50 = 150, length_m95 = 180, 
             catch_UR_on = 0, catch_surf_case = 1, catch_LL1_case = 1, 
             scenarios_surf = scenarios_surface, scenarios_LL1 = scenarios_LL1,
             sel_min_age_f = c(2, 2, 2, 8, 6, 0), 
             sel_max_age_f = c(17, 9, 17, 22, 25, 7),
             sel_end_f = c(1, 0, 1, 1, 1, 0),
             sel_change_sd_fy = t(as.matrix(sel_change_sd[,-1])),
             sel_smooth_sd_f = data_labrep1$sel.smooth.sd,
             removal_switch = c(0, 0, 0, 0, 0, 0), # 0=standard removals, 1=direct removals
             pop_switch = 1, 
             hsp_switch = 1, hsp_false_negative = 0.7467647, 
             gt_switch = 1,
             cpue_switch = 1, cpue_a1 = 5, cpue_a2 = 17,
             aerial_switch = 4, aerial_tau = 0.3, 
             troll_switch = 1, 
             af_switch = 3, # 0=multinomial, 1=Dirichlet, 2=Dirichlet-multinomial, 3=old
             lf_switch = 3, lf_minbin = c(1, 1, 1, 11),
             tag_switch = 1, tag_var_factor = 1.82
)

Data <- get_data(data_in = Data)

# Create parameter list ----

Params <- list(par_log_B0 = data_par1$ln_B0,
               par_log_psi = log(data_par1$psi),
               par_log_m0 = log(data_par1$m0), 
               par_log_m4 = log(data_par1$m4),
               par_log_m10 = log(data_par1$m10), 
               par_log_m30 = log(data_par1$m30),
               par_log_h = log(data_par1$steep), 
               par_log_sigma_r = log(data_par1$sigma_r), 
               par_log_cpue_q = data_par1$lnq,
               par_log_cpue_sigma = log(data_par1$sigma_cpue),
               par_log_cpue_omega = log(data_par1$cpue_omega),
               par_log_aerial_tau = log(data_par1$tau_aerial),
               par_log_aerial_sel = data_par1$ln_sel_aerial,
               par_log_troll_tau = log(data_par1$tau_troll),
               par_log_hsp_q = data_par1$lnqhsp, 
               # par_logit_hstar_i = qlogis(exp(data_par1$par_log_hstar_i)),
               par_log_tag_H_factor = log(data_par1$tag_H_factor),
               par_log_af_alpha = c(1, 1),
               par_log_lf_alpha = rep(1, 4),
               par_rdev_y = data_par1$Reps,
               par_sels_init_i = data_par1$par_sels_init_i,
               par_sels_change_i = data_par1$par_sels_change_i
)

# Use TMB's Map option to turn parameters on/off ----

Map <- list()
# Map[["par_log_B0"]] <- factor(NA) # est
Map[["par_log_psi"]] <- factor(NA)
Map[["par_log_m0"]] <- factor(NA)
# Map[["par_log_m4"]] <- factor(NA) # est
Map[["par_log_m10"]] <- factor(NA)
# Map[["par_log_m30"]] <- factor(NA) # est
Map[["par_log_h"]] <- factor(NA)
Map[["par_log_sigma_r"]] <- factor(NA)
# Map[["par_rdev_y"]] <- factor(rep(NA, length(Params$par_rdev_y))) # est
# Map[["par_sels_init_i"]] <- factor(rep(NA, length(Params$par_sels_init_i))) # est
# Map[["par_sels_change_i"]] <- factor(rep(NA, length(Params$par_sels_change_i))) # est
# Map[["par_log_cpue_q"]] <- factor(NA) # est
Map[["par_log_cpue_sigma"]] <- factor(NA)
Map[["par_log_cpue_omega"]] <- factor(NA)
Map[["par_log_aerial_tau"]] <- factor(NA)
Map[["par_log_aerial_sel"]] <- factor(rep(NA, 2))
Map[["par_log_troll_tau"]] <- factor(NA)
Map[["par_log_hsp_q"]] <- factor(NA)
# Map[["par_logit_hstar_i"]] <- factor(rep(NA, 17)) # est
Map[["par_log_tag_H_factor"]] <- factor(NA)

# Specify the random effects
Random <- c()

# Create the AD object ----

obj <- MakeADFun(data = Data, parameters = Params, map = Map, random = Random, 
                 hessian = TRUE, inner.control = list(maxit = 1000), DLL = "sbt")

# Set up estimation ----

unique(names(obj$par)) # List of parameters that are "on"
bnd <- get_bounds(obj = obj)
check_bounds(opt = obj, lower = bnd$lower, upper = bnd$upper)

# Optimize ----

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, 
              lower = bnd$lower, upper = bnd$upper,
              control = list(eval.max = 1000, iter.max = 1000))

# Run grid ----

Grid <- get_grid(par = Params)
control <- list(eval.max = 10000, iter.max = 10000)
grd <- run_grid(data = Data, grid = Grid, bounds = bnd, map = Map, 
                control = control)

# Next step ----

attach(data_csv1)

Data <- list(last_yr = 2022, age_increase_M = 25, 
             length_m50 = 150, length_m95 = 180, 
             catch_UR_on = 0, catch_surf_case = 1, catch_LL1_case = 1, 
             scenarios_surf = scenarios_surface, scenarios_LL1 = scenarios_LL1,
             sel_min_age_f = c(2, 2, 2, 8, 6, 0), 
             sel_max_age_f = c(17, 9, 17, 22, 25, 7),
             sel_end_f = c(1, 0, 1, 1, 1, 0),
             sel_change_sd_fy = t(as.matrix(sel_change_sd[,-1])), 
             sel_smooth_sd_f = data_labrep1$sel.smooth.sd,
             removal_switch = c(0, 0, 0, 0, 0, 0), # 0=standard, 1=direct
             pop_switch = 1, 
             hsp_switch = 1, hsp_false_negative = 0.7467647, 
             gt_switch = 1, 
             cpue_switch = 1, cpue_a1 = 5, cpue_a2 = 17,
             aerial_switch = 4, aerial_tau = data_labrep1$tau.aerial, 
             troll_switch = 0,
             af_switch = 3,
             lf_switch = 3, lf_minbin = c(1, 1, 1, 11),
             tag_switch = 1, tag_var_factor = 1.82
)

Data <- get_data(data_in = Data)

Params <- list(par_log_B0 = data_par1$ln_B0, 
               par_log_psi = log(1.5),
               par_log_m0 = log(data_par1$m0), 
               par_log_m4 = log(data_par1$m4),
               par_log_m10 = log(data_par1$m10), 
               par_log_m30 = log(data_par1$m30),
               par_log_h = log(0.72), 
               par_log_sigma_r = log(data_labrep1$sigma.r), 
               par_log_cpue_q = data_par1$lnq,
               par_log_cpue_sigma = log(data_par1$sigma_cpue),
               par_log_cpue_omega = log(data_par1$cpue_omega),
               par_log_aerial_tau = log(data_par1$tau_aerial),
               par_log_aerial_sel = data_par1$ln_sel_aerial,
               par_log_troll_tau = log(data_par1$tau_troll),
               par_log_hsp_q = data_par1$lnqhsp, 
               par_log_tag_H_factor = log(data_par1$tag_H_factor),
               par_log_af_alpha = c(1, 1),
               par_log_lf_alpha = c(1, 1, 1, 1),
               par_rdev_y = data_par1$Reps,
               par_sels_init_i = data_par1$par_sels_init_i, 
               par_sels_change_i = data_par1$par_sels_change_i
)

Map <- list()
Map[["par_log_psi"]] <- factor(NA)
Map[["par_log_h"]] <- factor(NA)
Map[["par_log_sigma_r"]] <- factor(NA)
Map[["par_log_cpue_sigma"]] <- factor(NA)
Map[["par_log_cpue_omega"]] <- factor(NA)
Map[["par_log_aerial_tau"]] <- factor(NA)
Map[["par_log_aerial_sel"]] <- factor(rep(NA, 2))
Map[["par_log_troll_tau"]] <- factor(NA)
Map[["par_log_hsp_q"]] <- factor(NA)
Map[["par_log_tag_H_factor"]] <- factor(NA)
Map[["par_log_af_alpha"]] <- factor(rep(NA, 2))
Map[["par_log_lf_alpha"]] <- factor(rep(NA, 4))

obj <- MakeADFun(data = Data, parameters = Params, map = Map, random = c(), 
                 hessian = TRUE, inner.control = list(maxit = 50), DLL = "sbt")

n_grid <- length(grd)
gdf <- matrix(NA, nrow = n_grid, ncol = 5)
colnames(gdf) <- c("nll", "m0", "m10", "h", "psi")
for (i in 1:n_grid) {
  gdf[i, 1] <- grd[[i]]$report()$nll
  gdf[i, 2] <- grd[[i]]$report()$par_m0
  gdf[i, 3] <- grd[[i]]$report()$par_m10
  gdf[i, 4] <- grd[[i]]$report()$par_h
  gdf[i, 5] <- grd[[i]]$report()$par_psi
}

dfstrat <- data.frame(gdf) %>% 
  mutate(ll = -(nll - min(nll))) %>%
  group_by(h, psi) %>% 
  summarise(sum_ll = sum(exp(ll)))
head(dfstrat)

df <- data.frame(gdf) %>%
  mutate(ll = -(nll - min(nll))) %>%
  full_join(dfstrat) %>%
  mutate(prob = exp(ll) / sum_ll)
head(df)
range(df$prob)

df %>% group_by(h, psi) %>% summarise(sum_prob = sum(prob))

n_samps <- 2000
grid_ints <- sample(x = 1:108, size = n_samps, replace = TRUE, prob = df$prob)
sbio <- matrix(NA, n_samps, 93)
for (i in 1:n_samps) {
  ii <- grid_ints[i]
  sbio[i,] <- grd[[ii]]$report()$spawning_biomass_y
}
df_sbio <- reshape2::melt(sbio) %>% rename(cell = Var1, year = Var2)
y1 <- df_sbio %>% filter(year == 1) %>% select(-year) %>% rename(y1 = value)
df_sbio <- df_sbio %>% full_join(y1) %>% mutate(value = value / y1)

for (i in 1:12) {
  load(paste0("/home/darcy/Projects/CCSBT/sbt/", "mcmc", i, ".rda"))
}
mc_grid <- list()
mc_grid[[1]] <- mcmc1
mc_grid[[2]] <- mcmc2
mc_grid[[3]] <- mcmc3
mc_grid[[4]] <- mcmc4
mc_grid[[5]] <- mcmc5
mc_grid[[6]] <- mcmc6
mc_grid[[7]] <- mcmc7
mc_grid[[8]] <- mcmc8
mc_grid[[9]] <- mcmc9
mc_grid[[10]] <- mcmc10
mc_grid[[11]] <- mcmc11
mc_grid[[12]] <- mcmc12

post <- list()
for (i in 1:12) {
  post[[i]] <- rstan::extract(mc_grid[[i]], pars = "lp__", permuted = FALSE, include = FALSE)
}
pp <- expand.grid(h = c(0.55, 0.63, 0.72, 0.8), psi = c(1.5, 1.75, 2))
n_iter = 2000
ex <- expand.grid(run = 1:12, chain = 1:2, iter = 1:1000)
mc_ints <- sample(x = 1:nrow(ex), size = n_samps, replace = TRUE)
mc_sbio <- matrix(NA, n_samps, 93)
for (i in 1:n_samps) {
  ii <- mc_ints[i]
  irun <- ex[ii,]$run
  ichain <- ex[ii,]$chain
  iter <- ex[ii,]$iter
  Params$par_log_h <- log(pp$h[irun])
  Params$par_log_psi <- log(pp$psi[irun])
  objj <- MakeADFun(data = Data, parameters = Params, map = Map, random = c(), DLL = "sbt")
  # unique(names(objj$par))
  # length(names(objj$par))
  # dim(post[[ii]])
  mc_sbio[i,] <- objj$report(par = post[[irun]][iter, ichain,])$spawning_biomass_y
}
df_sbio2 <- reshape2::melt(mc_sbio) %>% rename(cell = Var1, year = Var2)
y1 <- df_sbio2 %>% filter(year == 1) %>% select(-year) %>% rename(y1 = value)
df_sbio2 <- df_sbio2 %>% full_join(y1) %>% mutate(value = value / y1)

ggplot(data = df_sbio) +
  geom_line(data = df_sbio2, aes(x = factor(year), y = value, group = cell), color = "red", alpha = 0.15) +
  geom_line(aes(x = factor(year), y = value, group = cell), alpha = 0.15)

df <- bind_rows(df_sbio %>% mutate(Grid = "MPD"),
                df_sbio2 %>% mutate(Grid = "MCMC"))
ggplot(data = df %>% filter(year > 50)) +
  geom_boxplot(aes(x = factor(year + 1930), y = value, color = Grid, fill = Grid), alpha = 0.15) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "Year", y = "Relative TRO") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
# save(grd, file = "grd.rda")
# load("grd.rda")

# Plot ----

plot_biomass_spawning(data = Data, object = grd[[1]])
plot_biomass_spawning(data = Data, object = grd[[1]], grid = grd)
ggsave(filename = "biomass_spawning_grid.png", width = 7, height = 4)
