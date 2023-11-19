# Initial set up ----

rm(list = ls())

library(tidyverse)
library(sbt)
library(tmbstan)

theme_set(theme_bw())

options(mc.cores = parallel::detectCores()) # Specify the number of cores for MCMC

# Create data list ----

attach(data_csv1)

Data1 <- list(last_yr = 2022, age_increase_M = 25,
              length_m50 = 150, length_m95 = 180, 
              length_mean = length_mean, length_sd = length_sd, 
              catch = catch, catch_UA = catch_UA, 
              catch_UR_on = 0, catch_surf_case = 1, catch_LL1_case = 1, 
              scenarios_surf = scenarios_surface, scenarios_LL1 = scenarios_LL1,
              sel_min_age_f = c(2, 2, 2, 8, 6, 0), 
              sel_max_age_f = c(17, 9, 17, 22, 25, 7),
              sel_end_f = c(1, 0, 1, 1, 1, 0),
              sel_change_sd_fy = t(as.matrix(sel_change_sd[,-1])), 
              sel_smooth_sd_f = data_labrep1$sel.smooth.sd,
              hsp_switch = 1, HSPs = HSPs, hsp_false_negative = 0.7467647, 
              pop_switch = 1, POPs = POPs, 
              gt_switch = 1, GTs = GTs,
              cpue_switch = 1, cpue = cpue, cpue_a1 = 5, cpue_a2 = 17,
              aerial_switch = 4, aerial_surv = aerial_surv, 
              aerial_cov = aerial_cov, aerial_tau = data_labrep1$tau.aerial, 
              troll_switch = 0, troll = troll, 
              af = af, 
              lf = lf, lf_minbin = c(1, 1, 1, 11),
              tag_var_factor = 1.82, tag_switch = 1
)

Data <- get_data(data_in = Data1)

# Create parameter list ----

Params <- get_parameters()

# Use TMB's Map option to turn parameters on/off ----

Map <- list()
# Map[["par_log_B0"]] <- factor(NA) # est
Map[["par_log_psi"]] <- factor(NA)
# Map[["par_log_m0"]] <- factor(NA)
# Map[["par_log_m4"]] <- factor(NA) # est
# Map[["par_log_m10"]] <- factor(NA)
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

# Create the AD object ----

obj <- MakeADFun(data = Data, parameters = Params, map = Map, random = c(), 
                 hessian = TRUE, inner.control = list(maxit = 1000), DLL = "sbt")

# Specify parameter bounds ----

bnd <- get_bounds(obj = obj)

# Set up estimation ----

unique(names(obj$par)) # List of parameters that are "on"
check_bounds(opt = obj, lower = bnd$lower, upper = bnd$upper)

# Optimize ----

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, 
              lower = bnd$lower, upper = bnd$upper)

# Run MCMC ----

# control <- list(max_treedepth = 12, adapt_delta = 0.9)
control <- list(max_treedepth = 12)

if (FALSE) {
  Params$par_log_h <- log(0.55)
  Params$par_log_psi <- log(1.5)
  mcmc_g1 <- tmbstan(obj = obj, lower = bnd$lower, upper = bnd$upper,
                     init = rep(list(Params), 2), chains = 2, control = control)
  save(mcmc_g1, file = "mcmc_g1.rda")
  
  Params$par_log_h <- log(0.8)
  Params$par_log_psi <- log(2)
  mcmc_g2 <- tmbstan(obj = obj, lower = bnd$lower, upper = bnd$upper,
                     init = rep(list(Params), 2), chains = 2, control = control)
  save(mcmc_g2, file = "mcmc_g2.rda")
  
  Params$par_log_h <- log(0.55)
  Params$par_log_psi <- log(2)
  mcmc_g3 <- tmbstan(obj = obj, lower = bnd$lower, upper = bnd$upper,
                     init = rep(list(Params), 2), chains = 2, control = control)
  save(mcmc_g3, file = "mcmc_g3.rda")
  
  Params$par_log_h <- log(0.8)
  Params$par_log_psi <- log(1.5)
  mcmc_g4 <- tmbstan(obj = obj, lower = bnd$lower, upper = bnd$upper,
                     init = rep(list(Params), 2), chains = 2, control = control)
  save(mcmc_g4, file = "mcmc_g4.rda")
} else {
  load("mcmc_g1.rda")
  load("mcmc_g2.rda")
  load("mcmc_g3.rda")
  load("mcmc_g4.rda")
}

mcmc_grd <- list(mcmc_g1, mcmc_g2, mcmc_g3, mcmc_g4)

# Run grid ----

# M0 and M10 were free in MCMC, but here I fix for grid
Map[["par_log_m0"]] <- factor(NA)
Map[["par_log_m10"]] <- factor(NA)
names(Map)

Grid <- get_grid(par = Params)
grd <- run_grid(data = Data, grid = Grid, bounds = bnd, map = Map)

# Plots ----

pars <- c("lp__", "par_log_B0", 
          "par_log_m0", "par_log_m10",
          "par_log_m4", "par_log_m30", "par_log_cpue_q", 
          "par_rdev_y[1]", "par_rdev_y[45]", "par_rdev_y[92]", 
          "par_logit_hstar_i[1]", "par_logit_hstar_i[17]",
          "par_sels_init_i[1]", "par_sels_init_i[78]", 
          "par_sels_change_i[1]", "par_sels_change_i[1132]")

stan_trace(object = mcmc_grd[[1]], pars = pars)
stan_hist(object = mcmc_grd[[1]], pars = pars)
stan_dens(object = mcmc_grd[[1]], pars = pars)
stan_plot(object = mcmc_grd[[1]], pars = "par_logit_hstar_i")

plot_natural_mortality(data = Data, object = obj, posterior = mcmc_grd[[1]])

# only this function has posterior and grid so far
plot_biomass_spawning(data = Data, object = obj, grid = grd, posterior = mcmc_grd, relative = FALSE)
ggsave(filename = "biomass_spawning_grid_mcmc_grid.png", width = 7, height = 4)

# plot_biomass_spawning(data = Data, object = obj, grid = grd, posterior = mcmc_grd[[1]], relative = FALSE)
# plot_biomass_spawning(data = Data, object = obj, grid = grd, posterior = mcmc_grd[[1]], relative = FALSE)
# ggsave(filename = "biomass_spawning_grid_mcmc.png", width = 7, height = 4)

# Example of model averaging ----

if (FALSE) {
  loo_grid1 <- get_loo(data = Data, object = obj, posterior = mcmc_grd[[1]])
  loo_grid2 <- get_loo(data = Data, object = obj, posterior = mcmc_grd[[2]])
  loo_grid3 <- get_loo(data = Data, object = obj, posterior = mcmc_grd[[3]])
  loo_grid4 <- get_loo(data = Data, object = obj, posterior = mcmc_grd[[4]])
  save(loo_grid1, loo_grid2, loo_grid3, loo_grid4, file = "loo_grid.rda")
} else {
  load("loo_grid.rda")
}

print(loo_grid1)
print(loo_grid4)
plot_loo(x = loo_grid1)
plot_loo(x = loo_grid4)
loo::loo_model_weights(x = list(loo_grid1, loo_grid2, loo_grid3, loo_grid4), method = "pseudobma")
# loo::loo_model_weights(x = list(loo_grid1, loo_grid2, loo_grid3, loo_grid4), method = "stacking")
# mcmc2 <- sflist2stanfit(sflist = list(mcmc1, mcmc1))
