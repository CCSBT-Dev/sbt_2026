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

Params <- list(par_log_B0 = data_par1$ln_B0, 
               par_log_psi = log(data_par1$psi),
               par_log_m0 = log(data_par1$m0), 
               par_log_m4 = log(data_par1$m4),
               par_log_m10 = log(data_par1$m10), 
               par_log_m30 = log(data_par1$m30),
               par_log_h = log(data_par1$steep), 
               par_log_sigma_r = log(data_labrep1$sigma.r),
               par_rdev_y = data_par1$Reps,
               par_sels_init_i = data_par1$par_sels_init_i,
               par_sels_change_i = data_par1$par_sels_change_i,
               par_log_cpue_q = data_par1$lnq,
               par_log_cpue_sigma = log(data_par1$sigma_cpue),
               par_log_cpue_omega = log(data_par1$cpue_omega),
               par_log_aerial_tau = log(data_par1$tau_aerial),
               par_log_aerial_sel = data_par1$ln_sel_aerial,
               par_log_troll_tau = log(data_par1$tau_troll),
               par_log_hsp_q = data_par1$lnqhsp, 
               par_logit_hstar_i = qlogis(exp(data_par1$par_log_hstar_i)),
               par_log_tag_H_factor = log(data_par1$tag_H_factor)
)

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

# Set up estimation ----

unique(names(obj$par)) # List of parameters that are "on"
bnd <- get_bounds(obj = obj)
check_bounds(opt = obj, lb = bnd$lb, ub = bnd$ub)

# Optimize ----

opt <- nlminb(start = obj$par, objective = obj$fn, gr = obj$gr, lower = bnd$lb, upper = bnd$ub)

# Run MCMC ----

mcmc1 <- tmbstan(obj = obj, lower = bnd$lb, upper = bnd$ub, 
                 init = rep(list(Params), 2), chains = 2, control = list(max_treedepth = 12))
# mcmc2 <- tmbstan(obj = obj, lower = Lwr, upper = Upr, init = list(Params), chains = 1, laplace = TRUE)
# get_stancode(mcmc1)
save(obj, mcmc1, file = "mcmc1.rda")
load("mcmc1.rda")

# Run grid ----

# M0 and M10 were free in MCMC, but fix for grid
Map[["par_log_m0"]] <- NULL
Map[["par_log_m10"]] <- NULL
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

stan_trace(object = mcmc1, pars = pars)
stan_hist(object = mcmc1, pars = pars)
stan_dens(object = mcmc1, pars = pars)
stan_plot(object = mcmc1, pars = "par_logit_hstar_i")

plot_natural_mortality(data = Data, object = obj, posterior = mcmc1)

# only this function has posterior and grid so far
plot_biomass_spawning(data = Data, object = obj, grid = grd, posterior = mcmc1)
