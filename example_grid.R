# Initial set up ----

rm(list = ls())

library(tidyverse)
library(sbt)

theme_set(theme_bw())

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

unique(names(obj$par)) # List of parameters that are "on"

bnd <- get_bounds(obj = obj)
check_bounds(opt = obj, lb = bnd$lb, ub = bnd$ub)

# Optimize ----

opt <- nlminb(start = obj$par, objective = obj$fn, gr = obj$gr, lower = bnd$lb, upper = bnd$ub)

# Run grid ----

g_pars <- get_grid(par = Params)

library(foreach)
library(doParallel)
n_cores <- parallel::detectCores() - 1
# create the cluster
my.cluster <- parallel::makeCluster(
  n_cores, 
  type = "PSOCK"
)
# check cluster definition (optional)
print(my.cluster)
# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
# check if it is registered (optional)
foreach::getDoParRegistered()
# how many workers are available? (optional)
foreach::getDoParWorkers()

g_obj <- foreach(i = 1:length(g_pars), .packages = "sbt") %dopar% {
  obj1 <- MakeADFun(data = Data, parameters = g_pars[[i]], map = Map, random = Random, 
                    hessian = TRUE, inner.control = list(maxit = 1000), DLL = "sbt")
  opt1 <- nlminb(start = obj1$par, objective = obj1$fn, gr = obj1$gr, lower = bnd$lb, upper = bnd$ub)
  obj1
}

parallel::stopCluster(cl = my.cluster)

# save(g_obj, file = "g_obj.rda")
# load("g_obj.rda")

# Plot ----

plot_biomass_spawning(data = Data, object = g_obj[[1]])
plot_biomass_spawning(data = Data, object = g_obj[[1]], grid = g_obj)
