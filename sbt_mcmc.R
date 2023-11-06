# Load required R libraries ----

rm(list = ls())

library(TMB)
library(tmbstan)
library(tidyverse)
library(sbt)

theme_set(theme_bw())

options(mc.cores = parallel::detectCores())

# Compile the model ----

compile("src/sbt_v100.cpp")
dyn.load(dynlib("src/sbt_v100"))

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
             HSPs = HSPs, hsp_false_negative = 0.7467647, 
             POPs = POPs, GTs = GTs,
             aerial_surv = aerial_surv, aerial_cov = aerial_cov, 
             aerial_tau = data_labrep1$tau.aerial, aerial_switch = 4,
             cpue = cpue, cpue_a1 = 4, cpue_a2 = 18,
             cpue_adjust = c(1, 1.005, 1.01, 1.015, 1.02, 1.025, 1.03, 1.035, 
                             1.04, 1.045, 1.05, 1.055, 1.06, 1.065, 1.07, 1.075, 
                             1.08, 1.085, 1.09, 1.095, 1.1, 1.105, 1.11, 1.115, 
                             1.12, 1.125, 1.13, 1.135, 1.14, 1.145, 1.15, 1.155, 
                             1.16, 1.165, 1.17, 1.175, 1.18, 1.185, 1.19, 1.195,
                             1.2, 1.205, 1.21, 1.215, 1.22, 1.225, 1.23, 1.235, 
                             1.24, 1.245, 1.25, 1.255, 1.26, 1.265),
             af = af, lf = lf, lf_minbin = c(1, 1, 1, 11),
             tag_rel_min_age = c(2, 1, 1, 1, 1, 1),
             tag_rel_max_age = c(3, 3, 3, 3, 3, 3),
             tag_recap_max_age = c(7, 7, 6, 5, 4, 3),
             tag_shed_immediate = c(0.9737, 0.9608, 1, 1, 0.9342, 0.9666),
             tag_shed_continuous = c(0.0391, 0.0492, 0.0672, 0.0925, 0.0885, 0.1601),
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
               par_log_cpue_q = data_par1$`lnq:`,
               par_log_cpue_sigma = log(data_par1$sigma_cpue),
               par_log_cpue_omega = log(data_par1$cpue_omega),
               par_log_aerial_tau = log(data_par1$tau_aerial),
               par_log_aerial_sel = data_par1$ln_sel_aerial,
               par_log_hsp_q = data_par1$lnqhsp, 
               par_log_hstar_i = data_par1$par_log_hstar_i
)

# Use TMB's Map option to turn parameters on/off ----

Map <- list()
# Map[["par_rdev_y"]] <- factor(rep(NA, 92))
# Map[["par_log_B0"]] <- factor(NA)
Map[["par_log_psi"]] <- factor(NA)
Map[["par_log_m0"]] <- factor(NA)
# Map[["par_log_m4"]] <- factor(NA)
Map[["par_log_m10"]] <- factor(NA)
# Map[["par_log_m30"]] <- factor(NA)
Map[["par_log_h"]] <- factor(NA)
Map[["par_log_sigma_r"]] <- factor(NA)
Map[["par_sels_init_i"]] <- factor(rep(NA, 83))
Map[["par_sels_change_i"]] <- factor(rep(NA, 1132))
Map[["par_log_cpue_sigma"]] <- factor(NA)
Map[["par_log_cpue_omega"]] <- factor(NA)
Map[["par_log_aerial_tau"]] <- factor(NA)
Map[["par_log_aerial_sel"]] <- factor(rep(NA, 2))
Map[["par_log_hsp_q"]] <- factor(NA)
Map[["par_log_hstar_i"]] <- factor(rep(NA, 17))

# Specify the random effects ----

Random <- c()

# Create the AD object ----

obj <- MakeADFun(data = Data, parameters = Params, map = Map, random = Random, 
                 hessian = TRUE, inner.control = list(maxit = 50), DLL = "sbt_v100")

# List of parameters that are "on"
unique(names(obj$par))

# Set up estimation
# newtonOption(smartsearch = TRUE)
obj$fn(obj$par)
obj$gr(obj$par)
obj$control <- list(trace = 100)
ConvergeTol <- 2 # 1:Normal; 2:Strong
#obj$env$inner.control$step.tol <- c(1e-12,1e-15)[ConvergeTol] # Default : 1e-8 # Change in parameters limit inner optimization
#obj$env$inner.control$tol10 <- c(1e-8,1e-12)[ConvergeTol]  # Default : 1e-3 # Change in pen.like limit inner optimization
#obj$env$inner.control$grad.tol <- c(1e-12,1e-15)[ConvergeTol] # # Default : 1e-8 # Maximum gradient limit inner optimization
summary(obj)

Lwr <- rep(-Inf, length(obj$par))
Upr <- rep(Inf, length(obj$par))
Lwr[grep("par_rdev_y", names(obj$par))] <- rep(-5, length(Params$par_rdev_y))
Upr[grep("par_rdev_y", names(obj$par))] <- rep(5, length(Params$par_rdev_y))
Lwr[grep("par_log_psi", names(obj$par))] <- log(1.499999)
Upr[grep("par_log_psi", names(obj$par))] <- log(2.000001)
Lwr[grep("par_log_h", names(obj$par))] <- log(0.21)
Upr[grep("par_log_h", names(obj$par))] <- log(1.0)
Lwr[grep("par_log_sigma_r", names(obj$par))] <- log(0.4)
Upr[grep("par_log_sigma_r", names(obj$par))] <- log(2.0)
Lwr[grep("par_log_m0", names(obj$par))] <- log(0.2)
Upr[grep("par_log_m0", names(obj$par))] <- log(0.55)
Lwr[grep("par_log_m10", names(obj$par))] <- log(0.029)
Upr[grep("par_log_m10", names(obj$par))] <- log(0.21)
Lwr[grep("par_log_m30", names(obj$par))] <- log(0.2)
Upr[grep("par_log_m30", names(obj$par))] <- log(0.7)
cbind(Lwr, obj$par, Upr)

# Optimize ----

opt <- nlminb(start = obj$par, objective = obj$fn, gr = obj$gr, upper = Upr, lower = Lwr, control = list(eval.max = 1e4, iter.max = 1e4, rel.tol = c(1e-10, 1e-8)[ConvergeTol], trace = 1))
opt[["final_gradient"]] <- obj$gr(opt$par)
Diag <- obj$report()
Report <- sdreport(obj)
print(Report$pdHess) # Is the fit positive definite Hessian?

# MCMC ----

mcmc1 <- tmbstan(obj = obj, lower = Lwr, upper = Upr, init = list(Params), chains = 2)

traceplot(mcmc1, pars = c("par_log_B0", "par_log_m4", "par_log_m30", "par_log_cpue_q", "par_rdev_y[1]", "lp__"), inc_warmup = FALSE)
