# Make sure R is clean
rm(list = ls())

# Load required R libraries
library(TMB)
# library(tmbstan)
library(tidyverse)
library(sbt)

theme_set(theme_bw())

# Compile the model
# Load the model
# dyn.unload(dynlib("src/sbt"))
compile("src/sbt_v100.cpp")
dyn.load(dynlib("src/sbt_v100"))
# compile(file = "src/sbt_model.cpp", flags = "-O0 -g")
# gdbsource(file = "sbt_tmp.R")

# Create data list ----

attach(data_csv1)

Data1 <- list(last_yr = 2022, age_increase_M = 25,
              length_m50 = 150, length_m95 = 180, 
              length_mean = length_mean, length_sd = length_sd, 
              catch = catch, catch_UA = catch_UA, 
              catch_UR_on = 0, catch_surf_case = 1, catch_LL1_case = 1, 
              scenarios_surf = scenarios_surface, scenarios_LL1 = scenarios_LL1,
              sel_min_age_f = c(2, 2, 2, 8, 6, 0), 
              # sel_max_age_f = c(17, 9, 17, 17, 25, 7),
              # sel_end_f = c(0, 0, 0, 0, 1, 0),
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
              troll_switch = 1, troll = troll, 
              af = af, 
              lf = lf, lf_minbin = c(1, 1, 1, 11),
              tag_var_factor = 1.82, tag_switch = 1
)

# library(testthat)
# source("R/get-data.R")
Data <- get_data(data_in = Data1)

# Data$catch_obs_ysf[,2,1]
# alk <- reshape2::melt(Data$alk_ysal) %>%
#   rename(year = Var1, season = Var2, age = Var3, length = Var4)
# range(alk$value)
# ggplot(data = alk %>% filter(year %in% c(1, 46, 92), season == 1)) +
#   geom_line(aes(x = length, y = value, colour = factor(year))) +
#   facet_wrap(age ~ .)

# Create parameter list ----

# ADD A FUNCTION THAT DOES SOME PARAM CHECKS TOO (e.g., checks there are enough par sels, sets init vals)

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
               # par_sels_init_i = data_par1$par_sels_init_i[c(1:64, 70:83)], 
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

# n_sel_init <- sum(Data$sel_max_age_f - Data$sel_min_age_f + 1)
# n_sel_init
# length(Params$par_sels_init_i)
# 
# x <- Data$sel_change_sd_fy
# x[x > 0] <- 1
# n_sel_change <- sum(rowSums(x) * (Data$sel_max_age_f - Data$sel_min_age_f + 1))
# n_sel_change
# length(Params$par_sels_change_i)

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
# Map[["par_log_troll_tau"]] <- factor(NA) # est
Map[["par_log_hsp_q"]] <- factor(NA)
# Map[["par_logit_hstar_i"]] <- factor(rep(NA, 17)) # est
Map[["par_log_tag_H_factor"]] <- factor(NA)

# Specify the random effects
# Random <- c("par_rdev_y", "par_sels_change_i")
# Random <- c("par_rdev_y")
Random <- c()

# Create the AD object
obj <- MakeADFun(data = Data, parameters = Params, map = Map, random = Random, 
                 hessian = TRUE, inner.control = list(maxit = 1000), DLL = "sbt_v100")

obj$report()$lp_lf
obj$report()$lp_af
obj$report()$lp_cpue
sum(obj$report()$lp_aerial)
obj$report()$lp_troll
obj$report()$lp_tags #
sum(obj$report()$lp_pop)
sum(obj$report()$lp_hsp)
sum(obj$report()$lp_gt)

data_labrep1$lnlike
c(obj$report()$lp_lf[1:4], obj$report()$lp_af[5:6])
data_labrep1$lnlike[1:6]

# warning("ending now")
# stop()

df1 <- data_labrep1$tag.pred
df2 <- obj$report()$tag_pred %>%
  reshape2::melt() %>%
  pivot_wider(names_from = Var4) %>%
  arrange(Var1, Var2, Var3) %>%
  select(-Var1, -Var2, -Var3) %>%
  as.matrix()
head(df1)
head(df2)
tail(df1)
tail(df2)
df1 - df2

data_labrep1$tag.pred
data_labrep1$tag.obs

# List of parameters that are "on"
unique(names(obj$par))

# Set up estimation
newtonOption(obj = obj, smartsearch = TRUE)
obj$fn(obj$par)
obj$gr(obj$par)
obj$control <- list(trace = 100)
ConvergeTol <- 2 # 1:Normal; 2:Strong
# obj$env$inner.control$step.tol <- c(1e-12, 1e-15)[ConvergeTol] # Default : 1e-8 # Change in parameters limit inner optimization
# obj$env$inner.control$tol10 <- c(1e-8, 1e-12)[ConvergeTol]  # Default : 1e-3 # Change in pen.like limit inner optimization
# obj$env$inner.control$grad.tol <- c(1e-12, 1e-15)[ConvergeTol] # # Default : 1e-8 # Maximum gradient limit inner optimization
summary(obj)

bnd <- get_bounds(obj = obj)
check_bounds(opt = obj, lb = lb$lb, ub = lb$ub)

# Optimize
# Hess <- optimHess(par = opt$par, fn = obj$fn)
# opt <- nlminb(start = obj$par, objective = obj$fn, gr = obj$gr, upper = Upr, lower = Lwr, control = list(eval.max = 1e4, iter.max = 1e4, rel.tol = c(1e-10, 1e-8)[ConvergeTol], trace = 1))
opt <- nlminb(start = obj$par, objective = obj$fn, gr = obj$gr, upper = Upr, lower = Lwr)
              # control = list(eval.max = 2e4, iter.max = 1e4, rel.tol = 1e-7, trace = 1))

opt[["final_gradient"]] <- obj$gr(opt$par)
Diag <- obj$report()
Report <- sdreport(obj)
print(Report$pdHess) # Is the fit positive definite Hessian?
check_bounds(opt = opt, lb = Lwr, ub = Upr)
cbind(1:length(obj$par), Lwr, obj$par, Upr)
plot_selectivity(data = Data, object = obj, posterior = NULL, years = 1931:1958)
get_sel_list(data = Data, opt = opt)

# Set up a grid - now need to set up a way to run grids - pref in parallel, also this func
# needs to ensure that the Map below is fixed for grid elements
g_pars <- get_grid(par = Params)
library(foreach)
library(doParallel)
n_cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n_cores, 
  type = "PSOCK"
)
#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#check if it is registered (optional)
foreach::getDoParRegistered()
#how many workers are available? (optional)
foreach::getDoParWorkers()
# g_opt <- list()
# for (i in 1:nrow(g_pars)) {
# POOP
g_opt <- foreach(i = 1:5, .packages = "TMB") %dopar% {
  dyn.load(dynlib("src/sbt_v100"))
  g_obj <- MakeADFun(data = Data, parameters = g_pars[[i]], map = Map, random = Random, 
                   hessian = TRUE, inner.control = list(maxit = 1000), DLL = "sbt_v100")
  nlminb(start = g_obj$par, objective = g_obj$fn, gr = g_obj$gr, upper = Upr, lower = Lwr)
  g_obj
}
plot_recruitment(data = Data, object = g_opt[[5]], posterior = NULL)
# need to make these funcs work for grid as well as posterior
parallel::stopCluster(cl = my.cluster)


data <- Data
ages <- data$min_age:data$max_age
yrs <- data$first_yr:data$last_yr
fsh <- c("LL1", "LL2", "LL3", "LL4", "Indonesian", "Australian surface")
df_sel <- get_array(obj$report()$sel_fya) %>%
  reshape::melt() %>%
  mutate(fishery = fsh[X1], year = yrs[X2], age = ages[X3]) %>%
  filter(year >= data$first_yr_catch, year == 1952, value != 0, value < exp(-20))
df_sel

he <- obj$he()
he_inv <- solve(he)
he_ch <- chol(he)
ev <- eigen(he)
range(ev$values)

# Check generated quantities from fixed parameter run ----

plot_recruitment(data = Data, object = obj, posterior = NULL)
plot_biomass_spawning(data = Data, object = obj, posterior = NULL)
plot_hrate(data = Data, object = obj, posterior = NULL, years = 1990:2010)
plot_selectivity(data = Data, object = obj, posterior = NULL, years = 1990:2010, fisheries = c("LL1", "Indonesian"))
plot_natural_mortality(data = Data, object = obj, posterior = NULL)
plot_initial_numbers(data = Data, object = obj, posterior = NULL)

plot_compare <- function(x, ADMB, sbt, xlab, ylab) {
  df <- data.frame(x = x, ADMB = ADMB, sbt = sbt) %>%
    pivot_longer(cols = ADMB:sbt, names_to = "Model")
  
  ggplot(data = df, aes(x = x, y = value, colour = Model, linetype = Model)) +
    geom_line() + 
    labs(x = xlab, y = ylab) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
}

plot_compare(x = 1931:2023, 
             ADMB = data_labrep1$Sbio, sbt = obj$report()$spawning_biomass_y, 
             xlab = "Year", ylab = "Spawning biomass (tonnes)")

plot_compare(x = 1931:2023, 
             ADMB = data_labrep1$Recruitment, sbt = obj$report()$number_ysa[,1,1], 
             xlab = "Year", ylab = "Spawning biomass (tonnes)")

plot(obj$report()$number_ysa[1,1,])
lines(obj$report()$number_ysa[1,2,])
lines(obj$report()$number_ysa[2,1,], col = 2)
lines(obj$report()$number_ysa[2,2,], col = 2)

# Checks (to be moved to testthat directory eventually) ----

library(testthat)
expect_equal(obj$report()$par_B0, data_labrep1$B0, tolerance = 1e-06)
expect_equal(obj$report()$alpha, data_labrep1$alpha[1], tolerance = 1e-06)
expect_equal(obj$report()$lp_sel, as.numeric(data_labrep1$penal[c(2, 3, 13)]), tolerance = 1e-06)

# Save outputs
# capture.output(Report, file = paste(folder, "Report.txt", sep = ""))
# save(obj, file = paste(folder, "obj.RData", sep = ""))
# save(opt, file = paste(folder, "opt.RData", sep = ""))
# save(Report, file = paste(folder, "Report.RData", sep = ""))
# write.csv(data.frame(names(Report$value), Report$value), file = paste(folder, "Pars.csv", sep = ""), row.names = TRUE)

######################################################################################################
# Inspect results
######################################################################################################

Delta <- rep(0, length(opt$par))
Delta[2] <- 1e-5
(obj$fn(opt$par - Delta/2) - obj$fn(opt$par + Delta/2))  / abs(max(Delta))

# REs_zi <- Report$par.random[names(Report$par.random) %in% "z1"]
# plot(REs_zi)
# head(Report$value, 15)
# t(t(tapply(X = Report$value, INDEX = names(Report$value), FUN = length)))

# Likelihood profile ----

prof_sd_weight <- tmbprofile(obj = obj, name = "par_log_m4", ystep = 0.01, parm.range = c(-5, 0))
plot(x = prof_sd_weight)
# conf_sd_weight <- confint(object = prof_sd_weight)
# plot(x = conf_sd_weight)
# prof_sd_length <- tmbprofile(obj = obj, name = "log_sd_length", parm.range = c(-7, -1))
# conf_sd_length <- confint(object = prof_sd_length)
# plot(x = prof_sd_length)
# plot(x = conf_sd_length)

# Bayesian ----

library(tmbstan)

# Specify the number of cores for MCMC
options(mc.cores = parallel::detectCores())

mcmc1 <- tmbstan(obj = obj, lower = Lwr, upper = Upr, init = list(Params), chains = 2)
save(mcmc1, file = "mcmc1.rda")
load("mcmc1.rda")
# mcmc2 <- tmbstan(obj = obj, lower = Lwr, upper = Upr, init = list(Params), chains = 1, laplace = TRUE)
# get_stancode(mcmc1)

# traceplot(mcmc1, pars = c("log_alpha", "log_beta", "log_sd_length", "log_sd_weight", "lp__"), inc_warmup = FALSE)
pars <- c("lp__", "par_log_B0", "par_log_m4", "par_log_m30", "par_log_cpue_q", 
          "par_rdev_y[1]", "par_rdev_y[92]", 
          "par_sels_init_i[1]", "par_sels_init_i[78]", 
          "par_sels_change_i[1]", "par_sels_change_i[1132]")
traceplot(object = mcmc1, pars = pars, inc_warmup = FALSE)

plot_natural_mortality(data = Data, object = obj, posterior = mcmc1)

pp <- get_posterior(object = obj, posterior = mcmc1, pars = "par_B0")
hist(pp$value)

# Dynamically unload the model
dyn.unload(dynlib("sbt_model"))



# length_mean <- read_csv("data-raw/base22/csv/mean_length.csv")
# length_sd <- read_csv("data-raw/base22/csv/sd_length.csv")
# catch <- read_csv("data-raw/base22/csv/catch.csv")
# catch_UA <- read_csv("data-raw/base22/csv/catch_UA.csv")
# scenarios_surface <- read_csv("data-raw/base22/csv/scenarios_surface.csv")
# scenarios_LL1 <- read_csv("data-raw/base22/csv/scenarios_LL1.csv")
# sel_change_sd <- read_csv("data-raw/base22/csv/sd_sel_change.csv")
# POPs <- read_csv("data-raw/base22/csv/POPs.csv")
# HSPs <- read_csv("data-raw/base22/csv/HSPs.csv")
# GTs <- read_csv("data-raw/base22/csv/GTs.csv")
# aerial_surv <- read_csv("data-raw/base22/csv/aerial_survey.csv")
# aerial_cov <- read_csv("data-raw/base22/csv/aerial_cov.csv", col_names = FALSE) %>% as.matrix()
# troll_ind <- read_csv("data-raw/base22/csv/trolling_index.csv")
# af <- read_csv("data-raw/base22/csv/age_freq.csv")
# lf <- read_csv("data-raw/base22/csv/length_freq.csv")
