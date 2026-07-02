# Initial set up ----

rm(list = ls())

# remotes::install_github("janoleko/RTMBdist")
# remotes::install_github("andrjohns/StanEstimators")
# remotes::install_github("noaa-afsc/SparseNUTS")
# remotes::install_github(repo = "quantifish/sbt")

if (basename(getwd()) != "OMMP16") setwd(file.path(getwd(), "OMMP16"))

library(tidyverse)
library(sbt)
library(SparseNUTS)
library(bayesplot)

theme_set(theme_bw())

# run_suffix <- "troll_switch0_legacy9_free_m0_m10_psi1.25_h0.63_sparsenuts"
# pdf(paste0("Rplots_", run_suffix, ".pdf"))
# on.exit(dev.off(), add = TRUE)

# Read in the data ----

data_loc <- "csv_2026"

length_mean <- read_csv(file.path(data_loc, "mean_length.csv"))
length_sd <- read_csv(file.path(data_loc, "sd_length.csv"))
catch <- read_csv(file.path(data_loc, "catch.csv"))
catch_UA <- read_csv(file.path(data_loc, "catch_UA.csv"))
scenarios_surface <- read_csv(file.path(data_loc, "scenarios_surface.csv"))
scenarios_LL1 <- read_csv(file.path(data_loc, "scenarios_LL1.csv"))
POPs <- read_csv(file.path(data_loc, "POPs.csv"))
HSPs <- read_csv(file.path(data_loc, "HSPs.csv"))
GTs <- read_csv(file.path(data_loc, "GTs.csv"))
troll <- read_csv(file.path(data_loc, "trolling_index.csv"))
cpue <- read_csv(file.path(data_loc, "cpue.csv"))
age_freq <- read_csv(file.path(data_loc, "age_freq.csv"))
length_freq <- read_csv(file.path(data_loc, "lf_assessment.csv")) # this is the new version from Jim with the CPUE LFs

data <- list(
  last_yr = 2025, age_increase_M = 25,
  length_m50 = 150, length_m95 = 180, 
  catch_UR_on = 0, catch_surf_case = 1, catch_LL1_case = 1, 
  scenarios_surf = scenarios_surface, scenarios_LL1 = scenarios_LL1,
  removal_switch_f = c(0, 0, 0, 1, 0, 0), # 0=harvest rate, 1=direct removals
  sel_min_age_f = c(2, 2, 2, 8, 6, 0, 4),
  # sel_max_age_f = c(17, 9, 17, 22, 25, 7, 17),
  sel_max_age_f = c(17, 9, 17, 21, 25, 7, 17),
  sel_end_f = c(1, 0, 1, 1, 1, 0, 1), # 0=zero, 1=constant
  sel_LL1_yrs = c(1952, 1957, 1961, 1965, 1969, 1973, 1977, 1981, 1985, 1989, 1993, 1997, 2001, 2006, 2007, 2008, 2011, 2014, 2017, 2020, 2023),
  sel_LL2_yrs = c(1969, 2001, 2005, 2008, 2011, 2014, 2017, 2020, 2023),
  sel_LL3_yrs = c(1954, 1961, 1965, 1969, 1970, 1971, 2005, 2006, 2007),
  sel_LL4_yrs = c(1953),
  sel_Ind_yrs = c(1976, 1997, 1999, 2002, 2004, 2006, 2008, 2010, 2012:2021),
  sel_Aus_yrs = c(1952, 1969, 1973, 1977, 1981, 1985, 1989, 1993, 1997:2025),
  sel_CPUE_yrs = c(1969, 1973, 1977, 1981, 1985, 1989, 1993, 1997, 2001, 2006, 2007, 2008, 2011, 2014, 2017, 2020),
  af_switch = 1, # 1=multinomial, 2=Dirichlet, 3=Dirichlet-multinomial
  lf_switch = 1, lf_minbin = c(1, 1, 1, 11, 6), # seq(87.5, 184, 4)
  cpue_switch = 1, cpue_a1 = 5, cpue_a2 = 17,
  aerial_switch = 4, aerial_tau = 0.3, 
  troll_switch = 0, 
  pop_switch = 1, 
  hsp_switch = 1, 
  hsp_false_negative = 0.6840729,
  gt_switch = 1,
  tag_switch = 1, tag_var_factor = 1.82
)

data_in <- data
data <- get_data(data_in = data_in)

plot(seq(87.5, 184, 4), data$cpue_lfs[57,], type = "b")
lines(as.numeric(names(length_freq)[4:113]), length_freq[270, 4:113], col = 2, type = "b")
abline(h = 0)

# Parameters ----

parameters <- get_parameters(data = data)

exp(parameters$par_log_h)
parameters$par_log_h <- log(0.72)
parameters$par_log_psi <- log(1.75)
# parameters$par_log_psi <- log(1.25)

parameters$par_sel_rho_a[2] <- sel_rho_to_par(0.5)
parameters$par_sel_rho_y[2] <- sel_rho_to_par(0.7)
parameters$par_log_sel_sigma[2] <- log(0.5)

parameters$par_sel_rho_a[5] <- sel_rho_to_par(0.95)
parameters$par_sel_rho_y[5] <- sel_rho_to_par(0.6)
parameters$par_log_sel_sigma[5] <- log(0.2)

parameters$par_sel_rho_a[7] <- sel_rho_to_par(0.9)
parameters$par_sel_rho_y[7] <- sel_rho_to_par(0.9)
parameters$par_log_sel_sigma[7] <- log(0.1)

tibble(
  fishery = c("LL1", "LL2", "LL3", "LL4", "Indonesia", "Australia", "CPUE"),
  rho_a = sel_rho_from_par(parameters$par_sel_rho_a),
  rho_y = sel_rho_from_par(parameters$par_sel_rho_y),
  sigma = exp(parameters$par_log_sel_sigma)
)

map <- get_map(parameters = parameters)
map$par_log_m0 <- NULL
map$par_log_m10 <- NULL

data$priors <- get_priors(parameters = parameters)
evaluate_priors(parameters = parameters, priors = data$priors)

obj <- MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)

bounds <- get_bounds(obj = obj, parameters = parameters)

unique(names(obj$par))
obj$fn()

# Optimise a single grid cell ----

control <- list(eval.max = 10000, iter.max = 10000)

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
              lower = bounds$lower, upper = bounds$upper, control = control)
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
              lower = bounds$lower, upper = bounds$upper, control = control)

check_estimability(obj = obj)
TMBhelper::check_estimability(obj = obj)
# ce <- check_estimability(obj = obj)
# ce[[4]] %>% filter(Param_check != "OK"
obj$env$last.par.best[1:13]
exp(obj$env$last.par.best[1:13])
plot_selectivity(data = data, object = obj, years = 1:2021, fisheries = "Indonesian")
p1 <- plot_selectivity(data = data, object = obj, years = 1:3000, fisheries = "LL1")
p2 <- plot_selectivity(data = data, object = obj, years = 1:3000, fisheries = "CPUE")
p1 + p2

plot_cpue_residuals(data, obj)
plot_gt_residuals(data, obj)
plot_hsp_residuals(data, obj)
plot_pop_residuals(data, obj)
plot_tag_residuals(data, obj)
plot_troll_residuals(data, obj)
plot_aerial_residuals(data, obj)

plot_af_residuals(data, obj, fishery = "Indonesian")
plot_af_residuals(data, obj, fishery = "Australian")
plot_lf_residuals(data, obj, fishery = "LL1")
plot_lf_residuals(data, obj, fishery = "LL2")
plot_lf_residuals(data, obj, fishery = "LL3")
plot_lf_residuals(data, obj, fishery = "LL4")
plot_lf_residuals(data, obj, fishery = "CPUE")

plot_lf(data, obj, fishery = "LL3")
plot_lf(data, obj, fishery = "CPUE")
plot_lf(data = data, object = obj, fishery = "CPUE")

estimated_parameters <- make_parameter_table(obj = obj, data = data)
print(estimated_parameters, n = Inf)

# Sensitivity 1: change q from 2008 ----

seed_sensitivity_parameters <- function(parameters_sens, reference_obj) {
  reference_parameters <- reference_obj$env$parList(reference_obj$env$last.par.best)
  for (name in intersect(names(parameters_sens), names(reference_parameters))) {
    fitted <- as.numeric(reference_parameters[[name]])
    if (length(fitted) == 1L && length(parameters_sens[[name]]) > 1L) {
      parameters_sens[[name]][] <- fitted
    } else {
      n <- min(length(parameters_sens[[name]]), length(fitted))
      parameters_sens[[name]][seq_len(n)] <- fitted[seq_len(n)]
    }
  }
  parameters_sens
}

fit_sensitivity <- function(data_sens, label, reference_obj = obj) {
  parameters_sens <- get_parameters(data = data_sens)
  parameters_sens <- seed_sensitivity_parameters(parameters_sens, reference_obj)
  map_sens <- get_map(parameters = parameters_sens)
  map_sens$par_log_m0 <- NULL
  map_sens$par_log_m10 <- NULL
  data_sens$priors <- get_priors(parameters = parameters_sens)
  obj_sens <- MakeADFun(func = cmb(sbt_model, data_sens), parameters = parameters_sens, map = map_sens)
  bounds_sens <- get_bounds(obj = obj_sens, parameters = parameters_sens)
  opt_sens <- nlminb(
    start = obj_sens$par, objective = obj_sens$fn, gradient = obj_sens$gr, hessian = obj_sens$he,
    lower = bounds_sens$lower, upper = bounds_sens$upper, control = control
  )
  opt_sens <- nlminb(
    start = opt_sens$par, objective = obj_sens$fn, gradient = obj_sens$gr, hessian = obj_sens$he,
    lower = bounds_sens$lower, upper = bounds_sens$upper, control = control
  )
  obj_sens$fn(opt_sens$par)
  gr_sens <- obj_sens$gr(opt_sens$par)
  message(
    label, ": convergence = ", opt_sens$convergence,
    " (", opt_sens$message, ")",
    ", objective = ", signif(opt_sens$objective, 10),
    ", max |gradient| = ", signif(max(abs(gr_sens)), 5),
    ", max positive gradient = ", signif(max(gr_sens), 5)
  )
  list(data = data_sens, parameters = parameters_sens, map = map_sens,
       obj = obj_sens, bounds = bounds_sens, opt = opt_sens)
}

data_q2008_in <- data_in
data_q2008_in$cpue_q_yrs <- c(2008)
data_q2008 <- get_data(data_in = data_q2008_in)
sens_q2008 <- fit_sensitivity(data_sens = data_q2008, label = "Sensitivity 1: CPUE q split from 2008")
data_q2008 <- sens_q2008$data
parameters_q2008 <- sens_q2008$parameters
map_q2008 <- sens_q2008$map
obj_q2008 <- sens_q2008$obj
bounds_q2008 <- sens_q2008$bounds
opt_q2008 <- sens_q2008$opt

# Sensitivity 2: time varying CPUE CV ----

get_data_with_cpue <- function(data_in, cpue_sens) {
  cpue_base <- cpue
  on.exit(cpue <<- cpue_base, add = TRUE)
  cpue <<- cpue_sens
  get_data(data_in = data_in)
}

cpue_tv_cv <- cpue %>% mutate(CV = seq(0.30, 0.18, length.out = n()))
data_cpue_tv_cv <- get_data_with_cpue(data_in = data_in, cpue_sens = cpue_tv_cv)
sens_cpue_tv_cv <- fit_sensitivity(data_sens = data_cpue_tv_cv, label = "Sensitivity 2: time-varying CPUE CV")
data_cpue_tv_cv <- sens_cpue_tv_cv$data
parameters_cpue_tv_cv <- sens_cpue_tv_cv$parameters
map_cpue_tv_cv <- sens_cpue_tv_cv$map
obj_cpue_tv_cv <- sens_cpue_tv_cv$obj
bounds_cpue_tv_cv <- sens_cpue_tv_cv$bounds
opt_cpue_tv_cv <- sens_cpue_tv_cv$opt

plot_biomass_spawning(data_list = list(data, data_q2008, data_cpue_tv_cv), object_list = list(obj, obj_q2008, obj_cpue_tv_cv))

# Profiling psi ----

# map_psi <- map
# map_psi$par_log_psi <- NULL
# par_psi <- obj$env$parList(obj$env$last.par.best)
# obj_psi <- MakeADFun(func = cmb(sbt_model, data), parameters = par_psi, map = map_psi)
# prof_psi1 <- sbtprofile(obj = obj_psi, name = "par_log_psi")
# plot_profile(obj = obj_psi, x = prof_psi1, xlab = "Psi")
# prof_psi2 <- TMB::tmbprofile(obj_psi, name = "par_log_psi", parm.range = c(0, 0.7))
# prof_psi2$par_log_psi <- exp(prof_psi2$par_log_psi)
# plot(prof_psi2)

# Try turning on random effects ----

map2 <- map
map2$par_sel_rho_y <- factor(c(NA, NA, NA, NA, 1, NA, NA))
map2$par_sel_rho_a <- factor(c(NA, NA, NA, NA, 1, NA, NA))
map2$par_log_sel_sigma <- factor(c(NA, NA, NA, NA, 1, NA, NA))
# map2$par_sel_rho_y <- factor(c(1, NA, NA, NA, NA, NA, NA))
# map2$par_sel_rho_a <- factor(c(1, NA, NA, NA, NA, NA, NA))
# map2$par_log_sel_sigma <- factor(c(1, NA, NA, NA, NA, NA, NA))
parameters2 <- obj$env$parList(obj$env$last.par.best)
# parameters2$par_sel_rho_y[5] <- sel_rho_to_par(0.95)
# parameters2$par_sel_rho_a[5] <- sel_rho_to_par(0.95)
# parameters2$par_log_sel_sigma[5] <- log(0.15)

tibble(
  fishery = c("LL1", "LL2", "LL3", "LL4", "Indonesia", "Australia", "CPUE"),
  rho_a = sel_rho_from_par(parameters2$par_sel_rho_a),
  rho_y = sel_rho_from_par(parameters2$par_sel_rho_y),
  sigma = exp(parameters2$par_log_sel_sigma)
)

obj2 <- MakeADFun(func = cmb(sbt_model, data), parameters = parameters2, map = map2, random = c("par_log_sel_5"))
bounds2 <- get_bounds(obj = obj2, parameters = parameters2)
opt2 <- nlminb(start = obj2$par, objective = obj2$fn, gradient = obj2$gr,
               lower = bounds2$lower, upper = bounds2$upper, control = control)
opt2 <- nlminb(start = opt2$par, objective = obj2$fn, gradient = obj2$gr,
               lower = bounds2$lower, upper = bounds2$upper, control = control)
opt2$par[7:9]

plot_lf(data = data, object = obj2, fishery = "CPUE")
plot_selectivity(data = data, object = obj2, years = 1:3000, fisheries = "CPUE")
plot_selectivity(data = data, object = obj2, years = 1:3000, fisheries = "Indonesian")
plot_selectivity(data = data, object = obj2, years = 1:3000, fisheries = "LL1")
plot_af(data = data, object = obj2, fishery = "Indonesian")
plot_hsps(data, obj2)
obj2$env$last.par.best[1:10]

# Inspect model outputs ----

plot_selectivity(data = data, object = obj, years = 1:3000, fisheries = "LL1")
plot_selectivity(data = data, object = obj, years = 1:3000, fisheries = "LL2")
plot_selectivity(data = data, object = obj, years = 1:3000, fisheries = "LL3")
# plot_selectivity(data = data, object = obj, years = 1:3000, fisheries = "LL4")
plot_selectivity(data = data, object = obj, years = 1:3000, fisheries = "Indonesian")
plot_selectivity(data = data, object = obj, years = 1960:3000, fisheries = "Australian")

plot_cpue_residuals(data, obj)
plot_gt_residuals(data, obj)
plot_hsp_residuals(data, obj)
plot_pop_residuals(data, obj)
plot_tag_residuals(data, obj)

plot_hsps(data, obj)

p1 <- plot_hsps(data, obj1)
p2 <- plot_hsps(data, obj2)
p1 / p2

p1 <- plot_hsp_residuals(data, obj1)
p2 <- plot_hsp_residuals(data, obj2)
p1 + p2

exp(obj1$env$last.par.best[4])
exp(obj2$env$last.par.best[4])

p1 <- plot_selectivity(data = data, object = obj1, years = 1900:2026, fisheries = "Indonesian")
p2 <- plot_selectivity(data = data, object = obj2, years = 1900:2026, fisheries = "Indonesian")
p1 + p2

plot_af_residuals(data = data, obj = obj, fishery = 5)
plot_af_residuals(data = data, obj = obj, fishery = 6)

plot_cpue(data = data, object = obj, nsim = 1)
plot_cpue_residuals(data = data, obj = obj, type = "OSA")

plot_aerial_survey(data = data, object = obj, nsim = 1)

plot_lf(data = data, object = obj, fishery = "CPUE")

plot_af(data = data, object = obj2, fishery = "Indonesian")
plot_selectivity(data = data, object = obj, years = 2013:2026, fisheries = "Indonesian")
plot_selectivity(data = data, object = obj, years = 1900:2026, fisheries = "Indonesian")

plot_af(data = data, object = obj2, fishery = "Australian")
plot_selectivity(data = data, object = obj, years = 1900:2026, fisheries = "Australian")

plot_lf(data = data, object = obj, fishery = "LL1")
plot_lf(data = data, object = obj, fishery = "LL2")
plot_lf(data = data, object = obj, fishery = "LL3")
plot_lf(data = data, object = obj, fishery = "LL4")

p1 <- plot_selectivity(data = data, object = obj, years = 1969:2025, fisheries = "LL1")
p2 <- plot_selectivity(data = data, object = obj, years = 1969:2025, fisheries = "CPUE")
p1 + p2
plot_selectivity(data = data, object = obj, years = 2013:2026)
# plot_selectivity(data = data, object = obj, years = 1969:2022, fisheries = "LL3")
# plot_selectivity(data = data, object = obj, years = 1969:2022, fisheries = "LL4")

plot_biomass_spawning(data_list = list(data, data), object_list = list(obj1, obj2))

# Run MLE grid ----

do_run <- FALSE

if (do_run) {
  # grid_pars <- get_grid(parameters = parameters, m0 = c(0.4), m10 = c(0.065), h = c(0.8), psi = c(1.5, 1.75, 2))
  # grid_pars <- get_grid(parameters = parameters, h = c(0.63, 0.72), psi = c(1.75))
  grid_pars <- get_grid(parameters = parameters)
  grid_list <- run_grid(data = data, grid_parameters = grid_pars, bounds = bounds, map = map, control = control, parallel = FALSE)
  grid_check <- check_grid(grid = grid_list)
  # if some grid cells have not converged then you can run them again using rerun_grid
  # idc <- which(grid_check$grid_summary$Check != "All parameters are estimable")
  # grid_list <- rerun_grid(grid = grid_list, bounds = bounds, cells = idc, control = control) # this crashes my machine
  grid_cells <- sample_grid(grid = grid_list, seed = 42, prior_psi = c(1))
  grid_tmbfit <- grid_to_tmbfit(data = data, parameters = parameters, grid = grid_list, grid_parameters = grid_pars, grid_cells = grid_cells)
  save(grid_pars, grid_list, grid_check, grid_cells, grid_tmbfit, file = paste0("grid_run1_", run_suffix, ".rda"))
  # save_grid(grid = grid_list, dir = "inst/extdata/grid_list", overwrite = TRUE, compress = "gzip")
  # save(grid_pars, file = "inst/extdata/grid_pars.rda")
  # save(grid_list, file = "inst/extdata/grid_list.rda")
  # save(grid_check, file = "inst/extdata/grid_check.rda")
  # save(grid_cells, file = "inst/extdata/grid_cells.rda")
  # save(grid_tmbfit, file = "inst/extdata/grid_tmbfit.rda")
} else {
  # load(system.file("extdata", "grid_pars.rda", package = "sbt"))
  # # grid_list <- load_grid(dir = "inst/extdata/grid_list")
  # load(system.file("extdata", "grid_list.rda", package = "sbt"))
  # load(system.file("extdata", "grid_check.rda", package = "sbt"))
  # load(system.file("extdata", "grid_cells.rda", package = "sbt"))
  # load(system.file("extdata", "grid_tmbfit.rda", package = "sbt"))
}

if (do_run) {
  length(grid_pars)
  sapply(grid_list, function(x) !is.null(x$opt) && x$opt$convergence != 0)
  data.frame(grid_check$grid_summary)
  check_estimability(grid_list[[6]])
  left_join(grid_cells$grid_freq, grid_check$grid_summary, by = join_by(m0, m10, h, psi)) %>% 
    filter(Freq > 0) %>%
    data.frame()
  # psi_values <- sapply(grid_list, function(x) exp(x$report()$par_log_psi))
  # unique(psi_values)
  table(exp(grid_tmbfit$samples[,,'par_log_psi']))
  table(exp(grid_tmbfit$samples[,,'par_log_m0']))
  
  table(grid_check$grid_summary$m10[grid_cells$grid_cells])
  table(exp(grid_tmbfit$samples[,,'par_log_m10']))
  data$priors$par_log_m10 # the m10 prior could be doing this?
  
  plot_biomass_spawning(data_list = rep(list(data), times = length(grid_list)), object_list = grid_list)
  
  # extract_samples(grid_tmbfit)
  post <- as.data.frame(grid_tmbfit)
  pars <- grid_tmbfit$par_names[1:8]
  mcmc_trace(x = post, pars = pars)
}

run_mcmc <- TRUE
if (!run_mcmc) quit(save = "no")

# MCMC for single grid cell ----

mcmc <- sample_snuts(
  obj = obj, metric = "auto", num_samples = 3, num_warmup = 3,
  # iter = 1000, warmup = 750, chains = 4, cores = 4,
  # obj = obj, metric = "auto", iter = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.99), init = "last.par.best",
  # lower = bounds$lower, upper = bounds$upper, # these bounds dont seem to work
  globals = sbt_globals())

save(data, parameters, obj, opt, mcmc, file = paste0("mcmc_", run_suffix, ".rda"))

quit(save = "no")

# save(data, parameters, obj, opt, mcmc, file = "mcmc_0divergences.rda")
# save(data, parameters, obj, opt, mcmc, file = "mcmc_new_rec.rda")
# save(data, parameters, obj, opt, mcmc, file = "mcmc_old_rec.rda")
# save(data, parameters, obj, opt, mcmc, file = "mcmc_new2_rec.rda")
# 
# load("mcmc_old_rec.rda")
# pold <- plot_recruitment(data = data, object = obj, posterior = mcmc)
# plot_rec_devs(data = data, object = obj, posterior = mcmc)
# load("mcmc_new_rec.rda")
# pnew <- plot_recruitment(data = data, object = obj, posterior = mcmc)
# plot_rec_devs(data = data, object = obj, posterior = mcmc)
# load("mcmc_new2_rec.rda")
# pnew2 <- plot_recruitment(data = data, object = obj, posterior = mcmc)
# plot_rec_devs(data = data, object = obj, posterior = mcmc)
# pold + pnew2

plot_sampler_params(fit = mcmc, plot = TRUE)
plot_uncertainties(fit = mcmc, log = TRUE, plot = TRUE)
decamod::pairs_rtmb(fit = mcmc, order = "slow", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "mismatch", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "fast", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "divergent", pars = 1:5)

# Do grid of MCMCs ----

map$par_log_h <- NULL
obj <- MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
lincomb <- numeric(length(obj$par))
lincomb[1] <- 1
lincomb[1] <- 6
# x2 <- TMB::tmbprofile(obj, lincomb = lincomb)
x2 <- TMB::tmbprofile(obj, name = "par_log_h", parm.range = c(-0.15, -0.1), ystep = 0.001)
plot(x2)

hh <- seq(0.8, 0.55, length.out = 8)
nll <- numeric(length(hh))
for (i in 1:length(hh)) {
  parameters$par_log_h <- log(hh[i])
  obj <- MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he, 
                lower = bounds$lower, upper = bounds$upper, control = control)
  opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he, 
                lower = bounds$lower, upper = bounds$upper, control = control)
  nll[i] <- obj$fn(opt$par)
}
plot(hh, nll, type = "l")

# mcmc_list <- list()
# obj_list <- list()
run_suffix <- "test_grid"
grd <- expand.grid(h = c(0.55, 0.63, 0.72, 0.8), psi = c(1.5, 1.75, 2))
grid_dir <- paste0("grid1_", run_suffix)
dir.create(grid_dir)
for (i in 1:nrow(grd)) {
  parameters$par_log_h <- log(grd$h[i])
  parameters$par_log_psi <- log(grd$psi[i])
  obj <- MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he, 
                lower = bounds$lower, upper = bounds$upper, control = control)
  opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he, 
                lower = bounds$lower, upper = bounds$upper, control = control)
  
  # Q <- obj$he(opt$par) # Qinv and use metric dense
  # sdreport(obj, hessian.fixed = Q)
  # x <- TMB::tmbprofile(obj, "par_log_B0")
  # x2 <- sbtprofile(obj, name = "par_log_B0")
  # plot_profile(obj = obj, x = x2, xlab = "B0")
  
  # he <- obj$he()
  # he_inv <- solve(he)
  # he_ch <- chol(he)
  # ev <- eigen(he)
  # range(ev$values)

  mcmc <- sample_snuts(
    obj = obj, 
    # metric = "auto", # set this to dense?
    metric = "dense", # set this to dense?
    num_samples = 500, num_warmup = 150,
    # iter = 1000, warmup = 750, chains = 4, cores = 4,
    # obj = obj, metric = "auto", iter = 2000, chains = 4, cores = 4,
    control = list(adapt_delta = 0.99), init = "last.par.best",
    # lower = bnd$lower, upper = bnd$upper, # these bounds dont seem to work
    skip_optimization = TRUE,
    Qinv = obj$he(obj$env$last.par.best),
    globals = sbt_globals())
  
  save(opt, mcmc, file = paste0(grid_dir, "/grid", i, ".rda"))
  # obj_list[[i]] <- obj
  # mcmc_list[[i]] <- mcmc
}

load("grid1_test_grid/grid5.rda")
plot_marginals(fit = mcmc, pars = 1:6)
pairs(mcmc, pars = 1:8)
decamod::pairs_rtmb(fit = mcmc, order = "slow", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "mismatch", pars = 1:5)
decamod::pairs_rtmb(fit = mcmc, order = "fast", pars = 1:5)
plot_sampler_params(fit = mcmc, plot = TRUE)

load("grid1_test_grid/grid2.rda")
mcmc1 <- mcmc
load("grid1_test_grid/grid3.rda")
mcmc2 <- mcmc

grid_mcmc_to_tmbfit <- function(data, parameters, grid, grid_parameters, grid_cells) {
  map <- get_map(parameters = parameters)
  map$par_log_psi <- NULL
  map$par_log_m0 <- NULL
  map$par_log_m10 <- NULL
  map$par_log_h <- NULL
  obj <- MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
  parnames <- names(obj$par)
  n_samples <- length(grid_cells$grid_cells)
  samples <- array(dim = c(n_samples, 1, length(obj$par)), dimnames = list(NULL, NULL, parnames))
  for (i in 1:n_samples) {
    j <- grid_cells$grid_cells[i]
    # obj$env$parList() alternatively could use this below?
    samples[i, 1,] <- c(grid[[j]]$par[1], # par_log_B0
                        grid_parameters[[j]]$par_log_psi, 
                        grid_parameters[[j]]$par_log_m0, 
                        grid[[j]]$par[2], # par_log_m4
                        grid_parameters[[j]]$par_log_m10,
                        grid[[j]]$par[3], # par_log_m30
                        grid_parameters[[j]]$par_log_h,
                        grid[[j]]$par[-c(1:3)])
  }
  # getS3method("print", "tmbfit")
  # Create tmbfit object
  x <- list(samples = samples, sampler_params = list(), mle = NULL,
            monitor = NULL, model = "RTMB", metric = "", par_names = parnames, 
            max_treedepth = NA, warmup = 0, iter = n_samples, thin = 1, 
            time.warmup = 0, time.sampling = NA, time.total = NA, algorithm = "grid")
  class(x) <- c("tmbfit", "list")
  return(x)
}

# save(mcmc_list, obj_list, file = "mcmc_grid46.rda")

# Sensitivities ----

# To do at OMMP16

# Projections ----

# To do at OMMP16

# proj_fix <- repair_af_slices_weighted(
#   catch_ysf = proj_catch_ysf,
#   af_ysfa   = proj_af_sliced_ysfa,
#   N = 3, p = 1
# )
# proj_af_sliced_ysfa <- proj_fix$af_ysfa
