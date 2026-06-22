# Initial set up ----

rm(list = ls())

# remotes::install_github("janoleko/RTMBdist")
# remotes::install_github("andrjohns/StanEstimators")
# remotes::install_github("noaa-afsc/SparseNUTS")
# remotes::install_github(repo = "quantifish/sbt")

# setwd(file.path(getwd(), "OMMP16"))

library(tidyverse)
library(sbt)
library(SparseNUTS)
library(bayesplot)

theme_set(theme_bw())

# Read in the data ----

data_loc <- "csv_2026"

length_mean <- read_csv(file.path(data_loc, "mean_length.csv"))
length_sd <- read_csv(file.path(data_loc, "sd_length.csv"))
catch <- read_csv(file.path(data_loc, "catch.csv"))
catch_UA <- read_csv(file.path(data_loc, "catch_UA.csv"))
scenarios_surface <- read_csv(file.path(data_loc, "scenarios_surface.csv"))
scenarios_LL1 <- read_csv(file.path(data_loc, "scenarios_LL1.csv"))
# POPs <- read_csv(file.path(data_loc, "POPs.csv"))
POPs <- read_csv("POPs_test.csv")
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
  # removal_switch_f = c(0, 0, 0, 0, 0, 0), # 0=harvest rate, 1=direct removals
  sel_min_age_f = c(2, 2, 2, 8, 6, 0, 2),
  # sel_max_age_f = c(17, 9, 17, 22, 25, 7, 17),
  sel_max_age_f = c(17, 9, 17, 21, 25, 7, 17),
  sel_end_f = c(1, 0, 1, 1, 1, 0, 1), # 0=zero, 1=constant
  sel_LL1_yrs = c(1952, 1957, 1961, 1965, 1969, 1973, 1977, 1981, 1985, 1989, 1993, 1997, 2001, 2006, 2007, 2008, 2011, 2014, 2017, 2020, 2023),
  sel_LL2_yrs = c(1969, 2001, 2005, 2008, 2011, 2014, 2017, 2020),
  sel_LL3_yrs = c(1954, 1961, 1965, 1969, 1970, 1971, 2005, 2006, 2007),
  sel_LL4_yrs = c(1953),
  sel_Ind_yrs = c(1976, 1995, 1997, 1999, 2002, 2004, 2006, 2008, 2010, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023),
  sel_Aus_yrs = c(1952, 1969, 1973, 1977, 1981, 1985, 1989, 1993, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023),
  sel_CPUE_yrs = c(1969, 1973, 1977, 1981, 1985, 1989, 1993, 1997, 2001, 2006, 2007, 2008, 2011, 2014, 2017, 2020),
  # af_switch = 9,
  # af_switch = 1, # CAUSES ISSUES WHY?
  af_switch = 9, # 1=multinomial, 2=Dirichlet, 3=Dirichlet-multinomial, 9=old
  lf_switch = 9, lf_minbin = c(1, 1, 1, 11),
  cpue_switch = 1, cpue_a1 = 5, cpue_a2 = 17,
  aerial_switch = 4, aerial_tau = 0.3, 
  troll_switch = 0, 
  pop_switch = 1, 
  hsp_switch = 1, 
  hsp_false_negative = 0.6840729, # If I set the hsp false negative to 1 and af_switch = 9 then model OK
  gt_switch = 1,
  tag_switch = 1, tag_var_factor = 1.82
)

data <- get_data(data_in = data)

# Need to fix LFs ----

repair_af_slices_weighted <- function(catch_ysf,
                                      af_ysfa,
                                      N = 3L,
                                      p = 1,
                                      eps = 1e-12,
                                      expand_if_empty = TRUE,
                                      verbose = TRUE) {
  # af_ysfa is assumed [year, season, fishery, age]
  d_af <- dim(af_ysfa)
  stopifnot(length(d_af) == 4)
  n_y <- d_af[1]; n_s <- d_af[2]; n_f <- d_af[3]; n_a <- d_af[4]
  
  d_c <- dim(catch_ysf)
  stopifnot(length(d_c) == 3)
  
  # map catch dims to y/s/f by matching sizes in af
  map_dim <- function(target, used) {
    idx <- which(d_c == target & !(seq_along(d_c) %in% used))
    if (length(idx) != 1) return(NA_integer_)
    idx
  }
  used <- integer(0)
  iy <- map_dim(n_y, used); used <- c(used, iy)
  is <- map_dim(n_s, used); used <- c(used, is)
  iff <- map_dim(n_f, used)
  
  if (anyNA(c(iy, is, iff))) {
    stop(sprintf(
      "Could not align catch_ysf dims %s to af_ysfa y/s/f dims %s. ",
      paste(d_c, collapse = "x"),
      paste(c(n_y, n_s, n_f), collapse = "x")
    ))
  }
  
  # helper to read catch regardless of dim order
  get_catch <- function(y, s, f) {
    idx <- vector("list", 3)
    idx[[iy]] <- y; idx[[is]] <- s; idx[[iff]] <- f
    do.call(`[`, c(list(catch_ysf), idx, list(drop = TRUE)))
  }
  
  repaired <- vector("list", 0L)
  
  for (f in seq_len(n_f)) {
    for (s in seq_len(n_s)) {
      af_sum_y <- vapply(seq_len(n_y), function(yi) sum(af_ysfa[yi, s, f, ], na.rm = TRUE), numeric(1))
      usable_y <- which(af_sum_y > eps)
      
      for (y in seq_len(n_y)) {
        cval <- get_catch(y, s, f)
        asum <- af_sum_y[y]
        
        if (!is.na(cval) && cval > 0 && (is.na(asum) || asum <= eps)) {
          lo <- max(1L, y - N); hi <- min(n_y, y + N)
          cand <- setdiff(seq.int(lo, hi), y)
          donors <- intersect(cand, usable_y)
          if (length(donors) == 0L && expand_if_empty) donors <- setdiff(usable_y, y)
          if (length(donors) == 0L) {
            stop(sprintf("Cannot repair AF at y=%d s=%d f=%d: no donor years with nonzero AF.", y, s, f))
          }
          
          dist <- abs(donors - y); dist[dist == 0] <- 1e-9
          w <- (1 / (dist^p)); w <- w / sum(w)
          
          donor_mat <- vapply(donors, function(yy) af_ysfa[yy, s, f, ], numeric(n_a))
          af_new <- as.numeric(donor_mat %*% w)
          af_new[is.na(af_new)] <- 0
          z <- sum(af_new)
          if (z <= eps) stop(sprintf("Repair produced zero AF at y=%d s=%d f=%d.", y, s, f))
          af_ysfa[y, s, f, ] <- af_new / z
          af_sum_y[y] <- 1.0
          
          repaired[[length(repaired) + 1L]] <- data.frame(
            year = y, season = s, fishery = f, n_donors = length(donors)
          )
        }
      }
    }
  }
  
  repaired_df <- if (length(repaired)) do.call(rbind, repaired) else data.frame()
  if (verbose && nrow(repaired_df) > 0) {
    message(sprintf("Repaired %d year-season-fishery AF slice(s).", nrow(repaired_df)))
  }
  list(af_ysfa = af_ysfa, repaired = repaired_df, catch_dim_map = c(year = iy, season = is, fishery = iff))
}

for (f in seq_len(data$n_fishery)) {
  if (data$removal_switch_f[f] > 0) {
    x1 <- data$catch_obs_ysf[,,f]
    x2 <- apply(data$af_sliced_ysfa[,,f,], 1:2, sum)
    # x1 <- proj_catch_ysf[,,f]
    # x2 <- apply(proj_af_sliced_ysfa[,,f,], 1:2, sum)
    df1 <- as.data.frame(x1) %>%
      rownames_to_column("Year") %>%
      pivot_longer(cols = c(`1`, `2`), names_to = "Season", values_to = "catch")
    df2 <- as.data.frame(x2) %>%
      rownames_to_column("Year") %>%
      pivot_longer(cols = c(`1`, `2`), names_to = "Season", values_to = "AF")
    df <- left_join(df1, df2, by = c("Year", "Season")) %>%
      mutate(Year = as.numeric(Year)) %>%
      mutate(violation = (catch > 0 & AF == 0) | (AF > 0 & catch == 0)) %>%
      filter(violation)
    if (nrow(df) > 0) warning("You have specified catch for a fishery using direct removals but no associated sliced AFs.")
    # data$catch_obs_ysf[39:44,1,f]
    # rowSums(data$af_sliced_ysfa[60:65,1,f,])
  }
}
df

# obs_fix <- repair_af_slices_weighted(
#   catch_ysf = data$catch_obs_ysf,
#   af_ysfa   = data$af_sliced_ysfa,
#   N = 3, p = 1
# )
# data$af_sliced_ysfa <- obs_fix$af_ysfa

# proj_fix <- repair_af_slices_weighted(
#   catch_ysf = proj_catch_ysf,
#   af_ysfa   = proj_af_sliced_ysfa,
#   N = 3, p = 1
# )
# proj_af_sliced_ysfa <- proj_fix$af_ysfa

# Parameters ----

parameters <- get_parameters(data = data)
parameters$par_log_h <- log(0.8)

exp(parameters$par_log_h)
exp(parameters$par_log_af_alpha)
exp(parameters$par_log_lf_alpha)

map <- get_map(parameters = parameters)
# map$par_log_af_alpha <- NULL
# map$par_log_lf_alpha <- factor(c(1, 2, 3, NA, 4))

data$priors <- get_priors(parameters = parameters)
evaluate_priors(parameters = parameters, priors = data$priors)

obj <- MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
bounds <- get_bounds(obj = obj, parameters = parameters)

unique(names(obj$par))
obj$report()$lp_lf
obj$report()$lp_af
obj$fn()

# Optimise a single grid cell ----

control <- list(eval.max = 10000, iter.max = 10000)

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
              lower = bounds$lower, upper = bounds$upper, control = control)
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he,
              lower = bounds$lower, upper = bounds$upper, control = control)

check_estimability(obj = obj)
# ce <- check_estimability(obj = obj)
# ce[[4]] %>% filter(Param_check != "OK")

obj$env$last.par.best[1:13]

# Inspect model outputs ----

plot_cpue(data = data, object = obj, nsim = 10)
plot_cpue_residuals(data = data, obj = obj, type = "OSA")

plot_aerial_survey(data = data, object = obj, nsim = 10)

plot_cpue_lf(data = data, object = obj)
plot_af(data = data, object = obj, fishery = "Indonesian")
plot_af(data = data, object = obj, fishery = "Australian")
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

plot_biomass_spawning(data_list = list(data), object_list = list(obj))

# Run MLE grid ----

do_run <- FALSE

if (do_run) {
  # grid_pars <- get_grid(parameters = parameters, m0 = c(0.4), m10 = c(0.065), h = c(0.8), psi = c(1.5, 1.75, 2))
  grid_pars <- get_grid(parameters = parameters, h = c(0.8), psi = c(1.75))
  grid_list <- run_grid(data = data, grid_parameters = grid_pars, bounds = bounds, map = map, control = control, parallel = FALSE)
  grid_check <- check_grid(grid = grid_list)
  # if some grid cells have not converged then you can run them again using rerun_grid
  # idc <- which(grid_check$grid_summary$Check != "All parameters are estimable")
  # grid_list <- rerun_grid(grid = grid_list, bounds = bounds, cells = idc, control = control) # this crashes my machine
  grid_cells <- sample_grid(grid = grid_list, seed = 42, prior_psi = c(1))
  grid_tmbfit <- grid_to_tmbfit(data = data, parameters = parameters, grid = grid_list, grid_parameters = grid_pars, grid_cells = grid_cells)
  save(grid_pars, grid_list, grid_check, grid_cells, grid_tmbfit, file = "grid_run4.rda")
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

length(grid_pars) # 108
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

plot_biomass_spawning(data_list = rep(list(data), times = 108), object_list = grid_list)

# extract_samples(grid_tmbfit)
post <- as.data.frame(grid_tmbfit)
pars <- grid_tmbfit$par_names[1:8]
mcmc_trace(x = post, pars = pars)

# MCMC for single grid cell ----

mcmc <- sample_snuts(
  obj = obj, metric = "auto", num_samples = 250, num_warmup = 750,
  # iter = 1000, warmup = 750, chains = 4, cores = 4,
  # obj = obj, metric = "auto", iter = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.99), init = "last.par.best",
  # lower = bnd$lower, upper = bnd$upper, # these bounds dont seem to work
  globals = sbt_globals())

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

# mcmc_list <- list()
# obj_list <- list()
grd <- expand.grid(h = c(0.55, 0.63, 0.72, 0.8), psi = c(1.5, 1.75, 2))

grid_dir <- "grid1"
dir.create(grid_dir)
for (i in 1:nrow(grd)) {
  parameters$par_log_h <- log(grd$h[i])
  parameters$par_log_psi <- log(grd$psi[i])
  obj <- MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he, lower = bnd$lower, upper = bnd$upper, control = control)
  mcmc <- sample_sparse_tmb(obj = obj, metric = "auto", iter = 1000, warmup = 750, chains = 4, cores = 4,
                            # lower = Lwr, upper = Upr, # these bounds dont seem to work
                            control = list(adapt_delta = 0.99), init = "last.par.best", globals = sbt_globals())
  save(opt, mcmc, file = paste0(grid_dir, "/grid", i, ".rda"))
  # obj_list[[i]] <- obj
  # mcmc_list[[i]] <- mcmc
}

# save(mcmc_list, obj_list, file = "mcmc_grid46.rda")

# Sensitivities ----

# To do at OMMP16

# Projections ----

# To do at OMMP16
