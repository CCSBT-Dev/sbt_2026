# Initial set up ----

rm(list = ls())

# remotes::install_github("janoleko/RTMBdist")
# remotes::install_github("andrjohns/StanEstimators")
# remotes::install_github("noaa-afsc/SparseNUTS")
# remotes::install_github(repo = "quantifish/sbt")

if (basename(getwd()) != "OMMP16") {
  setwd(file.path(getwd(), "OMMP16"))
}

library(tidyverse)
library(sbt)
library(SparseNUTS)
library(bayesplot)

attach(data_csv1)
lr <- data_labrep1

theme_set(theme_bw())

sbt_model_oldsel <- function(parameters, data) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(data, parameters, warn = FALSE)
  
  # Natural mortality ----
  
  par_m0 <- exp(par_log_m0)
  par_m4 <- exp(par_log_m4)
  par_m10 <- exp(par_log_m10)
  par_m30 <- exp(par_log_m30)
  M_a <- get_M(min_age, max_age, age_increase_M, par_m0, par_m4, par_m10, par_m30)
  REPORT(par_m0)
  REPORT(par_m4)
  REPORT(par_m10)
  REPORT(par_m30)
  REPORT(M_a)
  # M_a2 <- get_M_length(min_age, max_age, age_increase_M, par_m0, par_m4, par_m10, par_m30, length_mu_ysa)
  # plot(M_a, ylim = c(0, 0.45))
  # lines(M_a2, col = 2)
  
  # Selectivity ----
  
  # par_log_sel_fya <- list(par_log_sel_1, par_log_sel_2, par_log_sel_3, par_log_sel_4, par_log_sel_5, par_log_sel_6, par_log_sel_7)
  # lp_sel <- get_selectivity_prior(par_sel_rho_y, par_sel_rho_a, par_log_sel_sigma, par_log_sel_fya)
  # sel_fya <- get_selectivity(n_age, max_age, first_yr, first_yr_catch, sel_min_age_f, sel_max_age_f, sel_end_f, sel_change_year_fy, par_log_sel_fya)
  sel_min_age_f <- sel_min_age_f[1:6]
  sel_max_age_f <- sel_max_age_f[1:6]
  sel_end_f <- sel_end_f[1:6]
  sel_change_year_fy <- sel_change_year_fy[1:6, , drop = FALSE]
  sel_change_sd_fy <- sel_change_sd_fy[1:6, , drop = FALSE]
  sel_smooth_sd_f <- sel_smooth_sd_f[1:6]
  first_yr_catch_f <- first_yr_catch_f[1:6]
  
  sel_fya <- get_selectivity_v1(n_age, max_age, first_yr, first_yr_catch, 
                                sel_min_age_f, sel_max_age_f, sel_end_f, 
                                sel_change_year_fy, par_sels_init_i, par_sels_change_i)
  
  lp_sel <- get_sel_like_v1(first_yr, first_yr_catch_f, 
                            sel_min_age_f = sel_min_age_f, 
                            sel_max_age_f = sel_max_age_f, 
                            sel_change_year_fy, 
                            sel_change_sd_fy, 
                            sel_smooth_sd_f, 
                            par_sels_init_i, 
                            par_sels_change_i, 
                            sel_fya = sel_fya)
  
  new_fya <- array(NA, dim = c(7, 92, 31))
  new_fya[1:6,,] <- sel_fya
  new_fya[7,,] <- sel_fya[1,,]
  sel_fya <- new_fya
  
  REPORT(sel_fya)
  
  # Recruitment ----
  
  sigma_r <- exp(par_log_sigma_r)
  tau_ac2 <- get_rho(first_yr, last_yr, par_rdev_y)
  lp_rec <- get_recruitment_prior(par_rdev_y, sigma_r, tau_ac2)
  rdev_y <- par_rdev_y
  # for (y in (n_year - 2):n_year) rdev_y[y] <- tau_ac2 * rdev_y[y - 1] + par_rdev_y[y]
  
  # Spawning output per recruit ----
  
  phi_ya <- get_phi(par_log_psi, length_m50, length_m95, length_mu_ysa, length_sd_a, dl_yal)
  # for (i in 2:93) phi_ya[i,] <- phi_ya[1,]
  REPORT(phi_ya)
  
  # Main population loop ----
  
  B0 <- exp(par_log_B0)
  par_h <- exp(par_log_h)
  init <- get_initial_numbers(B0 = B0, h = par_h, M_a = M_a, phi_ya = phi_ya)
  R0 <- init$R0
  alpha <- init$alpha
  beta <- init$beta
  
  dyn <- do_dynamics(first_yr, first_yr_catch, 
                     B0 = B0, R0 = R0, alpha, beta, h = par_h, sigma_r, rdev_y, M_a, phi_ya,
                     init_number_a = init$Ninit,
                     removal_switch_f, catch_obs_ysf, sel_fya, weight_fya, af_sliced_ysfa)
  
  hrate_ysfa  <- dyn$hrate_ysfa
  hrate_ysa  <- dyn$hrate_ysa
  catch_pred_fya <- dyn$catch_pred_fya
  catch_pred_ysf <- dyn$catch_pred_ysf
  number_ysa <- dyn$number_ysa
  spawning_biomass_y <- dyn$spawning_biomass_y
  # recruitment_y <- dyn$recruitment_y
  lp_penalty <- dyn$lp_penalty
  
  # Likelihoods and priors ----
  
  lp_af <- get_age_like(af_switch, removal_switch_f, af_year, af_fishery, af_min_age, af_max_age, af_obs, af_n, par_log_af_alpha, catch_pred_fya)
  lp_lf <- get_length_like(lf_switch, removal_switch_f, lf_year, lf_season, lf_fishery, lf_minbin, lf_obs, lf_n, par_log_lf_alpha, catch_pred_fya, alk_ysal)
  lp_cpue_lf <- get_cpue_length_like(lf_switch, cpue_years, cpue_lfs, cpue_n, par_log_lf_alpha, number_ysa, sel_fya, alk_ysal)
  
  lp_troll <- get_troll_like(troll_switch, troll_years, troll_obs, troll_sd, par_log_troll_tau, number_ysa)
  lp_tags <- get_tag_like(tag_switch, min_K + 1, n_K, n_T, n_I, n_J, first_yr, M_a, hrate_ysa,
                          tag_release_cta, tag_recap_ctaa,
                          minI = tag_rel_min_age, maxI = tag_rel_max_age, maxJ = tag_recap_max_age,
                          shed1 = tag_shed_immediate, shed2 = tag_shed_continuous,
                          tag_rep_rates_ya, tag_H_factor = exp(par_log_tag_H_factor), tag_var_factor)
  
  x <- get_aerial_survey_like(aerial_switch, aerial_survey, aerial_cov, first_yr, par_log_aerial_tau, par_log_aerial_sel, number_ysa, weight_fya)
  lp_aerial <- x$lp
  lp_aerial_tau <- x$lp_aerial_tau
  
  # The likelihoods below return lists as they are also used in simulation
  x <- get_cpue_like(cpue_switch, cpue_years, cpue_obs, cpue_sd, cpue_a1, cpue_a2, par_log_cpue_q, par_log_cpue_sigma, par_log_cpue_omega, par_cpue_creep, creep_init = 1, number_ysa, sel_fya)
  lp_cpue <- x$lp
  x <- get_POP_like(pop_switch, pop_obs, paly, phi_ya, spawning_biomass_y)
  lp_pop <- x$lp
  x <- get_HSP_like(hsp_switch, hsp_obs, hsp_false_negative, first_yr, par_log_hsp_q, number_ysa, phi_ya, M_a, spawning_biomass_y, hrate_ysa)
  lp_hsp <- x$lp
  x <- get_GT_like(gt_switch, gt_obs, first_yr, par_log_gt_q, number_ysa)
  lp_gt <- x$lp
  
  lp_prior <- evaluate_priors(parameters, priors)
  
  nll <- lp_prior + sum(lp_sel) + lp_rec + lp_penalty +
    sum(lp_af) + sum(lp_lf) + sum(lp_cpue_lf) + 
    sum(lp_cpue) + sum(lp_aerial) + lp_aerial_tau + sum(lp_troll) + 
    sum(lp_tags) + sum(lp_pop) + sum(lp_hsp) + sum(lp_gt)
  
  # Reporting ----
  
  REPORT(B0)
  REPORT(R0)
  REPORT(alpha)
  REPORT(beta)
  REPORT(par_h)
  REPORT(sigma_r)
  # ADREPORT(sigma_r)
  REPORT(tau_ac2)
  REPORT(par_rdev_y)
  # REPORT(rec_dev_y)
  REPORT(rdev_y)
  
  REPORT(par_log_psi)
  REPORT(spawning_biomass_y)
  REPORT(number_ysa)
  REPORT(hrate_ysa)
  REPORT(catch_pred_ysf)
  REPORT(catch_pred_fya)
  
  REPORT(lp_sel)
  REPORT(lp_rec)
  REPORT(lp_prior)
  REPORT(lp_penalty)
  REPORT(lp_af)
  REPORT(lp_lf)
  REPORT(lp_cpue_lf)
  REPORT(lp_cpue)
  REPORT(lp_aerial)
  REPORT(lp_aerial_tau)
  REPORT(lp_troll)
  REPORT(lp_tags)
  REPORT(lp_pop)
  REPORT(lp_hsp)
  REPORT(lp_gt)
  
  return(nll)
}

data <- list(
  last_yr = 2022, age_increase_M = 25,
  length_m50 = 150, length_m95 = 180, 
  catch_UR_on = 0, catch_surf_case = 1, catch_LL1_case = 1, 
  scenarios_surf = scenarios_surface, scenarios_LL1 = scenarios_LL1,
  # removal_switch_f = c(0, 0, 0, 1, 0, 0), # 0=harvest rate, 1=direct removals
  removal_switch_f = c(0, 0, 0, 0, 0, 0), # 0=harvest rate, 1=direct removals
  sel_min_age_f = c(2, 2, 2, 8, 6, 0, 2),
  sel_max_age_f = c(17, 9, 17, 22, 25, 7, 17),
  # sel_max_age_f = c(17, 9, 17, 21, 25, 7, 17),
  sel_end_f = c(1, 0, 1, 1, 1, 0, 1), # 0=zero, 1=constant
  sel_LL1_yrs = c(1952, 1957, 1961, 1965, 1969, 1973, 1977, 1981, 1985, 1989, 1993, 1997, 2001, 2006, 2007, 2008, 2011, 2014, 2017, 2020, 2023),
  sel_LL2_yrs = c(1969, 2001, 2005, 2008, 2011, 2014, 2017, 2020),
  sel_LL3_yrs = c(1954, 1961, 1965, 1969, 1970, 1971, 2005, 2006, 2007),
  sel_LL4_yrs = c(1953),
  sel_Ind_yrs = c(1976, 1995, 1997, 1999, 2002, 2004, 2006, 2008, 2010, 2012:2022),
  sel_Aus_yrs = c(1952, 1969, 1973, 1977, 1981, 1985, 1989, 1993, 1997:2022),
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
data$sel_change_sd_fy <- t(as.matrix(sel_change_sd[,-1]))
data$sel_smooth_sd_f <- lr$sel.smooth.sd

# Parameters ----

parameters <- get_parameters(data = data)
parameters$par_sel_rho_y <- NULL
parameters$par_sel_rho_a <- NULL
parameters$par_log_sel_sigma <- NULL
parameters$par_log_sel_1 <- NULL
parameters$par_log_sel_2 <- NULL
parameters$par_log_sel_3 <- NULL
parameters$par_log_sel_4 <- NULL
parameters$par_log_sel_5 <- NULL
parameters$par_log_sel_6 <- NULL
parameters$par_log_sel_7 <- NULL
parameters$par_sels_init_i <- data_par1$par_sels_init_i
parameters$par_sels_change_i <- data_par1$par_sels_change_i

map <- list()
map[["par_log_psi"]] <- factor(NA)
map[["par_log_m0"]] <- factor(NA)
map[["par_log_m10"]] <- factor(NA)
map[["par_log_h"]] <- factor(NA)
map[["par_log_sigma_r"]] <- factor(NA)
map[["par_log_cpue_sigma"]] <- factor(NA)
map[["par_log_cpue_omega"]] <- factor(NA)
map[["par_cpue_creep"]] <- factor(NA)
map[["par_log_aerial_tau"]] <- factor(NA)
map[["par_log_aerial_sel"]] <- factor(rep(NA, 2))
map[["par_log_troll_tau"]] <- factor(NA)
map[["par_log_gt_q"]] <- factor(NA)
map[["par_log_hsp_q"]] <- factor(NA)

data$priors <- get_priors(parameters = parameters)
evaluate_priors(parameters = parameters, priors = data$priors)

obj <- MakeADFun(func = cmb(sbt_model_oldsel, data), parameters = parameters, map = map)
bounds <- get_bounds(obj = obj, parameters = parameters)

unique(names(obj$par))
obj$report()$lp_lf
obj$report()$lp_af
obj$fn()
