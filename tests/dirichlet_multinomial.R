# Clean up R
rm(list = ls())

# R libraries
library(tidyverse)
library(brms)
library(TMB)
library(extraDistr)
library(sbt)

# setwd("tests")

compile(file = "dm.cpp")
dyn.load(dynlib("dm"))

attach(data_csv1)

Data1 <- list(last_yr = 2022, age_increase_M = 25,
              length_m50 = 150, length_m95 = 180, 
              catch_UR_on = 0, catch_surf_case = 1, catch_LL1_case = 1, 
              scenarios_surf = scenarios_surface, scenarios_LL1 = scenarios_LL1,
              sel_min_age_f = c(2, 2, 2, 8, 6, 0),
              # sel_max_age_f = c(17, 9, 17, 22, 25, 7),
              sel_max_age_f = c(17, 9, 17, 16, 25, 7),
              # sel_end_f = c(1, 0, 1, 1, 1, 0),
              sel_end_f = c(0, 0, 0, 0, 1, 0),
              sel_change_sd_fy = t(as.matrix(sel_change_sd[,-1])), 
              sel_smooth_sd_f = data_labrep1$sel.smooth.sd,
              hsp_switch = 1, hsp_false_negative = 0.7467647, 
              pop_switch = 1, 
              gt_switch = 1,
              cpue_switch = 1, cpue = cpue, cpue_a1 = 5, cpue_a2 = 17,
              aerial_switch = 4, aerial_tau = data_labrep1$tau.aerial, 
              troll_switch = 0, 
              lf_minbin = c(1, 1, 1, 11),
              tag_switch = 1, tag_var_factor = 1.82
)

Data <- get_data(data_in = Data1)

dmultinom_r <- function(x, size = NULL, prob, log = FALSE) {
  K <- length(prob)
  if (length(x) != K) 
    stop("x[] and prob[] must be equal length vectors.")
  if (any(!is.finite(prob)) || any(prob < 0) || (s <- sum(prob)) == 0) 
    stop("probabilities must be finite, non-negative and not all 0")
  prob <- prob/s
  # x <- as.integer(x + 0.5)
  if (any(x < 0)) 
    stop("'x' must be non-negative")
  N <- sum(x)
  if (is.null(size)) 
    size <- N
  else if (size != N) 
    stop("size != sum(x), i.e. one is wrong")
  i0 <- prob == 0
  if (any(i0)) {
    if (any(x[i0] != 0)) 
      return(if (log) -Inf else 0)
    if (all(i0)) 
      return(if (log) 0 else 1)
    # x <- x[!i0]
    # prob <- prob[!i0]
  }
  r <- lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
  if (log) 
    r
  else exp(r)
}

# Test on simulated data ----

alpha <- matrix(runif(3, 1, 3), nrow = 3, ncol = 1)
p <- brms::rdirichlet(n = 1, alpha = alpha[,1])
p
sum(p)
x <- rmultinom(n = 1, size = 10, prob = p)
x
colSums(x)

obj <- MakeADFun(data = list(x = x, prob = t(p), alpha = alpha), 
                 parameters = list(par = 1), DLL = "dm")

obj$report()$multinomial_dens
dmultinom(x = x[,1], prob = p, log = TRUE)

obj$report()$dirichlet_dens
brms::ddirichlet(x = p, alpha = alpha[,1], log = TRUE)

obj$report()$dirmult_dens
obj$report()$dirmult_dens_wham
extraDistr::ddirmnom(x = x[,1], size = 10, alpha = alpha[,1], log = TRUE)

# Test on real data ----

x <- t(Data$lf_n * Data$lf_obs)
x <- x[,Data$lf_fishery == 4]
# x <- t(Data$af_n * Data$af_obs)
p <- rowSums(x)
p <- p / sum(p)
p <- matrix(p, nrow = 25, ncol = 38)
alpha <- matrix(1:nrow(p), nrow = 25, ncol = 38)

obj <- MakeADFun(data = list(x = x, prob = p, alpha = alpha), 
                 parameters = list(par = 1), DLL = "dm")

obj$report()$multinomial_dens
# dmultinom(x = x[,1], prob = p, log = TRUE)
apply(x, 2, FUN = function(x) dmultinom_r(x = x, prob = p[,1], log = TRUE))

obj$report()$dirmult_dens
obj$report()$dirmult_dens_wham
sum(obj$report()$dirmult_dens - obj$report()$dirmult_dens_wham)
extraDistr::ddirmnom(x = x[,1], size = sum(x[,1]), alpha = alpha, log = TRUE)

prob <- x + 1e-6
prob <- t(t(prob) / colSums(prob))
colSums(prob)
alpha <- t(t(p) * colSums(x))
colSums(p)
colSums(x)
colSums(alpha)

data <- list(x = x, prob = prob, alpha = alpha)

obj <- MakeADFun(data = data, parameters = list(par = 1), DLL = "dm")

res <- rep(NA, ncol(x))
for (i in 1:ncol(x)) {
  res[i] <- brms::ddirichlet(x = data$prob[,i], alpha = data$alpha[,i], log = TRUE)
}

obj$report()$dirichlet_dens
res
