# Clean up R
rm(list = ls())

# R libraries
library(tidyverse)
library(TMB)
library(sbt)

# setwd("tests")

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

# Likelihood vs. pdf  for ADMB multinomial distribution ----

# https://discourse.mc-stan.org/t/multinomial-with-non-integer-data/9220
multiF_stan <- function(x, p) {
  lp <- x * (p - log(sum(exp(p))))
  return(sum(lp))
}

# TMB
# Type dmultinom(vector<Type> x, vector<Type> p, int give_log=0) {
#   vector<Type> xp1 = x+Type(1);
#   Type logres = lgamma(x.sum() + Type(1)) - lgamma(xp1).sum() + (x*log(p)).sum();
#   if(give_log) return logres;
#   else return exp(logres);
# }
multiF_tmb <- function(x, p) {
  xp1 <- x + 1
  lp <- lgamma(sum(x) + 1) - sum(lgamma(xp1)) + sum(x * log(p))
  return(sum(lp))
}

multiF_casal <- function(x, p) {
  n <- sum(x)
  x <- x / sum(x)
  # lp <- log(factorial(n)) + sum(log(factorial(n * x)) - (n * x) * log(p)) # CASAL manual
  lp <- -lgamma(n + 1) + sum(lgamma((n * x) + 1) - (n * x * log(p)))
  # lp <- sum(lgamma((n * x) + 1) - (n * x * log(p)))
  return(-sum(lp))
}

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

# From CCSBT ADMB code
# lp(i) -= lf_n(i) * (x * log(pred)).sum(); // ln_like(iff) -= Nsamp(iff,iy)*((Nrobust+obs_len_freq_il(ii)(mbin,nbins)*log(Nrobust+pred_len_freq_il(ii)(mbin,nbins))));
# x += Type(1e-6);
# lp(i) += lf_n(i) * (x * log(x)).sum(); // mult_constant(iff) += Nsamp(iff,iy)*(1e-6+obs_len_freq_il(irec)(mbin,nbins))*log(1e-6+obs_len_freq_il(irec)(mbin,nbins));

multiF_age <- function(x, p) {
  Nsamp <- sum(x)
  obs_age_freq <- x / sum(x)
  obs_age_freq <- obs_age_freq + 1e-6
  pred_age_freq <- p + 1e-6
  lp <- -1 * sum(Nsamp * obs_age_freq * log(pred_age_freq))
  lp <- lp + sum(Nsamp * obs_age_freq * log(obs_age_freq))
  return(-lp)
}

multiF_len <- function(x, p) {
  Nsamp <- sum(x)
  obs_len_freq <- x / sum(x)
  pred_len_freq <- p + 1e-6
  lp <- -1 * sum(Nsamp * obs_len_freq * log(pred_len_freq)) # ln_like(iff) -= Nsamp(iff,iy)*((Nrobust+obs_len_freq_il(ii)(mbin,nbins)*log(Nrobust+pred_len_freq_il(ii)(mbin,nbins))));
  obs_len_freq <- obs_len_freq + 1e-6
  lp <- lp + sum(Nsamp * obs_len_freq * log(obs_len_freq)) # mult_constant(iff) += Nsamp(iff,iy)*(1e-6+obs_len_freq_il(irec)(mbin,nbins))*log(1e-6+obs_len_freq_il(irec)(mbin,nbins));
  return(-lp)
}

# Integer case

alpha <- 1:10
p <- 1:10
p <- p / sum(p)
x <- rmultinom(n = 100, size = 42, prob = p)
p <- rowSums(x)
p <- p / sum(p)
dim(x)

# compile(file = "dm.cpp")
# dyn.load(dynlib("dm"))
# 
# obj <- MakeADFun(data = list(x = x), parameters = list(prob = p, alpha = alpha), DLL = "dm")
# ll_tmb1 <- obj$report()$multinomial_dens
# sum(ll_tmb1 - ll_tmb)
# ll_stan <- apply(x, 2, FUN = function(x) multiF_stan(x, p))
ll_a <- apply(x, 2, FUN = function(x) multiF_age(x, p))
ll_l <- apply(x, 2, FUN = function(x) multiF_len(x, p))
ll_tmb <- apply(x, 2, FUN = function(x) multiF_tmb(x, p))
ll_casal <- apply(x, 2, FUN = function(x) multiF_casal(x, p))
ll_r <- apply(x, 2, FUN = function(x) dmultinom(x = x, prob = p, log = TRUE))
ll_r2 <- apply(x, 2, FUN = function(x) dmultinom_r(x = x, prob = p, log = TRUE))

par(mfrow = c(2, 3))
plot(ll_tmb, ll_r, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. R")
plot(ll_tmb, ll_r2, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. R (continuous)")
plot(ll_tmb, ll_casal, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. CASAL")
plot(ll_tmb, ll_a, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. CCSBT age")
plot(ll_tmb, ll_l, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. CCSBT length")
# plot(ll_tmb, ll_stan, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. Stan")

# Non-integer case

alpha <- 1:10
p <- 1:10
p <- p / sum(p)
x <- rmultinom(n = 100, size = 42, prob = p) + runif(n = 100, 0, 1)
p <- rowSums(x)
p <- p / sum(p)
dim(x)

# obj <- MakeADFun(data = list(x = x), parameters = list(prob = p, alpha = alpha), DLL = "dm")
# ll_tmb1 <- obj$report()$multinomial_dens
# sum(ll_tmb1 - ll_tmb)
# ll_stan <- apply(x, 2, FUN = function(x) multiF_stan(x, p))
ll_a <- apply(x, 2, FUN = function(x) multiF_age(x, p))
ll_l <- apply(x, 2, FUN = function(x) multiF_len(x, p))
ll_tmb <- apply(x, 2, FUN = function(x) multiF_tmb(x, p))
ll_casal <- apply(x, 2, FUN = function(x) multiF_casal(x, p))
ll_r <- apply(x, 2, FUN = function(x) dmultinom(x = x, prob = p, log = TRUE))
ll_r2 <- apply(x, 2, FUN = function(x) dmultinom_r(x = x, prob = p, log = TRUE))

par(mfrow = c(2, 3))
plot(ll_tmb, ll_r, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. R")
plot(ll_tmb, ll_r2, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. R (continuous)")
plot(ll_tmb, ll_casal, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. CASAL")
plot(ll_tmb, ll_a, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. CCSBT age")
plot(ll_tmb, ll_l, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. CCSBT length")
# plot(ll_tmb, ll_stan, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. Stan")

# Non-integer real-data case

x <- t(Data$lf_n * Data$lf_obs)
x <- x[,Data$lf_fishery == 4]
# x <- t(Data$af_n * Data$af_obs)
p <- rowSums(x)
p <- p / sum(p)
alpha <- 1:length(p)
# dim(x)
# plot(x[,1], type = "b")
# plot(p, type = "b")
# colSums(x) == Data$lf_n

# obj <- MakeADFun(data = list(x = x), parameters = list(prob = p, alpha = alpha), DLL = "dm")
# ll_tmb1 <- obj$report()$multinomial_dens
# sum(ll_tmb1 - ll_tmb)
# ll_stan <- apply(x, 2, FUN = function(x) multiF_stan(x, p))
ll_a <- apply(x, 2, FUN = function(x) multiF_age(x, p))
ll_l <- apply(x, 2, FUN = function(x) multiF_len(x, p))
ll_tmb <- apply(x, 2, FUN = function(x) multiF_tmb(x, p))
ll_casal <- apply(x, 2, FUN = function(x) multiF_casal(x, p))
ll_r <- apply(x, 2, FUN = function(x) dmultinom(x = x, prob = p, log = TRUE))
ll_r2 <- apply(x, 2, FUN = function(x) dmultinom_r(x = x, prob = p, log = TRUE))

par(mfrow = c(2, 3))
plot(ll_tmb, ll_r, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. R")
plot(ll_tmb, ll_r2, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. R (continuous)")
plot(ll_tmb, ll_casal, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. CASAL")
plot(ll_tmb, ll_a, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. CCSBT age")
plot(ll_tmb, ll_l, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. CCSBT length")
# plot(ll_tmb, ll_stan, col = 2); abline(a = 0, b = 1); title(main = "TMB vs. Stan")

# Likelihood vs. pdf  for normal distribution ----

normalF <- function(x, mu, sigma) { sum(-0.5 * log(sigma) - 0.5 * (x - mu)^2 / sigma) }

x <- matrix(rnorm(n = 100 * 10, 1, 1), nrow = 100, ncol = 10)

normalF(x[1,], 1, 1) # log likelihood function value for given x and mu=sd=1
sum(dnorm(x[1,], 1, 1, log = TRUE))

ll1 <- apply(x, 1, FUN = function(x) normalF(x, 1, 1))
ll2 <- apply(x, 1, FUN = function(x) sum(dnorm(x, 1, 1, log = TRUE)))

plot(ll1, ll2)

# Likelihood vs. pdf  for lognormal distribution ----

# This is the CPUE likelihood
# *cpue_resid = (log(cpue_obs) - cpue_log_pred) / cpue_sigma;
# vector<Type> lp = log(cpue_sigma) + Type(0.5) * square(*cpue_resid);
lnormalF1 <- function(x, mu, sigma) {
  res <- (log(x) - log(mu)) / sigma
  lp <- log(sigma) + 0.5 * res^2
  sum(lp)
}

# template<class Type>
#   Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log) {
#     Type resid = (log(x) - meanlog) / sdlog;
#     Type logans = sdlog + Type(0.5) * square(resid);
#     if (give_log) return logans;
#     else return exp(logans);
#   }
# VECTORIZE4_ttti(dlnorm);
# y = (log(x) - meanlog) / sdlog;
# return (give_log ?
#           -(M_LN_SQRT_2PI   + 0.5 * y * y + log(x * sdlog)) :
#           M_1_SQRT_2PI * exp(-0.5 * y * y)  /	 (x * sdlog));
lnormalF <- function(x, meanlog, sdlog, givelog = TRUE) {
  res <- (log(x) - meanlog) / sdlog
  if (givelog) {
    M_LN_SQRT_2PI <- log(sqrt(2 * pi))
    lp <- -(M_LN_SQRT_2PI + 0.5 * res^2 + log(x * sdlog))
  } else {
    M_1_SQRT_2PI <- 0.3989422804014327
    lp <- M_1_SQRT_2PI * exp(-0.5 * res^2) / (x * sdlog)
  }
  return(lp)
}

x <- matrix(rlnorm(n = 100 * 10, meanlog = log(2), sdlog = 0.1), nrow = 100, ncol = 10)
hist(x)

ll0 <- apply(x, 1, FUN = function(x) lnormalF1(x, 2, 0.1))
ll1 <- apply(x, 1, FUN = function(x) sum(lnormalF(x, log(2), 0.1, T)))
ll2 <- apply(x, 1, FUN = function(x) sum(dlnorm(x, meanlog = log(2), sdlog = 0.1, log = T)))

sum(ll0 - ll1)
sum(ll1 - ll2)

plot(ll0, ll1)
plot(ll1, ll2)
