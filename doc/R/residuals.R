#' Standard deviation of the normalised residual (SDNR)
#' 
#' Calculate the standard deviation of normalised residuals (SDNR).
#' 
#' @param x a \code{numeric} vector of residuals.
#' @return a \code{numeric} value representing the standard deviation.
#' @importFrom stats sd
#' @export
#' 
sdnr <- function(x) sd(x)

#' Median of the absolute residual (MAR)
#' 
#' Calculate the median of absolute residuals (MAR).
#' 
#' @param x a \code{numeric} vector of residuals.
#' @return a \code{numeric} value representing the median of absolute values.
#' @importFrom stats median
#' @export
#' 
mar <- function(x) median(abs(x))

#' Plot CPUE residuals
#' 
#' Obtain the mean size (cm) at recruitment.
#' 
#' @param data data list
#' @param obj obj
#' @param type type of residuals, Raw, PIT, Pearson, OSA
#' @return a \code{ggplot}.
#' @import dplyr
#' @import ggplot2
#' @importFrom scales breaks_pretty
#' @importFrom RTMB oneStepPredict
#' @export
#' 
plot_cpue_residuals <- function(data, obj, type = "OSA") {
  
  cpue_sigma <- exp(obj$env$parList(obj$env$last.par.best)$par_log_cpue_sigma)
  
  df <- data.frame(
    year = data$cpue_year + data$first_yr,
    obs = data$cpue_obs, 
    pred = obj$report(obj$env$last.par.best)$cpue_pred, 
    sigma = sqrt(data$cpue_sd^2 + cpue_sigma^2)
  )
  
  if (type == "Raw") {
    ylab <- "Raw residual"
    cap <- "Raw residuals."
    df$resid <- df$obs - df$pred
    p <- ggplot(data = df) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "#00BA38") +
      geom_segment(aes(x = .data$year, y = .data$resid, xend = .data$year, yend = 0), color = "#619CFF")
  }
  
  if (type == "OSA") {
    ylab <- "OSA residual"
    cap <- "One step ahead (OSA) residuals for lognormal distribution."
    osa_cpue <- oneStepPredict(obj = obj, observation.name = "cpue_log_obs", method = "oneStepGeneric", trace = FALSE)
    df$resid <- osa_cpue$residual
  }

  if (type != "Raw") {
    res <- df %>% summarise(SDNR = sdnr(.data$resid), MAR = mar(.data$resid))
    cat(sprintf("\nDiagnostic Summary:\n  SDNR: %.5f  (Target: ~1.0)\n  MAR:  %.5f  (Target: ~0.67)\n\n", res$SDNR, res$MAR))
    p <- ggplot(data = df) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "#00BA38") + 
      geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#F8766D") +
      geom_segment(aes(x = .data$year, y = .data$resid, xend = .data$year, yend = 0), color = "#619CFF")
  }
  
  p <- p +
    geom_point(aes(x = .data$year, y = .data$resid), color = "#619CFF") +
    # scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 5)]) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Year", y = ylab, caption = cap)
  
  return(p)
}

#' Plot HSP residuals
#' 
#' @param data data list
#' @param obj obj
#' @return a \code{ggplot}.
#' @import dplyr
#' @import ggplot2
#' @importFrom RTMB oneStepPredict
#' @export
#' 
plot_hsp_residuals <- function(data, obj) {
  
  osa_res <- oneStepPredict(
    obj = obj, 
    observation.name = "hsp_nK", 
    method = "oneStepGeneric", 
    trace = FALSE
  )

  df <- data$hsp_obs %>%
    mutate(resid = osa_res$residual)
  # cohort1 cmin: 
  # cohort2 cmax: year of birth of oldest animal

  res <- df %>% summarise(SDNR = sdnr(.data$resid), MAR = mar(.data$resid))
  cat(sprintf("\nDiagnostic Summary:\n  SDNR: %.5f  (Target: ~1.0)\n  MAR:  %.5f  (Target: ~0.67)\n\n", res$SDNR, res$MAR))
  
  p <- ggplot(data = df) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "#00BA38") + 
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#F8766D") +
    geom_segment(aes(x = cmax, y = resid, xend = cmax, yend = 0), color = "#619CFF") +
    geom_point(aes(x = cmax, y = resid), color = "#619CFF") +
    facet_wrap(cmin ~ .) +
    # scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 5)]) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Year", y = "OSA residual", caption = "One step ahead (OSA) residuals for binomial distribution.")
  
  return(p)
}

#' Plot GT residuals
#' 
#' @param data data list
#' @param obj obj
#' @return a \code{ggplot}.
#' @import dplyr
#' @import ggplot2
#' @importFrom RTMB oneStepPredict
#' @export
#' 
plot_gt_residuals <- function(data, obj) {
  
  osa_res <- oneStepPredict(
    obj = obj, 
    observation.name = "gt_nrec", 
    method = "oneStepGeneric", 
    trace = FALSE
  )
  
  df <- data$gt_obs %>%
    mutate(resid = osa_res$residual)

  res <- df %>% summarise(SDNR = sdnr(resid), MAR = mar(resid))
  cat(sprintf("\nDiagnostic Summary:\n  SDNR: %.5f  (Target: ~1.0)\n  MAR:  %.5f  (Target: ~0.67)\n\n", res$SDNR, res$MAR))
  
  p <- ggplot(data = df) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "#00BA38") + 
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#F8766D") +
    geom_segment(aes(x = RelYear, y = resid, xend = RelYear, yend = 0), color = "#619CFF") +
    geom_point(aes(x = RelYear, y = resid), color = "#619CFF") +
    # facet_wrap(cmin ~ .) +
    # scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 5)]) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Year", y = "OSA residual", caption = "One step ahead (OSA) residuals for binomial distribution.")
  
  return(p)
}

#' Plot AF residuals
#' 
#' @param data data list
#' @param obj obj
#' @param fishery fishery
#' @return a \code{ggplot}.
#' @import dplyr
#' @import ggplot2
#' @importFrom RTMB oneStepPredict
#' @export
#' 
# library(compResidual)
# library(gridExtra)
plot_af_residuals <- function(data, obj, fishery = 1) {
  res_rep <- obj$report()
  af_pred <- res_rep$af_pred
  af_obs     <- data$af_obs
  af_n       <- data$af_n
  af_fishery <- data$af_fishery
  af_year    <- data$af_year
  af_min_age <- data$af_min_age
  af_max_age <- data$af_max_age
  n_obs <- nrow(af_obs)
  n_age <- ncol(af_obs)
  plot_df <- data.frame()
  for (i in seq_len(n_obs)) {
    if (af_fishery[i] != fishery) next
    amin <- af_min_age[i] + 1
    amax <- af_max_age[i] + 1
    # observed and predicted are stored on the same full-age scale
    obs_vec <- af_obs[i, ]
    pred_vec <- af_pred[i, ]
    # only keep the modeled age range
    obs_use <- obs_vec[amin:amax]
    pred_use <- pred_vec[amin:amax]
    # skip rows with missing or zero-sum data
    if (sum(obs_use, na.rm = TRUE) <= 0 || sum(pred_use, na.rm = TRUE) <= 0) next
    obs_use <- obs_use / sum(obs_use)
    pred_use <- pred_use / sum(pred_use)
    resid <- (obs_use - pred_use) / sqrt(pmax(pred_use * (1 - pred_use), 1e-8))
    row_df <- data.frame(
      Year = af_year[i] + data$first_yr - 1,
      Age = (amin:amax) - 1,
      Residual = resid
    )
    plot_df <- rbind(plot_df, row_df)
  }
  
  if (nrow(plot_df) == 0) {
    stop(paste("No calculated residuals found matching fishery =", fishery))
  }
  
  plot_df$Sign <- ifelse(plot_df$Residual >= 0, "Positive", "Negative")
  fishery_sdnr <- sd(plot_df$Residual, na.rm = TRUE)
  
  cat("\n=============================================\n")
  cat(paste("DIAGNOSTICS FOR FISHERY:", fishery, "\n"))
  cat(paste("Number of OSA Residuals Evaluated:", nrow(plot_df), "\n"))
  cat(paste("Empirical OSA SDNR:", round(fishery_sdnr, 4), "\n"))
  cat("=============================================\n\n")
  
  p1 <- ggplot(plot_df, aes(x = Year, y = Age, size = abs(Residual), fill = Sign)) +
    geom_point(shape = 21, alpha = 0.7) +
    scale_fill_manual(values = c("Positive" = "white", "Negative" = "black")) +
    scale_size_continuous(range = c(1, 8), name = "|Residual|") +
    labs(
      title = paste("OSA Residuals - Fishery", fishery),
      subtitle = paste("SDNR =", round(fishery_sdnr, 3)),
      x = "Year",
      y = "Age Class",
      caption = "White = Under-predicting (Positive), Black = Over-predicting (Negative)"
    ) +
    theme(panel.border = element_rect(fill = NA, color = "grey50"))
  
  p2 <- ggplot(data = plot_df, aes(sample = Residual)) +
    stat_qq() +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(title = "Q-Q Normal Distribution Check", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme(panel.border = element_rect(fill = NA, color = "grey50"))
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
  invisible(plot_df)
}
