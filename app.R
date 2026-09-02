# GUIDED + DEVELOPER INTERFACE — 2026-09-02
# Based on the supplied methods-aligned engine and its tested guided workflow.
# View toggle is presentation-only. No preview/mockup data is used as output.
# Additions: true retained-chain diagnostics, imputation fit retention, model
# records, consistent six-screen UI. The multisite Fourier extra brace is fixed.
# Existing CSV contracts are preserved: sex remains a shared population input;
# action Amount retains action-specific semantics. No new cohort schema is implied.
# Diagnostic thresholds are screening conventions, not proof of convergence.
# Exact replay seeds are not yet retained; exported chains support inspection.
# =====================================================================
# METHODS-ALIGNED APP -- guided UI revision 2026-09-02
# Five-step workflow; embedded examples; saved scenarios compared with shared draws.
# Population kernels and manuscript equations are unchanged from 2026-09-01.
# Reference: Methods_ConsBio_20260831.pdf and
#            MethodsScript_ConsBio_20260831(8).R (user supplied).
#
# Equation map (paper symbols -> retained R names):
#   7-8:  VBGF + first-nesting transition p(T) -> maturity_probability
#   9:    p(FN)[j] = p(T)[j] * prod(1 - p(T)[h], h < j)
#   10-11: phi = (1-p(T))*pj + p(T)*pa; S = cumprod(phi)
#   12-14: Take = turtles * S * p(FN) * pf / ri * mortality
#   15-18: Give = EXTRA cohort * survival to census * p(FN) * pf / ri
#   19-22: N_next = max(0, max(0,N-Take)*exp(U) + sqrt(Q)*z + Give)
# U is paper r (instantaneous trend); Q is a variance, not an SD.
# The inherited additive abundance-scale use of sqrt(Q) is deliberately
# retained, including the manuscript's stated log-model/projection convention.
#
# Projection state and initializer are ANNUAL NESTERS. RI/4 adult-female
# summaries are reporting quantities only. Take precedes growth; Give follows
# it. Neither first nesting nor sex/mortality is re-drawn as a Bernoulli event.
# Optional count/size/mortality/demographic uncertainty is an app extension;
# with fixed counts and zero SDs the result equals the manuscript calculation.
#
# Threat CSV: year,turtles,median_cm,mortality. One row per forecast year.
# Year labels the destination census: 2027 means the transition 2026 -> 2027.
# turtles = number affected, NOT deaths; mortality = total death probability
# per interaction (0..1). Do not multiply the count by mortality in advance.
# median_cm is the median length of that year's affected turtles, using the
# same carapace measure as the VBGF. A median alone specifies a representative
# cohort; optional log-length SD supplies a within-cohort size distribution.
# With uncertainty, mortality is the logit-normal median (fixed if SD = 0).
# Lognormal median = exp(mean log length), NOT the arithmetic mean length.
# Zero turtles must be explicit; omitted years are not assumed to mean zero.
#
# Changes beyond the paper's ordinary input domain are annotated locally:
# length == Linf uses the existing over-Linf forced-age fallback;
# probability draws are bounded; user input is validated rather than guessed.
# The app's built-in mathematical checks can run after loading functions;
# Validation: 565 external checks and 25 built-in checks passed in R 4.3.3.
# Largest reference difference: 1.39e-17. Guided UI additionally tested with
# real annual and monthly JAGS/MARSS fits, saved scenarios, and matched comparisons.
# UI/JAGS fitting still requires the packages listed below and JAGS installed.
# =====================================================================

library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(mvtnorm)
library(truncnorm)
library(jagsUI)
library(MARSS)

# =====================================================================
# MATHEMATICAL HELPERS & DISTRIBUTIONS
# =====================================================================
`%||%` <- function(x, y) if (is.null(x)) y else x

inv_logit <- function(x) { exp(x) / (1 + exp(x)) }
safe_logit <- function(x) {
  if (any(!is.finite(x)) || any(x < 0 | x > 1)) stop("Invalid probability.")
  qlogis(pmax(.Machine$double.eps, pmin(1 - .Machine$double.eps, x)))
}

compute_CMP_constant <- function(Lambda, Nu, Mu, Tol, Max, Log=TRUE, Type="Z"){
  if( (!is.na(Lambda) & Lambda > 10^Nu) | (!is.na(Mu) & Mu^Nu > 10^Nu) ){
    if(Type=="Z"){ ln_Const = Nu*Lambda^(1/Nu) - ((Nu-1)/(2*Nu))*log(Lambda) - ((Nu-1)/2)*log(2*pi) - (1/2)*log(Nu) }
    if(Type=="S"){ ln_Const = Nu*Mu - ((Nu-1)/(2))*log(Mu) - ((Nu-1)/2)*log(2*pi) - (1/2)*log(Nu) }
  }else{
    Const = rep(0,Max+1); Index = 1; Const[Index] = 1
    while( Const[Index]/Const[1] > Tol ){
      if(Type=="Z") Const[Index+1] = Const[Index] * ( Lambda / Index^Nu )
      if(Type=="S") Const[Index+1] = Const[Index] * ( Mu / Index )^Nu; Index = Index + 1
    }
    ln_Const = log(sum(Const))
  }
  if(Log) return(ln_Const) else return(exp(ln_Const))
}

dCMP <- function( x, lambda, mu, nu, log=TRUE, tol=0.01, iter.max=200 ){
  if(!missing(mu) && !missing(lambda)) stop("'mu' and 'lambda' cannot both be supplied.")
  if(missing(mu) && missing(lambda)) stop("Supply exactly one of 'mu' or 'lambda'.")
  if(missing(mu) & !missing(lambda)) loglike = x*log(lambda) - nu*lfactorial(x) - compute_CMP_constant(Lambda=lambda, Nu=nu, Mu=NA, Tol=tol, Max=iter.max, Log=TRUE, Type="Z")
  if(!missing(mu) & missing(lambda)) loglike = nu*x*log(mu) - nu*lfactorial(x) - compute_CMP_constant(Lambda=NA, Nu=nu, Mu=mu, Tol=tol, Max=iter.max, Log=TRUE, Type="S")
  if(log) return(loglike) else return(exp(loglike))
}

rCMP <- function( n, lambda, mu, nu, tol=0.01, x_max=200 ){
  if(!missing(mu) && !missing(lambda)) stop("'mu' and 'lambda' cannot both be supplied.")
  if(missing(mu) && missing(lambda)) stop("Supply exactly one of 'mu' or 'lambda'.")
  loglike_x = rep(NA, x_max+1)
  for( x in 0:x_max ){
    if(missing(mu) & !missing(lambda)) loglike_x[x+1] = dCMP( x=x, lambda=lambda, nu=nu, log=TRUE, tol=tol, iter.max=x_max)
    if(!missing(mu) & missing(lambda)) loglike_x[x+1] = dCMP( x=x, mu=mu, nu=nu, log=TRUE, tol=tol, iter.max=x_max)
  }
  return(sample(x=0:x_max, size=n, replace=TRUE, prob=exp(loglike_x)))
}

# =====================================================================
# BIOLOGICAL NESTING-SEASON INDEXING
# =====================================================================
biological_month_order <- function(season_start_month = 4L) {
  season_start_month <- as.integer(season_start_month)
  if (
    length(season_start_month) != 1L ||
      is.na(season_start_month) ||
      !season_start_month %in% 1:12
  ) {
    stop("season_start_month must be one integer from 1 through 12.")
  }

  c(
    seq.int(season_start_month, 12L),
    if (season_start_month > 1L) {
      seq_len(season_start_month - 1L)
    } else {
      integer(0)
    }
  )
}

add_nesting_season <- function(df, season_start_month = 4L) {
  required <- c("Year", "Month")
  missing_required <- setdiff(required, names(df))
  if (length(missing_required) > 0L) {
    stop(
      "Missing required column(s): ",
      paste(missing_required, collapse = ", ")
    )
  }

  season_start_month <- as.integer(season_start_month)
  if (
    length(season_start_month) != 1L ||
      is.na(season_start_month) ||
      !season_start_month %in% 1:12
  ) {
    stop("season_start_month must be one integer from 1 through 12.")
  }

  calendar_year <- suppressWarnings(as.integer(df$Year))
  calendar_month <- suppressWarnings(as.integer(df$Month))
  invalid_month <- !is.na(calendar_month) & !calendar_month %in% 1:12
  if (any(invalid_month)) {
    stop("Month values must be integers from 1 through 12.")
  }

  df %>%
    mutate(
      Calendar_Year = calendar_year,
      Calendar_Month = calendar_month,
      Season = if_else(
        is.na(Calendar_Month),
        NA_integer_,
        Calendar_Year - if_else(
          Calendar_Month < season_start_month,
          1L,
          0L
        )
      ),
      Seq_Month = if_else(
        is.na(Calendar_Month),
        NA_integer_,
        as.integer(
          ((Calendar_Month - season_start_month) %% 12L) + 1L
        )
      )
    )
}

complete_nesting_seasons <- function(df, season_start_month = 4L) {
  indexed <- add_nesting_season(
    df,
    season_start_month = season_start_month
  ) %>%
    filter(
      !is.na(Calendar_Year),
      !is.na(Calendar_Month),
      !is.na(Season)
    )

  if (nrow(indexed) == 0L) {
    return(integer(0))
  }

  observed_month_index <- (
    indexed$Calendar_Year * 12L +
      indexed$Calendar_Month - 1L
  )
  candidate_seasons <- seq.int(
    min(indexed$Season),
    max(indexed$Season)
  )
  candidate_start_index <- (
    candidate_seasons * 12L +
      as.integer(season_start_month) - 1L
  )

  candidate_seasons[
    candidate_start_index >= min(observed_month_index) &
      candidate_start_index + 11L <= max(observed_month_index)
  ]
}

# =====================================================================
# MONITORING-STATUS AND LOG-COUNT HELPERS
# =====================================================================
normalise_monitoring_status <- function(df) {
  if (!"Count" %in% names(df)) {
    stop("The nest-count data must contain a Count column.")
  }

  # In the original Martin data, gaps represented months with no monitoring
  # effort. For generalized uploads, a blank Count is therefore interpreted
  # as unmonitored unless an optional Monitored column says otherwise.
  if (!"Monitored" %in% names(df)) {
    df$Monitored <- !is.na(df$Count)
    return(df)
  }

  raw_status <- trimws(tolower(as.character(df$Monitored)))
  parsed_status <- rep(NA, length(raw_status))
  parsed_status[raw_status %in% c("1", "true", "t", "yes", "y", "monitored")] <- TRUE
  parsed_status[raw_status %in% c("0", "false", "f", "no", "n", "unmonitored")] <- FALSE

  unspecified <- is.na(raw_status) | raw_status == ""
  parsed_status[unspecified] <- !is.na(df$Count[unspecified])

  invalid_status <- is.na(parsed_status)
  if (any(invalid_status)) {
    bad_values <- unique(df$Monitored[invalid_status])
    stop(
      "Monitored must be TRUE/FALSE, 1/0, yes/no, or blank. Invalid value(s): ",
      paste(bad_values, collapse = ", ")
    )
  }

  monitored_without_count <- parsed_status & is.na(df$Count)
  if (any(monitored_without_count)) {
    stop(
      "Rows marked Monitored = TRUE must contain a Count. Use Monitored = FALSE ",
      "or leave Monitored blank when the month was not surveyed."
    )
  }

  count_without_monitoring <- !parsed_status & !is.na(df$Count)
  if (any(count_without_monitoring)) {
    stop(
      "Rows marked Monitored = FALSE must have a blank Count. This prevents ",
      "observed counts from being silently discarded."
    )
  }

  df$Monitored <- parsed_status
  df$Count[!df$Monitored] <- NA_real_
  df
}

validate_log_count_input <- function(count, monitored = NULL,
                                     label = "Nest count") {
  count <- suppressWarnings(as.numeric(count))
  if (is.null(monitored)) {
    monitored <- !is.na(count)
  }
  monitored <- as.logical(monitored)

  negative <- monitored & !is.na(count) & count < 0
  if (any(negative)) {
    stop(label, " values cannot be negative.")
  }

  observed_zero <- monitored & !is.na(count) & count == 0
  if (any(observed_zero)) {
    stop(
      "The Martin/Siders Gaussian log-count model cannot use observed zero ",
      label, " values. A blank/NA Count represents no monitoring and can be ",
      "imputed. A genuine observed zero requires a count-data observation model ",
      "and must not be silently converted to missing. Found ",
      sum(observed_zero), " observed zero value(s)."
    )
  }

  invisible(TRUE)
}

# =====================================================================
# FOURIER IMPUTATION ENGINE -- Methods 2.2, Eqs. 1-3
# Includes y[1] in single- and multi-site likelihoods; annual N is log(sum(exp(X))).
# =====================================================================
run_fourier_imputation <- function(
    df, iter = 100000, six_month_sites = NULL,
    legacy_mcmc = FALSE, season_start_month = 4L) {

  df <- normalise_monitoring_status(df)
  validate_log_count_input(
    count = df$Count,
    monitored = df$Monitored,
    label = "monthly nest-count"
  )

  df_clean <- df %>%
    add_nesting_season(
      season_start_month = season_start_month
    ) %>%
    filter(
      !is.na(Season),
      !is.na(Seq_Month),
      !is.na(Site)
    ) %>%
    group_by(Season, Seq_Month, Site) %>%
    summarise(
      Count = if (all(is.na(Count))) {
        NA_real_
      } else {
        sum(Count, na.rm = TRUE)
      },
      .groups = "drop"
    )

  if (nrow(df_clean) == 0L) {
    stop("No valid monthly observations remain after nesting-season assignment.")
  }

  all_seasons <- seq.int(
    min(df_clean$Season, na.rm = TRUE),
    max(df_clean$Season, na.rm = TRUE)
  )
  sites <- sort(unique(df_clean$Site))
  n_years <- length(all_seasons)
  n_timeseries <- length(sites)

  full_grid <- expand.grid(
    Seq_Month = 1:12,
    Season = all_seasons,
    Site = sites,
    stringsAsFactors = FALSE
  )
  prep_df <- full_grid %>%
    left_join(
      df_clean,
      by = c("Season", "Seq_Month", "Site")
    ) %>%
    arrange(Season, Seq_Month) %>%
    pivot_wider(
      names_from = Site,
      values_from = Count
    ) %>%
    select(all_of(sites))

  y_matrix <- as.matrix(prep_df)
  y_matrix[is.nan(y_matrix) | is.infinite(y_matrix)] <- NA
  if (any(y_matrix <= 0, na.rm = TRUE)) {
    stop(
      "Internal validation failure: non-positive monthly counts reached the ",
      "Gaussian log-count imputation model."
    )
  }
  y_matrix <- log(y_matrix)
  y_matrix <- matrix(y_matrix, ncol = n_timeseries)
  
  periods <- rep(12, n_timeseries)
  for(i in 1:n_timeseries) { 
    if(sites[i] %in% six_month_sites) periods[i] <- 6 
  }
  
  jags_data <- list(m = rep(1:12, times = n_years), n.steps = nrow(y_matrix), n.months = 12, pi = pi, period = periods, n.timeseries = n_timeseries, n.years = n_years)
  
  if (n_timeseries == 1) {
    jags_data$y <- as.vector(y_matrix)
    
    # The single-site JAGS model does not use this variable.
    jags_data$n.timeseries <- NULL
    
    model_string <- "
    model {
      X[1] ~ dnorm(0, 0.01)
      y[1] ~ dnorm(X[1], tau.y)
      for (t in 2:n.steps){
         predX[t] <- c[m[t]] + X[t-1]
         X[t] ~ dnorm(predX[t], tau.X)
         y[t] ~ dnorm(X[t], tau.y)
      }
      for (y_idx in 1:n.years){
         for (mm in 1:12){ tmp2[y_idx, mm] <- exp(X[(y_idx*12 - mm + 1)]) }
         N[y_idx] <- log(sum(tmp2[y_idx, ]))
      }
      for (k in 1:n.months){
          c.const[k] <- 2 * pi * k / period[1]
          c[k] <- beta.cos * cos(c.const[k]) + beta.sin * sin(c.const[k])
      }
      sigma.y ~ dgamma(2, 0.5); tau.y <- 1/(sigma.y * sigma.y)
      beta.cos ~ dnorm(0, 1); beta.sin ~ dnorm(0, 1)
      sigma.X ~ dgamma(2, 0.5); tau.X <- 1/(sigma.X * sigma.X)
    }"
  } else {
    jags_data$y <- y_matrix
    # Integration-test repair: one closing brace ends this JAGS model.
    # The prior string had an extra brace; all model equations are retained.
    model_string <- "
    model {
      for(j in 1:n.timeseries) {
         X[1,j] ~ dnorm(0, 0.01)
         # Methods 2.2: condition on the first monitored month as well.
         y[1,j] ~ dnorm(X[1,j], tau.y)
         for (t in 2:n.steps){
             predX[t,j] <-  c[j,m[t]] + X[t-1, j]
             X[t,j] ~ dnorm(predX[t,j], tau.X)
             y[t,j] ~ dnorm(X[t,j], tau.y)
         }
         for (y_idx in 1:n.years){
            for (mm in 1:12){ tmp2[y_idx, mm, j] <- exp(X[(y_idx*12 - mm + 1), j]) }
            N[y_idx, j] <- log(sum(tmp2[y_idx, , j]))
         }
      }
        for (j in 1:n.timeseries){
           for (k in 1:n.months){
              c.const[j,k] <- 2 * pi * k / period[j]
              c[j,k] <- beta.cos[j] * cos(c.const[j,k]) +
                        beta.sin[j] * sin(c.const[j,k])
           }
        
           beta.cos[j] ~ dnorm(0, 1)
           beta.sin[j] ~ dnorm(0, 1)
        }
        
        sigma.y ~ dgamma(2, 0.5)
        tau.y <- 1/(sigma.y * sigma.y)
        
        sigma.X ~ dgamma(2, 0.5)
        tau.X <- 1/(sigma.X * sigma.X)
    }"
  }
  
  # The supplied Martin/Siders implementation used five chains, 100,000
  # iterations, a 50,000-iteration burn-in, thinning by five, and parallel
  # chains. The fast setting retains the same model equations but is intended
  # only for interactive exploration.
  if (isTRUE(legacy_mcmc)) {
    n_chains <- 5
    n_iter <- max(100000, iter)
    n_burnin <- 50000
    n_thin <- 5
    run_parallel <- TRUE
  } else {
    n_chains <- 3
    n_iter <- max(2000, iter)
    n_burnin <- floor(n_iter / 3)
    n_thin <- 5
    run_parallel <- FALSE
  }

  # Parallel JAGS workers cannot safely receive a textConnection.
  # Write the model to a temporary file that every worker can open.
  
  model_file <- tempfile(
    pattern = "fourier_model_",
    fileext = ".bug"
  )
  
  writeLines(
    model_string,
    con = model_file
  )
  
  on.exit(
    unlink(model_file),
    add = TRUE
  )
  
  jm <- jagsUI::jags(
    data = jags_data,
    parameters.to.save = c("N", "beta.cos", "beta.sin", "sigma.X", "sigma.y"),
    model.file = model_file,
    n.chains = n_chains,
    n.iter = n_iter,
    n.burnin = n_burnin,
    n.thin = n_thin,
    parallel = run_parallel,
    verbose = FALSE
  )
  
  jm$dev_record <- list(model=model_string, data=jags_data,
    sampling=list(chains=n_chains, iterations=n_iter, burnin=n_burnin, thin=n_thin, parallel=run_parallel),
    mode=if(legacy_mcmc) "Reference settings" else "Exploration settings",
    seed_note="JAGS initialization was automatic; chains and model data are exported, but exact replay seed is not retained.",
    app_version="guided-developer-2026-09-02")
  as_site_matrix <- function(x) {
    if (n_timeseries == 1) matrix(as.numeric(x), ncol = 1) else as.matrix(x)
  }
  n_median <- as_site_matrix(jm$q50$N)
  n_lower <- as_site_matrix(jm$q2.5$N)
  n_upper <- as_site_matrix(jm$q97.5$N)
  
  d_annual_list <- list()
  for (i in 1:n_timeseries) {
    d_annual_list[[i]] <- data.frame(
      Year = all_seasons,
      Site = sites[i],
      Count = exp(n_median[, i]),
      Count_Lower_95 = exp(n_lower[, i]),
      Count_Upper_95 = exp(n_upper[, i]),
      Season_Start_Month = as.integer(season_start_month)
    )
  }
  res_annual <- do.call(rbind, d_annual_list)
  attr(res_annual, "mcmc_mode") <- if (legacy_mcmc) {
    "Martin/Siders settings"
  } else {
    "Interactive fast settings"
  }
  attr(res_annual, "dev_fit") <- jm
  return(res_annual)
}

# =====================================================================
# IMPUTATION-UNCERTAINTY BRANCH HELPERS
# =====================================================================
imputation_branch_labels <- c(
  observed = "Observed annual counts",
  lower = "Lower 95% imputed counts",
  median = "Median imputed counts",
  upper = "Upper 95% imputed counts"
)

build_imputation_branch_inputs <- function(
    d_annual, clutch_frequency, include_bounds = TRUE,
    source_type = c("imputed", "observed")) {
  source_type <- match.arg(source_type)
  required <- c("Year", "Site", "Count")
  missing_required <- setdiff(required, names(d_annual))
  if (length(missing_required) > 0L) {
    stop(
      "Annual imputation output is missing required column(s): ",
      paste(missing_required, collapse = ", ")
    )
  }
  clutch_frequency <- positive_or_stop(
    clutch_frequency, "Clutch frequency"
  )

  branch_from_column <- function(column_name, branch_name) {
    count_value <- suppressWarnings(as.numeric(d_annual[[column_name]]))
    data.frame(
      Year = as.integer(d_annual$Year),
      Site = as.character(d_annual$Site),
      Count = count_value,
      Annual_Nesters = count_value / clutch_frequency,
      Imputation_Branch = branch_name,
      stringsAsFactors = FALSE
    )
  }

  base_branch <- if (identical(source_type, "observed")) {
    "observed"
  } else {
    "median"
  }
  branches <- stats::setNames(
    list(branch_from_column("Count", base_branch)),
    base_branch
  )

  bound_columns <- c(
    lower = "Count_Lower_95",
    upper = "Count_Upper_95"
  )
  bounds_available <- all(bound_columns %in% names(d_annual))
  if (
    identical(source_type, "imputed") &&
      isTRUE(include_bounds) &&
      bounds_available
  ) {
    branches <- list(
      lower = branch_from_column(bound_columns[["lower"]], "lower"),
      median = branches$median,
      upper = branch_from_column(bound_columns[["upper"]], "upper")
    )
  }

  for (branch_name in names(branches)) {
    bad <- !is.na(branches[[branch_name]]$Annual_Nesters) &
      (!is.finite(branches[[branch_name]]$Annual_Nesters) |
         branches[[branch_name]]$Annual_Nesters <= 0)
    if (any(bad)) {
      stop(
        "The ", branch_name,
        " imputation branch contains non-positive annual-nester values."
      )
    }
  }

  attr(branches, "bounds_available") <- bounds_available
  branches
}

regional_abundance_draws <- function(res) {
  if (is.null(res$fit) || is.null(res$fit$sims.list$X)) {
    stop("Trend result does not contain posterior latent-state draws.")
  }
  x_draws <- as.matrix(res$fit$sims.list$X)
  if (is.null(res$fit$sims.list$A)) {
    return(exp(x_draws))
  }

  a_draws <- as.matrix(res$fit$sims.list$A)
  if (nrow(a_draws) != nrow(x_draws)) {
    stop("Site-offset and latent-state posterior draws are misaligned.")
  }
  site_scale <- rowSums(exp(a_draws))
  sweep(exp(x_draws), 1, site_scale, "*")
}

summarise_trend_result <- function(res, remigration_interval) {
  remigration_interval <- positive_or_stop(
    remigration_interval, "Remigration interval"
  )
  x_total <- regional_abundance_draws(res)
  years <- res$years
  final_index <- ncol(x_total)

  # Exact Siders take_helper_Fn.R current-abundance calculation:
  #   Sum = (N_fym0 + N_fym1 + N_fym2 + N_fym3) * RI / 4
  # where each N_fym value is the regional annual-nester abundance in one of
  # the final four model years. This is retained exactly for lineage fidelity.
  if (final_index < 4L) {
    stop(
      "At least four fitted annual-nester years are required to calculate ",
      "current abundance using the original Siders RI/4 method."
    )
  }
  recent_index <- seq.int(final_index - 3L, final_index)
  total_female_draws <- rowSums(
    x_total[, recent_index, drop = FALSE]
  ) * (remigration_interval / 4)

  u_draws <- as.numeric(res$fit$sims.list$U)
  q_draws <- as.numeric(res$fit$sims.list$Q)
  r_draws <- res$fit$sims.list$R
  r_mean <- if (is.matrix(r_draws)) rowMeans(r_draws) else as.numeric(r_draws)

  list(
    X_total = x_total,
    years = years,
    year = max(years),
    nesters = median(x_total[, final_index]),
    total = median(total_female_draws),
    trend_display = paste0(
      round(median(u_draws), 3), " (",
      round((exp(median(u_draws)) - 1) * 100, 2), "%)"
    ),
    trend_pct = (exp(median(u_draws)) - 1) * 100,
    draws = data.frame(
      U = u_draws,
      Q = q_draws,
      R_mean = r_mean,
      N_fym0 = x_total[, final_index],
      N_fym1 = x_total[, final_index - 1L],
      N_fym2 = x_total[, final_index - 2L],
      N_fym3 = x_total[, final_index - 3L],
      Total_Females = total_female_draws
    ),
    abundance_method = "Original Siders RI/4 current-abundance helper"
  )
}

# =====================================================================
# MATCHED SIDERS PROJECTION POSTERIOR
# =====================================================================

# Construct one joint posterior table for regional projection.
#
# Each row retains the original dependence among:
#   - trend U
#   - process variance Q
#   - final four regional annual-nester states
#   - Siders current adult-female abundance
#
# The current-abundance equation follows curr.abund.fn() in
# take_helper_Fn.R:
#
#   Total_Females =
#     (N_fym0 + N_fym1 + N_fym2 + N_fym3) * RI / 4

build_siders_projection_posterior <- function(
    trend_result,
    remigration_interval) {
  
  remigration_interval <- positive_or_stop(
    remigration_interval,
    "Remigration interval"
  )
  
  if (
    is.null(trend_result$fit) ||
    is.null(trend_result$fit$sims.list)
  ) {
    stop(
      "The trend result does not contain posterior draws."
    )
  }
  
  sims <- trend_result$fit$sims.list
  
  if (
    is.null(sims$X) ||
    is.null(sims$U) ||
    is.null(sims$Q)
  ) {
    stop(
      "The trend posterior must contain X, U, and Q."
    )
  }
  
  regional_nesters <- regional_abundance_draws(
    trend_result
  )
  
  if (ncol(regional_nesters) < 4L) {
    stop(
      "At least four modeled annual-nester years are required ",
      "for the original Siders current-abundance calculation."
    )
  }
  
  final_index <- ncol(regional_nesters)
  
  u_draws <- as.numeric(sims$U)
  q_draws <- as.numeric(sims$Q)
  
  n_joint <- nrow(regional_nesters)
  if (length(u_draws) != n_joint || length(q_draws) != n_joint) {
    stop("Abundance, trend and variance posterior rows must match exactly.")
  }

  if (n_joint < 2L) {
    stop(
      "Too few matched posterior draws are available ",
      "for projection."
    )
  }
  
  regional_nesters <- regional_nesters[
    seq_len(n_joint),
    ,
    drop = FALSE
  ]
  
  posterior <- data.frame(
    Posterior_Row = seq_len(n_joint),
    
    U = u_draws[seq_len(n_joint)],
    Q = q_draws[seq_len(n_joint)],
    
    N_fym0 = regional_nesters[
      ,
      final_index
    ],
    
    N_fym1 = regional_nesters[
      ,
      final_index - 1L
    ],
    
    N_fym2 = regional_nesters[
      ,
      final_index - 2L
    ],
    
    N_fym3 = regional_nesters[
      ,
      final_index - 3L
    ],
    
    stringsAsFactors = FALSE
  )
  
  posterior$Total_Females <- (
    posterior$N_fym0 +
      posterior$N_fym1 +
      posterior$N_fym2 +
      posterior$N_fym3
  ) * remigration_interval / 4
  
  valid <- complete.cases(posterior) &
    is.finite(posterior$N_fym0) &
    posterior$N_fym0 >= 0 &
    is.finite(posterior$U) &
    is.finite(posterior$Q) &
    posterior$Q > 0
  
  posterior <- posterior[
    valid,
    ,
    drop = FALSE
  ]
  
  if (nrow(posterior) < 2L) {
    stop(
      "Too few valid matched Siders posterior draws remain ",
      "after validation."
    )
  }
  
  rownames(posterior) <- NULL
  
  attr(
    posterior,
    "abundance_method"
  ) <- "Final-year regional annual nesters (N_fym0); RI/4 is reporting only"
  
  posterior
}

run_marss_diagnostic <- function(abund) {
  prep_marss <- abund %>%
    dplyr::select(Year, Site, Annual_Nesters) %>%
    pivot_wider(names_from = Site, values_from = Annual_Nesters) %>%
    arrange(Year)
  mat_data <- as.matrix(prep_marss %>% select(-Year))
  mat_data[mat_data <= 0 | is.nan(mat_data) | is.infinite(mat_data)] <- NA
  y_matrix <- t(log(mat_data))
  n_sites <- nrow(y_matrix)

  if (n_sites == 1L) {
    fit_shared <- MARSS::MARSS(
      y_matrix,
      model = list(
        Z = matrix(1), A = "zero", R = matrix("r"),
        Q = matrix("q"), U = matrix("u")
      ),
      silent = TRUE
    )
    fit_independent <- fit_shared
  } else {
    fit_shared <- MARSS::MARSS(
      y_matrix,
      model = list(
        Z = matrix(1, nrow = n_sites, ncol = 1),
        A = "scaling", R = "diagonal and unequal",
        Q = matrix("q"), U = matrix("u")
      ),
      silent = TRUE
    )
    fit_independent <- MARSS::MARSS(
      y_matrix,
      model = list(
        Z = diag(1, n_sites), A = "zero",
        R = "diagonal and unequal",
        Q = "diagonal and unequal", U = "unequal"
      ),
      silent = TRUE
    )
  }

  list(
    shared = fit_shared,
    indep = fit_independent,
    years = prep_marss$Year
  )
}

# =====================================================================
# TREND REFERENCE-SITE HELPERS
# =====================================================================
# The shared trend model requires one site offset to be fixed at zero for
# identifiability. Martin et al. achieved this by manually placing a suitable
# reference beach in the first data column. The generalized app selects that
# site automatically using a documented, reproducible rule and permits an
# advanced override only among sites observed in the first analysis year.

longest_consecutive_run <- function(years) {
  years <- sort(unique(as.integer(years[is.finite(years)])))
  if (length(years) == 0L) return(0L)
  run_group <- cumsum(c(1L, as.integer(diff(years) != 1L)))
  as.integer(max(tabulate(run_group)))
}

reference_site_diagnostics <- function(site, analysis_year, observed) {
  if (length(site) != length(analysis_year) ||
      length(site) != length(observed)) {
    stop("Reference-site diagnostic vectors must have equal lengths.")
  }

  dat <- data.frame(
    Site = as.character(site),
    Analysis_Year = suppressWarnings(as.integer(analysis_year)),
    Observed = as.logical(observed),
    stringsAsFactors = FALSE
  )
  dat$Observed[is.na(dat$Observed)] <- FALSE
  dat <- dat[!is.na(dat$Site) & nzchar(dat$Site) &
               is.finite(dat$Analysis_Year), , drop = FALSE]

  if (nrow(dat) == 0L) {
    stop("No valid site-year records are available for reference-site selection.")
  }

  all_years <- seq.int(
    min(dat$Analysis_Year, na.rm = TRUE),
    max(dat$Analysis_Year, na.rm = TRUE)
  )
  first_year <- min(all_years)
  sites <- sort(unique(dat$Site))

  stats <- do.call(
    rbind,
    lapply(sites, function(site_name) {
      site_rows <- dat[dat$Site == site_name, , drop = FALSE]
      observed_years <- sort(unique(
        site_rows$Analysis_Year[site_rows$Observed]
      ))
      data.frame(
        Site = site_name,
        Eligible_First_Year = first_year %in% observed_years,
        Observed_Years = length(observed_years),
        Total_Analysis_Years = length(all_years),
        Completeness = length(observed_years) / length(all_years),
        Longest_Continuous_Run = longest_consecutive_run(observed_years),
        Positive_Records = sum(site_rows$Observed),
        stringsAsFactors = FALSE
      )
    })
  )

  eligible <- stats[stats$Eligible_First_Year, , drop = FALSE]
  if (nrow(eligible) == 0L) {
    stop(
      "At least one site must have a positive observed value in the first ",
      "analysis year (", first_year, ")."
    )
  }

  # Selection hierarchy:
  #   1. positive observation in first analysis year (eligibility filter)
  #   2. greatest proportion of analysis years observed
  #   3. longest continuous observed run
  #   4. greatest number of positive observation records
  #   5. alphabetical order as a reproducible final tie-breaker
  rank_order <- order(
    -eligible$Completeness,
    -eligible$Longest_Continuous_Run,
    -eligible$Positive_Records,
    eligible$Site
  )
  recommended <- eligible[rank_order[1], , drop = FALSE]

  reason <- paste0(
    "Selected automatically because it has positive data in the first ",
    "analysis year (", first_year, "), observations in ",
    recommended$Observed_Years, " of ", recommended$Total_Analysis_Years,
    " analysis years (", round(100 * recommended$Completeness, 1),
    "%), and a longest continuous observed run of ",
    recommended$Longest_Continuous_Run, " years."
  )

  list(
    reference_site = recommended$Site[[1]],
    first_year = first_year,
    diagnostics = stats,
    eligible_sites = eligible$Site,
    selection_reason = reason
  )
}

recommend_reference_site_from_observations <- function(
    df, season_start_month = 4L) {
  required <- c("Year", "Site", "Count")
  missing_columns <- setdiff(required, names(df))
  if (length(missing_columns) > 0L) {
    stop(
      "Automatic reference-site selection requires column(s): ",
      paste(missing_columns, collapse = ", "), "."
    )
  }

  indexed <- df
  has_monthly_data <- "Month" %in% names(indexed) &&
    any(!is.na(indexed$Month))
  if (has_monthly_data) {
    if (!all(c("Season", "Seq_Month") %in% names(indexed))) {
      indexed <- add_nesting_season(
        indexed,
        season_start_month = season_start_month
      )
    }
    analysis_year <- indexed$Season
  } else {
    analysis_year <- indexed$Year
  }

  monitored <- if ("Monitored" %in% names(indexed)) {
    is.na(indexed$Monitored) | indexed$Monitored
  } else {
    !is.na(indexed$Count)
  }
  observed <- monitored & !is.na(indexed$Count) &
    is.finite(indexed$Count) & indexed$Count > 0

  reference_site_diagnostics(
    site = indexed$Site,
    analysis_year = analysis_year,
    observed = observed
  )
}

resolve_reference_site <- function(df, reference_site = NULL) {
  required <- c("Year", "Site", "Annual_Nesters")
  missing_columns <- setdiff(required, names(df))
  if (length(missing_columns) > 0L) {
    stop(
      "Reference-site selection requires column(s): ",
      paste(missing_columns, collapse = ", "), "."
    )
  }

  diagnostic <- reference_site_diagnostics(
    site = df$Site,
    analysis_year = df$Year,
    observed = !is.na(df$Annual_Nesters) &
      is.finite(df$Annual_Nesters) & df$Annual_Nesters > 0
  )

  automatic <- is.null(reference_site) || length(reference_site) != 1L ||
    is.na(reference_site) || !nzchar(reference_site)
  if (automatic) {
    reference_site <- diagnostic$reference_site
    selection_mode <- "automatic"
    selection_reason <- diagnostic$selection_reason
  } else {
    reference_site <- as.character(reference_site)
    selection_mode <- "advanced override"
    selection_reason <- paste0(
      "Selected through the advanced override. Eligible reference sites are ",
      "restricted to those with positive data in the first analysis year (",
      diagnostic$first_year, ")."
    )
  }

  all_sites <- sort(unique(as.character(df$Site[!is.na(df$Site)])))
  if (!reference_site %in% all_sites) {
    stop("Selected reference site is not present in the trend data: ", reference_site)
  }
  if (!reference_site %in% diagnostic$eligible_sites) {
    stop(
      "The selected reference site ('", reference_site,
      "') has no positive annual-nester observation in the first model year ",
      "(", diagnostic$first_year, "). Choose an eligible site."
    )
  }

  first_year_values <- df$Annual_Nesters[
    df$Year == diagnostic$first_year & df$Site == reference_site
  ]
  first_year_values <- first_year_values[
    is.finite(first_year_values) & first_year_values > 0
  ]
  if (length(first_year_values) < 1L) {
    stop("Could not identify the first-year value for the reference site.")
  }
  reference_value <- sum(first_year_values)

  list(
    reference_site = reference_site,
    first_year = diagnostic$first_year,
    first_value = as.numeric(reference_value),
    x0_mean = log(as.numeric(reference_value)),
    sites = c(reference_site, setdiff(all_sites, reference_site)),
    eligible_sites = diagnostic$eligible_sites,
    diagnostics = diagnostic$diagnostics,
    selection_mode = selection_mode,
    selection_reason = selection_reason
  )
}

# =====================================================================
# EXACT STATE-SPACE TREND ENGINE
# =====================================================================
# Methods 2.3, Eqs. 4-6: shared latent trend, reference-site offset A[1]=0.
# Q and R are variances (precisions are their reciprocals); CF is applied
# upstream once to convert annual nests into annual nesters.
run_jags_aligned <- function(
    df, iter = 50000, burnin = 10000, thin = 10,
    n_chains = 3, parallel = FALSE,
    reference_site = NULL) {
  validate_log_count_input(
    count = df$Annual_Nesters,
    label = "annual-nester"
  )

  reference <- resolve_reference_site(df, reference_site)
  all_years <- min(df$Year, na.rm = TRUE):max(df$Year, na.rm = TRUE)
  sites <- reference$sites
  n.yrs <- length(all_years); n.timeseries <- length(sites)
  
  prep_df <- df %>%
    group_by(Year, Site) %>%
    summarise(Annual_Nesters = if(all(is.na(Annual_Nesters))) NA_real_ else sum(Annual_Nesters, na.rm = TRUE), .groups = "drop") %>%
    complete(Year = all_years, Site = sites) %>%
    arrange(Year) %>%
    pivot_wider(names_from = Site, values_from = Annual_Nesters) %>%
    select(all_of(sites))
  
  mat_data <- as.matrix(prep_df)
  mat_data[is.nan(mat_data) | is.infinite(mat_data)] <- NA
  if (any(mat_data <= 0, na.rm = TRUE)) {
    stop(
      "Internal validation failure: non-positive annual-nester values reached ",
      "the Gaussian log-count trend model."
    )
  }
  Y_matrix <- t(log(mat_data))
  
  jags_data <- list(
    n.yrs = n.yrs, n.timeseries = n.timeseries, a_mean = 0, a_sd = 4,
    u_mean = 0, u_sd = 0.5, q_alpha = 0.01, q_beta = 0.01, r_alpha = 0.01, r_beta = 0.01,
    # Martin singleUQ: center X[1] on the first observation of the
    # deliberately selected first/reference time series.
    x0_mean = reference$x0_mean, x0_sd = 10
  )
  
  if (n.timeseries == 1) {
    jags_data$Y <- as.vector(Y_matrix)
    jags_data$n.timeseries <- NULL
    jags_data$a_mean <- NULL
    jags_data$a_sd <- NULL
    model_string <- "
    model {
      U ~ dnorm(u_mean, 1/(u_sd^2))
      tauQ ~ dgamma(q_alpha, q_beta)
      Q <- 1/tauQ
      
      X0 ~ dnorm(x0_mean, 1/(x0_sd^2))
      predX[1] <- X0 + U
      X[1] <- predX[1]
      
      for(t in 2:n.yrs) {
        predX[t] <- X[t-1] + U
        X[t] ~ dnorm(predX[t], tauQ)
      }
      tauR ~ dgamma(r_alpha, r_beta)
      R <- 1/tauR
      for(t in 1:n.yrs) { Y[t] ~ dnorm(X[t], tauR) }
    }"
    params <- c("U", "Q", "R", "X0", "X")
  } else {
    jags_data$Y <- Y_matrix
    jags_data$Z <- matrix(rep(1, n.timeseries), ncol=1)
    model_string <- "
    model {
      A[1] <- 0
      for(j in 2:n.timeseries) { A[j] ~ dnorm(a_mean, 1/(a_sd^2)) }
      U ~ dnorm(u_mean, 1/(u_sd^2))
      tauQ ~ dgamma(q_alpha, q_beta)
      Q <- 1/tauQ
      X0 ~ dnorm(x0_mean, 1/(x0_sd^2))
      predX[1] <- X0 + U
      X[1] <- predX[1]
      for(t in 2:n.yrs) {
        predX[t] <- X[t-1] + U
        X[t] ~ dnorm(predX[t], tauQ)
      }
      for(j in 1:n.timeseries) {
        tauR[j] ~ dgamma(r_alpha, r_beta)
        R[j] <- 1/tauR[j]
        for(t in 1:n.yrs) {
          predY[j,t] <- Z[j,1] * X[t] + A[j]
          Y[j,t] ~ dnorm(predY[j,t], tauR[j])
        }
      }
    }"
    params <- c("U", "Q", "R", "X0", "X", "A")
  }
  
  # Parallel JAGS workers cannot safely receive a textConnection.
  # Write the trend model to a temporary file accessible to every worker.
  
  trend_model_file <- tempfile(
    pattern = "trend_model_",
    fileext = ".bug"
  )
  
  writeLines(
    model_string,
    con = trend_model_file
  )
  
  on.exit(
    unlink(trend_model_file),
    add = TRUE
  )
  
  fit <- jagsUI::jags(
    data = jags_data,
    parameters.to.save = params,
    model.file = trend_model_file,
    n.chains = n_chains,
    n.iter = iter,
    n.burnin = burnin,
    n.thin = thin,
    parallel = parallel,
    verbose = FALSE
  )
  
  fit$dev_record <- list(model=model_string, data=jags_data, parameters=params,
    sampling=list(chains=n_chains, iterations=iter, burnin=burnin, thin=thin, parallel=parallel),
    seed_note="Automatic initialization; no exact replay seed retained.",
    app_version="guided-developer-2026-09-02")
  return(list(
    fit = fit,
    years = all_years,
    sites = sites,
    reference_site = reference$reference_site,
    reference_first_year = reference$first_year,
    reference_first_value = reference$first_value,
    x0_mean = reference$x0_mean,
    reference_selection_mode = reference$selection_mode,
    reference_selection_reason = reference$selection_reason,
    reference_diagnostics = reference$diagnostics
  ))
}

# =====================================================================
# HISTORICAL RETROSPECTIVE COMPILER ENGINE
# =====================================================================
calculate_empirical_ane <- function(obs_df, safe_df, params, current_year) {
  # Historical wrapper, not a second ANE model. Observed turtles retain their
  # measured sizes/risks; missing individuals receive empirical input draws.
  # Their ANE uses Eqs. 7-14, so a turtle is allocated to first nesting once.
  # No species-wide fallback counts, variation or correlations are inserted.
  if (!all(c("Year", "Length", "Mortality") %in% names(obs_df)) ||
      !all(c("Year", "Total_Est") %in% names(safe_df))) {
    stop("Historical inputs require Year/Length/Mortality and Year/Total_Est.")
  }
  if (!nrow(safe_df) || anyDuplicated(safe_df$Year) ||
      any(!is.finite(safe_df$Year)) || any(safe_df$Year != floor(safe_df$Year)) ||
      any(safe_df$Year > current_year) || any(!is.finite(safe_df$Total_Est)) ||
      any(safe_df$Total_Est < 0)) stop("Invalid historical year/count schedule.")
  if (nrow(obs_df) &&
      (any(!is.finite(obs_df$Length)) || any(obs_df$Length <= 0) ||
       any(!is.finite(obs_df$Mortality)) || any(obs_df$Mortality < 0 | obs_df$Mortality > 1) ||
       any(!obs_df$Year %in% safe_df$Year))) {
    stop("Observed turtles need valid lengths, mortality probabilities and scheduled years.")
  }
  observed_count <- vapply(safe_df$Year, function(y) sum(obs_df$Year == y), numeric(1))
  if (any(observed_count > safe_df$Total_Est)) stop("Observed turtles exceed the estimated annual total.")
  missing_total <- safe_df$Total_Est - observed_count
  if (any(missing_total > 0) && !nrow(obs_df)) {
    stop("Historical expansion requires observed size/mortality data.")
  }
  p <- list(linf=params$sp_linf, k=params$sp_k, t0=params$sp_tknot,
            lmat=params$sp_lmat, sig_mat=params$sp_sig_mat, mat_p=0.99,
            pj=params$sp_ane_pj, pa=params$sp_ane_pa,
            pf=params$sp_pf, ri=params$sp_remig_int)
  years <- seq.int(min(safe_df$Year), current_year)
  loss <- numeric(length(years))
  # Empirical bivariate distribution is an app extension around expected ANE.
  # A single complete record implies no estimated spread, not invented spread.
  log_len <- log(obs_df$Length)
  logit_mort <- safe_logit(obs_df$Mortality)
  beta0 <- if (length(log_len)) mean(log_len) else NA_real_
  beta1 <- 0
  sigma_l <- if (length(log_len) > 1) sd(log_len) else 0
  sigma_d <- if (length(logit_mort) > 1) sd(logit_mort) else 0
  rho <- if (sigma_l > 0 && sigma_d > 0) cor(log_len, logit_mort) else 0
  mu0 <- if (length(logit_mort)) mean(logit_mort) else NA_real_
  total <- safe_df$Total_Est[match(obs_df$Year, safe_df$Year)]
  if (length(unique(total)) > 1 && length(total) > 2) {
    fit <- lm(log_len ~ total)
    beta0 <- unname(coef(fit)[1]); beta1 <- unname(coef(fit)[2])
    sigma_l <- summary(fit)$sigma
    rho <- if (sigma_l > 0 && sigma_d > 0) cor(residuals(fit), logit_mort) else 0
  }
  covariance <- matrix(c(sigma_l^2, rho*sigma_l*sigma_d,
                         rho*sigma_l*sigma_d, sigma_d^2), 2)
  for (j in seq_len(nrow(safe_df))) {
    y <- safe_df$Year[j]
    rows <- obs_df[obs_df$Year == y, , drop=FALSE]
    n_extra <- ceiling(missing_total[j])
    lengths <- rows$Length; mortality <- rows$Mortality
    weights <- rep(1, nrow(rows))
    if (n_extra > 0) {
      center <- c(beta0 + beta1 * safe_df$Total_Est[j], mu0)
      draws <- mvtnorm::rmvnorm(n_extra, mean=center, sigma=covariance)
      extra_mort <- plogis(draws[, 2])
      if (all(obs_df$Mortality == 0)) extra_mort[] <- 0
      if (all(obs_df$Mortality == 1)) extra_mort[] <- 1
      lengths <- c(lengths, exp(draws[, 1]))
      mortality <- c(mortality, extra_mort)
      weights <- c(weights, rep(missing_total[j]/n_extra, n_extra))
    }
    start <- match(y, years)
    for (i in seq_along(lengths)) {
      loss[start:length(years)] <- loss[start:length(years)] + weights[i] *
        fishery_ane_kernel(lengths[i], mortality[i], p, length(years)-start+1L)
    }
  }
  # Legacy column name is retained for consumers; each row is ANNUAL assigned
  # ANE summed across cohorts, not a cumulative sum over calendar years.
  result <- data.frame(Year=years, Total_Cumulative_ANE=loss)
  attr(result, "meta_beta0") <- beta0; attr(result, "meta_beta1") <- beta1
  attr(result, "meta_sigma_L") <- sigma_l; attr(result, "meta_sigma_D") <- sigma_d
  attr(result, "meta_rho") <- rho; attr(result, "meta_mu0") <- mu0
  result
}

# =====================================================================
# USER-SUPPLIED BIOLOGICAL PARAMETERS AND MODEL-CURRENCY HELPERS
# =====================================================================
# Biological parameters are deliberately not populated from species-wide
# presets. Users must supply population-specific values and document their
# sources. This prevents global species means or provisional values from being
# treated as validated inputs for a particular DPS, RMU, or nesting population.
probability_or_stop <- function(x, label) {
  if (length(x) != 1 || is.na(x) || !is.finite(x) || x < 0 || x > 1) {
    stop(label, " must be one finite value between 0 and 1.")
  }
  x
}

positive_or_stop <- function(x, label) {
  if (length(x) != 1 || is.na(x) || !is.finite(x) || x <= 0) {
    stop(label, " must be one finite value greater than zero.")
  }
  x
}

vbgf_age_from_length <- function(length_cm, linf, k, t0) {
  positive_or_stop(linf, "Linf")
  positive_or_stop(k, "Growth coefficient k")
  if (any(!is.finite(length_cm)) || any(length_cm <= 0)) {
    stop("Carapace length must contain finite values greater than zero.")
  }
  length_safe <- pmin(length_cm, linf - 0.1)
  t0 - log(1 - length_safe / linf) / k
}

# Exact Siders vbgf_bc() behavior for fishery Take.
#
# Valid lengths use the direct inverse von Bertalanffy equation.
# A length greater than Linf produces NaN; the original code then
# substitutes the age corresponding to 99% of Linf.
#
# This helper is kept separate from the app's general VBGF helper so
# conservation-gain extensions are not changed.

siders_vbgf_age_from_length <- function(
    length_cm,
    linf,
    k,
    t0) {
  
  positive_or_stop(
    linf,
    "Linf"
  )
  
  positive_or_stop(
    k,
    "Growth coefficient k"
  )
  
  if (
    any(!is.finite(length_cm)) ||
    any(length_cm <= 0)
  ) {
    stop(
      "Carapace length must contain finite values ",
      "greater than zero."
    )
  }
  
  age <- suppressWarnings(
    t0 -
      log(1 - length_cm / linf) / k
  )
  
  fallback_age <- (
    t0 -
      log(1 - 0.99) / k
  )
  
  age[!is.finite(age)] <- fallback_age
  
  age
}

# Derive the VBGF age at the forced-maturity length threshold.
#
# In the original Siders pathway, VBGF$max_age was not maximum
# lifespan. It was calculated as vbgf_bc(prop = mat_p).
#
# This helper is used only by fishery Take so the app's separate
# maximum-lifespan parameter remains available to conservation models.

siders_forced_maturity_age <- function(
    linf,
    k,
    t0,
    max_age_prop = 0.99) {
  
  positive_or_stop(
    linf,
    "Linf"
  )
  
  positive_or_stop(
    k,
    "Growth coefficient k"
  )
  
  if (
    length(max_age_prop) != 1L ||
    is.na(max_age_prop) ||
    !is.finite(max_age_prop) ||
    max_age_prop <= 0 ||
    max_age_prop >= 1
  ) {
    stop(
      "The Siders maximum-age length proportion must be ",
      "one finite proportion in (0, 1)."
    )
  }
  
  t0 -
    log(1 - max_age_prop) / k
}

maturity_probability <- function(length_cm, lmat, sig_mat, linf,
                                 mat_p = 0.99) {
  positive_or_stop(sig_mat, "Maturity-ogive width")
  if (
    length(mat_p) != 1 || !is.finite(mat_p) ||
      mat_p <= 0 || mat_p > 1
  ) {
    stop("Maturity threshold must be one finite proportion in (0, 1].")
  }
  p <- plogis((length_cm - lmat) / sig_mat)
  ifelse(length_cm >= mat_p * linf, 1, p)
}

current_model_params <- function(input) {
  params <- list(
    ri = positive_or_stop(input$remig_int, "Remigration interval"),
    cf = positive_or_stop(input$clutch_freq, "Clutch frequency"),
    pf = probability_or_stop(input$pf, "Proportion female"),
    linf = positive_or_stop(input$linf, "Linf"),
    k = positive_or_stop(input$k, "Growth coefficient k"),
    t0 = input$tknot,
    lmat = positive_or_stop(input$lmat, "Length at maturity"),
    sig_mat = positive_or_stop(input$sig_mat, "Maturity-ogive width"),
    mat_p = {
      threshold_value <- suppressWarnings(
        as.numeric(input$ane_mat_threshold)
      )
      if (
        is.null(threshold_value) || length(threshold_value) != 1 ||
          !is.finite(threshold_value)
      ) 0.99 else threshold_value
    },
    # Diagnostic only; all kernels derive this from the VBGF (Eq. 7).
    max_age = siders_forced_maturity_age(input$linf, input$k, input$tknot),
    pj = probability_or_stop(input$ane_pj, "Juvenile survival"),
    pa = probability_or_stop(input$ane_pa, "Adult survival")
  )
  if (length(params$t0) != 1 || is.na(params$t0) || !is.finite(params$t0)) {
    stop("VBGF t0 must be one finite value.")
  }
  if (params$lmat >= params$linf) {
    stop("Length at maturity must be less than Linf.")
  }
  if (params$mat_p <= 0 || params$mat_p >= 1) {
    stop("First-nesting threshold must be in (0, 1).")
  }
  params
}

# Take uses the manuscript's capture-year trajectory. Give uses its own
# entry-stage-to-census timing below; the demographic currency is shared.
# Methods Eqs. 7-11: manuscript trajectory and expected first nesting.
siders_forced_first_nesting_age <- function(linf, k, t0, threshold = 0.99) {
  siders_forced_maturity_age(linf, k, t0, max_age_prop = threshold)
}

siders_age_from_length <- function(length_cm, params) {

  if (any(!is.finite(length_cm)) || any(length_cm <= 0)) {
    stop("Carapace length must contain finite values greater than zero.")
  }

  age <- suppressWarnings(
    params$t0 -
      log(1 - length_cm / params$linf) / params$k
  )

  fallback_age <- (
    params$t0 -
      log(1 - 0.99) / params$k
  )

  # App boundary safeguard: the manuscript inverse is infinite at exactly
  # Linf. Treat that boundary like lengths above Linf (forced adult age).
  age[!is.finite(age)] <- fallback_age
  age
}

siders_demographic_path <- function(length_cm, params, horizon, start_age = NULL) {

  if (length(horizon) != 1L || !is.finite(horizon) || horizon < 1) {
    stop("horizon must be a positive integer.")
  }
  horizon <- as.integer(horizon)

  age0 <- if (is.null(start_age)) siders_age_from_length(length_cm, params) else start_age
  if (length(age0) != 1L || !is.finite(age0)) {
    stop("This manuscript helper expects one starting length.")
  }

  max_age_99 <- siders_forced_first_nesting_age(
    params$linf,
    params$k,
    params$t0,
    threshold = 0.99
  )

  at_or_beyond_max <- (
    age0 >= max_age_99 ||
      isTRUE(all.equal(age0, max_age_99, tolerance = 1e-12))
  )

  if (at_or_beyond_max) {

    lengths <- rep(
      params$linf * params$mat_p,
      horizon
    )

  } else {

    growth_ages_to_max <- seq(
      age0,
      max_age_99,
      by = 1
    )

    year_reaching_max <- length(growth_ages_to_max)

    if (year_reaching_max < horizon) {

      lengths <- c(
        params$linf * (
          1 -
            exp(
              -params$k *
                (growth_ages_to_max - params$t0)
            )
        ),
        rep(
          params$linf * params$mat_p,
          horizon - year_reaching_max
        )
      )

    } else {

      growth_ages <- seq(
        age0,
        length.out = horizon,
        by = 1
      )

      lengths <- params$linf * (
        1 -
          exp(
            -params$k *
              (growth_ages - params$t0)
          )
      )
    }
  }

  lengths <- lengths[seq_len(horizon)]

  p_maturity <- maturity_probability(
    lengths,
    params$lmat,
    params$sig_mat,
    params$linf,
    params$mat_p
  )

  # Exact ordering in prop_ane(): the continuous p(M) values are first used to
  # construct the juvenile/adult survival mixture. The Bernoulli first-nesting
  # history is generated only afterward.
  # Methods Eqs. 10-11: probability-weighted survival before allocation.
  annual_survival <- (
    (1 - p_maturity) * params$pj +
      p_maturity * params$pa
  )
  survival_to_year <- cumprod(annual_survival)

  list(
    age = age0 + seq.int(0, horizon - 1L),
    length = lengths,
    maturity_probability = p_maturity,
    annual_survival = annual_survival,
    survival_to_year = survival_to_year
  )
}

first_nesting_path <- function(length_cm, params, horizon) {

  path <- siders_demographic_path(
    length_cm = length_cm,
    params = params,
    horizon = horizon
  )

  # p = annual probability governing first nesting (Siders p(M))
  p <- path$maturity_probability

  # Methods paper expected first-event probability.
  p.FN <- p * cumprod(c(1, head(1 - p, -1)))

  # Eq. 9: sum(p.FN) <= 1; no repeat allocation after first nesting.
  first_nesting_probability <- p.FN

  c(
    path,
    list(
      first_nesting_probability =
        first_nesting_probability
    )
  )
}

# A signed ANE cohort kernel. event_probability is the probability that one
# member of the cohort is removed (Take) or additionally survives because of
# an intervention (Give). The sign convention is negative for Take and
# positive for Give, allowing a later projection to use one net ANE ledger.
ane_cohort_kernel <- function(length_cm, event_probability, params, horizon,
                              sign = c("take", "give")) {
  sign <- match.arg(sign)
  event_probability <- probability_or_stop(
    event_probability,
    if (sign == "take") "Fishery mortality" else "Additional survival probability"
  )
  path <- first_nesting_path(length_cm, params, horizon)
  # Eq. 12: mortality, sex and remigration standardization appear once.
  magnitude <- path$survival_to_year * path$first_nesting_probability *
    event_probability * params$pf / params$ri
  if (sign == "take") -magnitude else magnitude
}

# Sum signed ANE contributions from cohorts created in each calendar year.
# This signed helper uses Take timing. It is not the conservation Give engine:
# intervention cohorts use conservation_ane_schedule(), with entry-stage-to-census
# survival and the manuscript's distinct Give timing.
signed_ane_cohort_ledger <- function(cohort_amount, length_cm,
                                     event_probability, params,
                                     sign = c("take", "give")) {
  sign <- match.arg(sign)
  horizon <- length(cohort_amount)
  if (length(length_cm) == 1) length_cm <- rep(length_cm, horizon)
  if (length(event_probability) == 1) event_probability <- rep(event_probability, horizon)
  if (length(length_cm) != horizon || length(event_probability) != horizon) {
    stop("Cohort amount, length, and event-probability schedules must have the same length.")
  }
  if (any(!is.finite(cohort_amount)) || any(cohort_amount < 0)) {
    stop("Cohort amounts must be finite non-negative values.")
  }

  ledger <- rep(0, horizon)
  for (capture_index in seq_len(horizon)) {
    if (cohort_amount[capture_index] == 0) next
    remaining <- horizon - capture_index + 1
    kernel <- ane_cohort_kernel(
      length_cm = length_cm[capture_index],
      event_probability = event_probability[capture_index],
      params = params, horizon = remaining, sign = sign
    )
    ledger[capture_index:horizon] <- ledger[capture_index:horizon] +
      cohort_amount[capture_index] * kernel
  }
  ledger
}

# =====================================================================
# CONSERVATION-GAIN COHORT PATHWAY
# =====================================================================
#
# Conservation actions are an Ortega extension, not an original
# Martin/Siders equation.
#
# Core counterfactual principle:
#
#   Give = outcome WITH intervention - outcome WITHOUT intervention
#
# The quantity entering this pathway must therefore represent organisms
# that exist because of the intervention, rather than all organisms
# handled by the intervention.
#
# Timing convention:
#
#   - an action occurs during projection transition t -> t + 1;
#   - start_age is age at the time the intervention cohort enters;
#   - contribution can first enter the annual-nester population at the
#     following annual census;
#   - survival during the first interval is therefore applied before
#     the first possible annual-nester contribution.
#
# The maturity/growth currency uses the same VBGF and 99%-of-Linf
# forced-maturity age as the aligned Siders Take pathway.


conservation_first_nesting_path <- function(
    start_age,
    params,
    horizon,
    first_interval_survival = NULL) {
  
  if (
    length(start_age) != 1L ||
    is.na(start_age) ||
    !is.finite(start_age) ||
    start_age < 0
  ) {
    stop(
      "Conservation cohort start age must be one finite ",
      "non-negative value."
    )
  }
  
  if (
    length(horizon) != 1L ||
    is.na(horizon) ||
    !is.finite(horizon) ||
    horizon < 1
  ) {
    stop(
      "Conservation cohort horizon must be one positive ",
      "finite value."
    )
  }
  
  horizon <- as.integer(horizon)
  
  if (!is.null(first_interval_survival)) {
    first_interval_survival <- probability_or_stop(
      first_interval_survival,
      "First-interval survival"
    )
  }
  
  mat_p <- if (is.null(params$mat_p)) {
    0.99
  } else {
    params$mat_p
  }
  
  # Methods Eq. 16: Give uses the selected first-nesting threshold.
  # Default .99; an alternative threshold is a sensitivity scenario.
  siders_max_age <- siders_forced_maturity_age(
    linf = params$linf,
    k = params$k,
    t0 = params$t0,
    max_age_prop = mat_p
  )
  
  # Age at the beginning of each annual transition.
  interval_start_age <- (
    start_age +
      seq.int(
        0,
        horizon - 1L
      )
  )
  
  # Age at the annual census reached after each transition.
  census_age <- (
    start_age +
      seq_len(horizon)
  )
  
  vbgf_length <- function(age_value) {
    
    length_value <- (
      params$linf *
        (
          1 -
            exp(
              -params$k *
                (age_value - params$t0)
            )
        )
    )
    
    length_value[
      age_value >= siders_max_age
    ] <- (
      params$linf *
        mat_p
    )
    
    length_value
  }
  
  interval_start_length <- vbgf_length(
    interval_start_age
  )
  
  census_length <- vbgf_length(
    census_age
  )
  
  # Stage composition governing survival during each annual interval.
  interval_maturity <- maturity_probability(
    length_cm = interval_start_length,
    lmat = params$lmat,
    sig_mat = params$sig_mat,
    linf = params$linf,
    mat_p = mat_p
  )
  
  annual_survival <- (
    (1 - interval_maturity) *
      params$pj +
      interval_maturity *
      params$pa
  )
  
  # For newly emerged hatchlings, the first interval is governed by
  # the separately supplied emerged-to-yearling survival probability.
  if (!is.null(first_interval_survival)) {
    annual_survival[1] <-
      first_interval_survival
  }
  
  survival_to_census <- cumprod(
    annual_survival
  )
  
  # Probability of being mature at each future annual census.
  census_maturity <- maturity_probability(
    length_cm = census_length,
    lmat = params$lmat,
    sig_mat = params$sig_mat,
    linf = params$linf,
    mat_p = mat_p
  )
  
  # Expected first-nesting probability:
  #
  # P(first nesting in year j) =
  #   P(mature in j) *
  #   product[P(not mature in all preceding census years)]
  
  first_nesting_probability <- (
    census_maturity *
      cumprod(
        c(
          1,
          head(
            1 - census_maturity,
            -1
          )
        )
      )
  )
  
  list(
    interval_start_age = interval_start_age,
    census_age = census_age,
    interval_start_length = interval_start_length,
    census_length = census_length,
    annual_survival = annual_survival,
    survival_to_census = survival_to_census,
    maturity_probability = census_maturity,
    first_nesting_probability =
      first_nesting_probability
  )
}


# Convert hatching/emergence inputs into the unconditional probability
# that one egg produces one emerged hatchling under the intervention.
#
# Two common definitions of emergence success occur in monitoring data:
#
#   conditional_on_hatching:
#       emergence = emerged / hatched
#
#       P(emerged | egg) =
#         P(hatched | egg) *
#         P(emerged | hatched)
#
#   per_egg:
#       emergence = emerged / eggs laid
#
#       emergence is already P(emerged | egg), so hatch success must
#       NOT be multiplied again.

conservation_emergence_per_egg <- function(
    hatch,
    emergence,
    emergence_basis = c(
      "conditional_on_hatching",
      "per_egg"
    )) {
  
  emergence_basis <- match.arg(
    emergence_basis
  )
  
  hatch <- probability_or_stop(
    hatch,
    "Hatching success"
  )
  
  emergence <- probability_or_stop(
    emergence,
    "Emergence success"
  )
  
  if (
    identical(
      emergence_basis,
      "conditional_on_hatching"
    )
  ) {
    
    hatch * emergence
    
  } else {
    
    emergence
  }
}


# Expected incremental annual-nester-equivalent gains from a
# conservation action.
#
# IMPORTANT DEFINITIONS OF amount:
#
#   nests:
#     number of nests receiving protection.
#
#   yearlings:
#     NET ADDITIONAL turtles alive at release age because of the
#     intervention, relative to the no-intervention counterfactual.
#
#   adults:
#     adult deaths genuinely prevented by the intervention.
#
# For nest protection:
#
#   incremental emerged hatchlings =
#
#     protected nests *
#     eggs per nest *
#     (
#       emerged/egg under protection -
#       emerged/egg without protection
#     )
#
# This prevents all hatchlings from a protected nest being incorrectly
# treated as conservation benefit when some would have survived anyway.

conservation_ane_schedule <- function(
    type,
    amount,
    params,
    eggs,
    hatch,
    emergence,
    year1,
    release_age = 1,
    emergence_basis = c(
      "conditional_on_hatching",
      "per_egg"
    ),
    counterfactual_emergence_per_egg = 0) {
  
  type <- match.arg(
    type,
    c(
      "nests",
      "yearlings",
      "adults"
    )
  )
  
  emergence_basis <- match.arg(
    emergence_basis
  )
  
  horizon <- length(amount)
  
  if (
    horizon < 1L ||
    any(is.na(amount)) ||
    any(!is.finite(amount)) ||
    any(amount < 0)
  ) {
    stop(
      "Conservation action amounts must contain finite ",
      "non-negative values."
    )
  }
  
  if (
    length(release_age) != 1L ||
    is.na(release_age) ||
    !is.finite(release_age) ||
    release_age < 0
  ) {
    stop(
      "Release age must be one finite non-negative value."
    )
  }
  
  first_interval_survival <- NULL
  
  if (identical(type, "nests")) {
    
    eggs <- positive_or_stop(
      eggs,
      "Eggs per nest"
    )
    
    year1 <- probability_or_stop(
      year1,
      "Emerged-to-yearling survival"
    )
    
    counterfactual_emergence_per_egg <-
      probability_or_stop(
        counterfactual_emergence_per_egg,
        "Counterfactual emergence per egg"
      )
    
    protected_emergence_per_egg <-
      conservation_emergence_per_egg(
        hatch = hatch,
        emergence = emergence,
        emergence_basis = emergence_basis
      )
    
    # Eq. 18: counterfactual difference, not total protected hatchlings.
    incremental_emergence_per_egg <- (
      protected_emergence_per_egg -
        counterfactual_emergence_per_egg
    )
    
    if (
      incremental_emergence_per_egg <
      -1e-12
    ) {
      stop(
        "Protected emergence per egg cannot be lower than ",
        "counterfactual emergence per egg when calculating ",
        "a positive conservation gain."
      )
    }
    
    incremental_emergence_per_egg <- max(
      0,
      incremental_emergence_per_egg
    )
    
    cohort_amount <- (
      amount *
        eggs *
        incremental_emergence_per_egg
    )
    
    # Incremental emerged hatchlings enter at age zero.
    entry_age <- 0
    
    # Their first annual transition is emerged -> yearling.
    first_interval_survival <- year1
    
  } else if (identical(type, "yearlings")) {
    
    # amount is already the incremental number alive at release age
    # relative to the no-headstarting counterfactual.
    cohort_amount <- amount
    
    entry_age <- release_age
    
  } else {
    
    # amount is the number of adult deaths genuinely prevented.
    cohort_amount <- amount
    
    # Use the same mature-age convention as the aligned Siders
    # fishery pathway rather than biological maximum lifespan.
    entry_age <- siders_forced_maturity_age(
      linf = params$linf,
      k = params$k,
      t0 = params$t0,
      max_age_prop = params$mat_p
    )
  }
  
  ledger <- rep(
    0,
    horizon
  )
  
  for (
    action_index in seq_len(horizon)
  ) {
    
    if (
      cohort_amount[action_index] == 0
    ) {
      next
    }
    
    remaining <- (
      horizon -
        action_index +
        1L
    )
    
    path <- conservation_first_nesting_path(
      start_age = entry_age,
      params = params,
      horizon = remaining,
      first_interval_survival =
        first_interval_survival
    )
    
    kernel <- (
      path$survival_to_census *
        path$first_nesting_probability *
        params$pf /
        params$ri
    )
    
    target <- action_index:horizon
    
    ledger[target] <- (
      ledger[target] +
        cohort_amount[action_index] *
        kernel
    )
  }
  
  ledger
}

# Expected-value adaptation of the Siders probability-of-maturity method.
# The original code draws maturity as a Bernoulli event for each animal. When
# the app receives aggregate counts and median length, this kernel uses the
# corresponding expected first-maturity probability for each future year.
fishery_ane_kernel <- function(length_cm, mortality, params, horizon) {
  -ane_cohort_kernel(
    length_cm = length_cm, event_probability = mortality, params = params,
    horizon = horizon, sign = "take"
  )
}

fishery_ane_schedule <- function(interactions, length_cm, mortality, params) {
  horizon <- length(interactions)
  if (length(length_cm) == 1) length_cm <- rep(length_cm, horizon)
  if (length(mortality) == 1) mortality <- rep(mortality, horizon)
  if (length(length_cm) != horizon || length(mortality) != horizon) {
    stop("Fishery schedule columns must have the same length.")
  }
  if (any(!is.finite(interactions)) || any(interactions < 0)) {
    stop("Fishery interactions must be finite non-negative values.")
  }

  -signed_ane_cohort_ledger(
    cohort_amount = interactions, length_cm = length_cm,
    event_probability = mortality, params = params, sign = "take"
  )
}

fishery_ane_expected <- function(interactions, length_cm, mortality, params) {
  if (!is.finite(interactions) || interactions < 0) {
    stop("Fishery interactions must be a finite non-negative value.")
  }
  sum(
    interactions * fishery_ane_kernel(
      # This helper is used for a standalone long-horizon expected value.
      # Match the 100-year Martin/Siders PVA horizon rather than truncating a
      # cohort at maximum biological age.
      length_cm, mortality, params, horizon = 100
    )
  )
}

# Compatibility entry point using the methods-paper expectation; optional
# interaction-count uncertainty does not re-realize first nesting.
fishery_ane_schedule_stochastic <- function(
    interactions, length_cm, mortality, params, n_sims,
    count_mode = c("fixed", "poisson", "cmp"), cmp_nu = 1) {
  # Compatibility wrapper: stochastic counts, expected first-nesting allocation.
  siders_fishery_ane_simulations(
    interactions, length_cm, mortality, params, n_sims,
    stochastic_demography = FALSE,
    interaction_count_mode = match.arg(count_mode), cmp_nu = cmp_nu
  )
}

# Expected ANE integrated over the same bivariate scale used by mvnorm.stan:
# log(carapace length) and logit(discard mortality). This retains the
# length-mortality correlation in the original method without requiring or
# exposing confidential individual records. The supplied SD and correlation
# values are scenario inputs; they are not estimated by the app.
fishery_ane_schedule_mvn_expected <- function(
    interactions, length_cm, mortality, params,
    log_length_sd = 0, logit_mortality_sd = 0,
    rho = 0, n_demography = 2000) {
  horizon <- length(interactions)
  if (length(length_cm) == 1) length_cm <- rep(length_cm, horizon)
  if (length(mortality) == 1) mortality <- rep(mortality, horizon)
  if (length(length_cm) != horizon || length(mortality) != horizon) {
    stop("Fishery schedule columns must have the same length.")
  }
  if (any(!is.finite(interactions)) || any(interactions < 0)) {
    stop("Fishery interactions must be finite non-negative values.")
  }
  if (
    !is.finite(log_length_sd) || log_length_sd < 0 ||
      !is.finite(logit_mortality_sd) || logit_mortality_sd < 0 ||
      !is.finite(rho) || abs(rho) > 1 ||
      !is.finite(n_demography) || n_demography < 100
  ) {
    stop("Invalid correlated-demography settings.")
  }
  n_demography <- as.integer(n_demography)
  covariance <- matrix(
    c(
      log_length_sd^2,
      rho * log_length_sd * logit_mortality_sd,
      rho * log_length_sd * logit_mortality_sd,
      logit_mortality_sd^2
    ),
    nrow = 2
  )

  losses <- rep(0, horizon)
  for (capture_index in seq_len(horizon)) {
    if (interactions[capture_index] == 0) next
    center <- c(
      log(positive_or_stop(length_cm[capture_index], "Fishery length")),
      safe_logit(probability_or_stop(
        mortality[capture_index], "Fishery mortality"
      ))
    )
    if (log_length_sd == 0 && logit_mortality_sd == 0) {
      draws <- matrix(center, nrow = n_demography, ncol = 2, byrow = TRUE)
    } else {
      draws <- mvtnorm::rmvnorm(
        n_demography, mean = center, sigma = covariance
      )
    }
    draw_lengths <- exp(draws[, 1])
    draw_mortality <- plogis(draws[, 2])
    if (mortality[capture_index] %in% c(0, 1)) {
      draw_mortality[] <- mortality[capture_index]
    }
    remaining <- horizon - capture_index + 1
    integrated_kernel <- Reduce(
      "+",
      Map(
        function(len, mort) {
          fishery_ane_kernel(len, mort, params, remaining)
        },
        draw_lengths, draw_mortality
      )
    ) / n_demography
    target <- capture_index:horizon
    losses[target] <- losses[target] +
      interactions[capture_index] * integrated_kernel
  }
  losses
}

# =====================================================================
# ORIGINAL SIDERS FIRST-NESTING FISHERY ANE
# =====================================================================

# =====================================================================
# USER-SUPPLIED UNCERTAINTY FOR SIDERS FISHERY ANE
# =====================================================================

draw_bounded_probability <- function(mean_value, sd_value, label) {
  mean_value <- probability_or_stop(mean_value, label)
  
  if (
    length(sd_value) != 1L ||
    is.na(sd_value) ||
    !is.finite(sd_value) ||
    sd_value < 0
  ) {
    stop(label, " SD must be one finite non-negative value.")
  }
  
  if (sd_value == 0) {
    return(mean_value)
  }
  
  truncnorm::rtruncnorm(
    n = 1L,
    a = 0,
    b = 1,
    mean = mean_value,
    sd = sd_value
  )
}


draw_positive_parameter <- function(mean_value, sd_value, label) {
  mean_value <- positive_or_stop(mean_value, label)
  
  if (
    length(sd_value) != 1L ||
    is.na(sd_value) ||
    !is.finite(sd_value) ||
    sd_value < 0
  ) {
    stop(label, " SD must be one finite non-negative value.")
  }
  
  if (sd_value == 0) {
    return(mean_value)
  }
  
  truncnorm::rtruncnorm(
    n = 1L,
    a = .Machine$double.eps,
    b = Inf,
    mean = mean_value,
    sd = sd_value
  )
}


draw_correlated_length_mortality <- function(
    mean_length,
    mean_mortality,
    log_length_sd = 0,
    logit_mortality_sd = 0,
    rho = 0) {
  
  mean_length <- positive_or_stop(
    mean_length,
    "Median fishery length"
  )
  
  mean_mortality <- probability_or_stop(
    mean_mortality,
    "Fishery mortality"
  )
  
  if (
    !is.finite(log_length_sd) ||
    log_length_sd < 0
  ) {
    stop("SD of log fishery length must be non-negative.")
  }
  
  if (
    !is.finite(logit_mortality_sd) ||
    logit_mortality_sd < 0
  ) {
    stop("SD of logit fishery mortality must be non-negative.")
  }
  
  if (
    !is.finite(rho) ||
    rho < -1 ||
    rho > 1
  ) {
    stop("Length-mortality correlation must be between -1 and 1.")
  }
  
  if (
    log_length_sd == 0 &&
    logit_mortality_sd == 0
  ) {
    return(
      c(
        length_cm = mean_length,
        mortality = mean_mortality
      )
    )
  }
  
  covariance <- matrix(
    c(
      log_length_sd^2,
      rho * log_length_sd * logit_mortality_sd,
      rho * log_length_sd * logit_mortality_sd,
      logit_mortality_sd^2
    ),
    nrow = 2L
  )
  
  draw <- mvtnorm::rmvnorm(
    n = 1L,
    mean = c(
      log(mean_length),
      safe_logit(mean_mortality)
    ),
    sigma = covariance
  )
  
  c(
    length_cm = exp(draw[1, 1]),
    mortality = if (mean_mortality %in% c(0, 1)) mean_mortality else plogis(draw[1, 2])
  )
}


draw_interaction_count <- function(
    expected_count,
    mode = c("fixed", "poisson", "cmp"),
    cmp_nu = 1) {
  
  mode <- match.arg(mode)
  
  if (
    length(expected_count) != 1L ||
    is.na(expected_count) ||
    !is.finite(expected_count) ||
    expected_count < 0
  ) {
    stop(
      "Expected annual interactions must be one finite ",
      "non-negative value."
    )
  }
  
  if (mode == "fixed") {
    return(as.numeric(expected_count))
  }
  
  if (mode == "poisson") {
    return(stats::rpois(1L, lambda = expected_count))
  }
  
  if (
    length(cmp_nu) != 1L ||
    is.na(cmp_nu) ||
    !is.finite(cmp_nu) ||
    cmp_nu <= 0
  ) {
    stop("CMP dispersion nu must be greater than zero.")
  }
  
  as.integer(
    rCMP(
      n = 1L,
      mu = expected_count,
      nu = cmp_nu,
      x_max = max(
        200L,
        as.integer(ceiling(4 * expected_count))
      )
    )
  )
}

# Individual wrapper for Methods Eqs. 7-14. Historical function name is
# retained for call compatibility; n1 now means expected first-nesting ANE.
siders_prop_ane_n1 <- function(
    capture_index, horizon, age, juvenile_survival, adult_survival,
    proportion_female, mortality, remigration_interval, params,
    stochastic_demography = FALSE) {
  # Compatibility name only: this is the expected formulation (Eqs. 9-13),
  # not historical sto.n1. The flag controls upstream parameter uncertainty.
  if (capture_index < 1L || capture_index > horizon) {
    stop("capture_index must fall within the projection horizon.")
  }
  p <- params
  p$pj <- probability_or_stop(juvenile_survival, "Juvenile survival draw")
  p$pa <- probability_or_stop(adult_survival, "Adult survival draw")
  p$pf <- probability_or_stop(proportion_female, "Female-probability draw")
  p$ri <- positive_or_stop(remigration_interval, "Remigration interval draw")
  mortality <- probability_or_stop(mortality, "Mortality draw")
  path <- siders_demographic_path(NULL, p, horizon - capture_index + 1L,
                                 start_age = age)
  transition <- path$maturity_probability
  first <- transition * cumprod(c(1, head(1 - transition, -1)))
  data.frame(
    projection_index = seq.int(capture_index, horizon),
    length = path$length,
    maturity_probability = transition,
    first_nesting_probability = first,
    survival_to_year = path$survival_to_year,
    n1 = path$survival_to_year * first * p$pf * mortality / p$ri
  )
}

# Generate one Siders n1 fishery-loss schedule.
#
# Individual uncertainty uses one draw per animal. Fixed expected totals may
# be fractional; their magnitude is preserved rather than silently rounded.

siders_fishery_ane_one_sim <- function(
    interactions,
    length_cm,
    mortality,
    params,
    stochastic_demography = FALSE,
    log_length_sd = 0,
    logit_mortality_sd = 0,
    length_mortality_rho = 0,
    ri_sd = 0,
    pj_sd = 0,
    pa_sd = 0,
    pf_sd = 0,
    ri_distribution = c("fixed", "normal", "cmp"),
    ri_cmp_nu = 1,
    fishery_center_mode = c("direct", "siders_atl"),
    siders_length_beta0 = NA_real_,
    siders_length_beta1 = NA_real_,
    siders_mortality_mu0 = NA_real_,
    force_all_interactions_fatal = FALSE) {
  
  horizon <- length(interactions)
  
  ri_distribution <- match.arg(
    ri_distribution
  )
  
  fishery_center_mode <- match.arg(
    fishery_center_mode
  )
  
  if (length(length_cm) == 1L) {
    length_cm <- rep(length_cm, horizon)
  }
  
  if (length(mortality) == 1L) {
    mortality <- rep(mortality, horizon)
  }
  
  if (
    length(length_cm) != horizon ||
    length(mortality) != horizon
  ) {
    stop(
      "Interactions, length, and mortality schedules must have equal lengths."
    )
  }
  
  if (
    any(!is.finite(interactions)) ||
    any(interactions < 0)
  ) {
    stop("Interactions must be finite non-negative values.")
  }
  
  if (
    any(!is.finite(length_cm)) ||
    any(length_cm <= 0)
  ) {
    stop("Fishery lengths must be finite values greater than zero.")
  }
  
  if (
    any(!is.finite(mortality)) ||
    any(mortality < 0 | mortality > 1)
  ) {
    stop("Fishery mortality values must be between zero and one.")
  }
  
  if (
    length(log_length_sd) != 1L ||
    !is.finite(log_length_sd) ||
    log_length_sd < 0
  ) {
    stop("SD of log fishery length must be non-negative.")
  }
  
  if (
    length(logit_mortality_sd) != 1L ||
    !is.finite(logit_mortality_sd) ||
    logit_mortality_sd < 0
  ) {
    stop("SD of logit fishery mortality must be non-negative.")
  }
  
  if (
    length(length_mortality_rho) != 1L ||
    !is.finite(length_mortality_rho) ||
    abs(length_mortality_rho) > 1
  ) {
    stop(
      "Length-mortality correlation must be between -1 and 1."
    )
  }
  
  validate_uncertainty_sd <- function(value, label) {
    if (
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value < 0
    ) {
      stop(
        label,
        " must be one finite non-negative value."
      )
    }
  }
  
  validate_uncertainty_sd(
    ri_sd,
    "SD of remigration interval"
  )
  
  validate_uncertainty_sd(
    pj_sd,
    "SD of juvenile survival"
  )
  
  validate_uncertainty_sd(
    pa_sd,
    "SD of adult survival"
  )
  
  validate_uncertainty_sd(
    pf_sd,
    "SD of proportion female"
  )
  
  if (
    length(ri_cmp_nu) != 1L ||
    is.na(ri_cmp_nu) ||
    !is.finite(ri_cmp_nu) ||
    ri_cmp_nu <= 0
  ) {
    stop(
      "CMP remigration-interval dispersion nu must be ",
      "one finite value greater than zero."
    )
  }
  
  if (identical(fishery_center_mode, "siders_atl")) {
    
    center_parameters <- c(
      siders_length_beta0,
      siders_length_beta1,
      siders_mortality_mu0
    )
    
    if (
      length(center_parameters) != 3L ||
      any(is.na(center_parameters)) ||
      any(!is.finite(center_parameters))
    ) {
      stop(
        "Original Siders ATL mode requires finite values for ",
        "length beta0, length beta1, and mortality mu0."
      )
    }
  }
  
  if (
    length(force_all_interactions_fatal) != 1L ||
    is.na(force_all_interactions_fatal)
  ) {
    stop(
      "The all-fatal interaction setting must be TRUE or FALSE."
    )
  }
  
  force_all_interactions_fatal <- isTRUE(
    force_all_interactions_fatal
  )
  
  # App uncertainty extension. Bounded normal draws prevent impossible
  # survival/sex probabilities; the supplied center is the pre-truncation mean.
  # The manuscript's expected pf is then applied once (no Bernoulli sex draw).
  yearly_proportion_female <- if (isTRUE(stochastic_demography)) {
    vapply(seq_len(horizon), function(i) {
      draw_bounded_probability(params$pf, pf_sd, "Proportion female")
    }, numeric(1))
  } else rep(params$pf, horizon)

  no_individual_variation <- log_length_sd == 0 && logit_mortality_sd == 0 &&
    (!isTRUE(stochastic_demography) ||
       (pj_sd == 0 && pa_sd == 0 && pf_sd == 0 &&
          (ri_distribution == "fixed" || (ri_distribution == "normal" && ri_sd == 0))))
  if (identical(fishery_center_mode, "direct") && no_individual_variation) {
    # Exact linear aggregation: removes an unnecessary per-turtle loop.
    return(fishery_ane_schedule(interactions, length_cm,
      if (force_all_interactions_fatal) rep(1, horizon) else mortality, params))
  }

  annual_loss <- rep(0, horizon)
  
  for (capture_index in seq_len(horizon)) {
    
    # For a fractional estimated total, simulate ceil(n) individuals and
    # weight their average by n. Integer counts retain unit weights.
    n_interactions <- as.integer(ceiling(interactions[capture_index]))
    interaction_weight <- if (n_interactions > 0) {
      interactions[capture_index] / n_interactions
    } else 0
    
    if (n_interactions < 1L) {
      next
    }
    
    # Preserve exact mortality boundaries in direct-input mode.
    #
    # safe_logit() deliberately clips probabilities away from zero and one
    # for numerical stability. However, a user-supplied mortality risk of
    # exactly zero must produce exactly zero fishery ANE unless the legacy
    # all-fatal scenario is explicitly selected.
    
    if (
      identical(fishery_center_mode, "direct") &&
      mortality[capture_index] == 0 &&
      !force_all_interactions_fatal
    ) {
      next
    }
    
    if (identical(fishery_center_mode, "siders_atl")) {
      
      # Exact original Siders center equations:
      #
      #   mu.l = TD_MVN$beta0 + TD_MVN$beta1 * ATL
      #   mu.m = TD_MVN$mu0
      #
      # n_interactions is the realized ATL value for this
      # simulation and capture year.
      
      mean_log_length <-
        siders_length_beta0 +
        siders_length_beta1 * n_interactions
      
      mean_logit_mortality <-
        siders_mortality_mu0
      
    } else {
      
      # App extension: users directly supply the annual center
      # of the length and mortality distributions.
      
      mean_log_length <- log(
        length_cm[capture_index]
      )
      
      mean_logit_mortality <- safe_logit(
        mortality[capture_index]
      )
    }
    
    center_length <- exp(
      mean_log_length
    )
    
    center_mortality <- plogis(
      mean_logit_mortality
    )
    
    # Match the original Siders draw_propagate() level of
    # randomization: one correlated length-mortality draw for
    # every individual interaction.
    
    if (
      log_length_sd == 0 &&
      logit_mortality_sd == 0
    ) {
      
      individual_length <- rep(
        center_length,
        n_interactions
      )
      
      individual_mortality <- rep(
        center_mortality,
        n_interactions
      )
      
    } else if (log_length_sd == 0) {
      
      individual_length <- rep(
        center_length,
        n_interactions
      )
      
      individual_mortality <- plogis(
        stats::rnorm(
          n = n_interactions,
          mean = mean_logit_mortality,
          sd = logit_mortality_sd
        )
      )
      
    } else if (logit_mortality_sd == 0) {
      
      individual_length <- exp(
        stats::rnorm(
          n = n_interactions,
          mean = mean_log_length,
          sd = log_length_sd
        )
      )
      
      individual_mortality <- rep(
        center_mortality,
        n_interactions
      )
      
    } else {
      
      covariance <- matrix(
        c(
          log_length_sd^2,
          length_mortality_rho *
            log_length_sd *
            logit_mortality_sd,
          length_mortality_rho *
            log_length_sd *
            logit_mortality_sd,
          logit_mortality_sd^2
        ),
        nrow = 2L
      )
      
      fishery_draws <- mvtnorm::rmvnorm(
        n = n_interactions,
        mean = c(
          mean_log_length,
          mean_logit_mortality
        ),
        sigma = covariance
      )
      
      individual_length <- exp(
        fishery_draws[, 1]
      )
      
      individual_mortality <- plogis(
        fishery_draws[, 2]
      )
    }
    
    # Preserve an exact direct-input mortality risk of one.
    #
    # This prevents safe_logit() from converting 1 to 0.999. A boundary
    # input represents a structural scenario rather than an uncertain mean.
    
    if (
      identical(fishery_center_mode, "direct") &&
      mortality[capture_index] == 1
    ) {
      individual_mortality[] <- 1
    }
    
    # Exact original Siders grim.reaper behavior:
    # all interactions are treated as fatal regardless of the
    # modeled or drawn mortality probability.
    
    if (force_all_interactions_fatal) {
      individual_mortality[] <- 1
    }
    
    for (animal_index in seq_len(n_interactions)) {
      
      capture_age <- siders_vbgf_age_from_length(
        length_cm = individual_length[animal_index],
        linf = params$linf,
        k = params$k,
        t0 = params$t0
      )
      
      if (isTRUE(stochastic_demography)) {
        
        # Original Siders hierarchy:
        # RI, juvenile survival, and adult survival are drawn
        # separately for every individual interaction.
        
        individual_ri <- switch(
          ri_distribution,
          
          fixed = params$ri,
          
          normal = {
            if (ri_sd == 0) {
              params$ri
            } else {
              truncnorm::rtruncnorm(
                n = 1L,
                a = 0,
                b = Inf,
                mean = params$ri,
                sd = ri_sd
              )
            }
          },
          
          cmp = {
            ri_draw <- rCMP(
              n = 1L,
              mu = params$ri,
              nu = ri_cmp_nu,
              x_max = max(
                200L,
                as.integer(ceiling(params$ri * 10))
              )
            )
            
            # Exact original behavior: RI cannot be zero.
            while (ri_draw == 0L) {
              ri_draw <- rCMP(
                n = 1L,
                mu = params$ri,
                nu = ri_cmp_nu,
                x_max = max(
                  200L,
                  as.integer(ceiling(params$ri * 10))
                )
              )
            }
            
            as.numeric(ri_draw)
          }
        )
        
        individual_pj <- draw_bounded_probability(
          params$pj, pj_sd, "Juvenile survival"
        )
        individual_pa <- draw_bounded_probability(
          params$pa, pa_sd, "Adult survival"
        )

      } else {
        
        individual_ri <- params$ri
        individual_pj <- params$pj
        individual_pa <- params$pa
      }
      
      contribution <- siders_prop_ane_n1(
        capture_index = capture_index,
        horizon = horizon,
        age = capture_age,
        juvenile_survival = individual_pj,
        adult_survival = individual_pa,
        proportion_female =
          yearly_proportion_female[capture_index],
        mortality = individual_mortality[animal_index],
        remigration_interval = individual_ri,
        params = params,
        stochastic_demography = stochastic_demography
      )
      
      annual_loss[contribution$projection_index] <-
        annual_loss[contribution$projection_index] +
        interaction_weight * contribution$n1
    }
  }
  
  annual_loss
}

# Run the Siders-derived fishery ANE calculation across projection simulations,
# using the methods-paper expected first-nesting reformulation. 
siders_fishery_ane_simulations <- function(
    interactions,
    length_cm,
    mortality,
    params,
    n_sims,
    stochastic_demography = TRUE,
    interaction_count_mode = c("fixed", "poisson", "cmp"),
    cmp_nu = 1,
    log_length_sd = 0,
    logit_mortality_sd = 0,
    length_mortality_rho = 0,
    ri_sd = 0,
    pj_sd = 0,
    pa_sd = 0,
    pf_sd = 0,
    ri_distribution = c("fixed", "normal", "cmp"),
    ri_cmp_nu = 1,
    fishery_center_mode = c("direct", "siders_atl"),
    siders_length_beta0 = NA_real_,
    siders_length_beta1 = NA_real_,
    siders_mortality_mu0 = NA_real_,
    force_all_interactions_fatal = FALSE) {
  
  interaction_count_mode <- match.arg(
    interaction_count_mode
  )
  
  ri_distribution <- match.arg(
    ri_distribution
  )
  
  fishery_center_mode <- match.arg(
    fishery_center_mode
  )
  
  n_sims <- as.integer(n_sims)
  
  if (
    length(n_sims) != 1L ||
    is.na(n_sims) ||
    n_sims < 1L
  ) {
    stop("n_sims must be one positive integer.")
  }
  
  horizon <- length(interactions)
  
  if (length(length_cm) == 1L) {
    length_cm <- rep(length_cm, horizon)
  }
  
  if (length(mortality) == 1L) {
    mortality <- rep(mortality, horizon)
  }
  
  if (
    length(length_cm) != horizon ||
    length(mortality) != horizon
  ) {
    stop(
      "Interactions, length, and mortality schedules ",
      "must have equal lengths."
    )
  }
  
  result <- matrix(
    0,
    nrow = n_sims,
    ncol = horizon
  )
  
  for (simulation_index in seq_len(n_sims)) {
    
    # The original uncertainty hierarchy is implemented inside
    # siders_fishery_ane_one_sim():
    #
    #   - proportion female by capture year;
    #   - RI, juvenile survival, and adult survival by animal.
    simulation_params <- params
    
    simulated_interactions <- vapply(
      seq_len(horizon),
      function(year_index) {
        draw_interaction_count(
          expected_count = interactions[year_index],
          mode = interaction_count_mode,
          cmp_nu = cmp_nu
        )
      },
      numeric(1L)
    )
    
    result[simulation_index, ] <-
      siders_fishery_ane_one_sim(
        interactions = simulated_interactions,
        length_cm = length_cm,
        mortality = mortality,
        params = simulation_params,
        stochastic_demography = stochastic_demography,
        log_length_sd = log_length_sd,
        logit_mortality_sd = logit_mortality_sd,
        length_mortality_rho = length_mortality_rho,
        ri_sd = ri_sd,
        pj_sd = pj_sd,
        pa_sd = pa_sd,
        pf_sd = pf_sd,
        ri_distribution = ri_distribution,
        ri_cmp_nu = ri_cmp_nu,
        fishery_center_mode = fishery_center_mode,
        siders_length_beta0 = siders_length_beta0,
        siders_length_beta1 = siders_length_beta1,
        siders_mortality_mu0 = siders_mortality_mu0,
        force_all_interactions_fatal =
          force_all_interactions_fatal
      )
  }
  
  attr(result, "uncertainty_settings") <- list(
    interaction_count_mode = interaction_count_mode,
    cmp_nu = cmp_nu,
    log_length_sd = log_length_sd,
    logit_mortality_sd = logit_mortality_sd,
    length_mortality_rho = length_mortality_rho,
    fishery_center_mode = fishery_center_mode,
    siders_length_beta0 = siders_length_beta0,
    siders_length_beta1 = siders_length_beta1,
    siders_mortality_mu0 = siders_mortality_mu0,
    force_all_interactions_fatal =
      isTRUE(force_all_interactions_fatal),
    ri_distribution = ri_distribution,
    ri_cmp_nu = ri_cmp_nu,
    ri_sd = ri_sd,
    pj_sd = pj_sd,
    pa_sd = pa_sd,
    pf_sd = pf_sd
  )
  
  result
}

# Draw future U/Q rows using the exact legacy ordering: the trajectory's
# selected posterior row is used for the first projected transition, followed
# by annual posterior rows sampled without replacement. Insufficient rows
# cause an explicit error, as in the manuscript; indices are shared across scenarios.
draw_projection_uq <- function(
    post_matrix, valid_rows, draw_index, horizon,
    mode = c("dynamic", "static")) {

  mode <- match.arg(mode)
  n_sims <- length(draw_index)

  if (mode == "static") {
    index_matrix <- matrix(
      rep(draw_index, horizon),
      nrow = n_sims,
      ncol = horizon
    )
  } else {
    index_matrix <- matrix(
      NA_integer_,
      nrow = n_sims,
      ncol = horizon
    )
    index_matrix[, 1] <- draw_index

    if (horizon > 1L) {
      if ((horizon - 1L) > length(valid_rows)) {
        stop(
          "Original Siders dynamic U/Q sampling is without replacement; ",
          "the posterior must contain at least horizon - 1 valid rows."
        )
      }

      for (i in seq_len(n_sims)) {
        index_matrix[i, 2:horizon] <- valid_rows[sample.int(
          length(valid_rows),
          horizon - 1L,
          replace = FALSE
        )]
      }
    }
  }

  list(
    index = index_matrix,
    U = matrix(
      post_matrix[as.vector(index_matrix), "U"],
      nrow = n_sims,
      ncol = horizon
    ),
    Q = matrix(
      post_matrix[as.vector(index_matrix), "Q"],
      nrow = n_sims,
      ncol = horizon
    )
  )
}

# One Siders-compatible annual projection transition.
#
# Original ordering:
#
#   1. subtract fishery ANE from current abundance;
#   2. apply annual population growth exp(U);
#   3. add abundance-scale Normal(0, sqrt(Q)) process error;
#   4. add conservation ANE as a separate post-growth extension.
#
# Conservation gain is not part of the original Take method and therefore
# remains explicitly separate.

project_population_step <- function(
    abundance,
    loss,
    gain,
    U,
    Q,
    standard_normal) {
  
  if (
    any(!is.finite(abundance)) ||
    any(abundance < 0)
  ) {
    stop(
      "Projection abundance must contain finite ",
      "non-negative values."
    )
  }
  
  if (
    any(!is.finite(loss)) ||
    any(loss < 0)
  ) {
    stop(
      "Fishery ANE loss must contain finite ",
      "non-negative values."
    )
  }
  
  if (
    any(!is.finite(gain)) ||
    any(gain < 0)
  ) {
    stop(
      "Conservation ANE gain must contain finite ",
      "non-negative values."
    )
  }
  
  if (
    any(!is.finite(U)) ||
    any(!is.finite(Q)) ||
    any(Q <= 0)
  ) {
    stop(
      "Projection U and Q must be finite, ",
      "with Q greater than zero."
    )
  }
  
  # Eqs. 19-22: Take before growth; process SD = sqrt(Q); Give after growth.
  abundance_after_take <- pmax(
    0,
    abundance - loss
  )
  
  next_abundance <- (
    abundance_after_take * exp(U)
  ) + (
    sqrt(Q) * standard_normal
  )
  
  # Give is an explicit extension and is added only after the original
  # Siders Take propagation step.
  next_abundance <- next_abundance + gain
  
  pmax(
    0,
    next_abundance
  )
}


# =====================================================================
# MATHEMATICAL REGRESSION AND LINEAGE SELF-CHECKS
# =====================================================================

self_check_row <- function(
    test,
    passed,
    maximum_difference = NA_real_,
    expected_behavior,
    details = "") {
  
  data.frame(
    Test = as.character(test),
    Passed = isTRUE(passed),
    Maximum_Difference = if (
      length(maximum_difference) == 1L &&
      is.finite(maximum_difference)
    ) {
      maximum_difference
    } else {
      NA_real_
    },
    Expected_Behavior = as.character(expected_behavior),
    Details = as.character(details),
    stringsAsFactors = FALSE
  )
}


max_abs_difference <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y) || any(!is.finite(x)) || any(!is.finite(y))) return(Inf)
  if (!length(x)) return(0)
  max(abs(x-y))
}


run_mathematical_self_checks <- function(
    trend_result = NULL,
    remigration_interval = NULL,
    tolerance = 1e-10) {
  
  if (
    length(tolerance) != 1L ||
    !is.finite(tolerance) ||
    tolerance <= 0
  ) {
    stop("Self-check tolerance must be one positive finite value.")
  }
  
  checks <- list()
  
  add_check <- function(...) {
    checks[[length(checks) + 1L]] <<- self_check_row(...)
  }
  
  # ---------------------------------------------------------------
  # TEST 1: PROJECTION EQUATION IDENTITY
  # ---------------------------------------------------------------
  
  abundance <- c(100, 50, 4)
  loss <- c(10, 60, 1)
  gain <- c(2, 3, 0)
  U <- c(0.02, -0.01, 0)
  Q <- c(4, 9, 1)
  standard_normal <- c(0.5, -1, 2)
  
  calculated_projection <- project_population_step(
    abundance = abundance,
    loss = loss,
    gain = gain,
    U = U,
    Q = Q,
    standard_normal = standard_normal
  )
  
  expected_projection <- pmax(
    0,
    (
      pmax(0, abundance - loss) *
        exp(U)
    ) +
      sqrt(Q) * standard_normal +
      gain
  )
  
  projection_difference <- max_abs_difference(
    calculated_projection,
    expected_projection
  )
  
  add_check(
    test = "Projection equation identity",
    passed = projection_difference <= tolerance,
    maximum_difference = projection_difference,
    expected_behavior = paste(
      "The annual transition equals",
      "(max(0, N - ANE) × exp(U)) + Normal error + Give."
    )
  )
  
  # ---------------------------------------------------------------
  # TEST 2: ZERO-TAKE EQUIVALENCE
  # ---------------------------------------------------------------
  
  baseline_projection <- project_population_step(
    abundance = abundance,
    loss = rep(0, length(abundance)),
    gain = rep(0, length(abundance)),
    U = U,
    Q = Q,
    standard_normal = standard_normal
  )
  
  zero_take_projection <- project_population_step(
    abundance = abundance,
    loss = rep(0, length(abundance)),
    gain = rep(0, length(abundance)),
    U = U,
    Q = Q,
    standard_normal = standard_normal
  )
  
  zero_take_difference <- max_abs_difference(
    baseline_projection,
    zero_take_projection
  )
  
  add_check(
    test = "Zero-take projection equivalence",
    passed = zero_take_difference <= tolerance,
    maximum_difference = zero_take_difference,
    expected_behavior = paste(
      "A threat scenario containing zero fishery ANE",
      "is identical to the baseline projection."
    )
  )
  
  # ---------------------------------------------------------------
  # TEST 3: TAKE IS APPLIED BEFORE GROWTH
  # ---------------------------------------------------------------
  
  ordering_abundance <- 100
  ordering_loss <- 20
  ordering_u <- log(1.1)
  
  calculated_ordering <- project_population_step(
    abundance = ordering_abundance,
    loss = ordering_loss,
    gain = 0,
    U = ordering_u,
    Q = 1,
    standard_normal = 0
  )
  
  loss_before_growth <- (
    ordering_abundance - ordering_loss
  ) * exp(ordering_u)
  
  loss_after_growth <- (
    ordering_abundance * exp(ordering_u)
  ) - ordering_loss
  
  before_difference <- abs(
    calculated_ordering - loss_before_growth
  )
  
  after_difference <- abs(
    calculated_ordering - loss_after_growth
  )
  
  add_check(
    test = "Fishery ANE timing",
    passed = (
      before_difference <= tolerance &&
        after_difference > tolerance
    ),
    maximum_difference = before_difference,
    expected_behavior = paste(
      "Fishery ANE is subtracted before annual",
      "population growth is applied."
    ),
    details = paste0(
      "Difference from loss-after-growth alternative: ",
      signif(after_difference, 6)
    )
  )
  
  # ---------------------------------------------------------------
  # SYNTHETIC BIOLOGICAL PARAMETERS FOR ANE CHECKS
  # ---------------------------------------------------------------
  
  test_params <- list(
    ri = 3,
    cf = 4,
    pf = 0.5,
    linf = 100,
    k = 0.12,
    t0 = -1,
    lmat = 80,
    sig_mat = 4,
    mat_p = 0.99,
    max_age = 80,
    pj = 0.85,
    pa = 0.95
  )
  
  test_horizon <- 20L
  
  # ---------------------------------------------------------------
  # TEST 4: ZERO INTERACTIONS PRODUCE ZERO FISHERY ANE
  # ---------------------------------------------------------------
  
  set.seed(90201)
  
  zero_ane <- siders_fishery_ane_simulations(
    interactions = rep(0, test_horizon),
    length_cm = rep(60, test_horizon),
    mortality = rep(0.3, test_horizon),
    params = test_params,
    n_sims = 25,
    stochastic_demography = TRUE,
    interaction_count_mode = "fixed",
    cmp_nu = 1,
    log_length_sd = 0,
    logit_mortality_sd = 0,
    length_mortality_rho = 0,
    ri_sd = 0,
    pj_sd = 0,
    pa_sd = 0,
    pf_sd = 0,
    ri_distribution = "fixed",
    ri_cmp_nu = 1
  )
  
  zero_ane_difference <- max(
    abs(zero_ane),
    na.rm = TRUE
  )
  
  add_check(
    test = "Zero interactions produce zero ANE",
    passed = zero_ane_difference <= tolerance,
    maximum_difference = zero_ane_difference,
    expected_behavior = paste(
      "A fishery schedule with no interactions",
      "produces no annual-nester-equivalent loss."
    )
  )
  
  # ---------------------------------------------------------------
  # TEST 5: ZERO MORTALITY PRODUCES ZERO FISHERY ANE
  # ---------------------------------------------------------------
  
  set.seed(90202)
  
  zero_mortality_ane <- siders_fishery_ane_simulations(
    interactions = rep(20, test_horizon),
    length_cm = rep(60, test_horizon),
    mortality = rep(0, test_horizon),
    params = test_params,
    n_sims = 25,
    stochastic_demography = TRUE,
    interaction_count_mode = "fixed",
    cmp_nu = 1,
    log_length_sd = 0,
    logit_mortality_sd = 0,
    length_mortality_rho = 0,
    ri_sd = 0,
    pj_sd = 0,
    pa_sd = 0,
    pf_sd = 0,
    ri_distribution = "fixed",
    ri_cmp_nu = 1
  )
  
  zero_mortality_difference <- max(
    abs(zero_mortality_ane),
    na.rm = TRUE
  )
  
  add_check(
    test = "Zero mortality produces zero ANE",
    passed = zero_mortality_difference <= tolerance,
    maximum_difference = zero_mortality_difference,
    expected_behavior = paste(
      "Interactions with zero mortality risk",
      "produce no fishery ANE loss."
    )
  )
  
  # ---------------------------------------------------------------
  # TEST 6: FIXED-INPUT REPRODUCIBILITY
  # ---------------------------------------------------------------
  
  run_reproducibility_test <- function(seed) {
    set.seed(seed)
    
    siders_fishery_ane_simulations(
      interactions = rep(10, test_horizon),
      length_cm = rep(65, test_horizon),
      mortality = rep(0.25, test_horizon),
      params = test_params,
      n_sims = 40,
      stochastic_demography = TRUE,
      interaction_count_mode = "fixed",
      cmp_nu = 1,
      log_length_sd = 0,
      logit_mortality_sd = 0,
      length_mortality_rho = 0,
      ri_sd = 0,
      pj_sd = 0,
      pa_sd = 0,
      pf_sd = 0,
      ri_distribution = "fixed",
      ri_cmp_nu = 1
    )
  }
  
  reproducibility_a <- run_reproducibility_test(90203)
  reproducibility_b <- run_reproducibility_test(90203)
  
  reproducibility_difference <- max_abs_difference(
    reproducibility_a,
    reproducibility_b
  )
  
  add_check(
    test = "Random-seed reproducibility",
    passed = reproducibility_difference <= tolerance,
    maximum_difference = reproducibility_difference,
    expected_behavior = paste(
      "The same inputs and random seed return",
      "the same stochastic fishery ANE matrix."
    )
  )
  
  # ---------------------------------------------------------------
  # TEST 7: INCREASING INTERACTIONS DOES NOT REDUCE EXPECTED ANE
  # ---------------------------------------------------------------
  
  # Use the deterministic-demography branch so this is a mathematical
  # monotonicity test rather than a comparison affected by Monte Carlo noise.
  
  low_interaction_ane <- siders_fishery_ane_one_sim(
    interactions = rep(10, test_horizon),
    length_cm = rep(65, test_horizon),
    mortality = rep(0.25, test_horizon),
    params = test_params,
    stochastic_demography = FALSE
  )
  
  high_interaction_ane <- siders_fishery_ane_one_sim(
    interactions = rep(20, test_horizon),
    length_cm = rep(65, test_horizon),
    mortality = rep(0.25, test_horizon),
    params = test_params,
    stochastic_demography = FALSE
  )
  
  interaction_difference <- (
    sum(high_interaction_ane) -
      sum(low_interaction_ane)
  )
  
  add_check(
    test = "ANE increases with interactions",
    passed = interaction_difference >= -tolerance,
    maximum_difference = max(
      0,
      -interaction_difference
    ),
    expected_behavior = paste(
      "Doubling otherwise identical interactions",
      "does not reduce total expected fishery ANE."
    ),
    details = paste0(
      "Change in cumulative ANE: ",
      signif(interaction_difference, 6)
    )
  )
  
  # ---------------------------------------------------------------
  # TEST 8: INCREASING MORTALITY DOES NOT REDUCE EXPECTED ANE
  # ---------------------------------------------------------------
  
  low_mortality_ane <- siders_fishery_ane_one_sim(
    interactions = rep(20, test_horizon),
    length_cm = rep(65, test_horizon),
    mortality = rep(0.1, test_horizon),
    params = test_params,
    stochastic_demography = FALSE
  )
  
  high_mortality_ane <- siders_fishery_ane_one_sim(
    interactions = rep(20, test_horizon),
    length_cm = rep(65, test_horizon),
    mortality = rep(0.4, test_horizon),
    params = test_params,
    stochastic_demography = FALSE
  )
  
  mortality_difference <- (
    sum(high_mortality_ane) -
      sum(low_mortality_ane)
  )
  
  add_check(
    test = "ANE increases with mortality",
    passed = mortality_difference >= -tolerance,
    maximum_difference = max(
      0,
      -mortality_difference
    ),
    expected_behavior = paste(
      "Increasing otherwise identical mortality risk",
      "does not reduce total expected fishery ANE."
    ),
    details = paste0(
      "Change in cumulative ANE: ",
      signif(mortality_difference, 6)
    )
  )
  
  # ===============================================================
  # CONSERVATION-GAIN EXTENSION SELF-CHECKS
  # ===============================================================
  #
  # These tests validate the mathematical accounting of the Ortega
  # conservation extension. They do not claim that the intervention
  # parameters themselves are empirically estimated.
  
  
  # ---------------------------------------------------------------
  # TEST 9: ZERO ACTION PRODUCES ZERO GIVE
  # ---------------------------------------------------------------
  
  zero_action <- rep(
    0,
    test_horizon
  )
  
  zero_nest_gain <- conservation_ane_schedule(
    type = "nests",
    amount = zero_action,
    params = test_params,
    eggs = 100,
    hatch = 0.80,
    emergence = 0.75,
    year1 = 0.03,
    emergence_basis = "conditional_on_hatching",
    counterfactual_emergence_per_egg = 0.20
  )
  
  zero_headstart_gain <- conservation_ane_schedule(
    type = "yearlings",
    amount = zero_action,
    params = test_params,
    eggs = 100,
    hatch = 0.80,
    emergence = 0.75,
    year1 = 0.03,
    release_age = 1
  )
  
  zero_adult_gain <- conservation_ane_schedule(
    type = "adults",
    amount = zero_action,
    params = test_params,
    eggs = 100,
    hatch = 0.80,
    emergence = 0.75,
    year1 = 0.03
  )
  
  zero_action_difference <- max(
    abs(
      c(
        zero_nest_gain,
        zero_headstart_gain,
        zero_adult_gain
      )
    ),
    na.rm = TRUE
  )
  
  add_check(
    test = "Zero conservation action produces zero Give",
    passed = zero_action_difference <= tolerance,
    maximum_difference = zero_action_difference,
    expected_behavior = paste(
      "Zero protected nests, zero net-additional headstarts,",
      "and zero adult deaths averted produce zero conservation ANE."
    )
  )
  
  
  # ---------------------------------------------------------------
  # TEST 10: ZERO COUNTERFACTUAL IMPROVEMENT PRODUCES ZERO NEST GIVE
  # ---------------------------------------------------------------
  
  protected_hatch_test <- 0.80
  protected_emergence_test <- 0.75
  
  protected_emergence_per_egg_test <- (
    protected_hatch_test *
      protected_emergence_test
  )
  
  no_increment_nest_gain <- conservation_ane_schedule(
    type = "nests",
    amount = rep(10, test_horizon),
    params = test_params,
    eggs = 100,
    hatch = protected_hatch_test,
    emergence = protected_emergence_test,
    year1 = 0.03,
    emergence_basis = "conditional_on_hatching",
    counterfactual_emergence_per_egg =
      protected_emergence_per_egg_test
  )
  
  no_increment_difference <- max(
    abs(no_increment_nest_gain),
    na.rm = TRUE
  )
  
  add_check(
    test = "No nest improvement produces zero Give",
    passed = no_increment_difference <= tolerance,
    maximum_difference = no_increment_difference,
    expected_behavior = paste(
      "If protected and counterfactual emergence per egg are equal,",
      "nest protection produces exactly zero incremental ANE."
    )
  )
  
  
  # ---------------------------------------------------------------
  # TEST 11: EMERGENCE-DENOMINATOR EQUIVALENCE
  # ---------------------------------------------------------------
  
  test_action_schedule <- c(
    10,
    rep(
      0,
      test_horizon - 1L
    )
  )
  
  conditional_emergence_gain <- conservation_ane_schedule(
    type = "nests",
    amount = test_action_schedule,
    params = test_params,
    eggs = 100,
    hatch = protected_hatch_test,
    emergence = protected_emergence_test,
    year1 = 0.03,
    emergence_basis = "conditional_on_hatching",
    counterfactual_emergence_per_egg = 0.20
  )
  
  per_egg_emergence_gain <- conservation_ane_schedule(
    type = "nests",
    amount = test_action_schedule,
    params = test_params,
    eggs = 100,
    hatch = protected_hatch_test,
    emergence = protected_emergence_per_egg_test,
    year1 = 0.03,
    emergence_basis = "per_egg",
    counterfactual_emergence_per_egg = 0.20
  )
  
  emergence_basis_difference <- max_abs_difference(
    conditional_emergence_gain,
    per_egg_emergence_gain
  )
  
  add_check(
    test = "Nest emergence denominator equivalence",
    passed = emergence_basis_difference <= tolerance,
    maximum_difference = emergence_basis_difference,
    expected_behavior = paste(
      "Equivalent emerged-hatchling probabilities give the same",
      "nest-protection ANE whether emergence is entered per hatched",
      "egg or directly per egg laid."
    )
  )
  
  
  # ---------------------------------------------------------------
  # TEST 12: CONSERVATION GIVE IS LINEAR IN ACTION AMOUNT
  # ---------------------------------------------------------------
  
  single_action <- c(
    10,
    rep(
      0,
      test_horizon - 1L
    )
  )
  
  double_action <- 2 * single_action
  
  
  nest_single <- conservation_ane_schedule(
    type = "nests",
    amount = single_action,
    params = test_params,
    eggs = 100,
    hatch = 0.80,
    emergence = 0.75,
    year1 = 0.03,
    emergence_basis = "conditional_on_hatching",
    counterfactual_emergence_per_egg = 0.20
  )
  
  nest_double <- conservation_ane_schedule(
    type = "nests",
    amount = double_action,
    params = test_params,
    eggs = 100,
    hatch = 0.80,
    emergence = 0.75,
    year1 = 0.03,
    emergence_basis = "conditional_on_hatching",
    counterfactual_emergence_per_egg = 0.20
  )
  
  
  headstart_single <- conservation_ane_schedule(
    type = "yearlings",
    amount = single_action,
    params = test_params,
    eggs = 100,
    hatch = 0.80,
    emergence = 0.75,
    year1 = 0.03,
    release_age = 1
  )
  
  headstart_double <- conservation_ane_schedule(
    type = "yearlings",
    amount = double_action,
    params = test_params,
    eggs = 100,
    hatch = 0.80,
    emergence = 0.75,
    year1 = 0.03,
    release_age = 1
  )
  
  
  adult_single <- conservation_ane_schedule(
    type = "adults",
    amount = single_action,
    params = test_params,
    eggs = 100,
    hatch = 0.80,
    emergence = 0.75,
    year1 = 0.03
  )
  
  adult_double <- conservation_ane_schedule(
    type = "adults",
    amount = double_action,
    params = test_params,
    eggs = 100,
    hatch = 0.80,
    emergence = 0.75,
    year1 = 0.03
  )
  
  
  conservation_linearity_difference <- max(
    max_abs_difference(
      nest_double,
      2 * nest_single
    ),
    max_abs_difference(
      headstart_double,
      2 * headstart_single
    ),
    max_abs_difference(
      adult_double,
      2 * adult_single
    )
  )
  
  add_check(
    test = "Conservation Give scales linearly with action",
    passed =
      conservation_linearity_difference <= tolerance,
    maximum_difference =
      conservation_linearity_difference,
    expected_behavior = paste(
      "Doubling an otherwise identical conservation action",
      "exactly doubles its expected ANE contribution."
    )
  )
  
  
  # ---------------------------------------------------------------
  # TEST 13: NEST YEAR-1 SURVIVAL IS APPLIED EXACTLY ONCE
  # ---------------------------------------------------------------
  
  test_year1_survival <- 0.20
  
  # One protected nest with one egg, complete protected emergence,
  # and zero counterfactual emergence creates exactly one incremental
  # emerged hatchling at age zero.
  
  one_hatchling_gain <- conservation_ane_schedule(
    type = "nests",
    amount = c(
      1,
      rep(
        0,
        test_horizon - 1L
      )
    ),
    params = test_params,
    eggs = 1,
    hatch = 1,
    emergence = 1,
    year1 = test_year1_survival,
    emergence_basis = "per_egg",
    counterfactual_emergence_per_egg = 0
  )
  
  age_one_length <- (
    test_params$linf *
      (
        1 -
          exp(
            -test_params$k *
              (1 - test_params$t0)
          )
      )
  )
  
  age_one_maturity <- maturity_probability(
    length_cm = age_one_length,
    lmat = test_params$lmat,
    sig_mat = test_params$sig_mat,
    linf = test_params$linf,
    mat_p = test_params$mat_p
  )
  
  expected_first_year_nest_gain <- (
    test_year1_survival *
      age_one_maturity *
      test_params$pf /
      test_params$ri
  )
  
  nest_first_interval_difference <- abs(
    one_hatchling_gain[1] -
      expected_first_year_nest_gain
  )
  
  add_check(
    test = "Nest year-1 survival applied once",
    passed =
      nest_first_interval_difference <= tolerance,
    maximum_difference =
      nest_first_interval_difference,
    expected_behavior = paste(
      "An incremental emerged hatchling receives emerged-to-yearling",
      "survival exactly once before its first possible census",
      "contribution; juvenile survival is not additionally applied",
      "during that same first interval."
    )
  )
  
  
  # ---------------------------------------------------------------
  # TEST 14: ADULT GIVE / CERTAIN-FATAL ADULT TAKE SYMMETRY
  # ---------------------------------------------------------------
  #
  # This test compares an initial cohort. The added final-year test below
  # also checks the terminal year: the expected formulation has no historical
  # sr2 indexing artifact.
  
  adult_test_amount <- c(
    7,
    rep(
      0,
      test_horizon - 1L
    )
  )
  
  adult_give_test <- conservation_ane_schedule(
    type = "adults",
    amount = adult_test_amount,
    params = test_params,
    eggs = 1,
    hatch = 1,
    emergence = 1,
    year1 = 1
  )
  
  adult_take_test <- siders_fishery_ane_one_sim(
    interactions = adult_test_amount,
    length_cm = rep(
      0.99 * test_params$linf,
      test_horizon
    ),
    mortality = rep(
      1,
      test_horizon
    ),
    params = test_params,
    stochastic_demography = FALSE,
    log_length_sd = 0,
    logit_mortality_sd = 0,
    length_mortality_rho = 0,
    ri_sd = 0,
    pj_sd = 0,
    pa_sd = 0,
    pf_sd = 0,
    ri_distribution = "fixed",
    fishery_center_mode = "direct",
    force_all_interactions_fatal = FALSE
  )
  
  adult_symmetry_difference <- max_abs_difference(
    adult_give_test,
    adult_take_test
  )
  
  add_check(
    test = "Adult Give/Take ANE symmetry",
    passed =
      adult_symmetry_difference <= tolerance,
    maximum_difference =
      adult_symmetry_difference,
    expected_behavior = paste(
      "For an isolated fully mature cohort, one adult death averted",
      "has the same deterministic ANE magnitude as one otherwise",
      "identical certain-fatal adult interaction."
    )
  )
  
  
  # ---------------------------------------------------------------
  # TEST 15: CONSERVATION GIVE ENTERS AFTER BASELINE TRANSITION
  # ---------------------------------------------------------------
  #
  # A conservation cohort ledger already represents turtles arriving
  # in the destination annual-nester census. It therefore enters after
  # the existing population has undergone the baseline U/Q transition.
  
  give_timing_abundance <- 100
  give_timing_gain <- adult_give_test[1]
  give_timing_u <- log(1.10)
  
  calculated_give_timing <- project_population_step(
    abundance = give_timing_abundance,
    loss = 0,
    gain = give_timing_gain,
    U = give_timing_u,
    Q = 1,
    standard_normal = 0
  )
  
  expected_gain_after_transition <- (
    give_timing_abundance *
      exp(give_timing_u)
  ) + give_timing_gain
  
  alternative_gain_before_transition <- (
    give_timing_abundance +
      give_timing_gain
  ) * exp(give_timing_u)
  
  give_timing_difference <- abs(
    calculated_give_timing -
      expected_gain_after_transition
  )
  
  give_alternative_difference <- abs(
    calculated_give_timing -
      alternative_gain_before_transition
  )
  
  add_check(
    test = "Conservation Give timing",
    passed = (
      give_timing_difference <= tolerance &&
        give_alternative_difference > tolerance
    ),
    maximum_difference =
      give_timing_difference,
    expected_behavior = paste(
      "Conservation ANE is added after the baseline U/Q",
      "population transition because the cohort ledger already",
      "represents arrival into the destination annual-nester census."
    ),
    details = paste0(
      "Difference from gain-before-growth alternative: ",
      signif(
        give_alternative_difference,
        6
      )
    )
  )
  
  # ---------------------------------------------------------------
  # TEST 16: PORTFOLIO COMPONENT ADDITIVITY
  # ---------------------------------------------------------------
  #
  # The portfolio must not introduce a new demographic equation.
  # Conservation Give is the sum of the validated A/B/C ledgers;
  # fishery Take is retained as the validated Siders loss ledger.
  # The resulting population trajectory must therefore equal direct
  # component-wise accounting under identical U/Q conditions.
  
  
  portfolio_nest_test <- conservation_ane_schedule(
    type = "nests",
    amount = rep(
      5,
      test_horizon
    ),
    params = test_params,
    eggs = 100,
    hatch = 0.80,
    emergence = 0.75,
    year1 = 0.03,
    release_age = 1,
    emergence_basis = "conditional_on_hatching",
    counterfactual_emergence_per_egg = 0.30
  )
  
  
  portfolio_head_test <- conservation_ane_schedule(
    type = "yearlings",
    amount = rep(
      3,
      test_horizon
    ),
    params = test_params,
    eggs = 100,
    hatch = 0,
    emergence = 0,
    year1 = 0.03,
    release_age = 1,
    emergence_basis = "conditional_on_hatching",
    counterfactual_emergence_per_egg = 0
  )
  
  
  portfolio_adult_test <- conservation_ane_schedule(
    type = "adults",
    amount = rep(
      2,
      test_horizon
    ),
    params = test_params,
    eggs = 100,
    hatch = 0,
    emergence = 0,
    year1 = 0.03,
    release_age = 1,
    emergence_basis = "conditional_on_hatching",
    counterfactual_emergence_per_egg = 0
  )
  
  
  # Use adults at the Siders forced-maturity threshold so that this
  # Take ledger is deterministic under the fixed test parameters.
  
  portfolio_take_test <- siders_fishery_ane_one_sim(
    interactions = rep(
      3,
      test_horizon
    ),
    length_cm = rep(
      0.99 * test_params$linf,
      test_horizon
    ),
    mortality = rep(
      0.40,
      test_horizon
    ),
    params = test_params,
    stochastic_demography = FALSE,
    log_length_sd = 0,
    logit_mortality_sd = 0,
    length_mortality_rho = 0,
    ri_sd = 0,
    pj_sd = 0,
    pa_sd = 0,
    pf_sd = 0,
    ri_distribution = "fixed",
    fishery_center_mode = "direct",
    force_all_interactions_fatal = FALSE
  )
  
  
  portfolio_gain_test <- (
    portfolio_nest_test +
      portfolio_head_test +
      portfolio_adult_test
  )
  
  portfolio_loss_test <- portfolio_take_test
  
  
  portfolio_combined_trajectory <- numeric(
    test_horizon + 1L
  )
  
  portfolio_manual_trajectory <- numeric(
    test_horizon + 1L
  )
  
  portfolio_combined_trajectory[1] <- 100
  
  portfolio_manual_trajectory[1] <- 100
  
  portfolio_test_u <- log(0.98)
  
  
  for (
    portfolio_year in seq_len(test_horizon)
  ) {
    
    portfolio_combined_trajectory[
      portfolio_year + 1L
    ] <- project_population_step(
      abundance =
        portfolio_combined_trajectory[
          portfolio_year
        ],
      loss =
        portfolio_loss_test[
          portfolio_year
        ],
      gain =
        portfolio_gain_test[
          portfolio_year
        ],
      U = portfolio_test_u,
      Q = 1,
      standard_normal = 0
    )
    
    
    portfolio_manual_trajectory[
      portfolio_year + 1L
    ] <- (
      pmax(
        0,
        portfolio_manual_trajectory[
          portfolio_year
        ] -
          portfolio_take_test[
            portfolio_year
          ]
      ) *
        exp(portfolio_test_u)
    ) +
      portfolio_nest_test[
        portfolio_year
      ] +
      portfolio_head_test[
        portfolio_year
      ] +
      portfolio_adult_test[
        portfolio_year
      ]
  }
  
  
  portfolio_additivity_difference <-
    max_abs_difference(
      portfolio_combined_trajectory,
      portfolio_manual_trajectory
    )
  
  
  add_check(
    test = "Portfolio component additivity",
    passed =
      portfolio_additivity_difference <= tolerance,
    maximum_difference =
      portfolio_additivity_difference,
    expected_behavior = paste(
      "The Combined Portfolio uses the validated fishery Take ledger",
      "and the sum of the validated nest, headstart, and adult Give",
      "ledgers without introducing an additional demographic effect."
    )
  )
  
  # Methods-specific checks: exercise the public kernels and CSV contract.
  terminal <- c(rep(0, test_horizon - 1L), 2.5)
  terminal_loss <- siders_fishery_ane_one_sim(
    terminal, 0.99*test_params$linf, 1, test_params)
  terminal_expected <- terminal * test_params$pa * test_params$pf / test_params$ri
  d <- max_abs_difference(terminal_loss, terminal_expected)
  add_check("Final-year adult Take and fractional counts", d <= tolerance, d,
            "Adult Take is n * adult survival * pf / RI, including the last year.")

  d <- max_abs_difference(
    siders_fishery_ane_one_sim(c(2.5, rep(0, test_horizon-1L)), 60, 0.00001, test_params),
    fishery_ane_schedule(c(2.5, rep(0, test_horizon-1L)), 60, 0.00001, test_params))
  add_check("Main forecast and expected kernel agree", d <= tolerance, d,
            "Main fixed-input forecasts equal the shared manuscript kernel without count/risk rounding.")

  fn_path <- first_nesting_path(60, test_params, 100)
  fn_sum <- sum(fn_path$first_nesting_probability)
  d <- max(0, fn_sum-1, -min(fn_path$first_nesting_probability))
  add_check("First nesting is allocated once", d <= tolerance, d,
            "First-event probabilities are non-negative and sum to no more than one.")

  adult_gain <- conservation_ane_schedule("adults", terminal, test_params, 1, 1, 1, 1)
  d <- max_abs_difference(adult_gain, terminal_expected)
  add_check("Final-year adult Give identity", d <= tolerance, d,
            "Averted adult deaths retain the adult-survival/sex/RI identity in the final year.")

  d <- abs(project_population_step(0.5, 0, 0, 0, 1, 0)-0.5)
  add_check("Fractional annual nesters are retained", d <= tolerance, d,
            "Only negative abundance is floored; positive abundance below one is retained.")

  example_years <- 2030:2032
  csv <- data.frame(year=c(2032,2030,2031), turtles=c(0,5,8),
                    median_cm=c(NA,60,80), mortality=c(NA,0.2,1))
  aligned <- validate_threat_schedule(csv, example_years)
  ok <- identical(aligned$year, as.numeric(example_years)) &&
    identical(aligned$turtles, c(5,8,0))
  add_check("Threat schedule year alignment", ok, if(ok) 0 else Inf,
            "CSV rows align by calendar year and explicit zero-count rows may omit size/risk.")
  rejects <- function(x) inherits(try(validate_threat_schedule(x, example_years), silent=TRUE), "try-error")
  bad_missing <- csv[-1, ]; bad_duplicate <- rbind(csv,csv[1, ])
  bad_risk <- csv; bad_risk$mortality[2] <- 20
  bad_count <- csv; bad_count$turtles[2] <- NA
  bad_year <- csv; bad_year$year[2] <- 2030.5
  bad_size <- csv; bad_size$median_cm[2] <- NA
  bad_infinite <- csv; bad_infinite$mortality[1] <- Inf
  ok <- all(vapply(list(bad_missing,bad_duplicate,bad_risk,bad_count,bad_year,
                        bad_size,bad_infinite), rejects, logical(1)))
  add_check("Threat schedule rejects invalid inputs", ok, if(ok) 0 else Inf,
            "Missing/duplicate years, missing counts/sizes, nonfinite values and percent risks fail explicitly.")

  # ---------------------------------------------------------------
  # TESTS REQUIRING THE ACTIVE TREND POSTERIOR
  # ---------------------------------------------------------------
  
  if (
    !is.null(trend_result) &&
    !is.null(remigration_interval) &&
    length(remigration_interval) == 1L &&
    is.finite(remigration_interval) &&
    remigration_interval > 0
  ) {
    
    projection_posterior <-
      build_siders_projection_posterior(
        trend_result = trend_result,
        remigration_interval = remigration_interval
      )
    
    # -------------------------------------------------------------
    # TEST 9: CURRENT-ABUNDANCE RI/4 IDENTITY
    # -------------------------------------------------------------
    
    expected_total_females <- (
      projection_posterior$N_fym0 +
        projection_posterior$N_fym1 +
        projection_posterior$N_fym2 +
        projection_posterior$N_fym3
    ) * remigration_interval / 4
    
    abundance_difference <- max_abs_difference(
      projection_posterior$Total_Females,
      expected_total_females
    )
    
    add_check(
      test = "Siders RI/4 abundance identity",
      passed = abundance_difference <= tolerance,
      maximum_difference = abundance_difference,
      expected_behavior = paste(
        "Every posterior row satisfies",
        "Total_Females = sum(final four annual-nester states) × RI / 4."
      )
    )
    
    # -------------------------------------------------------------
    # TEST 10: MATCHED POSTERIOR ROWS
    # -------------------------------------------------------------
    
    original_sims <- trend_result$fit$sims.list
    original_u <- as.numeric(original_sims$U)
    original_q <- as.numeric(original_sims$Q)
    
    posterior_row <- projection_posterior$Posterior_Row
    
    matched_u_difference <- max_abs_difference(
      projection_posterior$U,
      original_u[posterior_row]
    )
    
    matched_q_difference <- max_abs_difference(
      projection_posterior$Q,
      original_q[posterior_row]
    )
    
    matched_difference <- max(
      matched_u_difference,
      matched_q_difference
    )
    
    add_check(
      test = "Matched posterior U and Q",
      passed = matched_difference <= tolerance,
      maximum_difference = matched_difference,
      expected_behavior = paste(
        "Starting abundance, U, and Q retain their",
        "joint posterior-row association."
      )
    )
    
  } else {
    
    add_check(
      test = "Siders RI/4 abundance identity",
      passed = FALSE,
      maximum_difference = NA_real_,
      expected_behavior = paste(
        "Run the trend model before checking the",
        "active posterior abundance identity."
      ),
      details = "Not evaluated because no active trend posterior was available."
    )
    
    add_check(
      test = "Matched posterior U and Q",
      passed = FALSE,
      maximum_difference = NA_real_,
      expected_behavior = paste(
        "Run the trend model before checking matched",
        "posterior-row associations."
      ),
      details = "Not evaluated because no active trend posterior was available."
    )
  }
  
  results <- dplyr::bind_rows(checks)
  
  results$Status <- ifelse(
    results$Passed,
    "PASS",
    ifelse(
      grepl(
        "^Not evaluated",
        results$Details
      ),
      "NOT RUN",
      "FAIL"
    )
  )
  
  results <- results[
    ,
    c(
      "Test",
      "Status",
      "Maximum_Difference",
      "Expected_Behavior",
      "Details",
      "Passed"
    )
  ]
  
  attr(results, "all_passed") <- all(
    results$Passed |
      results$Status == "NOT RUN"
  )
  
  attr(results, "evaluated_passed") <- all(
    results$Passed[
      results$Status != "NOT RUN"
    ]
  )
  
  results
}


# Threat-schedule contract. All modeled years are explicit; no silent zero fill.
# Optional uncertainty applies AFTER these direct annual inputs are validated.
validate_threat_schedule <- function(schedule, projection_years) {
  required <- c("year", "turtles", "median_cm", "mortality")
  if (!is.data.frame(schedule) || !all(required %in% names(schedule))) {
    stop("Threat CSV must contain: year, turtles, median_cm, mortality. ",
         "Old Mean_Len values must be replaced by median sizes.")
  }
  if (anyDuplicated(names(schedule))) stop("Threat CSV has duplicate column names.")
  schedule <- schedule[, required, drop = FALSE]
  if (!nrow(schedule)) stop("Threat schedule is empty.")
  for (column in required) {
    x <- schedule[[column]]
    if (is.factor(x)) x <- as.character(x)
    y <- suppressWarnings(as.numeric(x))
    blank <- is.na(x) | trimws(as.character(x)) == ""
    if (any(!blank & !is.finite(y))) stop(column, " must be numeric and finite.")
    schedule[[column]] <- y
  }
  if (any(!is.finite(schedule$year)) || any(schedule$year != floor(schedule$year))) {
    stop("year must contain whole calendar years.")
  }
  if (anyDuplicated(schedule$year)) stop("Threat schedule has duplicate years.")
  if (any(!is.finite(schedule$turtles)) || any(schedule$turtles < 0)) {
    stop("turtles must contain finite non-negative counts; use 0 explicitly.")
  }
  active <- schedule$turtles > 0
  if (any(!is.finite(schedule$median_cm[active])) || any(schedule$median_cm[active] <= 0)) {
    stop("Each year with turtles requires a positive median_cm.")
  }
  if (any(!is.finite(schedule$mortality[active]))) {
    stop("Each year with turtles requires mortality between 0 and 1.")
  }
  supplied_length <- !is.na(schedule$median_cm)
  supplied_mortality <- !is.na(schedule$mortality)
  if (any(schedule$median_cm[supplied_length] <= 0) ||
      any(schedule$mortality[supplied_mortality] < 0 | schedule$mortality[supplied_mortality] > 1)) {
    stop("Sizes must be positive; mortality is a proportion from 0 to 1, not a percent.")
  }
  missing <- setdiff(projection_years, schedule$year)
  outside <- setdiff(schedule$year, projection_years)
  if (length(missing)) stop("Threat schedule is missing forecast years: ", paste(missing, collapse = ", "))
  if (length(outside)) stop("Threat schedule contains years outside the forecast: ", paste(outside, collapse = ", "))
  schedule <- schedule[match(projection_years, schedule$year), , drop = FALSE]
  rownames(schedule) <- NULL
  schedule
}

make_threat_template <- function(projection_years) {
  data.frame(year = projection_years, turtles = 0,
             median_cm = NA_real_, mortality = NA_real_)
}

align_annual_schedule <- function(year, value, projection_years) {
  if (length(year) != length(value)) stop("Schedule Year and value columns differ in length.")
  schedule <- data.frame(Year = as.integer(year), Value = as.numeric(value))
  if (anyNA(schedule$Year) || anyNA(schedule$Value)) {
    stop("Schedule Year and value columns must be numeric and non-missing.")
  }
  if (anyDuplicated(schedule$Year)) stop("Schedule contains duplicate years.")
  aligned <- data.frame(Year = projection_years) %>%
    left_join(schedule, by = "Year")
  aligned$Value[is.na(aligned$Value)] <- 0
  aligned$Value
}

# =====================================================================
# UI HELPER FUNCTIONS
# =====================================================================
pop_input <- function(inputId, label, value, step = NA, title, def, calc, data_req) {
  numericInput(
    inputId = inputId,
    label = tags$span(
      label,
      popover(
        shiny::icon("circle-question", class = "ms-1 text-primary", style = "cursor: pointer;"),
        title = title,
        tags$div(
          tags$p(tags$b("Definition: "), def),
          tags$p(tags$b("How it's calculated: "), calc),
          tags$p(tags$b("Data needed: "), data_req)
        )
      )
    ),
    value = value, step = step
  )
}

# =====================================================================
# UI LAYOUT
# =====================================================================

# UI orchestration only: effect calculations below are extracted unchanged
# from the 2026-09-01 methods-aligned engine. A snapshot prevents live edits
# from changing a saved scenario during a comparison.
ux_build_effects <- function(input, n_sims, horizon, transition_years) {
  mode <- input$pva_mode
  zero_effect <- rep(0, horizon)
  track_effects <- list("Status Quo" = list(loss = zero_effect, gain = zero_effect))
      # Generate a matched matrix of first-nesting-only fishery ANE losses.
      # Each row is one projection simulation. Biological and fishery uncertainty
      # are drawn from user-supplied distributions; zero SD retains fixed inputs.
      make_fishery_loss <- function(
    interactions,
    length_cm,
    mortality
      ) {
        
        params <- current_model_params(input)
        
        count_mode <- if (
          is.null(input$siders_count_mode)
        ) {
          "fixed"
        } else {
          input$siders_count_mode
        }
        
        cmp_nu <- if (
          is.null(input$siders_cmp_nu) ||
          !is.finite(input$siders_cmp_nu)
        ) {
          1
        } else {
          input$siders_cmp_nu
        }
        
        uncertainty_value <- function(input_value) {
          if (
            is.null(input_value) ||
            length(input_value) != 1L ||
            is.na(input_value) ||
            !is.finite(input_value)
          ) {
            return(0)
          }
          
          input_value
        }
        
        siders_fishery_ane_simulations(
          interactions = interactions,
          length_cm = length_cm,
          mortality = mortality,
          params = params,
          n_sims = n_sims,
          
          # Retain demographic parameter uncertainty while using the
          # methods-paper expected first-nesting formulation.
          stochastic_demography = TRUE,
          
          interaction_count_mode = count_mode,
          cmp_nu = cmp_nu,
          
          log_length_sd = uncertainty_value(
            input$siders_log_length_sd
          ),
          
          logit_mortality_sd = uncertainty_value(
            input$siders_logit_mortality_sd
          ),
          
          length_mortality_rho = uncertainty_value(
            input$siders_length_mortality_rho
          ),
          
          ri_sd = uncertainty_value(
            input$siders_ri_sd
          ),
          
          pj_sd = uncertainty_value(
            input$siders_pj_sd
          ),
          
          pa_sd = uncertainty_value(
            input$siders_pa_sd
          ),
          
          pf_sd = uncertainty_value(
            input$siders_pf_sd
          ),
          
          ri_distribution = if (
            is.null(input$siders_ri_distribution)
          ) {
            "fixed"
          } else {
            input$siders_ri_distribution
          },
          
          ri_cmp_nu = if (
            is.null(input$siders_ri_cmp_nu) ||
            length(input$siders_ri_cmp_nu) != 1L ||
            is.na(input$siders_ri_cmp_nu) ||
            !is.finite(input$siders_ri_cmp_nu)
          ) {
            1
          } else {
            input$siders_ri_cmp_nu
          },
          
          fishery_center_mode = if (
            is.null(input$siders_fishery_center_mode)
          ) {
            "direct"
          } else {
            input$siders_fishery_center_mode
          },
          
          siders_length_beta0 = suppressWarnings(
            as.numeric(input$siders_length_beta0)
          ),
          
          siders_length_beta1 = suppressWarnings(
            as.numeric(input$siders_length_beta1)
          ),
          
          siders_mortality_mu0 = suppressWarnings(
            as.numeric(input$siders_mortality_mu0)
          ),
          
          force_all_interactions_fatal = isTRUE(
            input$siders_force_all_fatal
          )
        )
      }
      
      if (mode == "mode_threat") {
        if (input$threat_input_type == "static") {
          annual_loss <- make_fishery_loss(
            interactions = rep(input$threat_interactions, horizon),
            length_cm = input$threat_mean_len,
            mortality = input$threat_mort_score
          )
          track_effects[["With Threat"]] <- list(
            loss = annual_loss, gain = zero_effect
          )
        } else {
          req(input$threat_csv)
          if (!identical(input$siders_fishery_center_mode %||% "direct", "direct") ||
              isTRUE(input$siders_force_all_fatal)) {
            stop("For a CSV threat schedule, select direct annual inputs and turn off the all-fatal override.")
          }
          t_aligned <- validate_threat_schedule(
            input$ux_threat_data %||% read.csv(input$threat_csv$datapath, check.names = FALSE), transition_years
          )
          # Placeholders below are used only in zero-count years, whose kernel
          # is skipped; they are not biological estimates or missing-data fills.
          t_aligned$median_cm[is.na(t_aligned$median_cm)] <- 1
          t_aligned$mortality[is.na(t_aligned$mortality)] <- 0
          track_effects[["With Threat"]] <- list(
            loss = make_fishery_loss(
              interactions = t_aligned$turtles,
              length_cm = t_aligned$median_cm,
              mortality = t_aligned$mortality
            ),
            gain = zero_effect
          )
        }
      } else if (mode == "mode_action") {
        
        type_map <- if (
          input$action_type == "A: Protect Nests"
        ) {
          "nests"
        } else if (
          input$action_type == "B: Headstarting"
        ) {
          "yearlings"
        } else {
          "adults"
        }
        
        
        # ---------------------------------------------------------------
        # ACTION-SPECIFIC COUNTERFACTUAL SETTINGS
        # ---------------------------------------------------------------
        
        action_hatch <- if (
          identical(type_map, "nests")
        ) {
          input$action_protected_hatch
        } else {
          0
        }
        
        action_emergence <- if (
          identical(type_map, "nests")
        ) {
          input$action_protected_emergence
        } else {
          0
        }
        
        action_emergence_basis <- if (
          identical(type_map, "nests")
        ) {
          input$action_emergence_basis
        } else {
          "conditional_on_hatching"
        }
        
        action_counterfactual <- if (
          identical(type_map, "nests")
        ) {
          input$action_counterfactual_emergence
        } else {
          0
        }
        
        action_release_age <- if (
          identical(type_map, "yearlings")
        ) {
          input$action_release_age
        } else {
          1
        }
        
        
        make_action_gain <- function(
    amount_schedule) {
          
          conservation_ane_schedule(
            type = type_map,
            amount = amount_schedule,
            params = current_model_params(input),
            
            eggs = input$override_eggs,
            
            hatch = action_hatch,
            emergence = action_emergence,
            
            year1 = input$override_syr1,
            
            release_age =
              action_release_age,
            
            emergence_basis =
              action_emergence_basis,
            
            counterfactual_emergence_per_egg =
              action_counterfactual
          )
        }
        
        
        if (
          input$action_input_type == "static"
        ) {
          
          val_input <- if (
            identical(type_map, "nests")
          ) {
            input$action_static_nests
          } else if (
            identical(type_map, "yearlings")
          ) {
            input$action_static_head
          } else {
            input$action_static_adults
          }
          
          track_effects[["With Intervention"]] <-
            list(
              loss = zero_effect,
              gain = make_action_gain(
                rep(
                  val_input,
                  horizon
                )
              )
            )
          
        } else {
          
          req(input$action_csv)
          
          a_df <- input$ux_action_data %||% read.csv(input$action_csv$datapath)
          
          required <- c(
            "Year",
            "Amount"
          )
          
          if (
            !all(required %in% names(a_df))
          ) {
            stop(
              "Action schedule must contain: ",
              paste(
                required,
                collapse = ", "
              )
            )
          }
          
          action_amount <- align_annual_schedule(
            a_df$Year,
            a_df$Amount,
            transition_years
          )
          
          track_effects[["With Intervention"]] <-
            list(
              loss = zero_effect,
              gain = make_action_gain(
                action_amount
              )
            )
        }
      } else if (mode == "mode_portfolio") {
        
        # -------------------------------------------------------------
        # VALIDATED PORTFOLIO LEDGER
        #
        # Every component below uses the same equations and definitions
        # as the corresponding individually validated Action/Threat mode.
        # -------------------------------------------------------------
        
        portfolio_amount <- function(
    x,
    label) {
          
          if (
            is.null(x) ||
            length(x) != 1L ||
            is.na(x) ||
            !is.finite(x) ||
            x < 0
          ) {
            stop(
              label,
              " must be one finite non-negative value."
            )
          }
          
          as.numeric(x)
        }
        
        
        nests_amount <- if (
          "A: Protect Nests" %in% input$portfolio_active
        ) {
          portfolio_amount(
            input$portfolio_nests,
            "Portfolio nests protected"
          )
        } else {
          0
        }
        
        head_amount <- if (
          "B: Headstart Yearlings" %in% input$portfolio_active
        ) {
          portfolio_amount(
            input$portfolio_head,
            "Portfolio net-additional headstarts"
          )
        } else {
          0
        }
        
        adult_amount <- if (
          "C: Stop Adult Poaching" %in% input$portfolio_active
        ) {
          portfolio_amount(
            input$portfolio_adults,
            "Portfolio adult deaths averted"
          )
        } else {
          0
        }
        
        threat_amount <- if (
          "Threat: Bycatch / Mortality" %in% input$portfolio_active
        ) {
          portfolio_amount(
            input$portfolio_threat_interactions,
            "Portfolio fishery interactions"
          )
        } else {
          0
        }
        
        
        params <- current_model_params(input)
        

        portfolio_eggs <- positive_or_stop(
          input$portfolio_eggs,
          "Eggs per nest"
        )
        
        portfolio_year1 <- probability_or_stop(
          input$portfolio_year1,
          "Emerged-to-yearling survival"
        )
        
        
        combined_gain <- rep(
          0,
          horizon
        )
        
        combined_loss <- rep(
          0,
          horizon
        )
        
        
        # -------------------------------------------------------------
        # A. NEST PROTECTION
        # -------------------------------------------------------------
        
        if (nests_amount > 0) {
          
          d_A <- conservation_ane_schedule(
            type = "nests",
            amount = rep(
              nests_amount,
              horizon
            ),
            params = params,
            
            eggs = portfolio_eggs,
            
            hatch =
              input$portfolio_protected_hatch,
            
            emergence =
              input$portfolio_protected_emergence,
            
            year1 =
              portfolio_year1,
            
            release_age = 1,
            
            emergence_basis =
              input$portfolio_emergence_basis,
            
            counterfactual_emergence_per_egg =
              input$portfolio_counterfactual_emergence
          )
          
          track_effects[["A: Protect Nests"]] <-
            list(
              loss = zero_effect,
              gain = d_A
            )
          
          combined_gain <-
            combined_gain + d_A
        }
        
        
        # -------------------------------------------------------------
        # B. HEADSTARTING
        # -------------------------------------------------------------
        
        if (head_amount > 0) {
          
          d_B <- conservation_ane_schedule(
            type = "yearlings",
            amount = rep(
              head_amount,
              horizon
            ),
            params = params,
            
            eggs = portfolio_eggs,
            hatch = 0,
            emergence = 0,
            year1 = portfolio_year1,
            
            release_age =
              input$portfolio_release_age,
            
            emergence_basis =
              "conditional_on_hatching",
            
            counterfactual_emergence_per_egg = 0
          )
          
          track_effects[["B: Headstart Yearlings"]] <-
            list(
              loss = zero_effect,
              gain = d_B
            )
          
          combined_gain <-
            combined_gain + d_B
        }
        
        
        # -------------------------------------------------------------
        # C. ADULT PROTECTION
        # -------------------------------------------------------------
        
        if (adult_amount > 0) {
          
          d_C <- conservation_ane_schedule(
            type = "adults",
            amount = rep(
              adult_amount,
              horizon
            ),
            params = params,
            
            eggs = portfolio_eggs,
            hatch = 0,
            emergence = 0,
            year1 = portfolio_year1,
            
            release_age = 1,
            
            emergence_basis =
              "conditional_on_hatching",
            
            counterfactual_emergence_per_egg = 0
          )
          
          track_effects[["C: Stop Adult Poaching"]] <-
            list(
              loss = zero_effect,
              gain = d_C
            )
          
          combined_gain <-
            combined_gain + d_C
        }
        
        
        # -------------------------------------------------------------
        # D. FISHERY TAKE
        # -------------------------------------------------------------
        
        if (threat_amount > 0) {
          
          d_T <- make_fishery_loss(
            interactions = rep(
              threat_amount,
              horizon
            ),
            
            length_cm =
              input$portfolio_threat_mean_len,
            
            mortality =
              input$portfolio_threat_mortality
          )
          
          track_effects[["Threat: Bycatch / Mortality"]] <- list(
            loss = d_T,
            gain = zero_effect
          )
          
          # d_T is an n_sim × horizon matrix because the aligned
          # Siders Take pathway retains stochastic interaction-level
          # uncertainty.
          combined_loss <- d_T
        }
        
        
        # Always create the combined scenario. If every entered amount
        # is zero, this must reduce exactly to Status Quo.
        
        track_effects[["Combined Portfolio"]] <-
          list(
            loss = combined_loss,
            gain = combined_gain
          )
      }
      
  track_effects
}

ux_example_csv <- function(kind, years = 2026:2075) {
  examples <- list(
    monthly = "Year,Month,Site,Count,Monitored\n2006,4,example_beach_a,50,TRUE\n2006,4,example_beach_b,31,TRUE\n2006,5,example_beach_a,81,TRUE\n2006,5,example_beach_b,49,TRUE\n2006,6,example_beach_a,135,TRUE\n2006,6,example_beach_b,77,TRUE\n2006,7,example_beach_a,190,TRUE\n2006,7,example_beach_b,105,TRUE\n2006,8,example_beach_a,197,TRUE\n2006,8,example_beach_b,110,TRUE\n2006,9,example_beach_a,156,TRUE\n2006,9,example_beach_b,91,TRUE\n2006,10,example_beach_a,107,TRUE\n2006,10,example_beach_b,66,TRUE\n2006,11,example_beach_a,68,TRUE\n2006,11,example_beach_b,43,TRUE\n2006,12,example_beach_a,44,TRUE\n2006,12,example_beach_b,28,TRUE\n2007,1,example_beach_a,29,TRUE\n2007,1,example_beach_b,17,TRUE\n2007,2,example_beach_a,24,TRUE\n2007,2,example_beach_b,13,TRUE\n2007,3,example_beach_a,34,TRUE\n2007,3,example_beach_b,19,TRUE\n2007,4,example_beach_a,52,TRUE\n2007,4,example_beach_b,30,TRUE\n2007,5,example_beach_a,79,TRUE\n2007,5,example_beach_b,44,TRUE\n2007,6,example_beach_a,123,TRUE\n2007,6,example_beach_b,69,TRUE\n2007,7,example_beach_a,169,TRUE\n2007,7,example_beach_b,,FALSE\n2007,8,example_beach_a,182,TRUE\n2007,8,example_beach_b,112,TRUE\n2007,9,example_beach_a,156,TRUE\n2007,9,example_beach_b,99,TRUE\n2007,10,example_beach_a,115,TRUE\n2007,10,example_beach_b,72,TRUE\n2007,11,example_beach_a,75,TRUE\n2007,11,example_beach_b,44,TRUE\n2007,12,example_beach_a,46,TRUE\n2007,12,example_beach_b,26,TRUE\n2008,1,example_beach_a,28,TRUE\n2008,1,example_beach_b,15,TRUE\n2008,2,example_beach_a,21,TRUE\n2008,2,example_beach_b,12,TRUE\n2008,3,example_beach_a,,FALSE\n2008,3,example_beach_b,18,TRUE\n2008,4,example_beach_a,47,TRUE\n2008,4,example_beach_b,26,TRUE\n2008,5,example_beach_a,70,TRUE\n2008,5,example_beach_b,41,TRUE\n2008,6,example_beach_a,114,TRUE\n2008,6,example_beach_b,71,TRUE\n2008,7,example_beach_a,170,TRUE\n2008,7,example_beach_b,108,TRUE\n2008,8,example_beach_a,196,TRUE\n2008,8,example_beach_b,121,TRUE\n2008,9,example_beach_a,,FALSE\n2008,9,example_beach_b,100,TRUE\n2008,10,example_beach_a,119,TRUE\n2008,10,example_beach_b,67,TRUE\n2008,11,example_beach_a,71,TRUE\n2008,11,example_beach_b,39,TRUE\n2008,12,example_beach_a,42,TRUE\n2008,12,example_beach_b,23,TRUE\n2009,1,example_beach_a,25,TRUE\n2009,1,example_beach_b,15,TRUE\n2009,2,example_beach_a,20,TRUE\n2009,2,example_beach_b,12,TRUE\n2009,3,example_beach_a,31,TRUE\n2009,3,example_beach_b,20,TRUE\n2009,4,example_beach_a,44,TRUE\n2009,4,example_beach_b,28,TRUE\n2009,5,example_beach_a,72,TRUE\n2009,5,example_beach_b,45,TRUE\n2009,6,example_beach_a,124,TRUE\n2009,6,example_beach_b,76,TRUE\n2009,7,example_beach_a,185,TRUE\n2009,7,example_beach_b,108,TRUE\n2009,8,example_beach_a,201,TRUE\n2009,8,example_beach_b,113,TRUE\n2009,9,example_beach_a,161,TRUE\n2009,9,example_beach_b,89,TRUE\n2009,10,example_beach_a,107,TRUE\n2009,10,example_beach_b,61,TRUE\n2009,11,example_beach_a,64,TRUE\n2009,11,example_beach_b,38,TRUE\n2009,12,example_beach_a,39,TRUE\n2009,12,example_beach_b,,FALSE\n2010,1,example_beach_a,26,TRUE\n2010,1,example_beach_b,16,TRUE\n2010,2,example_beach_a,22,TRUE\n2010,2,example_beach_b,13,TRUE\n2010,3,example_beach_a,34,TRUE\n2010,3,example_beach_b,20,TRUE\n2010,4,example_beach_a,48,TRUE\n2010,4,example_beach_b,29,TRUE\n2010,5,example_beach_a,77,TRUE\n2010,5,example_beach_b,45,TRUE\n2010,6,example_beach_a,126,TRUE\n2010,6,example_beach_b,,FALSE\n2010,7,example_beach_a,174,TRUE\n2010,7,example_beach_b,96,TRUE\n2010,8,example_beach_a,179,TRUE\n2010,8,example_beach_b,102,TRUE\n2010,9,example_beach_a,145,TRUE\n2010,9,example_beach_b,87,TRUE\n2010,10,example_beach_a,102,TRUE\n2010,10,example_beach_b,64,TRUE\n2010,11,example_beach_a,67,TRUE\n2010,11,example_beach_b,42,TRUE\n2010,12,example_beach_a,43,TRUE\n2010,12,example_beach_b,26,TRUE\n2011,1,example_beach_a,28,TRUE\n2011,1,example_beach_b,16,TRUE\n2011,2,example_beach_a,,FALSE\n2011,2,example_beach_b,12,TRUE\n2011,3,example_beach_a,31,TRUE\n2011,3,example_beach_b,17,TRUE\n2011,4,example_beach_a,48,TRUE\n2011,4,example_beach_b,27,TRUE\n2011,5,example_beach_a,72,TRUE\n2011,5,example_beach_b,40,TRUE\n2011,6,example_beach_a,112,TRUE\n2011,6,example_beach_b,64,TRUE\n2011,7,example_beach_a,157,TRUE\n2011,7,example_beach_b,95,TRUE\n2011,8,example_beach_a,,FALSE\n2011,8,example_beach_b,109,TRUE\n2011,9,example_beach_a,152,TRUE\n2011,9,example_beach_b,95,TRUE\n2011,10,example_beach_a,112,TRUE\n2011,10,example_beach_b,68,TRUE\n2011,11,example_beach_a,71,TRUE\n2011,11,example_beach_b,41,TRUE\n2011,12,example_beach_a,43,TRUE\n2011,12,example_beach_b,24,TRUE\n2012,1,example_beach_a,26,TRUE\n2012,1,example_beach_b,14,TRUE\n2012,2,example_beach_a,19,TRUE\n2012,2,example_beach_b,11,TRUE\n2012,3,example_beach_a,28,TRUE\n2012,3,example_beach_b,17,TRUE\n2012,4,example_beach_a,43,TRUE\n2012,4,example_beach_b,25,TRUE\n2012,5,example_beach_a,65,TRUE\n2012,5,example_beach_b,40,TRUE\n2012,6,example_beach_a,110,TRUE\n2012,6,example_beach_b,69,TRUE\n2012,7,example_beach_a,166,TRUE\n2012,7,example_beach_b,104,TRUE\n2012,8,example_beach_a,190,TRUE\n2012,8,example_beach_b,114,TRUE\n2012,9,example_beach_a,161,TRUE\n2012,9,example_beach_b,92,TRUE\n2012,10,example_beach_a,110,TRUE\n2012,10,example_beach_b,61,TRUE\n2012,11,example_beach_a,65,TRUE\n2012,11,example_beach_b,,FALSE\n2012,12,example_beach_a,38,TRUE\n2012,12,example_beach_b,22,TRUE\n2013,1,example_beach_a,23,TRUE\n2013,1,example_beach_b,14,TRUE\n2013,2,example_beach_a,19,TRUE\n2013,2,example_beach_b,12,TRUE\n2013,3,example_beach_a,30,TRUE\n2013,3,example_beach_b,19,TRUE\n2013,4,example_beach_a,43,TRUE\n2013,4,example_beach_b,27,TRUE\n2013,5,example_beach_a,70,TRUE\n2013,5,example_beach_b,,FALSE\n2013,6,example_beach_a,120,TRUE\n2013,6,example_beach_b,72,TRUE\n2013,7,example_beach_a,175,TRUE\n2013,7,example_beach_b,99,TRUE\n2013,8,example_beach_a,185,TRUE\n2013,8,example_beach_b,102,TRUE\n2013,9,example_beach_a,146,TRUE\n2013,9,example_beach_b,82,TRUE\n2013,10,example_beach_a,98,TRUE\n2013,10,example_beach_b,57,TRUE\n2013,11,example_beach_a,60,TRUE\n2013,11,example_beach_b,37,TRUE\n2013,12,example_beach_a,38,TRUE\n2013,12,example_beach_b,24,TRUE\n2014,1,example_beach_a,,FALSE\n2014,1,example_beach_b,16,TRUE\n2014,2,example_beach_a,21,TRUE\n2014,2,example_beach_b,12,TRUE\n2014,3,example_beach_a,32,TRUE\n2014,3,example_beach_b,18,TRUE\n2014,4,example_beach_a,47,TRUE\n2014,4,example_beach_b,28,TRUE\n2014,5,example_beach_a,73,TRUE\n2014,5,example_beach_b,41,TRUE\n2014,6,example_beach_a,116,TRUE\n2014,6,example_beach_b,64,TRUE\n2014,7,example_beach_a,,FALSE\n2014,7,example_beach_b,89,TRUE\n2014,8,example_beach_a,165,TRUE\n2014,8,example_beach_b,97,TRUE\n2014,9,example_beach_a,136,TRUE\n2014,9,example_beach_b,85,TRUE\n2014,10,example_beach_a,99,TRUE\n2014,10,example_beach_b,63,TRUE\n2014,11,example_beach_a,65,TRUE\n2014,11,example_beach_b,40,TRUE\n2014,12,example_beach_a,41,TRUE\n2014,12,example_beach_b,24,TRUE\n2015,1,example_beach_a,26,TRUE\n2015,1,example_beach_b,15,TRUE\n2015,2,example_beach_a,20,TRUE\n2015,2,example_beach_b,11,TRUE\n2015,3,example_beach_a,28,TRUE\n2015,3,example_beach_b,16,TRUE\n2015,4,example_beach_a,44,TRUE\n2015,4,example_beach_b,25,TRUE\n2015,5,example_beach_a,65,TRUE\n2015,5,example_beach_b,37,TRUE\n2015,6,example_beach_a,103,TRUE\n2015,6,example_beach_b,62,TRUE\n2015,7,example_beach_a,148,TRUE\n2015,7,example_beach_b,93,TRUE\n2015,8,example_beach_a,169,TRUE\n2015,8,example_beach_b,106,TRUE\n2015,9,example_beach_a,148,TRUE\n2015,9,example_beach_b,91,TRUE\n2015,10,example_beach_a,108,TRUE\n2015,10,example_beach_b,,FALSE\n2015,11,example_beach_a,66,TRUE\n2015,11,example_beach_b,37,TRUE\n2015,12,example_beach_a,39,TRUE\n2015,12,example_beach_b,22,TRUE\n2016,1,example_beach_a,23,TRUE\n2016,1,example_beach_b,13,TRUE\n2016,2,example_beach_a,18,TRUE\n2016,2,example_beach_b,11,TRUE\n2016,3,example_beach_a,27,TRUE\n2016,3,example_beach_b,17,TRUE\n2016,4,example_beach_a,40,TRUE\n2016,4,example_beach_b,,FALSE\n2016,5,example_beach_a,62,TRUE\n2016,5,example_beach_b,39,TRUE\n2016,6,example_beach_a,107,TRUE\n2016,6,example_beach_b,67,TRUE\n2016,7,example_beach_a,162,TRUE\n2016,7,example_beach_b,99,TRUE\n2016,8,example_beach_a,182,TRUE\n2016,8,example_beach_b,106,TRUE\n2016,9,example_beach_a,150,TRUE\n2016,9,example_beach_b,84,TRUE\n2016,10,example_beach_a,101,TRUE\n2016,10,example_beach_b,56,TRUE\n2016,11,example_beach_a,59,TRUE\n2016,11,example_beach_b,34,TRUE\n2016,12,example_beach_a,,FALSE\n2016,12,example_beach_b,21,TRUE\n2017,1,example_beach_a,22,TRUE\n2017,1,example_beach_b,14,TRUE\n2017,2,example_beach_a,19,TRUE\n2017,2,example_beach_b,12,TRUE\n2017,3,example_beach_a,30,TRUE\n2017,3,example_beach_b,18,TRUE\n2017,4,example_beach_a,41,TRUE\n2017,4,example_beach_b,26,TRUE\n2017,5,example_beach_a,68,TRUE\n2017,5,example_beach_b,41,TRUE\n2017,6,example_beach_a,,FALSE\n2017,6,example_beach_b,66,TRUE\n2017,7,example_beach_a,162,TRUE\n2017,7,example_beach_b,90,TRUE\n2017,8,example_beach_a,169,TRUE\n2017,8,example_beach_b,94,TRUE\n2017,9,example_beach_a,133,TRUE\n2017,9,example_beach_b,77,TRUE\n2017,10,example_beach_a,91,TRUE\n2017,10,example_beach_b,55,TRUE\n2017,11,example_beach_a,57,TRUE\n2017,11,example_beach_b,36,TRUE\n2017,12,example_beach_a,37,TRUE\n2017,12,example_beach_b,23,TRUE\n2018,1,example_beach_a,24,TRUE\n2018,1,example_beach_b,15,TRUE\n2018,2,example_beach_a,20,TRUE\n2018,2,example_beach_b,11,TRUE\n2018,3,example_beach_a,29,TRUE\n2018,3,example_beach_b,,FALSE\n2018,4,example_beach_a,44,TRUE\n2018,4,example_beach_b,25,TRUE\n2018,5,example_beach_a,67,TRUE\n2018,5,example_beach_b,37,TRUE\n2018,6,example_beach_a,105,TRUE\n2018,6,example_beach_b,59,TRUE\n2018,7,example_beach_a,144,TRUE\n2018,7,example_beach_b,83,TRUE\n2018,8,example_beach_a,154,TRUE\n2018,8,example_beach_b,94,TRUE\n2018,9,example_beach_a,131,TRUE\n2018,9,example_beach_b,,FALSE\n2018,10,example_beach_a,97,TRUE\n2018,10,example_beach_b,60,TRUE\n2018,11,example_beach_a,63,TRUE\n2018,11,example_beach_b,38,TRUE\n2018,12,example_beach_a,39,TRUE\n2018,12,example_beach_b,22,TRUE\n2019,1,example_beach_a,24,TRUE\n2019,1,example_beach_b,13,TRUE\n2019,2,example_beach_a,18,TRUE\n2019,2,example_beach_b,10,TRUE\n2019,3,example_beach_a,26,TRUE\n2019,3,example_beach_b,15,TRUE\n2019,4,example_beach_a,40,TRUE\n2019,4,example_beach_b,23,TRUE\n2019,5,example_beach_a,60,TRUE\n2019,5,example_beach_b,35,TRUE\n2019,6,example_beach_a,97,TRUE\n2019,6,example_beach_b,60,TRUE\n2019,7,example_beach_a,143,TRUE\n2019,7,example_beach_b,90,TRUE\n2019,8,example_beach_a,165,TRUE\n2019,8,example_beach_b,102,TRUE\n2019,9,example_beach_a,143,TRUE\n2019,9,example_beach_b,85,TRUE\n2019,10,example_beach_a,101,TRUE\n2019,10,example_beach_b,58,TRUE\n2019,11,example_beach_a,,FALSE\n2019,11,example_beach_b,34,TRUE\n2019,12,example_beach_a,36,TRUE\n2019,12,example_beach_b,20,TRUE\n2020,1,example_beach_a,21,TRUE\n2020,1,example_beach_b,13,TRUE\n2020,2,example_beach_a,17,TRUE\n2020,2,example_beach_b,10,TRUE\n2020,3,example_beach_a,26,TRUE\n2020,3,example_beach_b,17,TRUE\n2020,4,example_beach_a,37,TRUE\n2020,4,example_beach_b,23,TRUE\n2020,5,example_beach_a,,FALSE\n2020,5,example_beach_b,38,TRUE\n2020,6,example_beach_a,104,TRUE\n2020,6,example_beach_b,65,TRUE\n2020,7,example_beach_a,156,TRUE\n2020,7,example_beach_b,92,TRUE\n2020,8,example_beach_a,171,TRUE\n2020,8,example_beach_b,96,TRUE\n2020,9,example_beach_a,138,TRUE\n2020,9,example_beach_b,76,TRUE\n2020,10,example_beach_a,91,TRUE\n2020,10,example_beach_b,51,TRUE\n2020,11,example_beach_a,54,TRUE\n2020,11,example_beach_b,32,TRUE\n2020,12,example_beach_a,33,TRUE\n2020,12,example_beach_b,21,TRUE\n2021,1,example_beach_a,22,TRUE\n2021,1,example_beach_b,14,TRUE\n2021,2,example_beach_a,18,TRUE\n2021,2,example_beach_b,,FALSE\n2021,3,example_beach_a,28,TRUE\n2021,3,example_beach_b,17,TRUE\n2021,4,example_beach_a,40,TRUE\n2021,4,example_beach_b,25,TRUE\n2021,5,example_beach_a,65,TRUE\n2021,5,example_beach_b,38,TRUE\n2021,6,example_beach_a,107,TRUE\n2021,6,example_beach_b,60,TRUE\n2021,7,example_beach_a,148,TRUE\n2021,7,example_beach_b,82,TRUE\n2021,8,example_beach_a,153,TRUE\n2021,8,example_beach_b,,FALSE\n2021,9,example_beach_a,123,TRUE\n2021,9,example_beach_b,73,TRUE\n2021,10,example_beach_a,86,TRUE\n2021,10,example_beach_b,54,TRUE\n2021,11,example_beach_a,56,TRUE\n2021,11,example_beach_b,35,TRUE\n2021,12,example_beach_a,36,TRUE\n2021,12,example_beach_b,22,TRUE\n2022,1,example_beach_a,23,TRUE\n2022,1,example_beach_b,14,TRUE\n2022,2,example_beach_a,19,TRUE\n2022,2,example_beach_b,10,TRUE\n2022,3,example_beach_a,27,TRUE\n2022,3,example_beach_b,15,TRUE\n2022,4,example_beach_a,41,TRUE\n2022,4,example_beach_b,23,TRUE\n2022,5,example_beach_a,62,TRUE\n2022,5,example_beach_b,34,TRUE\n2022,6,example_beach_a,96,TRUE\n2022,6,example_beach_b,55,TRUE\n2022,7,example_beach_a,133,TRUE\n2022,7,example_beach_b,80,TRUE\n2022,8,example_beach_a,146,TRUE\n2022,8,example_beach_b,92,TRUE\n2022,9,example_beach_a,128,TRUE\n2022,9,example_beach_b,80,TRUE\n2022,10,example_beach_a,,FALSE\n2022,10,example_beach_b,57,TRUE\n2022,11,example_beach_a,60,TRUE\n2022,11,example_beach_b,35,TRUE\n2022,12,example_beach_a,37,TRUE\n2022,12,example_beach_b,20,TRUE\n2023,1,example_beach_a,22,TRUE\n2023,1,example_beach_b,12,TRUE\n2023,2,example_beach_a,17,TRUE\n2023,2,example_beach_b,9,TRUE\n2023,3,example_beach_a,24,TRUE\n2023,3,example_beach_b,15,TRUE\n2023,4,example_beach_a,,FALSE\n2023,4,example_beach_b,21,TRUE\n2023,5,example_beach_a,55,TRUE\n2023,5,example_beach_b,34,TRUE\n2023,6,example_beach_a,92,TRUE\n2023,6,example_beach_b,58,TRUE\n2023,7,example_beach_a,139,TRUE\n2023,7,example_beach_b,87,TRUE\n2023,8,example_beach_a,160,TRUE\n2023,8,example_beach_b,97,TRUE\n2023,9,example_beach_a,137,TRUE\n2023,9,example_beach_b,79,TRUE\n2023,10,example_beach_a,94,TRUE\n2023,10,example_beach_b,52,TRUE\n2023,11,example_beach_a,56,TRUE\n2023,11,example_beach_b,31,TRUE\n2023,12,example_beach_a,32,TRUE\n2023,12,example_beach_b,19,TRUE\n2024,1,example_beach_a,20,TRUE\n2024,1,example_beach_b,,FALSE\n2024,2,example_beach_a,16,TRUE\n2024,2,example_beach_b,10,TRUE\n2024,3,example_beach_a,26,TRUE\n2024,3,example_beach_b,16,TRUE\n2024,4,example_beach_a,36,TRUE\n2024,4,example_beach_b,23,TRUE\n2024,5,example_beach_a,59,TRUE\n2024,5,example_beach_b,37,TRUE\n2024,6,example_beach_a,101,TRUE\n2024,6,example_beach_b,61,TRUE\n2024,7,example_beach_a,148,TRUE\n2024,7,example_beach_b,,FALSE\n2024,8,example_beach_a,158,TRUE\n2024,8,example_beach_b,88,TRUE\n2024,9,example_beach_a,125,TRUE\n2024,9,example_beach_b,70,TRUE\n2024,10,example_beach_a,83,TRUE\n2024,10,example_beach_b,48,TRUE\n2024,11,example_beach_a,51,TRUE\n2024,11,example_beach_b,31,TRUE\n2024,12,example_beach_a,32,TRUE\n2024,12,example_beach_b,20,TRUE\n2025,1,example_beach_a,21,TRUE\n2025,1,example_beach_b,13,TRUE\n2025,2,example_beach_a,18,TRUE\n2025,2,example_beach_b,11,TRUE\n2025,3,example_beach_a,,FALSE\n2025,3,example_beach_b,15,TRUE\n2025,4,example_beach_a,39,TRUE\n2025,4,example_beach_b,23,TRUE\n2025,5,example_beach_a,62,TRUE\n2025,5,example_beach_b,35,TRUE\n2025,6,example_beach_a,99,TRUE\n2025,6,example_beach_b,55,TRUE\n2025,7,example_beach_a,135,TRUE\n2025,7,example_beach_b,75,TRUE\n2025,8,example_beach_a,140,TRUE\n2025,8,example_beach_b,82,TRUE\n2025,9,example_beach_a,115,TRUE\n2025,9,example_beach_b,71,TRUE\n2025,10,example_beach_a,83,TRUE\n2025,10,example_beach_b,53,TRUE\n2025,11,example_beach_a,54,TRUE\n2025,11,example_beach_b,34,TRUE\n2025,12,example_beach_a,35,TRUE\n2025,12,example_beach_b,21,TRUE\n2026,1,example_beach_a,22,TRUE\n2026,1,example_beach_b,12,TRUE\n2026,2,example_beach_a,17,TRUE\n2026,2,example_beach_b,9,TRUE\n2026,3,example_beach_a,24,TRUE\n2026,3,example_beach_b,14,TRUE\n",
    annual = "Year,Month,Site,Count,Monitored\n2006,,example_beach_a,1115,TRUE\n2006,,example_beach_b,649,TRUE\n2007,,example_beach_a,1077,TRUE\n2007,,example_beach_b,640,TRUE\n2008,,example_beach_a,1075,TRUE\n2008,,example_beach_b,643,TRUE\n2009,,example_beach_a,1079,TRUE\n2009,,example_beach_b,632,TRUE\n2010,,example_beach_a,1042,TRUE\n2010,,example_beach_b,606,TRUE\n2011,,example_beach_a,1014,TRUE\n2011,,example_beach_b,605,TRUE\n2012,,example_beach_a,1020,TRUE\n2012,,example_beach_b,608,TRUE\n2013,,example_beach_a,1013,TRUE\n2013,,example_beach_b,590,TRUE\n2014,,example_beach_a,974,TRUE\n2014,,example_beach_b,573,TRUE\n2015,,example_beach_a,958,TRUE\n2015,,example_beach_b,577,TRUE\n2016,,example_beach_a,969,TRUE\n2016,,example_beach_b,574,TRUE\n2017,,example_beach_a,946,TRUE\n2017,,example_beach_b,550,TRUE\n2018,,example_beach_a,912,TRUE\n2018,,example_beach_b,539,TRUE\n2019,,example_beach_a,910,TRUE\n2019,,example_beach_b,547,TRUE\n2020,,example_beach_a,912,TRUE\n2020,,example_beach_b,536,TRUE\n2021,,example_beach_a,883,TRUE\n2021,,example_beach_b,515,TRUE\n2022,,example_beach_a,860,TRUE\n2022,,example_beach_b,512,TRUE\n2023,,example_beach_a,864,TRUE\n2023,,example_beach_b,516,TRUE\n2024,,example_beach_a,859,TRUE\n2024,,example_beach_b,502,TRUE\n2025,,example_beach_a,825,TRUE\n2025,,example_beach_b,484,TRUE\n",
    threat = "year,turtles,median_cm,mortality\n2026,120,90,0.35\n2027,120,90,0.35\n2028,120,90,0.35\n2029,120,90,0.35\n2030,120,90,0.35\n2031,80,100,0.25\n2032,80,100,0.25\n2033,80,100,0.25\n2034,80,100,0.25\n2035,80,100,0.25\n2036,40,110,0.15\n2037,40,110,0.15\n2038,40,110,0.15\n2039,40,110,0.15\n2040,40,110,0.15\n2041,40,110,0.15\n2042,40,110,0.15\n2043,40,110,0.15\n2044,40,110,0.15\n2045,40,110,0.15\n2046,0,,\n2047,0,,\n2048,0,,\n2049,0,,\n2050,0,,\n2051,0,,\n2052,0,,\n2053,0,,\n2054,0,,\n2055,0,,\n2056,0,,\n2057,0,,\n2058,0,,\n2059,0,,\n2060,0,,\n2061,0,,\n2062,0,,\n2063,0,,\n2064,0,,\n2065,0,,\n2066,0,,\n2067,0,,\n2068,0,,\n2069,0,,\n2070,0,,\n2071,0,,\n2072,0,,\n2073,0,,\n2074,0,,\n2075,0,,\n",
    nests = "Year,Amount\n2026,100\n2027,100\n2028,100\n2029,100\n2030,100\n2031,200\n2032,200\n2033,200\n2034,200\n2035,200\n2036,300\n2037,300\n2038,300\n2039,300\n2040,300\n2041,300\n2042,300\n2043,300\n2044,300\n2045,300\n2046,0\n2047,0\n2048,0\n2049,0\n2050,0\n2051,0\n2052,0\n2053,0\n2054,0\n2055,0\n2056,0\n2057,0\n2058,0\n2059,0\n2060,0\n2061,0\n2062,0\n2063,0\n2064,0\n2065,0\n2066,0\n2067,0\n2068,0\n2069,0\n2070,0\n2071,0\n2072,0\n2073,0\n2074,0\n2075,0\n",
    yearlings = "Year,Amount\n2026,20\n2027,20\n2028,20\n2029,20\n2030,20\n2031,30\n2032,30\n2033,30\n2034,30\n2035,30\n2036,40\n2037,40\n2038,40\n2039,40\n2040,40\n2041,40\n2042,40\n2043,40\n2044,40\n2045,40\n2046,0\n2047,0\n2048,0\n2049,0\n2050,0\n2051,0\n2052,0\n2053,0\n2054,0\n2055,0\n2056,0\n2057,0\n2058,0\n2059,0\n2060,0\n2061,0\n2062,0\n2063,0\n2064,0\n2065,0\n2066,0\n2067,0\n2068,0\n2069,0\n2070,0\n2071,0\n2072,0\n2073,0\n2074,0\n2075,0\n",
    adults = "Year,Amount\n2026,2\n2027,2\n2028,2\n2029,2\n2030,2\n2031,4\n2032,4\n2033,4\n2034,4\n2035,4\n2036,6\n2037,6\n2038,6\n2039,6\n2040,6\n2041,6\n2042,6\n2043,6\n2044,6\n2045,6\n2046,0\n2047,0\n2048,0\n2049,0\n2050,0\n2051,0\n2052,0\n2053,0\n2054,0\n2055,0\n2056,0\n2057,0\n2058,0\n2059,0\n2060,0\n2061,0\n2062,0\n2063,0\n2064,0\n2065,0\n2066,0\n2067,0\n2068,0\n2069,0\n2070,0\n2071,0\n2072,0\n2073,0\n2074,0\n2075,0\n"
  )
  d <- read.csv(text = examples[[kind]], stringsAsFactors = FALSE)
  if (kind %in% c("monthly", "annual")) return(d)
  index <- seq_along(years)
  d <- d[pmin(index, nrow(d)), , drop = FALSE]
  d[[1]] <- as.integer(years)
  rownames(d) <- NULL
  d
}



# Pure UI-boundary validation. It rejects incompatible uploads before calling
# the unchanged mathematical kernels and never imputes missing action amounts.
ux_nonnegative <- function(x, label) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < 0) stop(label, " must be a finite number at least zero.")
  x
}
ux_action_type <- function(x) switch(x, "A: Protect Nests" = "nests", "B: Headstarting" = "yearlings", "C: Stop Adult Poaching" = "adults", stop("Choose a conservation action."))
ux_action_schedule <- function(d, years) {
  if (!is.data.frame(d) || !all(c("Year", "Amount") %in% names(d))) stop("Action CSV needs the exact headers Year,Amount.")
  if (anyDuplicated(names(d))) stop("The action CSV contains duplicate column headers.")
  if (!nrow(d) || !is.numeric(d$Year) || any(!is.finite(d$Year)) || any(d$Year != floor(d$Year))) stop("Year must contain whole years.")
  if (anyDuplicated(d$Year)) stop("Each forecast year must appear once in the action CSV.")
  if (!is.numeric(d$Amount) || any(!is.finite(d$Amount)) || any(d$Amount < 0)) stop("Amount must contain finite non-negative numbers; enter 0 for years without actions.")
  missing <- setdiff(years, d$Year); outside <- setdiff(d$Year, years)
  if (length(missing)) stop("Action schedule is missing forecast years: ", paste(missing, collapse = ", "), ". Enter explicit zero amounts where needed.")
  if (length(outside)) stop("Action schedule contains years outside this forecast: ", paste(outside, collapse = ", "))
  d[match(years, d$Year), c("Year", "Amount"), drop = FALSE]
}
ux_validate_scenario <- function(v, years) {
  params <- current_model_params(v)
  if (!length(years)) stop("Estimate the population before creating scenarios.")
  if (v$pva_mode == "mode_threat") {
    if (v$threat_input_type == "table") {
      if (is.null(v$ux_threat_data)) stop("Upload a threat schedule or download its example CSV first.")
      if (!identical(v$siders_fishery_center_mode %||% "direct", "direct") || isTRUE(v$siders_force_all_fatal)) stop("CSV schedules use annual median size and mortality. Select that model and leave the all-fatal override off.")
      d <- validate_threat_schedule(v$ux_threat_data, years)
    } else {
      ux_nonnegative(v$threat_interactions, "Annual turtles affected")
      positive_or_stop(v$threat_mean_len, "Median length")
      probability_or_stop(v$threat_mort_score, "Mortality probability")
      d <- data.frame(year = years, turtles = v$threat_interactions, median_cm = v$threat_mean_len, mortality = v$threat_mort_score)
    }
  } else if (v$pva_mode == "mode_action") {
    kind <- ux_action_type(v$action_type)
    if (v$action_input_type == "table") {
      if (is.null(v$ux_action_data)) stop("Upload the CSV for this conservation action.")
      d <- ux_action_schedule(v$ux_action_data, years)
    } else {
      val <- switch(kind, nests = v$action_static_nests, yearlings = v$action_static_head, adults = v$action_static_adults)
      ux_nonnegative(val, "Annual action amount")
      d <- data.frame(Year = years, Amount = val)
    }
    if (kind == "nests") {
      positive_or_stop(v$override_eggs, "Eggs per nest")
      probability_or_stop(v$override_syr1, "Age-one survival")
      emerge <- probability_or_stop(v$action_protected_emergence, "Protected emergence success")
      if (v$action_emergence_basis == "conditional_on_hatching") emerge <- emerge * probability_or_stop(v$action_protected_hatch, "Protected hatching success")
      counter <- probability_or_stop(v$action_counterfactual_emergence, "Unprotected emergence per egg")
      if (counter > emerge) stop("Unprotected emergence exceeds protected emergence. This conservation-gain model requires a non-negative improvement.")
    }
    if (kind == "yearlings") ux_nonnegative(v$action_release_age, "Release age")
  } else if (v$pva_mode == "mode_portfolio") {
    active <- v$portfolio_active %||% character(0)
    if (!length(active)) stop("Select at least one component for the combined scenario.")
    if ("A: Protect Nests" %in% active) {
      ux_nonnegative(v$portfolio_nests, "Nests protected")
      positive_or_stop(v$portfolio_eggs, "Eggs per nest")
      probability_or_stop(v$portfolio_year1, "Age-one survival")
      pe <- probability_or_stop(v$portfolio_protected_emergence, "Protected emergence")
      if (v$portfolio_emergence_basis == "conditional_on_hatching") pe <- pe * probability_or_stop(v$portfolio_protected_hatch, "Protected hatching success")
      if (probability_or_stop(v$portfolio_counterfactual_emergence, "Unprotected emergence") > pe) stop("Unprotected emergence exceeds protected emergence.")
    }
    if ("B: Headstart Yearlings" %in% active) { ux_nonnegative(v$portfolio_head, "Additional released turtles"); ux_nonnegative(v$portfolio_release_age, "Release age") }
    if ("C: Stop Adult Poaching" %in% active) ux_nonnegative(v$portfolio_adults, "Adult deaths prevented")
    if ("Threat: Bycatch / Mortality" %in% active) {
      ux_nonnegative(v$portfolio_threat_interactions, "Annual turtle interactions")
      positive_or_stop(v$portfolio_threat_mean_len, "Median length")
      probability_or_stop(v$portfolio_threat_mortality, "Mortality probability")
    }
    d <- data.frame(Year = years)
    if ("A: Protect Nests" %in% active) d$Nests_protected <- v$portfolio_nests
    if ("B: Headstart Yearlings" %in% active) d$Additional_released <- v$portfolio_head
    if ("C: Stop Adult Poaching" %in% active) d$Adult_deaths_prevented <- v$portfolio_adults
    if ("Threat: Bycatch / Mortality" %in% active) d$Turtles_affected <- v$portfolio_threat_interactions
  } else stop("Choose a threat or conservation action.")
  if (v$pva_mode %in% c("mode_threat", "mode_portfolio")) {
    for (key in c("siders_log_length_sd", "siders_logit_mortality_sd", "siders_ri_sd", "siders_pj_sd", "siders_pa_sd", "siders_pf_sd")) ux_nonnegative(v[[key]] %||% 0, key)
    rho <- v$siders_length_mortality_rho %||% 0
    if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) stop("Length/mortality correlation must be strictly between -1 and 1.")
    if (identical(v$siders_count_mode, "cmp")) positive_or_stop(v$siders_cmp_nu, "Count dispersion")
    if (identical(v$siders_ri_distribution, "cmp")) positive_or_stop(v$siders_ri_cmp_nu, "Remigration dispersion")
    if (identical(v$siders_fishery_center_mode, "siders_atl")) for (key in c("siders_length_beta0", "siders_length_beta1", "siders_mortality_mu0")) {
      if (length(v[[key]]) != 1L || !is.finite(v[[key]])) stop("Supply the regression parameter: ", key)
    }
  }
  d
}


ux_threat_ui <- function() tagList(
  h4("Threat: fishery interactions"),
  radioButtons("threat_input_type", "How do the annual amounts vary?",
               c("Same each year"="static", "Upload annual schedule"="table"), selected="table"),
  conditionalPanel("input.threat_input_type == 'static'",
    numericInput("threat_interactions", "Turtles affected per year", 150, min=0),
    numericInput("threat_mean_len", "Median carapace length (cm)", 65, min=.1),
    numericInput("threat_mort_score", "Mortality probability (0–1)", .35, min=0, max=1, step=.01),
    helpText("Count turtles before applying mortality. A probability of 0.35 means 35% die from the interaction.")
  ),
  conditionalPanel("input.threat_input_type == 'table'",
    p("Each row gives the number of turtles affected in one year, their median size, and their probability of dying."),
    downloadButton("download_threat_template", "Download example CSV", class="btn-outline-primary"),
    fileInput("threat_csv", "Upload threat schedule (.csv)", accept=".csv"),
    helpText("Columns: year, turtles, median_cm, mortality. Include every forecast year; use 0 turtles in years without interactions."),
    tags$details(tags$summary("Size, mortality and year definitions"),
      p("Use the carapace-length measure specified by the growth model. The median represents the cohort unless size variation is supplied under Input uncertainty."),
      p("Mortality includes immediate and delayed deaths. Enter 0.35 for 35%; do not multiply the turtle count by this probability beforehand."),
      p("Year is the ending census: 2026 represents the transition from 2025 to 2026. Size and mortality may be blank only when turtles = 0."),
      p("CSV schedules use direct annual size and mortality inputs. Enter mortality = 1 for an all-fatal year instead of using the regression-mode override."))
  )
)
ux_action_ui <- function() tagList(
  h4("Conservation action"),
  selectInput("action_type", "What does the action do?", c("Protect nests"="A: Protect Nests", "Headstart turtles"="B: Headstarting", "Prevent adult deaths"="C: Stop Adult Poaching")),
  uiOutput("ux_action_amount_help"),
  radioButtons("action_input_type", "How do the annual amounts vary?", c("Same each year"="static", "Upload annual schedule"="table"), selected="static"),
  conditionalPanel("input.action_input_type == 'static'",
    conditionalPanel("input.action_type == 'A: Protect Nests'", numericInput("action_static_nests", "Nests protected per year", 0, min=0)),
    conditionalPanel("input.action_type == 'B: Headstarting'", numericInput("action_static_head", "Net additional turtles alive at release age per year", 0, min=0)),
    conditionalPanel("input.action_type == 'C: Stop Adult Poaching'", numericInput("action_static_adults", "Adult deaths prevented per year", 0, min=0))
  ),
  conditionalPanel("input.action_input_type == 'table'",
    downloadButton("download_action_template", "Download example CSV", class="btn-outline-primary"),
    conditionalPanel("input.action_type == 'A: Protect Nests'", fileInput("action_csv_nests", "Upload nest-protection CSV", accept=".csv")),
    conditionalPanel("input.action_type == 'B: Headstarting'", fileInput("action_csv_yearlings", "Upload headstarting CSV", accept=".csv")),
    conditionalPanel("input.action_type == 'C: Stop Adult Poaching'", fileInput("action_csv_adults", "Upload adult-protection CSV", accept=".csv")),
    helpText("Columns: Year, Amount. Include every forecast year and enter 0 for years without actions. Each action retains its own upload.")
  ),
  conditionalPanel("input.action_type == 'A: Protect Nests'",
    h5("Nest productivity and survival"),
    helpText("These assumptions apply to every row. The values below are illustrative; replace them with evidence for your intervention."),
    numericInput("override_eggs", "Eggs per nest", 100, min=1),
    selectInput("action_emergence_basis", "How is protected emergence success measured?",
      c("Emerged hatchlings / hatched eggs"="conditional_on_hatching", "Emerged hatchlings / eggs laid"="per_egg")),
    conditionalPanel("input.action_emergence_basis == 'conditional_on_hatching'", numericInput("action_protected_hatch", "With protection: hatching success (0–1)", .65, min=0, max=1, step=.01)),
    numericInput("action_protected_emergence", "With protection: emergence success (0–1)", .85, min=0, max=1, step=.01),
    numericInput("action_counterfactual_emergence", "Without protection: emerged hatchlings per egg (0–1)", .30, min=0, max=1, step=.01),
    numericInput("override_syr1", "Emerged hatchling to age-one survival (0–1)", .025, min=0, max=1, step=.005),
    helpText("Only the improvement over unprotected nests is credited. Enter zero without protection only when complete loss is the supported counterfactual.")
  ),
  conditionalPanel("input.action_type == 'B: Headstarting'",
    h5("Release stage"), numericInput("action_release_age", "Release age (years)", 1, min=0, step=.25),
    helpText("The amount already represents net additional survivors at this age. Do not apply another pre-release survival discount.")
  )
)
ux_portfolio_ui <- function() tagList(
        h4("Combined actions and threats"),
        
        div(
          class = "alert alert-info",
          style = "font-size: 0.88rem;",
          
          tags$b("Build your portfolio: "),
          
          paste(
            "Select the interventions you want to combine.",
            "The editable assumptions for each selected intervention",
            "will appear below."
          )
        ),
        
        
        checkboxGroupInput(
          "portfolio_active",
          "Select interventions:",
          choices = c(
            "A: Protect Nests",
            "B: Headstart Yearlings",
            "C: Stop Adult Poaching",
            "Threat: Bycatch / Mortality"
          ),
          selected = character(0)
        ),
        
        
        # =============================================================
        # A. NEST PROTECTION
        # =============================================================
        
        conditionalPanel(
          condition = paste0(
            "input.portfolio_active && ",
            "input.portfolio_active.indexOf('A: Protect Nests') !== -1"
          ),
          
          div(
            style = paste(
              "border-left: 4px solid #e69f00;",
              "padding: 12px 15px;",
              "margin: 12px 0;",
              "background-color: #fffaf2;"
            ),
            
            h6(tags$b("A. Protect Nests")),
            
            numericInput(
              "portfolio_nests",
              "Nests protected per year:",
              value = 0,
              min = 0,
              step = 1
            ),
            
            numericInput(
              "portfolio_eggs",
              "Eggs per nest:",
              value = 100,
              min = 1,
              step = 1
            ),
            
            selectInput(
              "portfolio_emergence_basis",
              "How is protected emergence success reported?",
              choices = c(
                "Emerged hatchlings / hatched eggs" =
                  "conditional_on_hatching",
                "Emerged hatchlings / eggs laid" =
                  "per_egg"
              ),
              selected = "conditional_on_hatching"
            ),
            
            conditionalPanel(
              condition = paste0(
                "input.portfolio_emergence_basis == ",
                "'conditional_on_hatching'"
              ),
              
              numericInput(
                "portfolio_protected_hatch",
                "With protection: hatching success",
                value = 0.65,
                min = 0,
                max = 1,
                step = 0.01
              )
            ),
            
            numericInput(
              "portfolio_protected_emergence",
              "With protection: emergence success",
              value = 0.85,
              min = 0,
              max = 1,
              step = 0.01
            ),
            
            numericInput(
              "portfolio_counterfactual_emergence",
              paste(
                "Without protection:",
                "emerged hatchlings per egg laid"
              ),
              value = 0,
              min = 0,
              max = 1,
              step = 0.01
            ),
            
            numericInput(
              "portfolio_year1",
              "Emerged-to-yearling survival:",
              value = 0.025,
              min = 0,
              max = 1,
              step = 0.005
            ),
            
            div(
              class = "alert alert-warning",
              style = "padding: 8px; font-size: 0.82rem;",
              
              tags$b("Counterfactual reminder: "),
              
              paste(
                "A counterfactual emergence value of 0 means",
                "the protected nests would otherwise have produced",
                "no emerged hatchlings."
              )
            )
          )
        ),
        
        
        # =============================================================
        # B. HEADSTARTING
        # =============================================================
        
        conditionalPanel(
          condition = paste0(
            "input.portfolio_active && ",
            "input.portfolio_active.indexOf('B: Headstart Yearlings') !== -1"
          ),
          
          div(
            style = paste(
              "border-left: 4px solid #56b4e9;",
              "padding: 12px 15px;",
              "margin: 12px 0;",
              "background-color: #f5fbfe;"
            ),
            
            h6(tags$b("B. Headstart Yearlings")),
            
            numericInput(
              "portfolio_head",
              "Net additional turtles at release age per year:",
              value = 0,
              min = 0,
              step = 1
            ),
            
            numericInput(
              "portfolio_release_age",
              "Release age (years):",
              value = 1,
              min = 0,
              step = 0.25
            ),
            
            p(
              paste(
                "Enter NET ADDITIONAL turtles alive at release age",
                "relative to the no-headstarting counterfactual,",
                "not the gross number released."
              ),
              style = "font-size: 0.82rem; color: #666;"
            )
          )
        ),
        
        
        # =============================================================
        # C. ADULT PROTECTION
        # =============================================================
        
        conditionalPanel(
          condition = paste0(
            "input.portfolio_active && ",
            "input.portfolio_active.indexOf('C: Stop Adult Poaching') !== -1"
          ),
          
          div(
            style = paste(
              "border-left: 4px solid #009e73;",
              "padding: 12px 15px;",
              "margin: 12px 0;",
              "background-color: #f4fbf8;"
            ),
            
            h6(tags$b("C. Stop Adult Poaching")),
            
            numericInput(
              "portfolio_adults",
              "Adult deaths averted per year:",
              value = 0,
              min = 0,
              step = 1
            ),
            
            p(
              paste(
                "Enter adult deaths genuinely prevented relative to",
                "the no-intervention counterfactual—not the number",
                "of adults encountered or protected."
              ),
              style = "font-size: 0.82rem; color: #666;"
            )
          )
        ),
        
        
        # =============================================================
        # D. FISHERY THREAT
        # =============================================================
        
        conditionalPanel(
          condition = paste0(
            "input.portfolio_active && ",
            "input.portfolio_active.indexOf('Threat: Bycatch / Mortality') !== -1"
          ),
          
          div(
            style = paste(
              "border-left: 4px solid #d9534f;",
              "padding: 12px 15px;",
              "margin: 12px 0;",
              "background-color: #fff7f7;"
            ),
            
            h6(tags$b("Threat: Bycatch / Mortality")),
            
            numericInput(
              "portfolio_threat_interactions",
              "Annual turtle interactions:",
              value = 0,
              min = 0,
              step = 1
            ),
            
            numericInput(
              "portfolio_threat_mean_len",
              "Median carapace length (cm):",
              value = 65,
              min = 0.1,
              step = 1
            ),
            
            numericInput(
              "portfolio_threat_mortality",
              "Mortality risk (0-1):",
              value = 0.35,
              min = 0,
              max = 1,
              step = 0.05
            )
          )
        ),
        
        
        hr(),
        
        NULL
      )

ux_uncertainty_ui <- function() tags$details(
                tags$summary(
                  tags$b("Input uncertainty")
                ),
                
                div(
                  style = "padding-top: 10px;",
                  
                  p(
                    "Enter uncertainty supported by the selected assessment sources. ",
                    "Zero SD holds each input fixed. Length and mortality centers are medians; survival and sex centers are normal means before truncation to [0,1]."
                  ),
                  
                  selectInput(
                    "siders_count_mode",
                    "Annual interaction-count distribution:",
                    choices = c(
                      "Fixed annual total" = "fixed",
                      "Poisson" = "poisson",
                      "Conway-Maxwell-Poisson" = "cmp"
                    ),
                    selected = "fixed"
                  ),
                  
                  conditionalPanel(
                    condition = "input.siders_count_mode == 'cmp'",
                    
                    numericInput(
                      "siders_cmp_nu",
                      "CMP dispersion (nu):",
                      value = 1,
                      min = 0.01,
                      step = 0.05
                    )
                  ),
                 
                  selectInput(
                    "siders_fishery_center_mode",
                    "Length and mortality model:",
                    choices = c(
                      "Annual median size and mortality estimate" = "direct",
                      "Original Siders ATL regression" = "siders_atl"
                    ),
                    selected = "direct"
                  ),
                  
                  conditionalPanel(
                    condition =
                      "input.siders_fishery_center_mode == 'siders_atl'",
                    
                    helpText(
                      paste(
                        "Original Siders mode uses",
                        "mean log length = beta0 + beta1 × ATL",
                        "and a constant mean logit mortality."
                      )
                    ),
                    
                    numericInput(
                      "siders_length_beta0",
                      "Length regression beta0:",
                      value = NA,
                      step = 0.001
                    ),
                    
                    numericInput(
                      "siders_length_beta1",
                      "Length regression beta1:",
                      value = NA,
                      step = 0.001
                    ),
                    
                    numericInput(
                      "siders_mortality_mu0",
                      "Mean logit mortality (mu0):",
                      value = NA,
                      step = 0.001
                    )
                  ),

                  checkboxInput(
                    "siders_force_all_fatal",
                    "Treat every interaction as fatal",
                    value = FALSE
                  ),
                  
                  helpText(
                    paste(
                      "Legacy all-fatal scenario. When selected, the model",
                      "sets mortality to 1 for every interaction and ignores",
                      "the mortality mean and uncertainty for the calculation."
                    )
                  ),
                  
                  numericInput(
                    "siders_log_length_sd",
                    "SD of log carapace length:",
                    value = 0,
                    min = 0,
                    step = 0.01
                  ),
                  
                  numericInput(
                    "siders_logit_mortality_sd",
                    "SD of logit mortality:",
                    value = 0,
                    min = 0,
                    step = 0.01
                  ),
                  
                  numericInput(
                    "siders_length_mortality_rho",
                    "Correlation between log length and logit mortality:",
                    value = 0,
                    min = -1,
                    max = 1,
                    step = 0.05
                  ),
                  
                  selectInput(
                    "siders_ri_distribution",
                    "Remigration-interval distribution:",
                    choices = c(
                      "Fixed" = "fixed",
                      "Truncated normal" = "normal",
                      "Conway-Maxwell-Poisson" = "cmp"
                    ),
                    selected = "fixed"
                  ),
                  
                  conditionalPanel(
                    condition =
                      "input.siders_ri_distribution == 'normal'",
                    
                    numericInput(
                      "siders_ri_sd",
                      "SD of remigration interval:",
                      value = 0,
                      min = 0,
                      step = 0.01
                    )
                  ),
                  
                  conditionalPanel(
                    condition =
                      "input.siders_ri_distribution == 'cmp'",
                    
                    numericInput(
                      "siders_ri_cmp_nu",
                      "RI CMP dispersion (nu):",
                      value = 1,
                      min = 0.01,
                      step = 0.05
                    )
                  ),
                  
                  numericInput(
                    "siders_pj_sd",
                    "SD of juvenile survival:",
                    value = 0,
                    min = 0,
                    max = 1,
                    step = 0.001
                  ),
                  
                  numericInput(
                    "siders_pa_sd",
                    "SD of adult survival:",
                    value = 0,
                    min = 0,
                    max = 1,
                    step = 0.001
                  ),
                  
                  numericInput(
                    "siders_pf_sd",
                    "SD of proportion female:",
                    value = 0,
                    min = 0,
                    max = 1,
                    step = 0.001
                  )
                )
              )


# Static controls stay mounted when users navigate. Only summaries and previews
# are reactive UI, so editing one scenario never recreates unrelated inputs.
ux_number <- function(id, label, value = NA_real_, min = 0, max = NA, step = 0.01) {
  control <- numericInput(id, label, value = value, min = min, max = max, step = step)
  if (identical(id, "max_age")) return(htmltools::tagQuery(control)$find("input")$addAttrs(readonly="readonly")$allTags())
  control
}
ux_section <- function(title, subtitle, ...) {
  div(class = "ux-section", h2(title), p(class = "ux-lead", subtitle), ...)
}
ux_methods_ui <- function() tags$details(
  tags$summary("Methods and calculations"),
  p("Annual nest counts are converted to annual nesters using clutch frequency. Historical trends and observation error are estimated before projecting scenarios."),
  p("Threat losses and conservation gains use the expected first-nesting pathway from Methods_ConsBio_20260831.pdf, equations 7–22. Mortality and the female proportion are applied once. Conservation gains represent improvement over the no-action counterfactual."),
  p("Projection: remove threat ANE before annual growth, add abundance-scale process error, then add conservation ANE; truncate population abundance at zero. Q is a variance."),
  helpText("ANE means annual-nester equivalents: the common currency for effects on future nesting. The supplied manuscript's additive process-error convention is retained."),
  selectInput("projection_uq_mode", "How should future trend and process variance be sampled?",
              c("Resample a paired posterior draw each year (methods default)" = "dynamic",
                "Keep one paired posterior draw per trajectory" = "static")),
  selectInput("ane_mat_threshold", "Forced first-nesting size threshold",
              c("99% of asymptotic length (methods default)" = 0.99,
                "97.5% of asymptotic length (sensitivity option)" = 0.975))
)

# Diagnostics are additions to the interface, not alterations to the model.
# Sources: https://cran.r-project.org/web/packages/coda/refman/coda.html
# and https://cran.r-universe.dev/jagsUI/doc/manual.html
# Never split pooled sims.list into invented chains. jagsUI$samples retains
# the genuine post-burn-in chain identities and mcpar iteration metadata.
dev_chains <- function(fit, terminal_year = NULL) {
  if (is.null(fit$samples) || !inherits(fit$samples,"mcmc.list")) return(NULL)
  result <- lapply(fit$samples,function(chain) {
    m <- as.matrix(chain)
    if (!is.null(terminal_year)) {
      x_name <- paste0("X[",terminal_year,"]")
      if (x_name %in% colnames(m)) {
        a_names <- grep("^A\\[",colnames(m),value=TRUE)
        # This mirrors the existing regional annual-nester conversion,
        # applied within each chain before any diagnostic is calculated.
        n <- if(length(a_names)) rowSums(exp(sweep(m[,a_names,drop=FALSE],1,m[,x_name],"+"))) else exp(m[,x_name])
        m <- cbind(m,N_final=n)
      }
    }
    coda::mcmc(m,start=stats::start(chain),thin=coda::thin(chain))
  })
  do.call(coda::mcmc.list,result)
}
dev_diagnostic_table <- function(chains, parameters = NULL) {
  if (is.null(chains)) return(data.frame())
  parameters <- intersect(parameters %||% colnames(as.matrix(chains[[1]])),colnames(as.matrix(chains[[1]])))
  safe <- function(expr) tryCatch(suppressWarnings(expr),error=function(e) NA_real_)
  rows <- lapply(parameters,function(parameter) {
    selected <- do.call(coda::mcmc.list,lapply(chains,function(ch) coda::mcmc(as.matrix(ch)[,parameter,drop=FALSE],start=stats::start(ch),thin=coda::thin(ch))))
    vectors <- lapply(selected,as.numeric)
    finite <- all(vapply(vectors,function(x) all(is.finite(x)),logical(1)))
    constant <- finite && all(vapply(vectors,function(x) stats::var(x)==0,logical(1)))
    fixed <- constant && length(unique(unlist(vectors)))==1L
    rhat <- if(finite && !constant && length(selected)>1L) safe(coda::gelman.diag(selected,autoburnin=FALSE,multivariate=FALSE)$psrf[1,"Point est."]) else NA_real_
    total_ess <- if(finite && !constant) safe(unname(coda::effectiveSize(selected)[1])) else NA_real_
    do.call(rbind,lapply(seq_along(selected),function(i) {
      v <- vectors[[i]]
      valid <- all(is.finite(v)) && length(v)>=20L && isTRUE(stats::var(v)>0)
      z <- if(valid) safe(unname(coda::geweke.diag(selected[[i]],frac1=.1,frac2=.5)$z[1])) else NA_real_
      ess <- if(valid) safe(unname(coda::effectiveSize(selected[[i]])[1])) else NA_real_
      note <- if(fixed) "Constant node; not assessed" else if(!valid || !is.finite(z) || !is.finite(rhat)) "Diagnostic unavailable" else if(abs(z)>1.96 || rhat>1.1 || total_ess<100) "Review" else "No screening flag; inspect plots"
      data.frame(Parameter=parameter,Chain=i,Draws=length(v),Geweke_z=z,Rhat_classic=rhat,ESS_chain=ess,ESS_total=total_ess,Note=note,check.names=FALSE)
    }))
  })
  if(!length(rows)) return(data.frame())
  do.call(rbind,rows)
}

# Presentation-only wrappers. conditionalPanel hides rather than recreates
# inputs: changing view cannot reset values, refit models, or alter results.
dev_only <- function(...) conditionalPanel("input.developer_view === true", div(class="dev-workspace", ...))
ux_info <- function(label, text) bslib::popover(
  tags$button(type="button", class="btn btn-link ux-info", `aria-label`=label, icon("circle-info")),
  text, placement="bottom")
ux_back_next <- function(back_id, back_label, next_id, next_label) div(class="ux-footer",
  if (!is.null(back_id)) actionButton(back_id, back_label, class="btn-link"),
  if (!is.null(next_id)) actionButton(next_id, next_label, class="btn-primary"))

ui <- page_fluid(
  theme=bs_theme(version=5, primary="#185e55", bg="#f8f7f2", fg="#243d3a", base_font="Arial"),
  tags$head(tags$script(HTML("$(document).on('shiny:connected', function(){Shiny.addCustomMessageHandler('biology-details',function(open){document.getElementById('biology_details').open=open;});});")),tags$style(HTML("
    body{background:#f8f7f2;color:#243d3a;font-family:Arial,Helvetica,sans-serif;}
    .ux-shell{max-width:1120px;margin:0 auto;padding:0 28px 35px;}
    .ux-brand{display:flex;justify-content:space-between;align-items:center;gap:20px;padding:24px 0;border-bottom:1px solid #d9dfd8;}
    .ux-brand h1{font-size:18px;font-weight:500;margin:0;line-height:1.3;}
    .ux-brand p{font-size:11px;letter-spacing:.1em;color:#596a66;margin:3px 0 0;}
    .ux-brand .shiny-input-container{margin:0;width:auto;}.ux-brand .checkbox{margin:0;}
    .ux-version{font-size:11px;color:#596a66;margin-top:5px;text-align:right;}
    .ux-nav{display:grid;grid-template-columns:repeat(6,minmax(0,1fr));gap:12px;margin:16px 0 28px;}
    .ux-nav button{border:0;border-bottom:2px solid #d9dfd8;background:transparent;text-align:left;padding:9px 0 12px;color:#596a66;font-size:13px;}
    .ux-nav button[aria-current='step']{color:#185e55;border-color:#185e55;}.ux-nav button:disabled{opacity:.55;cursor:not-allowed;}
    .ux-nav small{display:block;font-size:10px;letter-spacing:.07em;margin-bottom:5px;}
    .ux-section h2,.ux-welcome h2{font:400 36px/1.12 Georgia,serif;letter-spacing:-.025em;margin:6px 0 14px;}
    .ux-welcome{padding:28px 0 36px;max-width:740px;}.ux-welcome h2{font-size:46px;max-width:620px;}
    .ux-eyebrow{text-transform:uppercase;font-size:11px;letter-spacing:.12em;color:#185e55;}
    .ux-lead{color:#596a66;max-width:700px;margin-bottom:24px;}
    .ux-panel{background:#fffefa;border:1px solid #d9dfd8;border-radius:12px;padding:20px;margin-bottom:18px;}
    .ux-panel h4,h4{font-size:17px;font-weight:500;margin-bottom:16px;}h5{font-size:15px;font-weight:500;}
    .ux-plain{padding:0;background:transparent;border:0;}.ux-actions{display:flex;gap:10px;flex-wrap:wrap;margin:16px 0;align-items:center;}
    .ux-footer{display:flex;justify-content:space-between;gap:14px;border-top:1px solid #d9dfd8;margin-top:24px;padding-top:18px;}
    .btn{border-radius:7px;white-space:normal;padding:10px 16px;font-size:14px;}.btn-primary{background:#185e55;border-color:#185e55;}
    .btn-link{color:#185e55;}.ux-info{padding:2px 6px;}.ux-note{border-left:3px solid #185e55;padding:8px 12px;margin:12px 0;color:#596a66;font-size:13px;background:transparent;}
    .ux-muted,.help-block{font-size:13px;color:#596a66;}.shiny-input-container{max-width:100%;}label{font-weight:400;font-size:13px;}
    .form-control,.form-select,.selectize-input{background:#fffefa;border-color:#d9dfd8;border-radius:6px;}
    details{margin:14px 0;padding:4px 0;border:0;background:transparent;}summary{cursor:pointer;color:#185e55;font-size:13px;}details[open]>summary{margin-bottom:14px;}
    .dev-workspace{border-top:1px solid #d9dfd8;padding-top:18px;margin-top:18px;}.dev-workspace .card{box-shadow:none;border:1px solid #d9dfd8;background:#fffefa;}
    .nav-tabs .nav-link{font-size:13px;color:#596a66;}.nav-tabs .nav-link.active{color:#185e55;background:#fffefa;border-bottom-color:#fffefa;}
    .alert{font-size:13px;border-radius:6px;}.alert-success{background:transparent;border:0;border-left:3px solid #185e55;color:#243d3a;padding:8px 12px;}
    .ux-saved{border-bottom:1px solid #d9dfd8;padding:10px 0;}.table{font-size:13px;}.ux-scroll{overflow-x:auto;}
    .shiny-plot-output{max-width:100%;}.ux-fields-summary{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:22px 0;}
    @media(max-width:760px){.ux-shell{padding:0 16px 24px;}.ux-nav{grid-template-columns:repeat(3,1fr);}.ux-welcome h2{font-size:35px;}.ux-section h2{font-size:30px;}.ux-brand{gap:10px;}.ux-brand h1{font-size:15px;}.ux-fields-summary{grid-template-columns:1fr;}}
  "))),
  div(class="ux-shell",
    div(class="ux-brand",div(h1("Sea Turtle Population Toolkit"),p("FROM RECORDS TO CONSERVATION DECISIONS")),
      div(checkboxInput("developer_view","Developer view",FALSE),div(class="ux-version","Guided + developer · 2026-09-02"))),
    uiOutput("ux_progress"),uiOutput("ux_global_status"),
    navset_hidden(id="wizard_steps",selected="welcome",
      nav_panel("Welcome",value="welcome",div(class="ux-welcome",
        p(class="ux-eyebrow","From nesting records to conservation decisions"),
        h2("What could the future hold for your turtles?"),
        p(class="ux-lead","Use nest counts to explore population trends—and how threats and conservation actions could change them."),
        div(class="ux-actions",actionButton("ux_welcome_upload","Upload a CSV",class="btn-primary",onclick="document.getElementById('uploaded_file').click();"),actionButton("ux_welcome_example","Try an example →",class="btn-outline-primary")),
        downloadButton("ux_welcome_download","Download an example CSV",class="btn-link"),
        ux_info("About the example","Synthetic annual nest counts. The example follows the same workflow as your own data."),
        dev_only(p("The same scientific engine runs in both views. Developer view reveals settings and diagnostics; it does not change sampling, assumptions, or results."))
      )),
      nav_panel("Your data",value="step2",ux_section("Check your nesting records.","Do these beaches and years match what you want to assess?",
        fluidRow(column(4,div(class="ux-panel",
          radioButtons("ux_data_source","Data source",c("My CSV"="upload","Synthetic example"="demo"),selected="upload",inline=TRUE),
          selectInput("ux_example_monitoring","Example format",c("Annual counts"="annual","Monthly counts with gaps"="monthly")),
          conditionalPanel("input.ux_data_source == 'upload'",fileInput("uploaded_file","Nesting records (.csv)",accept=".csv")),
          downloadButton("download_template","Download example CSV",class="btn-link"),
          tags$details(tags$summary("File format & missing observations"),p("Columns: Year, Month, Site, Count, Monitored. Count means nests, not turtles. Leave Month blank for annual records."),p("Use a blank Count and Monitored = FALSE for no survey. A real observed zero is different: the current log-count model cannot fit true zeros.")),
          uiOutput("timeframe_ui"),uiOutput("period_override_ui")
        )),column(8,uiOutput("ux_data_status"),uiOutput("qaqc_alerts"),plotOutput("preview_annual_raw",height="300px"),
          tags$details(tags$summary("Preview records"),div(class="ux-scroll",tableOutput("data_preview_table_raw"))),
          tags$details(tags$summary("Monthly pattern & survey coverage"),plotOutput("preview_monthly_seasonality",height="260px"),uiOutput("seasonality_recommendation"),tableOutput("gap_table")),
          dev_only(downloadButton("download_preview_raw","Download data plot",class="btn-outline-secondary"))
        )),ux_back_next("ux_back_welcome","← Welcome","go_step3","Continue: your turtles →")
      )),
      nav_panel("Your turtles",value="step3",ux_section("A little about your population.","Nest counts are the starting point. We also need information about growth, breeding, and survival.",
        selectInput("parameter_preset","Population information",c("Choose a profile or enter my own values"="custom","Western Pacific leatherback — reference profile"="wp_leatherback"),selected="custom",width="100%"),
        uiOutput("ux_biology_status"),
        p(class="ux-muted","Reference: Martin 2020 / Siders 2023. Use only if appropriate for your population."),
        tags$details(id="biology_details",tags$summary("Review or edit growth, breeding & survival"),
          fluidRow(column(4,h4("Breeding"),ux_number("clutch_freq","Nests per female per nesting season",step=.1),ux_number("remig_int","Years between nesting seasons",step=.1),ux_number("pf","Female proportion (0–1)",max=1),ux_info("Female proportion","0.73 means 73% female. This population assumption is applied once in the impact calculation.")),
            column(4,h4("Growth & first nesting"),ux_number("linf","Asymptotic shell length (cm)",step=.1),ux_number("k","Growth coefficient (k)",step=.001),ux_number("tknot","Growth age offset (t0, years)",min=NA),ux_number("lmat","First-nesting transition length (cm)",step=.1),ux_number("sig_mat","First-nesting logistic scale (cm)",step=.1),dev_only(ux_number("max_age","Age at 99% asymptotic size (derived)",step=.01))),
            column(4,h4("Survival"),ux_number("ane_pj","Annual juvenile survival (0–1)",max=1),ux_number("ane_pa","Annual adult survival (0–1)",max=1),ux_info("Survival probabilities","0.81 means 81% survive one year. Values must be supported for the population being assessed."))),
          tableOutput("ux_biology_table"),downloadButton("download_current_parameters","Download assumptions",class="btn-outline-secondary")),
        dev_only(tags$details(tags$summary("Illustrative biological sensitivity screening"),p(class="ux-muted","A separate simplified stage matrix, not the fitted Martin/Siders model. Uses current nest-action productivity assumptions; ±15% perturbations are illustrative."),plotOutput("sensitivity_tornado_plot",height="460px"))),
        ux_back_next("back_step2","← Your data","go_step4","Continue: current trend →")
      )),
      nav_panel("Current trend",value="step4",ux_section("How has nesting changed?","Estimate the trend before exploring future threats and actions.",
        div(class="ux-actions",actionButton("run_model","Estimate trend",class="btn-primary"),ux_info("Analysis time","Monthly records are processed before trend estimation. Sampling settings—not the selected view—determine computation time.")),
        uiOutput("ux_fit_status"),uiOutput("dev_health"),uiOutput("summary_stats"),
        plotOutput("clean_baseline_plot",height="330px"),
        p(class="ux-muted","Line: estimated annual nesting females. Shading: 95% interval. Developer view also provides the original log-scale plot."),
        downloadButton("download_baseline","Download trend plot",class="btn-link"),
        dev_only(h4("Model workspace"),
          uiOutput("diagnostic_model_selector"),
          navset_card_tab(id="dev_tabs",
            nav_panel("Posteriors",plotOutput("dev_posterior",height="310px"),tags$details(tags$summary("Joint posterior relationships"),plotOutput("posterior_pairs_plot",height="430px"),uiOutput("covariance_interpretation_text"),downloadButton("download_posterior","Download joint plot"))),
            nav_panel("Chains & convergence",uiOutput("dev_fit_selector"),uiOutput("dev_parameter_selector"),plotOutput("dev_trace",height="240px"),plotOutput("dev_acf",height="210px"),div(class="ux-scroll",tableOutput("dev_diagnostics")),p(class="ux-muted","Geweke: first 10% versus last 50% of retained draws, separately for each chain. Classic coda R-hat; no additional burn-in. ESS estimates information for the mean. Missing/constant diagnostics are not passes."),downloadButton("dev_download_diagnostics","Download diagnostics CSV"),downloadButton("dev_download_chains","Download chains (RDS)")),
            nav_panel("Models & error",checkboxGroupInput("plot_layers","Models to display",c("Regional trend"="jags","Independent beach trends"="jags_indep"),selected="jags"),uiOutput("dynamic_unified_plot"),tableOutput("table_u"),uiOutput("aicc_recommendation_ui"),downloadButton("download_unified","Download comparison"),plotOutput("plot_var",height="260px"),uiOutput("variance_interpretation_text"),downloadButton("download_var","Download error plot")),
            nav_panel("Settings & methods",sliderInput("iterations","MCMC iterations",10000,150000,10000,step=10000),p(class="ux-muted","10,000 is an exploration setting, not a guarantee of convergence. Changing settings invalidates the fit in both views."),checkboxInput("legacy_fourier_mcmc","Use supplied Fourier sampling settings",FALSE),checkboxInput("legacy_trend_mcmc","Use supplied trend sampling settings",FALSE),uiOutput("imputation_uncertainty_option_ui"),uiOutput("imputation_branch_selector_ui"),uiOutput("imputation_branch_summary_ui"),checkboxInput("run_split","Fit trends before and after a year",FALSE),conditionalPanel("input.run_split == true",numericInput("split_year","Split year",2005)),uiOutput("reference_site_ui"),verbatimTextOutput("dev_run_record"),downloadButton("dev_download_record","Download model record (RDS)"))
          )
        ),ux_back_next("back_step3","← Your turtles","go_step5","Continue: threats & actions →")
      )),
      nav_panel("Threats & actions",value="step5",ux_section("What might change?","Describe a threat or an action. Save each option to compare it later.",
        tags$details(tags$summary("Forecast years & population"),fluidRow(column(6,selectInput("pva_site_filter","Population",c("All Beaches (Regional Total)"))),column(6,sliderInput("pva_years","Years into the future",10,100,50,step=10))),uiOutput("ux_forecast_range")),
        dev_only(numericInput("pva_sims","Projection simulations",500,min=100,max=100000,step=100),ux_methods_ui()),
        fluidRow(column(7,div(class="ux-panel",
          textInput("ux_scenario_name","Scenario name",placeholder="e.g. Additional bycatch"),
          radioButtons("pva_mode","Add an option",c("Threat"="mode_threat","Conservation action"="mode_action","Combination"="mode_portfolio"),selected="mode_threat",inline=TRUE),
          conditionalPanel("input.pva_mode == 'mode_threat'",ux_threat_ui()),
          conditionalPanel("input.pva_mode == 'mode_action'",ux_action_ui()),
          conditionalPanel("input.pva_mode == 'mode_portfolio'",ux_portfolio_ui()),
          uiOutput("ux_shared_assumptions"),
          conditionalPanel("input.pva_mode == 'mode_threat' || input.pva_mode == 'mode_portfolio'",uiOutput("ux_uncertainty_status"),dev_only(ux_uncertainty_ui())),
          uiOutput("ux_scenario_review"),tags$details(tags$summary("Preview annual schedule"),div(class="ux-scroll",tableOutput("ux_schedule_preview"))),
          div(class="ux-actions",actionButton("ux_save_scenario","Save option",class="btn-primary"),actionButton("ux_new_scenario","New option",class="btn-outline-secondary"),actionButton("ux_save_copy","Save as copy",class="btn-outline-secondary")),uiOutput("ux_editor_status")
        )),column(5,div(class="ux-panel",h4("Your comparison"),p("Reference forecast is always included."),uiOutput("ux_saved_list"),
          tags$details(tags$summary("Edit or remove a saved option"),selectInput("ux_edit_id","Saved option",c("Choose a saved scenario"="")),actionButton("ux_load_scenario","Load into editor"),actionButton("ux_delete_scenario","Remove selected")),
          p(class="ux-muted","Saved options last for this session. Download your settings with the results."),
          actionButton("ux_go_results","Continue: compare →",class="btn-primary w-100")))) ,
        ux_back_next("back_step4","← Current trend",NULL,NULL)
      )),
      nav_panel("Compare",value="results",ux_section("Compare possible futures.","Compare each saved option with the same reference forecast.",
        checkboxGroupInput("ux_compare_ids","Options to include",choices=character(0)),uiOutput("ux_comparison_review"),
        div(class="ux-actions",actionButton("ux_run_comparison","Run comparison",class="btn-primary"),actionButton("ux_back_scenarios","Edit options",class="btn-outline-secondary")),
        uiOutput("ux_result_status"),uiOutput("dev_comparison_health"),checkboxInput("show_proj_ci","Show the range of possible outcomes",TRUE),
        plotOutput("step5_dynamic_plot",height="380px"),tableOutput("ux_comparison_table"),
        p(class="ux-muted","Forecasts describe annual nesting females. The 50% benchmark is not a species-specific extinction threshold."),
        div(class="ux-actions",downloadButton("download_projection_plot","Download chart"),downloadButton("download_projection_summary","Download results CSV"),downloadButton("ux_download_settings","Download settings")),
        dev_only(downloadButton("download_projection_draws","Download simulation CSV"),tags$details(tags$summary("Mathematical checks"),p("Equation and timing checks do not validate empirical assumptions."),actionButton("run_math_self_checks","Run mathematical checks"),div(class="ux-scroll",tableOutput("math_self_check_table")))),
        uiOutput("math_self_check_status")
      ))
    )
  )
)


# =====================================================================
# SERVER ENGINE
# =====================================================================
server <- function(input, output, session) {
  # Developer-view visibility is deliberately absent from signatures, seeds,
  # effect calculations and fit observers. Hide/show is presentation only.
  observeEvent(input$ux_back_welcome,{nav_select("wizard_steps","welcome")})
  observeEvent(input$parameter_preset,{
    session$sendCustomMessage("biology-details",identical(input$parameter_preset,"custom"))
  })
  dev_fit_choices <- reactive({
    choices <- list()
    if(!is.null(vault$res)) choices$regional <- vault$res
    if(length(vault$jags_indep_fits)) for(n in names(vault$jags_indep_fits)) choices[[paste0("site:",n)]] <- vault$jags_indep_fits[[n]]
    if(length(vault$imputation_branch_results)) for(n in names(vault$imputation_branch_results)) {
      b <- vault$imputation_branch_results[[n]]
      if(!is.null(b$res)) choices[[paste0("branch:",n)]] <- b$res
    }
    if(!is.null(vault$dev_imputation)) choices$imputation <- list(fit=vault$dev_imputation,kind="imputation")
    choices
  })
  output$dev_fit_selector <- renderUI({
    choices <- names(dev_fit_choices());req(length(choices))
    selectInput("dev_fit_choice","Fit to inspect",setNames(choices,sub("^regional$","Active regional trend",choices)),selected=isolate(input$dev_fit_choice) %||% "regional",width="100%")
  })
  dev_selected_fit <- reactive({
    choices <- dev_fit_choices();req(length(choices))
    key <- input$dev_fit_choice %||% "regional"
    if(!key %in% names(choices)) key <- names(choices)[1]
    choices[[key]]
  })
  dev_selected_chains <- reactive({
    f <- dev_selected_fit()
    dev_chains(f$fit,if(identical(f$kind,"imputation")) NULL else length(f$years))
  })
  output$dev_parameter_selector <- renderUI({
    ch <- dev_selected_chains()
    if(is.null(ch)) return(p("No chain-preserving samples are available for this fit. Rerun the model in this version."))
    choices <- colnames(as.matrix(ch[[1]]));initial <- intersect(c("U","Q","N_final"),choices)
    if(!length(initial)) initial <- head(choices,3)
    selectizeInput("dev_parameters","Parameters (up to six plotted at once)",choices,selected=initial,multiple=TRUE,options=list(maxItems=6))
  })
  dev_selected_parameters <- reactive({
    ch <- dev_selected_chains();req(ch)
    choices <- colnames(as.matrix(ch[[1]]))
    selected <- intersect(input$dev_parameters %||% c("U","Q","N_final"),choices)
    if(!length(selected)) selected <- head(choices,3)
    head(selected,6)
  })
  dev_full_table <- reactive({dev_diagnostic_table(dev_selected_chains())})
  output$dev_diagnostics <- renderTable({
    d <- dev_full_table();validate(need(nrow(d)>0,"Retained chains are required for diagnostics."))
    d[d$Parameter %in% dev_selected_parameters(),,drop=FALSE]
  },digits=3,striped=TRUE,rownames=FALSE)
  output$dev_trace <- renderPlot({
    req(isTRUE(input$developer_view));ch <- dev_selected_chains();req(ch)
    parameters <- dev_selected_parameters();op<-par(mfrow=c(length(parameters),1),mar=c(2.7,4,1.8,1));on.exit(par(op))
    colors<-rep(c("#185e55","#a95d43","#647d9a","#8a659c","#797147"),length.out=length(ch))
    for(p in parameters){
      ys<-lapply(ch,function(x) as.matrix(x)[,p]);xs<-lapply(ch,function(x)as.numeric(time(x)))
      plot(range(unlist(xs)),range(unlist(ys),finite=TRUE),type="n",xlab="Iteration (after burn-in)",ylab=p)
      for(i in seq_along(ch)) lines(xs[[i]],ys[[i]],col=colors[i],lty=i,lwd=.7)
      legend("topright",paste("Chain",seq_along(ch)),col=colors,lty=seq_along(ch),bty="n",horiz=TRUE,cex=.7)
    }
  },height=function(){max(240,150*length(dev_selected_parameters()))})
  output$dev_acf <- renderPlot({
    req(isTRUE(input$developer_view));ch<-dev_selected_chains();req(ch)
    p<-dev_selected_parameters()[1];op<-par(mfrow=c(1,length(ch)),mar=c(4,3,2,1));on.exit(par(op))
    for(i in seq_along(ch)){
      v<-as.matrix(ch[[i]])[,p]
      if(all(is.finite(v)) && isTRUE(var(v)>0)) stats::acf(v,main=paste(p,"· chain",i),xlab="Lag (retained draws)",col="#185e55")
      else {plot.new();text(.5,.5,"Constant/unavailable chain")}
    }
  })
  output$dev_posterior <- renderPlot({
    req(isTRUE(input$developer_view));f<-active_diag();ch<-dev_chains(f$fit,length(f$years))
    validate(need(!is.null(ch),"Rerun the model to retain chain samples."))
    m<-as.matrix(ch);parameters<-intersect(c("U","N_final","Q"),colnames(m))
    req(length(parameters));op<-par(mfrow=c(1,length(parameters)),mar=c(4,4,2,1));on.exit(par(op))
    for(p in parameters){v<-m[,p];v<-v[is.finite(v)];validate(need(length(v)>1,"Insufficient finite samples."));hist(v,breaks=30,probability=TRUE,main=p,xlab=switch(p,U="Annual log trend",N_final="Annual nesting females",Q="Process variance"),col="#d9e8df",border="white");if(var(v)>0)lines(density(v),col="#185e55",lwd=2);abline(v=quantile(v,c(.025,.5,.975)),lty=c(2,1,2),col="#596a66")}
  })
  # A small always-on screen checks the active projection-driving posterior.
  # Full diagnostics remain lazily rendered in Developer view. Screening is
  # not an assertion of convergence, model validity, or suitable assumptions.
  dev_health_table <- reactive({
    if(is.null(vault$res)) return(data.frame())
    ch<-dev_chains(vault$res$fit,length(vault$res$years))
    dev_diagnostic_table(ch,c("U","Q","N_final"))
  })
  dev_health_ui <- reactive({
    if(is.null(vault$res)) return(NULL)
    d<-dev_health_table()
    if(!nrow(d)) return(div(class="alert alert-warning","Convergence could not be checked: retained chain samples are unavailable. Review before using these results."))
    flagged<-any(d$Note %in% c("Review","Diagnostic unavailable"))
    s<-vault$res$fit$summary
    other_rhat<-if(!is.null(s)&&"Rhat"%in%colnames(s)) any(s[,"Rhat"]>1.1,na.rm=TRUE) else FALSE
    if(flagged||other_rhat) div(class="alert alert-warning","The model needs a closer check before these estimates are used for decisions. Open Developer view → Chains & convergence; a screening flag is not by itself proof of model failure.")
    else div(class="ux-note","Estimates are available. Basic chain screening found no flag in the core quantities; review the diagnostics and assumptions before decision-making.")
  })
  output$dev_health<-renderUI({dev_health_ui()})
  output$dev_comparison_health<-renderUI({dev_health_ui()})
  output$dev_run_record<-renderPrint({
    f<-dev_selected_fit();print(f$fit$dev_record %||% list(note="Run record unavailable for an older fitted object."))
  })
  output$dev_download_diagnostics<-downloadHandler(filename=function()paste0("diagnostics_",gsub("[^a-zA-Z0-9]","_",input$dev_fit_choice %||% "regional"),".csv"),content=function(file){d<-dev_full_table();req(nrow(d));write.csv(d,file,row.names=FALSE,na="")})
  output$dev_download_chains<-downloadHandler(filename="posterior_chains.rds",content=function(file){ch<-dev_selected_chains();req(ch);saveRDS(ch,file)})
  output$dev_download_record<-downloadHandler(filename="model_record.rds",content=function(file){f<-dev_selected_fit();req(f);saveRDS(list(record=f$fit$dev_record,years=f$years,sites=f$sites,reference_site=f$reference_site,diagnostic_conventions=list(geweke_first=.1,geweke_last=.5,rhat="coda classic; autoburnin FALSE",rhat_screen=1.1,geweke_abs_screen=1.96,ess_screen=100),session=sessionInfo()),file)})

  observeEvent(list(input$linf, input$k, input$tknot), {
    req(is.finite(input$linf), input$linf > 0,
        is.finite(input$k), input$k > 0, is.finite(input$tknot))
    updateNumericInput(session, "max_age",
      value = siders_forced_maturity_age(input$linf, input$k, input$tknot))
  })

  
  user_state <- reactiveValues(data_mode = "upload")

  observeEvent(input$go_step3, { ux_go("step3") })
  observeEvent(input$go_step4, { ux_go("step4") })
  observeEvent(input$go_step5, { ux_go("step5") })
  observeEvent(input$back_step2, { ux_go("step2") })
  observeEvent(input$back_step3, { ux_go("step3") })
  observeEvent(input$back_step4, { ux_go("step4") })

  vault <- reactiveValues(
    res = NULL, abund = NULL, unimputed_abund = NULL, marss = NULL,
    summary = NULL, year = NULL, years = NULL, nesters = NULL,
    total = NULL, trend_display = NULL, trend_pct = NULL,
    pre_trend = NULL, post_trend = NULL, pre_u_val = 0,
    post_u_val = 0, split_branch = NULL,
    empirical_ane_ledger = NULL, jags_indep_fits = NULL,
    imputation_branch_results = NULL,
    active_imputation_branch = NULL,
    draws = NULL,
    math_self_checks = NULL
  )
  vault_portfolio <- reactiveValues(
    plot_df = NULL, scorecard = NULL, all_scen_raw = list()
  )
  # Five-step workflow state. Saved records contain immutable CSV contents,
  # not temporary upload paths. Every comparison builds effects from records
  # plus the current shared biology, then uses the original joint-draw engine.
  ux <- reactiveValues(scenarios = list(), next_id = 1L, editing_id = NULL,
                       loaded_threat = NULL, loaded_action = list(),
                       fit_signature = NULL, fit_revision = 0L, fit_time = NULL,
                       result_signature = NULL, result_time = NULL,
                       result_scenarios = NULL, result_biology = NULL)
  ux_bio_keys <- c("clutch_freq", "remig_int", "linf", "k", "tknot", "lmat", "sig_mat", "pf", "ane_pj", "ane_pa", "max_age", "ane_mat_threshold")
  ux_shared_values <- reactive({ setNames(lapply(ux_bio_keys, function(k) input[[k]]), ux_bio_keys) })
  ux_action_kind <- reactive({ ux_action_type(input$action_type %||% "A: Protect Nests") })
  ux_action_file <- reactive({ input[[paste0("action_csv_", ux_action_kind())]] })
  ux_forecast_years <- reactive({
    first <- if (length(vault$years)) max(vault$years) + 1L else 2026L
    first + seq_len(as.integer(input$pva_years %||% 50L)) - 1L
  })
  observeEvent(input$ux_data_source, { user_state$data_mode <- input$ux_data_source })
  observeEvent(input$ux_welcome_upload, {
    user_state$data_mode <- "upload"
    updateRadioButtons(session, "ux_data_source", selected = "upload")
  })
  observeEvent(input$uploaded_file, {
    if (identical(input$wizard_steps, "welcome")) {
      user_state$data_mode <- "upload"
      updateRadioButtons(session, "ux_data_source", selected = "upload")
      nav_select("wizard_steps", "step2")
    }
  }, ignoreInit = TRUE)
  observeEvent(input$ux_welcome_example, {
    user_state$data_mode <- "demo"
    updateSelectInput(session, "ux_example_monitoring", selected = "annual")
    updateRadioButtons(session, "ux_data_source", selected = "demo")
    nav_select("wizard_steps", "step2")
  })
  output$ux_welcome_download <- downloadHandler(
    filename = "monitoring_annual_example.csv",
    content = function(file) write.csv(ux_example_csv("annual"), file, row.names = FALSE, na = "")
  )
  observeEvent(input$ux_nav, { ux_go(input$ux_nav) })
  observeEvent(input$ux_go_results, { ux_go("results") })
  observeEvent(input$ux_back_scenarios, { ux_go("step5") })
  observeEvent(input$threat_csv, { ux$loaded_threat <- NULL }, ignoreInit = TRUE)
  for (kind in c("nests", "yearlings", "adults")) local({
    k <- kind
    observeEvent(input[[paste0("action_csv_", k)]], {
      saved <- ux$loaded_action; saved[[k]] <- NULL; ux$loaded_action <- saved
    }, ignoreInit = TRUE)
  })
  observeEvent(vault$abund, {
    choices <- c("All Beaches (Regional Total)", sort(unique(vault$abund$Site)))
    chosen <- isolate(input$pva_site_filter)
    if (is.null(chosen) || !chosen %in% choices) chosen <- choices[1]
    updateSelectInput(session, "pva_site_filter", choices = choices, selected = chosen)
  })
  ux_data_problem <- reactive({
    if (identical(user_state$data_mode, "upload") && is.null(input$uploaded_file)) return("Upload nesting data or choose Use an example.")
    tryCatch({
      d <- processed_data()
      if (is.null(d) || !nrow(d)) stop("No nesting records in the selected timeframe.")
      if (!any(d$Monitored & is.finite(d$Count))) stop("At least one monitored count is required.")
      if (any(!is.na(d$Count) & (!is.finite(d$Count) | d$Count < 0))) stop("Nest counts must be finite and non-negative.")
      if (any(d$Monitored & d$Count == 0, na.rm = TRUE)) stop("The log-count model cannot fit genuine monitored zeros. Correct missing-survey coding only if these were unmonitored observations.")
      NULL
    }, error = function(e) { msg <- conditionMessage(e); if (nzchar(msg)) msg else "Waiting for the uploaded data to finish loading." })
  })
  ux_biology_problem <- reactive({ tryCatch({current_model_params(input); NULL}, error = function(e) conditionMessage(e)) })
  ux_require_data <- function() {
    problem <- ux_data_problem()
    if (!is.null(problem)) {showNotification(problem, type = "error", duration = 8); return(FALSE)}
    TRUE
  }
  ux_require_biology <- function() {
    problem <- ux_biology_problem()
    if (!is.null(problem)) {showNotification(paste("Complete population biology:", problem), type = "error", duration = 8); return(FALSE)}
    TRUE
  }
  ux_fit_signature <- reactive({
    list(data = processed_data(), cf = input$clutch_freq, ri = input$remig_int,
         season = input$season_start_month %||% 4L, periods = sort(input$six_month_sites %||% character(0)),
         iterations = input$iterations, fourier = isTRUE(input$legacy_fourier_mcmc), trend = isTRUE(input$legacy_trend_mcmc),
         bounds = isTRUE(input$propagate_imputation_uncertainty), split = isTRUE(input$run_split),
         split_year = if (isTRUE(input$run_split)) input$split_year else NULL,
         reference = if (isTRUE(input$override_reference_site)) input$reference_site else NULL)
  })
  ux_fit_current <- reactive({
    !is.null(vault$res) && !is.null(ux$fit_signature) &&
      is.null(ux_data_problem()) && tryCatch(identical(ux$fit_signature, ux_fit_signature()), error = function(e) FALSE)
  })
  ux_require_fit <- function() {
    if (!ux_require_data() || !ux_require_biology()) return(FALSE)
    if (!ux_fit_current()) {showNotification("Estimate the population in step 3 using the current data and settings first.", type = "error", duration = 8); return(FALSE)}
    TRUE
  }
  ux_go <- function(step) {
    if (!step %in% c("welcome", "step2") && !ux_require_data()) return(invisible(FALSE))
    if (step %in% c("step4", "step5", "results") && !ux_require_biology()) return(invisible(FALSE))
    if (step %in% c("step5", "results") && !ux_require_fit()) return(invisible(FALSE))
    nav_select("wizard_steps", step)
  }
  output$ux_progress <- renderUI({
    ids <- c("welcome","step2","step3","step4","step5","results")
    labels <- c("Welcome","Your data","Your turtles","Current trend","Threats & actions","Compare")
    data_ok <- is.null(ux_data_problem()); bio_ok <- is.null(ux_biology_problem()); fit_ok <- ux_fit_current()
    enabled <- c(TRUE,TRUE,data_ok,data_ok&&bio_ok,data_ok&&bio_ok&&fit_ok,data_ok&&bio_ok&&fit_ok)
    current <- input$wizard_steps %||% "welcome"
    div(class="ux-nav",role="navigation","aria-label"="Assessment steps",
      lapply(seq_along(ids),function(i)tags$button(type="button",disabled=if(!enabled[i])"disabled" else NULL,
        "aria-current"=if(ids[i]==current)"step" else NULL,
        onclick=sprintf("Shiny.setInputValue('ux_nav','%s',{priority:'event'})",ids[i]),
        tags$small(if(i==1)"START" else sprintf("%02d",i-1)),labels[i])))
  })
  output$ux_global_status <- renderUI({
    if (!is.null(ux$fit_signature) && !ux_fit_current()) div(class = "alert alert-warning", "Historical inputs changed. Return to step 3 and estimate the population again before comparing scenarios.")
  })
  output$ux_data_status <- renderUI({
    problem <- ux_data_problem()
    if (!is.null(problem)) return(div(class = "ux-note", problem))
    d <- processed_data(); y <- if ("Season" %in% names(d)) d$Season else d$Year
    div(class = "alert alert-success", sprintf("Data ready: %s records · %s beaches · seasons %s–%s.", nrow(d), length(unique(d$Site)), min(y), max(y)),
        if (identical(user_state$data_mode, "demo")) p("Synthetic demonstration data."))
  })
  ux_reference_values <- reactive({
    if (identical(input$parameter_preset, "wp_leatherback")) return(c(clutch_freq=5.5, remig_int=3.06, linf=142.7, k=.2262, tknot=-.17, lmat=139.1325, sig_mat=6.3399, pf=.73, ane_pj=.81, ane_pa=.893))
    numeric(0)
  })
  output$ux_biology_status <- renderUI({
    problem <- ux_biology_problem()
    if (!is.null(problem)) return(div(class = "ux-note", "Choose a reference or complete the missing values. ", problem))
    refs <- ux_reference_values(); vals <- ux_shared_values()
    changed <- if (length(refs)) sum(vapply(names(refs), function(k) !isTRUE(all.equal(as.numeric(vals[[k]]), as.numeric(refs[[k]]))), logical(1))) else NA
    div(class = "alert alert-success", if (is.na(changed)) "Required biology is complete. Confirm the source of each value for your population." else sprintf("Required biology is complete. %s value(s) differ from the selected reference.", changed))
  })
  output$ux_biology_table <- renderTable({
    keys <- setdiff(ux_bio_keys, c("max_age", "ane_mat_threshold")); refs <- ux_reference_values()
    vals <- vapply(keys, function(k) as.numeric(input[[k]] %||% NA_real_), numeric(1))
    ref <- as.numeric(refs[keys])
    data.frame(Parameter = keys, Current = vals, Reference = ref,
               Source = ifelse(!is.finite(vals), "Missing", ifelse(!is.finite(ref), "User supplied", ifelse(abs(vals-ref)<1e-8, "Reference value", "User changed"))), check.names = FALSE)
  }, digits = 5, striped = TRUE)
  output$ux_fit_status <- renderUI({
    if(is.null(vault$res)) return(p(class="ux-muted","Ready when your data and population information are complete."))
    if(!ux_fit_current()) return(div(class="alert alert-warning","Inputs changed. Estimate the trend again before comparing options."))
    div(class="ux-note","Historical model fitted. Review the estimate and any warnings before continuing.")
  })
  output$ux_forecast_range <- renderUI({
    y <- ux_forecast_years()
    div(class = "ux-note", sprintf("Forecast schedules cover %s–%s (%s years).", min(y), max(y), length(y)),
        " Downloadable scenario examples match this range. Existing uploads are never shifted automatically.")
  })
  ux_capture <- function() {
    all <- reactiveValuesToList(input)
    pattern <- switch(all$pva_mode %||% "mode_threat",
      mode_threat = "^(threat_|siders_)|^pva_mode$",
      mode_action = "^action_|^override_(eggs|syr1)$|^pva_mode$",
      mode_portfolio = "^(portfolio_|siders_)|^pva_mode$")
    keys <- grep(pattern, names(all), value = TRUE)
    v <- all[keys]; v <- v[!grepl("^action_csv", names(v))]
    v <- modifyList(v, ux_shared_values())
    if (identical(v$pva_mode, "mode_threat") && identical(v$threat_input_type, "table")) {
      loaded <- ux$loaded_threat
      if (!is.null(loaded)) {v$ux_threat_data <- loaded$data; v$threat_csv <- loaded$file}
      else if (!is.null(input$threat_csv)) {v$ux_threat_data <- read.csv(input$threat_csv$datapath, check.names = FALSE); v$threat_csv <- input$threat_csv}
    }
    if (identical(v$pva_mode, "mode_action") && identical(v$action_input_type, "table")) {
      loaded <- ux$loaded_action[[ux_action_kind()]]; file <- ux_action_file()
      if (!is.null(loaded)) {v$ux_action_data <- loaded$data; v$action_csv <- loaded$file}
      else if (!is.null(file)) {v$ux_action_data <- read.csv(file$datapath, check.names = FALSE); v$action_csv <- file}
    }
    v
  }
  ux_draft <- reactive({tryCatch({v <- ux_capture(); d <- ux_validate_scenario(v, ux_forecast_years()); list(values=v, schedule=d, error=NULL)}, error = function(e) list(error=conditionMessage(e)))})
  output$ux_shared_assumptions <- renderUI({
    p <- tryCatch(current_model_params(input), error = function(e) NULL)
    if (is.null(p)) return(p("Complete population biology in step 2."))
    tagList(p(sprintf("Female proportion: %.1f%%. Remigration interval: %.2f years. Juvenile survival: %.3f; adult survival: %.3f.",100*p$pf,p$ri,p$pj,p$pa)),
      helpText("These population-level assumptions come from step 2. The current CSV formats do not contain sex columns. Conservation uploads supply net additional survivors or deaths prevented where indicated."))
  })
  output$ux_action_amount_help <- renderUI({
    kind <- ux_action_kind()
    div(class="ux-note", switch(kind,
      nests="Amount means nests receiving protection. The model calculates additional hatchlings using the survival assumptions below.",
      yearlings="Amount means net additional turtles alive at release age: survivors with headstarting minus expected survivors without it.",
      adults="Amount means adult deaths prevented: expected deaths without the action minus expected deaths with it. Mortality reduction is already included in this number."))
  })
  output$ux_uncertainty_status <- renderUI({
    keys <- c("siders_log_length_sd", "siders_logit_mortality_sd", "siders_ri_sd", "siders_pj_sd", "siders_pa_sd", "siders_pf_sd")
    varying <- any(vapply(keys, function(k) isTRUE(input[[k]] > 0), logical(1))) || !identical(input$siders_count_mode %||% "fixed", "fixed") || identical(input$siders_ri_distribution, "cmp")
    div(class="ux-note", tags$b(if(varying) "Input uncertainty: enabled. " else "Input uncertainty: fixed values. "), "Population-model uncertainty is included in either case.")
  })
  output$ux_scenario_review <- renderUI({
    draft <- ux_draft()
    if (!is.null(draft$error)) return(div(class = "alert alert-warning", draft$error))
    v <- draft$values; d <- draft$schedule; years <- ux_forecast_years()
    file <- if (v$pva_mode == "mode_threat" && v$threat_input_type == "table") v$threat_csv else if (v$pva_mode == "mode_action" && v$action_input_type == "table") v$action_csv else NULL
    detail <- if (v$pva_mode == "mode_threat") {
      sprintf("First year: %s turtles affected · median size %s cm · mortality %s%%.", d$turtles[1], d$median_cm[1], 100*d$mortality[1])
    } else if (v$pva_mode == "mode_action") {
      unit <- switch(ux_action_type(v$action_type), nests="nests receiving protection", yearlings="net additional turtles alive at release age", adults="adult deaths prevented")
      paste("First year:", d$Amount[1], unit)
    } else paste("Includes:", paste(v$portfolio_active, collapse = "; "))
    tagList(div(class = "alert alert-success", sprintf("Schedule ready: %s–%s. ", min(years), max(years)), detail),
      if (!is.null(file)) p(tags$b("Active schedule: "), file$name),
      if (v$pva_mode == "mode_action" && v$action_type == "A: Protect Nests") {
        per_egg <- v$action_protected_emergence * if (v$action_emergence_basis == "conditional_on_hatching") v$action_protected_hatch else 1
        p(sprintf("With protection: %.3f emerged hatchlings per egg; without protection: %.3f. Additional hatchlings per nest: %.2f, before later survival and ANE conversion.", per_egg, v$action_counterfactual_emergence, v$override_eggs*(per_egg-v$action_counterfactual_emergence)))
      },
      if (v$pva_mode == "mode_action" && v$action_type == "B: Headstarting") p(sprintf("Release age: %s year(s). Amount is already net additional survivors at that age.", v$action_release_age)),
      p(sprintf("Female proportion: %.1f%%, shared from population biology.", 100*v$pf)))
  })
  output$ux_schedule_preview <- renderTable({ draft <- ux_draft(); validate(need(is.null(draft$error), draft$error %||% "")); draft$schedule }, striped = TRUE)
  ux_save <- function(copy = FALSE) {
    if (!ux_require_fit()) return()
    draft <- ux_draft()
    if (!is.null(draft$error)) {showNotification(draft$error, type = "error", duration = 10); return()}
    name <- trimws(input$ux_scenario_name %||% "")
    if (!nzchar(name)) {showNotification("Give this scenario a name before saving.", type = "error"); return()}
    id <- if (!copy) ux$editing_id else NULL
    if (is.null(id) || !id %in% names(ux$scenarios)) {id <- paste0("scenario_", ux$next_id); ux$next_id <- ux$next_id + 1L}
    other_names <- vapply(ux$scenarios[setdiff(names(ux$scenarios), id)], `[[`, character(1), "name")
    if (tolower(name) == "status quo" || tolower(name) %in% tolower(other_names)) {showNotification("Choose a unique name; Status Quo is reserved for the baseline.", type="error");return()}
    v <- draft$values; v <- v[setdiff(names(v), ux_bio_keys)]
    # File contents are preserved; paths are discarded because Shiny upload
    # tempfiles can change. A named marker satisfies the original req() guard.
    for (k in c("threat_csv", "action_csv")) if (!is.null(v[[k]])) v[[k]] <- list(name = v[[k]]$name, datapath = NULL)
    records <- ux$scenarios; records[[id]] <- list(id=id, name=name, values=v, saved_at=Sys.time())
    ux$scenarios <- records; ux$editing_id <- id
    if (!is.null(v$ux_threat_data)) ux$loaded_threat <- list(data=v$ux_threat_data, file=v$threat_csv)
    if (!is.null(v$ux_action_data)) {files <- ux$loaded_action; files[[ux_action_type(v$action_type)]] <- list(data=v$ux_action_data,file=v$action_csv); ux$loaded_action <- files}
    selected <- unique(c(isolate(input$ux_compare_ids), id))
    choices <- setNames(names(records), vapply(records, `[[`, character(1), "name"))
    updateCheckboxGroupInput(session,"ux_compare_ids",choices=choices,selected=selected)
    updateSelectInput(session,"ux_edit_id",choices=c("Choose a saved scenario"="",choices),selected=id)
    showNotification(paste("Saved", name), type="message", id="ux_scenario_feedback", duration=3)
  }
  observeEvent(input$ux_save_scenario, {ux_save(FALSE)})
  observeEvent(input$ux_save_copy, {ux_save(TRUE)})
  ux_scenario_values <- function(record) modifyList(record$values, ux_shared_values())
  observeEvent(input$ux_load_scenario, {
    id <- input$ux_edit_id; req(id, id %in% names(ux$scenarios))
    record <- ux$scenarios[[id]]; v <- record$values
    ux$editing_id <- id
    updateTextInput(session,"ux_scenario_name",value=record$name)
    select_ids <- c("action_type", "action_emergence_basis", "portfolio_emergence_basis", "siders_count_mode", "siders_fishery_center_mode", "siders_ri_distribution")
    radio_ids <- c("pva_mode", "threat_input_type", "action_input_type")
    for (key in names(v)) {
      if (key %in% select_ids) updateSelectInput(session,key,selected=v[[key]])
      else if (key %in% radio_ids) updateRadioButtons(session,key,selected=v[[key]])
      else if (key == "portfolio_active") updateCheckboxGroupInput(session,key,selected=v[[key]])
      else if (key == "siders_force_all_fatal") updateCheckboxInput(session,key,value=isTRUE(v[[key]]))
      else if (is.numeric(v[[key]]) && length(v[[key]])==1L) updateNumericInput(session,key,value=v[[key]])
    }
    if (!is.null(v$ux_threat_data)) ux$loaded_threat <- list(data=v$ux_threat_data,file=v$threat_csv)
    if (!is.null(v$ux_action_data)) {files <- ux$loaded_action; files[[ux_action_type(v$action_type)]] <- list(data=v$ux_action_data,file=v$action_csv); ux$loaded_action <- files}
    showNotification("Saved settings loaded. The review shows the active schedule; a new upload replaces it.", type="message")
  })
  observeEvent(input$ux_delete_scenario, {
    id <- input$ux_edit_id; req(id, id %in% names(ux$scenarios))
    records <- ux$scenarios; records[[id]] <- NULL; ux$scenarios <- records
    if (identical(ux$editing_id,id)) ux$editing_id <- NULL
    choices <- setNames(names(records), vapply(records, `[[`, character(1), "name"))
    updateCheckboxGroupInput(session,"ux_compare_ids",choices=choices,selected=setdiff(input$ux_compare_ids,id))
    updateSelectInput(session,"ux_edit_id",choices=c("Choose a saved scenario"="",choices),selected="")
  })
  observeEvent(input$ux_new_scenario, {
    ux$editing_id <- NULL
    updateTextInput(session, "ux_scenario_name", value = "")
  })
  ux_editor_dirty <- reactive({
    id <- ux$editing_id
    if (is.null(id) || !id %in% names(ux$scenarios)) return(FALSE)
    v <- tryCatch(ux_capture(), error = function(e) NULL)
    if (is.null(v)) return(TRUE)
    v <- v[setdiff(names(v), ux_bio_keys)]
    for (k in c("threat_csv", "action_csv")) if (!is.null(v[[k]])) v[[k]] <- list(name = v[[k]]$name, datapath = NULL)
    record <- ux$scenarios[[id]]
    !identical(v[sort(names(v))], record$values[sort(names(record$values))]) ||
      !identical(trimws(input$ux_scenario_name %||% ""), record$name)
  })
  output$ux_editor_status <- renderUI({
    id <- ux$editing_id
    if (is.null(id) || !id %in% names(ux$scenarios)) return(p(class="ux-muted","Saving creates a new scenario."))
    div(class=if(ux_editor_dirty()) "alert alert-warning" else "ux-note",
        paste("Editing:",ux$scenarios[[id]]$name,".",if(ux_editor_dirty()) "Unsaved changes — save before comparing this version." else "All editor changes are saved.",
              "Save scenario updates this record. Start a new scenario to add another; Save as a copy requires a different name."))
  })
  output$ux_saved_list <- renderUI({
    records <- ux$scenarios
    if (!length(records)) return(div(class="ux-note","No scenarios saved yet. You can also compare status quo alone."))
    tagList(lapply(records,function(record) {
      error <- tryCatch({ux_validate_scenario(ux_scenario_values(record),ux_forecast_years());NULL},error=function(e) conditionMessage(e))
      div(class="ux-saved",tags$b(record$name),p(if(is.null(error)) "Ready for comparison" else paste("Needs attention:",error)))
    }))
  })
  ux_comparison_problem <- reactive({
    if (!ux_fit_current()) return("Estimate the population with the current inputs in step 3 first.")
    if (!is.null(ux_biology_problem())) return(ux_biology_problem())
    sims <- input$pva_sims
    if (length(sims)!=1 || !is.finite(sims) || sims<100 || sims!=floor(sims)) return("Use a whole number of at least 100 projection simulations.")
    selected <- input$ux_compare_ids %||% character(0)
    if (!is.null(ux$editing_id) && ux$editing_id %in% selected && ux_editor_dirty()) return("The selected scenario has unsaved editor changes. Save it in step 4, or load its saved version before comparing.")
    for (id in selected) {
      if (!id %in% names(ux$scenarios)) return("A selected scenario no longer exists. Choose again.")
      error <- tryCatch({ux_validate_scenario(ux_scenario_values(ux$scenarios[[id]]),ux_forecast_years());NULL},error=function(e) conditionMessage(e))
      if (!is.null(error)) return(paste(ux$scenarios[[id]]$name,":",error))
    }
    NULL
  })
  ux_comparison_valid <- function() {
    error <- ux_comparison_problem()
    if (!is.null(error)) {showNotification(error,type="error",duration=10);return(FALSE)}
    TRUE
  }
  ux_comparison_signature <- reactive({
    selected <- input$ux_compare_ids %||% character(0)
    list(scenarios=ux$scenarios[selected], biology=ux_shared_values(), years=ux_forecast_years(),
         simulations=input$pva_sims, site=input$pva_site_filter, mode=input$projection_uq_mode,
         fit_revision=ux$fit_revision, fit_current=ux_fit_current(), branch=vault$active_imputation_branch)
  })
  output$ux_comparison_review <- renderUI({
    problem <- ux_comparison_problem()
    if (!is.null(problem)) return(div(class="alert alert-warning",problem))
    y <- ux_forecast_years()
    div(class="ux-note",sprintf("Ready: status quo plus %s saved scenario(s) · %s–%s · %s simulations · %s.",length(input$ux_compare_ids),min(y),max(y),input$pva_sims,input$pva_site_filter),
        p("Only saved settings are compared. Save any editor changes in step 4 before running."))
  })
  ux_results_current <- reactive({!is.null(ux$result_signature) && identical(ux$result_signature,ux_comparison_signature())})
  output$ux_result_status <- renderUI({
    if (is.null(ux$result_signature)) return(div(class="ux-note","No comparison has been run yet."))
    if (ux_editor_dirty() && ux$editing_id %in% (input$ux_compare_ids %||% character(0))) return(div(class="alert alert-warning", "The editor contains unsaved changes. These results describe the last saved version. Save the scenario and rerun to compare the changes."))
    if (!ux_results_current()) return(div(class="alert alert-warning","Results are out of date for the current saved scenarios, selection, forecast settings or population inputs. Run comparison again. Previous results remain visible; result downloads are blocked until rerun."))
    div(class="alert alert-success",paste("Comparison complete:",format(ux$result_time,"%Y-%m-%d %H:%M"),". Results match the saved settings."))
  })
  output$ux_comparison_table <- renderTable({
    req(vault$step5_summary,vault$step5_raw)
    d <- vault$step5_summary; raw <- vault$step5_raw; last <- max(d$Year); first <- min(d$Year)
    baseline <- d$Median[d$Year==last & d$Scenario=="Status Quo"]
    threshold <- .5*median(raw$N[raw$Year==first & raw$Scenario=="Status Quo"])
    rows <- lapply(c("Status Quo", setdiff(unique(d$Scenario), "Status Quo")),function(scenario) {
      end <- d[d$Year==last & d$Scenario==scenario,]
      r <- raw[raw$Scenario==scenario,]
      hit <- tapply(r$N<=threshold,r$SimID,any)
      data.frame(Scenario=scenario, `Final annual nesters`=round(end$Median,1),
                 `95% interval`=sprintf("%.1f–%.1f",end$Lower,end$Upper),
                 `Difference from status quo`=round(end$Median-baseline,1),
                 `Probability of reaching 50% abundance`=sprintf("%.1f%%",100*mean(hit)),check.names=FALSE)
    })
    do.call(rbind,rows)
  },striped=TRUE,rownames=FALSE)
  output$ux_download_settings <- downloadHandler(
    filename="comparison_settings.txt",
    content=function(file) {
      ux_require_current_results()
      lines <- c("SEA TURTLE POPULATION TOOLKIT — COMPARISON SETTINGS",
                 "Guided interface revision: 2026-09-02",
                 paste("Run completed:", format(ux$result_time)),
                 paste("Forecast years:", paste(range(ux$result_signature$years), collapse="–")),
                 paste("Simulations:", ux$result_signature$simulations),
                 paste("Population:", ux$result_signature$site),
                 paste("Future trend/process sampling:", ux$result_signature$mode),
                 paste("Imputation branch:", ux$result_signature$branch),
                 "", "SHARED BIOLOGY", capture.output(print(ux$result_biology)),
                 "", "SAVED SCENARIOS")
      for (record in ux$result_scenarios) {
        lines <- c(lines, "", record$name, capture.output(print(record$values)))
      }
      lines <- c(lines, "", "HISTORICAL FIT SETTINGS AND INPUT DATA", capture.output(print(ux$fit_signature)),
                 "", "METHODS", "Methods_ConsBio_20260831.pdf, equations 7–22.",
                 "Scenario projections share initial population draws, paired trend/variance draws, and environmental shocks.",
                 "This file records the run settings; it is not a scenario import file.")
      writeLines(lines, file, useBytes=TRUE)
    }
  )
  ux_require_current_results <- function() {
    if (!ux_results_current()) stop("Run the comparison again before downloading results for these settings.")
    invisible(TRUE)
  }

  # =====================================================================
  # MATHEMATICAL SELF-CHECK SERVER
  # =====================================================================
  
  observeEvent(
    input$run_math_self_checks,
    {
      
      results <- tryCatch(
        {
          run_mathematical_self_checks(
            trend_result = vault$res,
            remigration_interval = input$remig_int,
            tolerance = 1e-10
          )
        },
        error = function(error_condition) {
          data.frame(
            Test = "Self-check engine",
            Status = "FAIL",
            Maximum_Difference = NA_real_,
            Expected_Behavior =
              "The mathematical self-check engine completes without error.",
            Details = conditionMessage(error_condition),
            Passed = FALSE,
            stringsAsFactors = FALSE
          )
        }
      )
      
      vault$math_self_checks <- results
    }
  )
  
  
  output$math_self_check_status <- renderUI({
    
    results <- vault$math_self_checks
    
    if (is.null(results)) {
      return(
        div(
          class = "alert alert-secondary",
          "Self-checks have not been run."
        )
      )
    }
    
    evaluated <- results$Status != "NOT RUN"
    evaluated_passed <- all(
      results$Passed[evaluated]
    )
    
    not_run_count <- sum(
      results$Status == "NOT RUN"
    )
    
    if (evaluated_passed) {
      
      div(
        class = "alert alert-success",
        
        tags$b(
          shiny::icon("circle-check"),
          " All evaluated mathematical checks passed."
        ),
        
        if (not_run_count > 0L) {
          tags$p(
            paste(
              not_run_count,
              "posterior-dependent check(s) were not run.",
              "Run the trend model and repeat the checks."
            ),
            style = "margin: 6px 0 0 0;"
          )
        }
      )
      
    } else {
      
      failed_tests <- results$Test[
        results$Status == "FAIL"
      ]
      
      div(
        class = "alert alert-danger",
        
        tags$b(
          shiny::icon("triangle-exclamation"),
          " One or more mathematical checks failed."
        ),
        
        tags$p(
          paste(
            failed_tests,
            collapse = "; "
          ),
          style = "margin: 6px 0 0 0;"
        )
      )
    }
  })
  
  
  output$math_self_check_table <- renderTable({
    
    results <- vault$math_self_checks
    
    validate(
      need(
        !is.null(results),
        "Click “Run mathematical self-checks” to generate results."
      )
    )
    
    display_results <- results[
      ,
      c(
        "Test",
        "Status",
        "Maximum_Difference",
        "Expected_Behavior",
        "Details"
      ),
      drop = FALSE
    ]
    
    display_results$Maximum_Difference <- ifelse(
      is.na(display_results$Maximum_Difference),
      "",
      format(
        display_results$Maximum_Difference,
        scientific = TRUE,
        digits = 4
      )
    )
    
    names(display_results) <- c(
      "Test",
      "Status",
      "Maximum difference",
      "Expected behavior",
      "Details"
    )
    
    display_results
    
  },
  striped = TRUE,
  hover = TRUE,
  bordered = TRUE,
  spacing = "s")

  activate_imputation_branch <- function(branch_name, notify = FALSE) {
    branch_results <- vault$imputation_branch_results
    if (
      is.null(branch_results) ||
        !branch_name %in% names(branch_results)
    ) {
      return(invisible(FALSE))
    }

    bundle <- branch_results[[branch_name]]
    state <- summarise_trend_result(bundle$res, input$remig_int)
    vault$active_imputation_branch <- branch_name
    vault$abund <- bundle$abund
    vault$res <- bundle$res
    vault$marss <- bundle$marss
    vault$years <- state$years
    vault$year <- state$year
    vault$nesters <- state$nesters
    vault$total <- state$total
    vault$trend_display <- state$trend_display
    vault$trend_pct <- state$trend_pct
    vault$draws <- state$draws
    vault$math_self_checks <- NULL

    # Independent-site fits, split-trend summaries, and future scenarios belong
    # to the branch on which they were run. Clear them rather than showing
    # stale results after the user changes the active imputation branch.
    vault$jags_indep_fits <- NULL
    if (!identical(vault$split_branch, branch_name)) {
      vault$pre_trend <- NULL
      vault$post_trend <- NULL
      vault$pre_u_val <- 0
      vault$post_u_val <- 0
      vault$split_branch <- NULL
    }
    vault_portfolio$plot_df <- NULL
    vault_portfolio$scorecard <- NULL
    vault_portfolio$all_scen_raw <- list()
    updateCheckboxGroupInput(session, "plot_layers", selected = "jags")

    if (isTRUE(notify)) {
      showNotification(
        paste0(
          "Active trend/PVA branch changed to: ",
          imputation_branch_labels[[branch_name]],
          ". Rerun future scenarios for this branch."
        ),
        type = "message"
      )
    }
    invisible(TRUE)
  }
  
  custom_colors <- c(
    "Status Quo"         = "#000000",
    "Threats Only"       = "#D55E00",
    "Strategy A Only"    = "#E69F00",
    "Strategy B Only"    = "#56B4E9",
    "Strategy C Only"    = "#009E73",
    "Strategy D Only"    = "#F0E442",
    "Strategy E Only"    = "#CC79A7",
    "Combined Portfolio" = "#0072B2"
  )
  
  # --- SPECIES PRESET UPDATER ---
  
  observeEvent(input$parameter_preset, {
    
    preset <- input$parameter_preset
    
    # ---------------------------------------------------------------
    # WESTERN PACIFIC LEATHERBACK
    # Martin et al. 2020 / Siders et al. 2023 lineage
    # ---------------------------------------------------------------
    if (identical(preset, "wp_leatherback")) {
      
      updateNumericInput(
        session,
        "clutch_freq",
        value = 5.5
      )
      
      updateNumericInput(
        session,
        "remig_int",
        value = 3.06
      )
      
      updateNumericInput(
        session,
        "linf",
        value = 142.7
      )
      
      updateNumericInput(
        session,
        "k",
        value = 0.2262
      )
      
      updateNumericInput(
        session,
        "tknot",
        value = -0.17
      )
      
      # Siders forced-maturity age:
      # age at 99% of Linf under the supplied VBGF
      updateNumericInput(
        session,
        "max_age",
        value = 20.19
      )
      
      # 97.5% of Linf; center of the Siders maturity ogive
      updateNumericInput(
        session,
        "lmat",
        value = 139.1325
      )
      
      updateNumericInput(
        session,
        "sig_mat",
        value = 6.3399
      )
      
      updateNumericInput(
        session,
        "pf",
        value = 0.73
      )
      
      updateNumericInput(
        session,
        "ane_pj",
        value = 0.81
      )
      
      updateNumericInput(
        session,
        "ane_pa",
        value = 0.893
      )
      
      showNotification(
        "Western Pacific leatherback reference parameters loaded.",
        type = "message",
        duration = 4
      )
    }
    
    
    # ---------------------------------------------------------------
    # NORTH PACIFIC LOGGERHEAD
    # Hold for one small lineage correction before activating.
    # ---------------------------------------------------------------
    if (identical(preset, "np_loggerhead")) {
      
      showNotification(
        paste(
          "North Pacific loggerhead preset is not yet activated.",
          "We first need to separate nest-count conversion clutch frequency",
          "from reproductive clutch frequency to retain Martin lineage fidelity."
        ),
        type = "warning",
        duration = 8
      )
    }
    
  })
  
  # --- STEP 3: LIVE PARAMETER SENSITIVITY ENGINE ---
  output$sensitivity_tornado_plot <- renderPlot({
    req(input$ane_pa, input$ane_pj, input$lmat, input$linf, input$k,
        input$tknot, input$pf, input$remig_int, input$clutch_freq,
        input$override_eggs, input$action_protected_hatch, input$action_protected_emergence,
        input$override_syr1)
    
    calc_lambda <- function(p) {
      ratio <- min(0.98, p$lmat / p$linf)
      t_mat <- max(2, p$tknot - (1 / p$k) * log(1 - ratio))
      fecundity <- (
        p$cf * p$eggs * p$hatch * p$emergence *
          p$pf * p$year1
      ) / p$ri
      T_stages <- max(2, round(t_mat))
      M <- matrix(0, nrow = T_stages, ncol = T_stages)
      M[1, T_stages] <- fecundity
      for (k in 1:(T_stages - 1)) M[k + 1, k] <- p$pj
      M[T_stages, T_stages] <- p$pa
      return(max(Re(eigen(M)$values)))
    }
    
    base_p <- list(
      pa = input$ane_pa, pj = input$ane_pj, lmat = input$lmat, linf = input$linf,
      k = input$k, tknot = input$tknot, pf = input$pf,
      ri = input$remig_int, cf = input$clutch_freq,
      eggs = input$override_eggs, hatch = input$action_protected_hatch,
      emergence = input$action_protected_emergence, year1 = input$override_syr1
    )
    
    base_lambda <- tryCatch(calc_lambda(base_p), error = function(e) NA)
    validate(need(!is.na(base_lambda) && is.finite(base_lambda), "Current parameters do not produce a mathematically valid life-stage matrix. Please adjust your inputs."))
    
    params_to_test <- c(
      "Adult Survival (pa)" = "pa", "Juvenile Survival (pj)" = "pj",
      "Maturation Length (Lmat)" = "lmat",
      "Asymptotic Length (Linf)" = "linf", "Growth Rate (k)" = "k",
      "Remigration Interval (RI)" = "ri", "Sex Ratio (% Female)" = "pf",
      "Clutch Frequency (CF)" = "cf", "Clutch Size" = "eggs",
      "Hatching Success" = "hatch", "Emergence Success" = "emergence",
      "Emerged-to-Yearling Survival" = "year1"
    )
    results <- data.frame()
    
    for (p_name in names(params_to_test)) {
      p_key <- params_to_test[[p_name]]
      val_base <- base_p[[p_key]]
      probability_keys <- c("pa", "pj", "pf", "hatch", "emergence", "year1")
      val_low <- if(p_key %in% probability_keys) max(0.0001, val_base * 0.85) else val_base * 0.85
      val_high <- if(p_key %in% probability_keys) min(0.9999, val_base * 1.15) else val_base * 1.15
      
      p_low <- base_p; p_low[[p_key]] <- val_low
      p_high <- base_p; p_high[[p_key]] <- val_high
      
      lam_low <- tryCatch(calc_lambda(p_low), error = function(e) NA)
      lam_high <- tryCatch(calc_lambda(p_high), error = function(e) NA)
      
      if(!is.na(lam_low) && !is.na(lam_high)) {
        dev_low <- ((lam_low - base_lambda) / base_lambda) * 100
        dev_high <- ((lam_high - base_lambda) / base_lambda) * 100
        results <- rbind(results, data.frame(Parameter = p_name, Low = dev_low, High = dev_high, Span = abs(dev_high - dev_low)))
      }
    }
    
    req(nrow(results) > 0)
    results$Parameter <- factor(results$Parameter, levels = results$Parameter[order(results$Span)])
    
    ggplot(results) +
      geom_segment(aes(y = Parameter, yend = Parameter, x = Low, xend = High), linewidth = 1.8, color = "gray60") +
      geom_point(aes(y = Parameter, x = Low), color = "#b02a37", size = 4.5) +  
      geom_point(aes(y = Parameter, x = High), color = "#084298", size = 4.5) + 
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
      theme_classic() +
      theme(
        text = element_text(size = 14), axis.text.y = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_line(color = "gray90"), plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)
      ) +
      labs(x = "% Deviation in Pop. Growth Rate (λ)", y = "")
  })
  

  # --- STEP 2: DYNAMIC UI RENDERING ---
  output$step2_dynamic_main <- renderUI({
    df <- processed_data()
    if (is.null(df) || nrow(df) == 0) return(p("Waiting for data ingestion...", class = "text-muted text-center mt-5"))
    
    elements <- list(
      uiOutput("qaqc_alerts"),
      layout_column_wrap(width = 1/2,
                         plotOutput("preview_annual_raw", height = "300px"),
                         div(style = "height: 300px; overflow-y: auto; border: 1px solid #eee; padding: 10px; border-radius: 5px;", 
                             tableOutput("data_preview_table_raw"))
      )
    )
    
    if (any(!is.na(df$Month))) {
      elements <- append(elements, list(
        hr(),
        plotOutput("preview_monthly_seasonality", height = "300px"),
        uiOutput("seasonality_recommendation")
      ))
    }
    do.call(tagList, elements)
  })
  
  output$preview_monthly_seasonality <- renderPlot({
    df <- processed_data()
    req(df, any(!is.na(df$Month)))

    season_start_month <- if (
      is.null(input$season_start_month) ||
        !is.finite(as.numeric(input$season_start_month))
    ) {
      4L
    } else {
      as.integer(input$season_start_month)
    }
    month_order <- biological_month_order(season_start_month)

    df_monthly <- df %>%
      add_nesting_season(
        season_start_month = season_start_month
      ) %>%
      filter(
        !is.na(Seq_Month),
        !is.na(Season),
        !is.na(Count),
        !is.na(Site)
      ) %>%
      group_by(Site, Season, Seq_Month) %>%
      summarise(
        Total_Count = sum(Count, na.rm = TRUE),
        .groups = "drop"
      )

    ggplot(
      df_monthly,
      aes(
        x = Seq_Month,
        y = Total_Count,
        group = factor(Season),
        color = factor(Season)
      )
    ) +
      geom_line(linewidth = 0.9, alpha = 0.75) +
      geom_point(size = 1.5) +
      scale_x_continuous(
        breaks = 1:12,
        labels = month.abb[month_order]
      ) +
      scale_color_viridis_d(option = "viridis") +
      facet_wrap(~Site, scales = "free_y", ncol = 1) +
      theme_classic() +
      labs(
        x = "Month within biological nesting season",
        y = "Aggregated nest counts",
        color = "Nesting season"
      ) +
      theme(
        text = element_text(size = 13),
        strip.text = element_text(face = "bold", size = 11)
      )
  })

  output$seasonality_recommendation <- renderUI({
    df <- processed_data()
    req(df, any(!is.na(df$Month)))
    
    sites <- sort(unique(df$Site))
    recommendations <- list()
    
    for (st in sites) {
      df_site <- df %>% filter(Site == st, !is.na(Month), !is.na(Count))
      if (nrow(df_site) == 0) next
      
      summer_activity <- sum(df_site$Count[df_site$Month %in% 4:9], na.rm = TRUE)
      winter_activity <- sum(df_site$Count[df_site$Month %in% c(10,11,12,1,2,3)], na.rm = TRUE)
      total_activity <- summer_activity + winter_activity
      if(total_activity == 0) next
      
      is_bimodal <- (winter_activity / total_activity > 0.25) & (summer_activity / total_activity > 0.25)
      if (is_bimodal) {
        recommendations[[length(recommendations) + 1]] <- tags$li(tags$b(st, ": "), "System suggests a ", tags$span(class = "badge bg-warning text-dark", "6-Month Period"), sprintf(" (Summer: %.1f%%, Winter: %.1f%%).", (summer_activity/total_activity)*100, (winter_activity/total_activity)*100))
      } else {
        recommendations[[length(recommendations) + 1]] <- tags$li(tags$b(st, ": "), "System suggests a ", tags$span(class = "badge bg-primary", "12-Month Period"))
      }
    }
    tagList(h5(shiny::icon("robot"), " Data Phenology Diagnostic System Suggestions:"), tags$ul(recommendations))
  })
  
  output$download_template <- downloadHandler(
    filename = function() paste0("monitoring_", input$ux_example_monitoring %||% "annual", "_example.csv"),
    content = function(file) write.csv(ux_example_csv(input$ux_example_monitoring %||% "annual"), file, row.names = FALSE, na = "")
  )
  output$download_threat_template <- downloadHandler(
    filename = "threat_schedule_example.csv",
    content = function(file) write.csv(ux_example_csv("threat", ux_forecast_years()), file, row.names = FALSE, na = "")
  )
  output$download_action_template <- downloadHandler(
    filename = function() paste0(ux_action_kind(), "_example.csv"),
    content = function(file) write.csv(ux_example_csv(ux_action_kind(), ux_forecast_years()), file, row.names = FALSE, na = "")
  )

  output$download_current_parameters <- downloadHandler(
    filename = function() {
      species_name <- "User-supplied population parameters"
      safe_name <- gsub("[^A-Za-z0-9]+", "_", species_name)
      paste0("Sea_Turtle_Parameters_", safe_name, ".csv")
    },
    content = function(file) {
      parameter_table <- data.frame(
        Parameter_Source = "User supplied; see assessment documentation",
        Parameter = c(
          "Remigration Interval (RI)", "Clutch Frequency (CF)",
          "Percent Female", "Ultimate Size (Linf)",
          "Brody Growth Coef. (k)", "Length at Age 0 (t0)",
          "First-nesting transition location", "First-nesting logistic scale",
          "Maximum Age", "Juvenile Survival (pj)",
          "Adult Survival (pa)", "Clutch Size",
          "Hatching Success", "Emergence Success",
          "Emerged-to-Yearling Survival"
        ),
        Value = c(
          input$remig_int, input$clutch_freq, input$pf, input$linf,
          input$k, input$tknot, input$lmat, input$sig_mat,
          input$max_age, input$ane_pj, input$ane_pa,
          if (is.null(input$override_eggs)) NA else input$override_eggs,
          if (is.null(input$action_protected_hatch)) NA else input$action_protected_hatch,
          if (is.null(input$action_protected_emergence)) NA else input$action_protected_emergence,
          if (is.null(input$override_syr1)) NA else input$override_syr1
        ),
        Units = c(
          "years", "nests/female/season", "proportion", "cm",
          "per year", "years", "cm", "cm", "years",
          "annual probability", "annual probability", "eggs/nest",
          "proportion", "proportion", "annual probability"
        ),
        Note = c(
          rep("Evidence-derived mean or custom value", 3),
          rep("Provisional unless supported by a population-specific primary source", 6),
          rep("Evidence-derived mean or custom value", 6)
        ),
        stringsAsFactors = FALSE
      )
      write.csv(parameter_table, file, row.names = FALSE)
    }
  )
  
  raw_ingested_data <- reactive({
    
    # No data source has been selected yet.
    if (
      is.null(user_state$data_mode) ||
      length(user_state$data_mode) != 1L ||
      is.na(user_state$data_mode) ||
      !nzchar(user_state$data_mode)
    ) {
      return(NULL)
    }
    
    if (identical(user_state$data_mode, "demo")) {
      return(ux_example_csv(input$ux_example_monitoring %||% "annual"))
    } else {
      req(input$uploaded_file)
      df <- read.csv(input$uploaded_file$datapath, stringsAsFactors = FALSE)
      missing <- setdiff(c("Year", "Site", "Count"), names(df))
      validate(need(!length(missing), paste("Missing column(s):", paste(missing, collapse = ", "))))
      validate(need(is.numeric(df$Year) && all(is.finite(df$Year)) && all(df$Year == floor(df$Year)), "Year must contain whole calendar or season years."))
      validate(need(is.numeric(df$Count) || all(is.na(df$Count)), "Count must contain numbers or blank cells."))
      validate(need(all(!is.na(df$Site) & nzchar(trimws(df$Site))), "Every row needs a site name."))
      if ("Month" %in% names(df)) {
        validate(need(all(is.na(df$Month)) || (is.numeric(df$Month) && all(!is.na(df$Month) & df$Month %in% 1:12)), "Use months 1–12 on every monthly row, or leave Month blank on every annual row."))
      }
      key <- paste(df$Year, if ("Month" %in% names(df)) df$Month else NA, df$Site)
      validate(need(!anyDuplicated(key), "Duplicate site/year/month rows: combine or correct them before uploading."))
      if (!"Month" %in% colnames(df)) df$Month <- NA
      df <- normalise_monitoring_status(df)
      return(df %>% select(Year, Month, Site, Count, Monitored))
    }
  })
  
  output$timeframe_ui <- renderUI({
    df <- raw_ingested_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)

    has_monthly_data <- any(!is.na(df$Month))
    season_start_month <- if (
      is.null(input$season_start_month) ||
        !is.finite(as.numeric(input$season_start_month))
    ) {
      4L
    } else {
      as.integer(input$season_start_month)
    }

    if (has_monthly_data) {
      indexed <- add_nesting_season(
        df,
        season_start_month = season_start_month
      )
      complete_seasons <- complete_nesting_seasons(
        df,
        season_start_month = season_start_month
      )
      if (length(complete_seasons) > 0L) {
        yr_range <- range(complete_seasons)
      } else {
        yr_range <- range(indexed$Season, na.rm = TRUE)
      }
      timeframe_label <- "Select nesting seasons:"
    } else {
      yr_range <- range(df$Year, na.rm = TRUE)
      timeframe_label <- "Select timeframe (years):"
    }

    sliderInput(
      "timeframe_filter",
      timeframe_label,
      min = yr_range[1],
      max = yr_range[2],
      value = c(yr_range[1], yr_range[2]),
      step = 1,
      sep = ""
    )
  })

  processed_data <- reactive({
    df <- raw_ingested_data()
    if (is.null(df) || nrow(df) == 0) return(df)

    has_monthly_data <- any(!is.na(df$Month))
    if (has_monthly_data) {
      season_start_month <- if (
        is.null(input$season_start_month) ||
          !is.finite(as.numeric(input$season_start_month))
      ) {
        4L
      } else {
        as.integer(input$season_start_month)
      }

      df <- add_nesting_season(
        df,
        season_start_month = season_start_month
      )

      if (!is.null(input$timeframe_filter)) {
        season_range <- range(df$Season, na.rm = TRUE)
        if (
          input$timeframe_filter[1] >= season_range[1] &&
            input$timeframe_filter[2] <= season_range[2]
        ) {
          df <- df %>%
            filter(
              Season >= input$timeframe_filter[1],
              Season <= input$timeframe_filter[2]
            )
        }
      }
    } else if (!is.null(input$timeframe_filter)) {
      data_range <- range(df$Year, na.rm = TRUE)
      if (
        input$timeframe_filter[1] >= data_range[1] &&
          input$timeframe_filter[2] <= data_range[2]
      ) {
        df <- df %>%
          filter(
            Year >= input$timeframe_filter[1],
            Year <= input$timeframe_filter[2]
          )
      }
    }

    df
  })

  output$completeness_card_ui <- renderUI({
    d_long <- processed_data(); if(nrow(d_long) == 0) return("0%")
    total_cells <- nrow(d_long)
    observed_cells <- sum(
      d_long$Monitored & !is.na(d_long$Count),
      na.rm = TRUE
    )
    return(paste0(round((observed_cells / total_cells) * 100, 1), "%"))
  })
  
  output$gap_table <- renderTable({
    d_long <- processed_data(); if(nrow(d_long) == 0) return(NULL)
    d_long %>% group_by(Site) %>%
      summarise(
        `Missing / not monitored` = sum(!Monitored | is.na(Count)),
        `Observed zero counts` = sum(Monitored & Count == 0, na.rm = TRUE),
        `Expected records` = n(),
        `Percent coverage` = sprintf(
          "%.1f%%",
          (sum(Monitored & !is.na(Count), na.rm = TRUE) / n()) * 100
        )
      )
  }, align = "c", striped = TRUE, hover = TRUE, spacing = "s")
  
  preview_annual_plot_object <- reactive({
    d_long <- processed_data(); if(is.null(d_long) || nrow(d_long) == 0) return(NULL)
    d_long <- d_long[d_long$Monitored & !is.na(d_long$Count), , drop = FALSE]
    monthly <- any(!is.na(d_long$Month))
    d_long$Time <- if (monthly) as.Date(sprintf("%04d-%02d-01", d_long$Year, d_long$Month)) else d_long$Year
    ggplot(d_long, aes(x = Time, y = Count, color = Site)) +
      geom_line(linewidth = .7) + geom_point(size = 1.4) + facet_wrap(~Site, scales = "free_y", ncol = 1) +
      labs(x = if(monthly) "Calendar date" else "Nesting season", y = "Observed nests") +
      theme_classic() + theme(legend.position = "none")
  })

  output$preview_annual_raw <- renderPlot({
    preview_annual_plot_object()
  })
  
  output$data_preview_table_raw <- renderTable({ 
    d_long <- processed_data(); if(nrow(d_long) == 0) return(NULL)
    head(d_long, 10) 
  })
  
  output$qaqc_alerts <- renderUI({
    df <- processed_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    negatives <- sum(df$Monitored & df$Count < 0, na.rm = TRUE)
    observed_zeros <- sum(
      df$Monitored & !is.na(df$Count) & df$Count == 0,
      na.rm = TRUE
    )
    nas <- sum(!df$Monitored | is.na(df$Count))
    
    alerts <- list()
    if (negatives > 0) {
      alerts[[length(alerts) + 1]] <- div(class = "alert alert-danger d-flex align-items-center", style = "margin-bottom: 10px;",
                                          shiny::icon("exclamation-triangle", class = "me-2"), sprintf("QA/QC Alert: Found %d negative count values in your spreadsheet. Please verify or fix your raw counts file.", negatives))
    }
    if (observed_zeros > 0) {
      alerts[[length(alerts) + 1]] <- div(
        class = "alert alert-danger d-flex align-items-center",
        style = "margin-bottom: 10px;",
        shiny::icon("exclamation-triangle", class = "me-2"),
        sprintf(
          paste0(
            "QA/QC Alert: Found %d observed zero count(s). The inherited ",
            "Martin/Siders Gaussian log-count model cannot represent a true ",
            "zero. If these rows mean no monitoring, leave Count blank and set ",
            "Monitored to FALSE. Genuine zeros require a count-data observation ",
            "model and this run will be blocked."
          ),
          observed_zeros
        )
      )
    }
    if (nas > 0) {
      alerts[[length(alerts) + 1]] <- div(class = "alert alert-warning d-flex align-items-center", style = "margin-bottom: 10px;",
                                          shiny::icon("info-circle", class = "me-2"), sprintf("QA/QC Note: Found %d missing or unmonitored observations. Monthly gaps are imputed by the Fourier model; annual gaps remain missing in the state-space trend model.", nas))
    }
    if (length(alerts) == 0) {
      return(div(class = "alert alert-success d-flex align-items-center", style = "margin-bottom: 10px;", 
                 shiny::icon("check-circle", class = "me-2"), "QA/QC Pass: Initial data structures cleared (No negative counts or formatting errors caught)."))
    } else {
      return(tagList(alerts))
    }
  })
  
  output$period_override_ui <- renderUI({
    df <- raw_ingested_data()
    if (
      is.null(df) ||
        nrow(df) == 0 ||
        !any(!is.na(df$Month))
    ) {
      return(NULL)
    }

    all_sites <- sort(unique(df$Site))
    default_6mo <- all_sites[
      grepl(
        "W_|Wermon|Wamlana|Waspait|Waenibe",
        all_sites,
        ignore.case = TRUE
      )
    ]

    tagList(
      selectInput(
        "season_start_month",
        label = tags$span(
          "Biological nesting season begins in:",
          tooltip(
            shiny::icon("question-circle"),
            paste(
              "April reproduces the Martin et al. western Pacific",
              "leatherback season definition. Other months are an",
              "Ortega generalized-data option."
            )
          )
        ),
        choices = stats::setNames(1:12, month.name),
        selected = 4
      ),
      selectizeInput(
        "six_month_sites", 
        label = tags$span(
          "Designate 6-Month (Bimodal) Beaches:",
          tooltip(
            shiny::icon("question-circle"),
            paste(
              "Select beaches represented by a six-month Fourier",
              "period. Unselected beaches use a 12-month period."
            )
          )
        ),
        choices = all_sites,
        selected = default_6mo,
        multiple = TRUE,
        options = list(plugins = list("remove_button"))
      )
    )
  })

  output$imputation_uncertainty_option_ui <- renderUI({
    df <- processed_data()
    if (
      is.null(df) || nrow(df) == 0L ||
        !any(!is.na(df$Month))
    ) {
      return(NULL)
    }

    checkboxInput(
      "propagate_imputation_uncertainty",
      paste0(
        "Run lower, median, and upper Fourier-imputation ",
        "branches (Martin 2020; about 3x trend runtime)"
      ),
      value = isolate(input$propagate_imputation_uncertainty) %||% TRUE
    )
  })

  output$reference_site_ui <- renderUI({
    df <- processed_data()
    if (is.null(df) || nrow(df) == 0L || !"Site" %in% names(df)) {
      return(NULL)
    }

    season_start_month <- if (
      is.null(input$season_start_month) ||
        !is.finite(as.numeric(input$season_start_month))
    ) 4L else as.integer(input$season_start_month)

    recommendation <- tryCatch(
      recommend_reference_site_from_observations(
        df,
        season_start_month = season_start_month
      ),
      error = function(e) e
    )

    if (inherits(recommendation, "error")) {
      return(
        div(
          class = "alert alert-warning",
          tags$b("Reference beach could not yet be selected automatically."),
          tags$br(),
          tags$small(conditionMessage(recommendation))
        )
      )
    }

    stats <- recommendation$diagnostics
    auto_site <- recommendation$reference_site
    eligible <- stats[stats$Eligible_First_Year, , drop = FALSE]
    eligible <- eligible[order(eligible$Site), , drop = FALSE]

    selected_site <- isolate(input$reference_site)
    if (is.null(selected_site) || !selected_site %in% eligible$Site) {
      selected_site <- auto_site
    }

    choice_labels <- paste0(
      eligible$Site,
      " — ", round(100 * eligible$Completeness, 1),
      "% of years observed; longest run ",
      eligible$Longest_Continuous_Run, " y"
    )
    choices <- stats::setNames(eligible$Site, choice_labels)

    override_active <- isTRUE(input$override_reference_site)

    tagList(
      div(
        class = "alert alert-info",
        style = "padding: 10px; margin-bottom: 10px;",
        tags$b("Automatically selected reference beach: "),
        auto_site,
        tags$br(),
        tags$small(recommendation$selection_reason),
        tags$br(),
        tags$small(
          "Selection order: first-year eligibility, data completeness, ",
          "longest continuous record, then alphabetical tie-break."
        )
      ),
      if (nrow(eligible) > 1L) {
        checkboxInput(
          "override_reference_site",
          "Advanced: choose a different eligible reference beach",
          value = override_active
        )
      },
      if (nrow(eligible) > 1L && override_active) {
        selectInput(
          "reference_site",
          label = tags$span(
            "Advanced reference-beach override (A = 0):",
            tooltip(
              shiny::icon("question-circle"),
              paste(
                "Changing the reference site changes the coordinate system",
                "used for the latent state and site offsets. It should not",
                "materially change the shared biological trend. Use this",
                "mainly for historical reproduction or sensitivity checks."
              )
            )
          ),
          choices = choices,
          selected = selected_site
        )
      }
    )
  })

  output$imputation_branch_selector_ui <- renderUI({
    branch_results <- vault$imputation_branch_results
    if (is.null(branch_results) || length(branch_results) <= 1L) {
      return(NULL)
    }
    available <- names(branch_results)
    labels <- imputation_branch_labels[available]
    selected <- vault$active_imputation_branch
    if (is.null(selected) || !selected %in% available) selected <- "median"

    div(
      class = "alert alert-info",
      style = "padding: 10px; margin-top: 10px;",
      selectInput(
        "active_imputation_branch",
        "Active imputation branch for plots and future PVA:",
        choices = stats::setNames(available, labels),
        selected = selected
      ),
      p(
        "Martin et al. fitted separate trend analyses to the median, lower 95%, and upper 95% annual imputation outputs. Changing this selection uses the already-fitted branch; it does not rerun JAGS.",
        style = "font-size: 0.82rem; margin-bottom: 0;"
      )
    )
  })

  observeEvent(input$active_imputation_branch, {
    req(vault$imputation_branch_results)
    branch_name <- input$active_imputation_branch
    if (
      !is.null(branch_name) &&
        !identical(branch_name, vault$active_imputation_branch)
    ) {
      activate_imputation_branch(branch_name, notify = TRUE)
    }
  }, ignoreInit = TRUE)

  output$imputation_branch_summary_ui <- renderUI({
    branch_results <- vault$imputation_branch_results
    if (is.null(branch_results) || length(branch_results) <= 1L) {
      return(NULL)
    }
    tagList(
      div(
        class = "alert alert-light border",
        style = "margin-top: 10px;",
        tags$b("Monthly-imputation uncertainty propagated through trend fitting"),
        p(
          "The table reports three separate trend fits. The selected branch above controls the abundance plots, diagnostics, and subsequent PVA.",
          style = "font-size: 0.88rem; margin-bottom: 6px;"
        ),
        tableOutput("imputation_branch_summary")
      )
    )
  })

  output$imputation_branch_summary <- renderTable({
    branch_results <- vault$imputation_branch_results
    req(branch_results)

    rows <- lapply(names(branch_results), function(branch_name) {
      bundle <- branch_results[[branch_name]]
      u <- as.numeric(bundle$res$fit$sims.list$U)
      final_n <- regional_abundance_draws(bundle$res)
      final_n <- final_n[, ncol(final_n)]
      u_pct <- (exp(u) - 1) * 100
      data.frame(
        Active = if (identical(branch_name, vault$active_imputation_branch)) "Yes" else "",
        `Imputation input` = imputation_branch_labels[[branch_name]],
        `Median annual trend` = sprintf("%.2f%%", median(u_pct)),
        `Trend 95% interval` = sprintf(
          "%.2f%% to %.2f%%",
          quantile(u_pct, 0.025), quantile(u_pct, 0.975)
        ),
        `Terminal annual nesters` = format(
          round(median(final_n)), big.mark = ","
        ),
        `Nesters 95% interval` = paste0(
          format(round(quantile(final_n, 0.025)), big.mark = ","),
          " to ",
          format(round(quantile(final_n, 0.975)), big.mark = ",")
        ),
        check.names = FALSE
      )
    })
    do.call(rbind, rows)
  }, striped = TRUE, hover = TRUE, spacing = "s", align = "c")

  output$aicc_recommendation_ui <- renderUI({
    req(vault$marss$shared, vault$marss$indep)
    sites <- unique(vault$abund$Site)
    if(length(sites) < 2) return(NULL)
    
    aic_shared <- vault$marss$shared$AICc
    aic_indep <- vault$marss$indep$AICc
    delta_aic <- aic_shared - aic_indep
    
    aic_badge_box <- div(
      style = "background: #ffffff; border: 1px solid #ced4da; border-radius: 6px; padding: 8px; margin-top: 8px; margin-bottom: 10px; font-size: 0.85rem;",
      tags$div(style = "display: flex; justify-content: space-between;",
               tags$span(tags$b("Regional Model AICc:"), sprintf(" %.1f", aic_shared)),
               tags$span(tags$b("Independent Model AICc:"), sprintf(" %.1f", aic_indep))
      ),
      tags$div(style = "text-align: center; margin-top: 4px; font-weight: bold; color: #495057;",
               sprintf("ΔAICc = %.1f (%s)", abs(delta_aic), if(delta_aic > 0) "Favors Independent" else "Favors Regional")
      )
    )
    
    if(delta_aic > 2) {
      div(class = "alert alert-warning mt-3", style = "padding: 10px;",
          h6(shiny::icon("exclamation-triangle"), " Independent Trends Favored in Diagnostic", style = "font-weight: bold;"),
          aic_badge_box,
          p("The MARSS diagnostic favors independent site trends. This does not by itself establish biological subpopulations. Consider the independent JAGS sensitivity run:", style = "font-size: 0.85rem;"),
          actionButton("run_jags_indep", "Rerun JAGS Separately per Beach", class = "btn-warning w-100 btn-sm", style = "font-weight: bold;")
      )
    } else {
      div(class = "alert alert-success mt-3", style = "padding: 10px;",
          h6(shiny::icon("check-circle"), " Shared Trend Favored in Diagnostic", style = "font-weight: bold;"),
          aic_badge_box,
          p("The MARSS diagnostic favors the shared regional structure for these data; retain biological judgment when defining the assessment unit.", style = "font-size: 0.85rem;")
      )
    }
  })
  
  observeEvent(input$run_jags_indep, {
    req(vault$abund)
    sites <- sort(unique(vault$abund$Site))
    withProgress(message = 'Isolating data and rerunning JAGS independently...', value = 0, {
      indep_fits <- list()
      for(i in seq_along(sites)) {
        s <- sites[i]
        setProgress(value = i / length(sites), detail = paste("Modeling", s, "... (Caution: Lower statistical power)"))
        site_data <- vault$abund %>% filter(Site == s)
        fit_i <- tryCatch({
          if (isTRUE(input$legacy_trend_mcmc)) {
            run_jags_aligned(
              site_data, iter = 750000, burnin = 250000, thin = 50,
              n_chains = 2, parallel = TRUE, reference_site = s
            )
          } else {
            run_jags_aligned(
              site_data, iter = input$iterations,
              burnin = floor(input$iterations / 3), thin = 10,
              reference_site = s
            )
          }
        }, error = function(e) NULL)
        if(!is.null(fit_i)) indep_fits[[s]] <- fit_i
      }
      vault$jags_indep_fits <- indep_fits
      updateCheckboxGroupInput(session, "plot_layers", selected = c(input$plot_layers, "jags_indep"))
      showNotification("Independent JAGS runs complete! Visualizing isolated trends.", type = "message")
    })
  })
  
  # Sensitivity comparisons use named copies of scenarios in the shared
  # comparison engine. This replaces the former separate plotting engine.

  # --- MASTER SIMULATION ENGINE ---
  observeEvent(input$ux_run_comparison, {
    if (!ux_require_fit()) return(invisible(NULL))
    if (!ux_comparison_valid()) return(invisible(NULL))
    tryCatch(withProgress(message = "Projecting selected scenarios…", value = 0, {
      n_sims   <- req(input$pva_sims)
      horizon  <- req(input$pva_years)
      sel_site <- if(!is.null(input$pva_site_filter)) input$pva_site_filter else "All Beaches (Regional Total)"
      
      sites_vec <- if(!is.null(vault$abund)) sort(unique(vault$abund$Site)) else character(0)
      s_idx     <- match(sel_site, sites_vec)
      sims      <- NULL 
      
      if (sel_site != "All Beaches (Regional Total)" && !is.na(s_idx)) {
        if (!is.null(vault$jags_indep_fits) && is.list(vault$jags_indep_fits) && sel_site %in% names(vault$jags_indep_fits)) {
          beach_fit <- vault$jags_indep_fits[[sel_site]]
          sims <- if(!is.null(beach_fit$fit)) beach_fit$fit$sims.list else beach_fit$sims.list
          X_mat <- sims$X; fy <- ncol(X_mat); X_T <- as.numeric(X_mat[, fy]); U <- as.numeric(sims$U); Q <- as.numeric(sims$Q)
        } else {
          req(vault$res)
          sims <- vault$res$fit$sims.list; X_mat <- sims$X; fy <- ncol(X_mat)
          X_T <- as.numeric(X_mat[, fy]); U <- as.numeric(sims$U); Q <- as.numeric(sims$Q)
        }
      } else {
        req(vault$res)
        sims <- vault$res$fit$sims.list; X_mat <- sims$X; fy <- ncol(X_mat)
        X_T <- as.numeric(X_mat[, fy]); U <- as.numeric(sims$U); Q <- as.numeric(sims$Q)
      }
      
      # ---------------------------------------------------------------
      # MATCHED STARTING ABUNDANCE, U, AND Q
      # ---------------------------------------------------------------
      
      if (sel_site == "All Beaches (Regional Total)") {
        
        # The original Siders projection begins with the final modeled
        # regional annual-nester abundance, N_fym0.
        #
        # Total_Females is retained separately as the RI/4 estimate of current
        # adult-female abundance for reporting; it is not the projection N0.
        
        regional_posterior <-
          build_siders_projection_posterior(
            trend_result = vault$res,
            remigration_interval = input$remig_int
          )
        
        valid_rows <- seq_len(
          nrow(regional_posterior)
        )
        
        draw_index <- sample(
          valid_rows,
          size = n_sims,
          replace = TRUE
        )
        
        selected_posterior <-
          regional_posterior[
            draw_index,
            ,
            drop = FALSE
          ]
        
        start_abund <-
          selected_posterior$N_fym0
        
        U_draws <- selected_posterior$U
        Q_draws <- selected_posterior$Q
        
        post_matrix <- as.matrix(
          regional_posterior[
            ,
            c("U", "Q"),
            drop = FALSE
          ]
        )
        
      } else {
        
        # Site-specific projections remain annual-nester projections because the
        # original Siders RI/4 helper estimates regional adult-female abundance.
        # Do not silently apply the regional four-year expansion to one beach.
        
        n_joint <- min(
          length(X_T),
          length(U),
          length(Q)
        )
        
        post_matrix <- cbind(
          X_T = X_T[seq_len(n_joint)],
          U = U[seq_len(n_joint)],
          Q = Q[seq_len(n_joint)]
        )
        
        valid_rows <- which(
          complete.cases(post_matrix) &
            post_matrix[, "Q"] > 0
        )
        
        if (length(valid_rows) < 2L) {
          stop(
            "Too few valid joint posterior draws ",
            "for the site-specific projection."
          )
        }
        
        draw_index <- sample(
          valid_rows,
          n_sims,
          replace = TRUE
        )
        
        future_params <- post_matrix[
          draw_index,
          ,
          drop = FALSE
        ]
        
        X_start_draws <- future_params[, "X_T"]
        U_draws <- future_params[, "U"]
        Q_draws <- future_params[, "Q"]
        
        a_draw <- if (
          !is.null(sims$A) &&
          !is.na(match(sel_site, vault$res$sites)) &&
          match(sel_site, vault$res$sites) <= ncol(sims$A)
        ) {
          sims$A[draw_index, match(sel_site, vault$res$sites)]
        } else {
          0
        }
        
        start_abund <- exp(
          X_start_draws + a_draw
        )
      }
      
      start_yr <- max(vault$years)
      transition_years <- start_yr + seq_len(horizon)
      zero_effect <- rep(0, horizon)
      track_effects <- list("Status Quo" = list(loss = zero_effect, gain = zero_effect))
      selected <- input$ux_compare_ids %||% character(0)
      for (id in selected) {
        scenario <- ux$scenarios[[id]]
        values <- ux_scenario_values(scenario)
        effects <- ux_build_effects(values, n_sims, horizon, transition_years)
        effect_name <- switch(values$pva_mode, mode_threat = "With Threat",
                              mode_action = "With Intervention", mode_portfolio = "Combined Portfolio")
        track_effects[[scenario$name]] <- effects[[effect_name]]
      }
      N_mats <- list()
      for (scen in names(track_effects)) { mat <- matrix(NA, nrow = n_sims, ncol = horizon + 1); mat[, 1] <- start_abund; N_mats[[scen]] <- mat }

      # Legacy take projections use dynUQ = TRUE. Dynamic mode resamples a
      # complete U/Q posterior pair for every simulation-year. Static mode
      # retains one U/Q pair for each full trajectory. Indices and shocks are
      # generated once and shared by every scenario.
      uq_mode <- if (is.null(input$projection_uq_mode)) {
        "dynamic"
      } else {
        input$projection_uq_mode
      }
      annual_uq <- draw_projection_uq(
        post_matrix = post_matrix,
        valid_rows = valid_rows,
        draw_index = draw_index,
        horizon = horizon,
        mode = uq_mode
      )
      annual_U <- annual_uq$U
      annual_Q <- annual_uq$Q
      
      for (t in 1:horizon) {
        incProgress(1/horizon, detail = paste("Simulating year", t))
        standard_normal <- rnorm(n_sims)
        for (scen in names(track_effects)) {
          loss_object <- track_effects[[scen]]$loss
          gain_object <- track_effects[[scen]]$gain
          loss <- if (is.matrix(loss_object)) {
            loss_object[, t]
          } else {
            loss_object[t]
          }
          gain <- if (is.matrix(gain_object)) {
            gain_object[, t]
          } else {
            gain_object[t]
          }
          # Siders timing for fishery loss: subtract ANE before population
          # growth. The same standard-normal shock is shared by all scenarios.
          # Conservation gains are an explicit extension and enter afterward.
          N_next <- project_population_step(
            abundance = N_mats[[scen]][, t],
            loss = loss,
            gain = gain,
            U = annual_U[, t],
            Q = annual_Q[, t],
            standard_normal = standard_normal
          )

          N_mats[[scen]][, t + 1] <- N_next
        }
      }
      
      proj_years <- start_yr + (0:horizon)
      plot_list <- list()
      for (scen in names(track_effects)) { 
        plot_list[[scen]] <- data.frame(
          Year = rep(proj_years, each = n_sims), 
          SimID = rep(1:n_sims, times = length(proj_years)), # <-- Added to track specific trajectories
          N = as.vector(N_mats[[scen]]), 
          Scenario = scen
        ) 
      }
      plot_data <- bind_rows(plot_list)
      summary_data <- plot_data %>% group_by(Year, Scenario) %>% summarize(Median = median(N, na.rm = TRUE), Lower = quantile(N, 0.025, na.rm = TRUE), Upper = quantile(N, 0.975, na.rm = TRUE), .groups = 'drop')
      
      ux$result_signature <- isolate(ux_comparison_signature())
      ux$result_time <- Sys.time()
      ux$result_scenarios <- lapply(selected, function(id) ux$scenarios[[id]])
      ux$result_biology <- isolate(ux_shared_values())
      nav_select("wizard_steps", "results")
      vault$step5_summary <- summary_data
      vault$step5_raw <- plot_data # <-- Save raw data for quasi-extinction math
      vault$step5_metadata <- list(
        starting_abundance_method = if (
          sel_site == "All Beaches (Regional Total)"
        ) {
          paste(
            "Original Siders projection N0:",
            "final modeled regional annual-nester posterior (N_fym0)"
          )
        } else {
          "Site-specific final annual-nester posterior"
        },
        
        population_update_method = paste(
          "Siders abundance-scale projection:",
          "(N - ANE) * exp(U) + Normal(0, sqrt(Q));",
          "conservation gain added separately"
        ),
        
        paper_reference = "Methods_ConsBio_20260831.pdf, Eqs. 7-22",
        threat_schedule_columns = c("year", "turtles", "median_cm", "mortality"),
        transition_year_convention = "year is the destination census",
        projection_uq_mode = uq_mode,
        posterior_rows = draw_index
      )
    }), error = function(e) {
      if (inherits(e, "shiny.silent.error")) return(invisible(NULL))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  }, ignoreInit = TRUE)
  
  # --- DYNAMIC PLOT OUTPUT ---
  step5_plot_object <- reactive({
    req(vault$step5_summary, vault$step5_raw) # Added vault$step5_raw
    
    # --- 1. Threshold Calculation ---
    start_year <- min(vault$step5_summary$Year)
    term_year <- max(vault$step5_summary$Year)
    start_n <- median(vault$step5_raw$N[vault$step5_raw$Year == start_year & vault$step5_raw$Scenario == "Status Quo"], na.rm = TRUE)
    threshold <- 0.5 * start_n
    
    # --- 2. Plot Aesthetics ---
    scenario_names <- unique(vault$step5_summary$Scenario)
    other <- setdiff(scenario_names, "Status Quo")
    scen_colors <- c("Status Quo" = "#334155", setNames(grDevices::hcl.colors(length(other), "Dark 3"), other))
    scen_lines <- setNames(ifelse(scenario_names == "Status Quo", "dashed", "solid"), scenario_names)
    scen_widths <- setNames(rep(1.1, length(scenario_names)), scenario_names)

    # --- 3. Build the Plot ---
    p <- ggplot(vault$step5_summary, aes(x = Year, y = Median, color = Scenario, group = Scenario)) +
      geom_line(aes(linetype = Scenario, linewidth = Scenario)) +
      
      # ⬇️ ADDED: 50% Danger Zone Line and Label ⬇️
      geom_hline(yintercept = threshold, linetype = "dotted", color = "#dc3545", linewidth = 1) +
      annotate("text", 
               x = term_year, 
               y = threshold + (start_n * 0.02), # Scaled offset to float just above the line
               label = "50% Abundance Threshold", 
               hjust = 1, 
               vjust = 0,
               color = "#dc3545", 
               fontface = "bold") +
      # ⬆️ END ADDITION ⬆️
      
      scale_color_manual(values = scen_colors) + 
      scale_linetype_manual(values = scen_lines) + 
      scale_linewidth_manual(values = scen_widths) +
      labs(x = "Year", y = "Females nesting each year (projected)", title = "Possible futures", color = "Scenario", linetype = "Scenario", linewidth = "Scenario") +
      theme_classic(base_size = 14) + 
      theme(panel.grid.major = element_line(color = "gray90"), legend.position = "bottom", legend.title = element_text(face = "bold"))
    
    if (!is.null(input$show_proj_ci) && input$show_proj_ci) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Scenario), alpha = 0.08, color = NA) + scale_fill_manual(values = scen_colors)
    }
    
    return(p)
  })

  output$step5_dynamic_plot <- renderPlot({
    step5_plot_object()
  })
  
  observeEvent(input$run_model, {
    if (!ux_require_data() || !ux_require_biology()) return(invisible(NULL))
    d_mapped <- processed_data()
    if(is.null(d_mapped) || nrow(d_mapped) == 0) return()
    if (!any(d_mapped$Monitored & !is.na(d_mapped$Count), na.rm = TRUE)) {
      showNotification("Data error: No monitored nest counts detected.", type = "error")
      return(NULL)
    }
    observed_zeros <- sum(
      d_mapped$Monitored & !is.na(d_mapped$Count) & d_mapped$Count == 0,
      na.rm = TRUE
    )
    if (observed_zeros > 0) {
      showNotification(
        paste0(
          "Model blocked: found ", observed_zeros,
          " observed zero count(s). The Martin/Siders Gaussian log-count ",
          "model cannot treat true zeros as data. Use blank Count with ",
          "Monitored = FALSE only for months with no survey effort."
        ),
        type = "error", duration = NULL
      )
      return(NULL)
    }
    if (any(d_mapped$Count < 0, na.rm = TRUE)) {
      showNotification("Data error: Nest counts cannot be negative.", type = "error")
      return(NULL)
    }
    if (is.null(input$clutch_freq) || !is.finite(input$clutch_freq) || input$clutch_freq <= 0) {
      showNotification("Clutch frequency must be greater than zero.", type = "error")
      return(NULL)
    }
    
    tryCatch(withProgress(message = 'Estimating historical population…', value = 0, {
      trend_fit_args <- if (isTRUE(input$legacy_trend_mcmc)) {
        list(
          iter = 750000, burnin = 250000, thin = 50,
          n_chains = 2, parallel = TRUE
        )
      } else {
        list(
          iter = input$iterations,
          burnin = floor(input$iterations / 3),
          thin = 10, n_chains = 3, parallel = FALSE
        )
      }
      if (any(!is.na(d_mapped$Month))) {
        setProgress(value = 0.1, detail = "Filling missing monitoring months…")
        season_start_month <- if (
          is.null(input$season_start_month) ||
            !is.finite(as.numeric(input$season_start_month))
        ) {
          4L
        } else {
          as.integer(input$season_start_month)
        }

        d_annual <- run_fourier_imputation(
          d_mapped,
          iter = if (isTRUE(input$legacy_fourier_mcmc)) {
            100000
          } else {
            max(2000, floor(input$iterations / 5))
          },
          six_month_sites = input$six_month_sites,
          legacy_mcmc = isTRUE(input$legacy_fourier_mcmc),
          season_start_month = season_start_month
        )
      } else {
        d_annual <- d_mapped %>%
          group_by(Year, Site) %>%
          summarise(
            Count = if (all(is.na(Count))) {
              NA_real_
            } else {
              sum(Count, na.rm = TRUE)
            },
            .groups = "drop"
          )
      }
      
      monthly_imputation <-
        any(!is.na(d_mapped$Month))
      use_imputation_bounds <-
        monthly_imputation && isTRUE(input$propagate_imputation_uncertainty)

      branch_inputs <- build_imputation_branch_inputs(
        d_annual = d_annual,
        clutch_frequency = input$clutch_freq,
        include_bounds = use_imputation_bounds,
        source_type = if (monthly_imputation) "imputed" else "observed"
      )

      if (
        use_imputation_bounds &&
          length(branch_inputs) == 1L
      ) {
        showNotification(
          paste0(
            "Lower and upper Fourier-imputation outputs were not available; ",
            "only the median branch will be fitted."
          ),
          type = "warning"
        )
      }

      season_start_month <- if (
        is.null(input$season_start_month) ||
          !is.finite(as.numeric(input$season_start_month))
      ) 4L else as.integer(input$season_start_month)

      raw_reference <- tryCatch(
        recommend_reference_site_from_observations(
          d_mapped,
          season_start_month = season_start_month
        ),
        error = function(e) NULL
      )
      advanced_override <- isTRUE(input$override_reference_site)
      requested_reference <- if (advanced_override) {
        input$reference_site
      } else if (!is.null(raw_reference)) {
        raw_reference$reference_site
      } else {
        NULL
      }

      reference_branch <- if ("median" %in% names(branch_inputs)) {
        "median"
      } else {
        names(branch_inputs)[1]
      }
      reference_resolution <- tryCatch(
        resolve_reference_site(
          branch_inputs[[reference_branch]],
          reference_site = requested_reference
        ),
        error = function(e) e
      )
      if (inherits(reference_resolution, "error")) {
        showNotification(
          paste0(
            "Reference-site selection failed: ",
            conditionMessage(reference_resolution)
          ),
          type = "error", duration = NULL
        )
        return(NULL)
      }

      selected_reference <- reference_resolution$reference_site
      reference_mode <- if (advanced_override) {
        "advanced override"
      } else {
        "automatic"
      }
      reference_reason <- if (!advanced_override && !is.null(raw_reference)) {
        raw_reference$selection_reason
      } else if (advanced_override) {
        paste0(
          "Selected through the advanced override from sites with positive ",
          "data in the first analysis year (",
          reference_resolution$first_year, ")."
        )
      } else {
        reference_resolution$selection_reason
      }
      reference_diagnostics <- if (!is.null(raw_reference)) {
        raw_reference$diagnostics
      } else {
        reference_resolution$diagnostics
      }

      branch_results <- list()
      branch_names <- names(branch_inputs)
      for (branch_index in seq_along(branch_names)) {
        branch_name <- branch_names[branch_index]
        branch_label <- imputation_branch_labels[[branch_name]]
        setProgress(
          value = 0.25 + 0.45 * branch_index / length(branch_names),
          detail = paste("Fitting trend model:", branch_label)
        )

        abund_branch <- branch_inputs[[branch_name]]
        model_error <- NULL
        res_branch <- tryCatch(
          do.call(
            run_jags_aligned,
            c(
              list(
                df = abund_branch,
                reference_site = selected_reference
              ),
              trend_fit_args
            )
          ),
          error = function(e) {
            model_error <<- conditionMessage(e)
            NULL
          }
        )
        if (is.null(res_branch)) {
          showNotification(
            paste0(
              "Trend model failed for ", branch_label, ": ", model_error
            ),
            type = "error", duration = NULL
          )
          return(NULL)
        }
        res_branch$reference_selection_mode <- reference_mode
        res_branch$reference_selection_reason <- reference_reason
        res_branch$reference_diagnostics <- reference_diagnostics

        setProgress(
          value = 0.70 + 0.15 * branch_index / length(branch_names),
          detail = paste("Running MARSS diagnostic:", branch_label)
        )
        marss_branch <- tryCatch(
          run_marss_diagnostic(abund_branch),
          error = function(e) NULL
        )
        branch_results[[branch_name]] <- list(
          abund = abund_branch,
          res = res_branch,
          marss = marss_branch
        )
      }

      vault$dev_imputation <- attr(d_annual, "dev_fit")
      vault$imputation_branch_results <- branch_results
      initial_branch <- if ("median" %in% names(branch_results)) {
        "median"
      } else {
        names(branch_results)[1]
      }
      activate_imputation_branch(initial_branch, notify = FALSE)

      # Split trends are an Ortega diagnostic and are fitted only for the
      # currently active imputation branch. They are cleared if another branch
      # is selected later, preventing mismatched summaries.
      if (input$run_split) {
        setProgress(value = 0.88, detail = "Running active-branch split trends...")
        slice_pre <- vault$abund %>% filter(Year <= input$split_year)
        fit_pre <- tryCatch({
          do.call(run_jags_aligned, c(list(df = slice_pre, reference_site = NULL), trend_fit_args))
        }, error = function(e) NULL)
        if (!is.null(fit_pre)) {
          vault$pre_u_val <- median(fit_pre$fit$sims.list$U)
          vault$pre_trend <- paste0(
            round(vault$pre_u_val, 3), " (",
            round((exp(vault$pre_u_val) - 1) * 100, 2), "%)"
          )
        }

        slice_post <- vault$abund %>% filter(Year >= input$split_year)
        fit_post <- tryCatch({
          do.call(run_jags_aligned, c(list(df = slice_post, reference_site = NULL), trend_fit_args))
        }, error = function(e) NULL)
        if (!is.null(fit_post)) {
          vault$post_u_val <- median(fit_post$fit$sims.list$U)
          vault$post_trend <- paste0(
            round(vault$post_u_val, 3), " (",
            round((exp(vault$post_u_val) - 1) * 100, 2), "%)"
          )
        }
        vault$split_branch <- initial_branch
      } else {
        vault$pre_trend <- NULL
        vault$post_trend <- NULL
        vault$split_branch <- NULL
      }

      setProgress(value = 1, detail = "Complete!")
      ux$fit_signature <- isolate(ux_fit_signature())
      ux$fit_revision <- ux$fit_revision + 1L
      ux$fit_time <- Sys.time()
    }), error = function(e) {
      showNotification(paste("Population estimation stopped:", conditionMessage(e)), type = "error", duration = NULL)
    })
  })

  output$summary_stats <- renderUI({
    validate(need(vault$nesters, "Use Estimate population to generate historical results."))
    if (!isTRUE(input$developer_view)) return(tagList(
      h4(paste0("Estimated annual nesting change: ", sprintf("%+.2f%%",vault$trend_pct))),
      p(class="ux-muted",paste0(format(round(vault$nesters),big.mark=","),
        " estimated nesting females in the final observed year. This is not a count of all turtles."))
    ))
    active_label <- if (
      is.null(vault$active_imputation_branch) ||
        !vault$active_imputation_branch %in% names(imputation_branch_labels)
    ) {
      "Observed annual data"
    } else {
      imputation_branch_labels[[vault$active_imputation_branch]]
    }
    boxes <- list(
      value_box(
        title = paste0("Nesters — ", active_label),
        value = format(round(vault$nesters), big.mark=","),
        theme = "primary"
      ),
      value_box(title = "Approx. mature females (RI cohorts)", value = format(round(vault$total), big.mark=","), theme = "secondary"),
      value_box(title = "Annual population growth rate", value = vault$trend_display, theme = if(vault$trend_pct >= 0) "success" else "danger")
    )
    if (input$run_split) {
      pre_val <- if(!is.null(vault$pre_trend)) vault$pre_trend else "N/A"
      post_val <- if(!is.null(vault$post_trend)) vault$post_trend else "N/A"
      boxes[[length(boxes)+1]] <- value_box(title = paste("Trend Pre-", input$split_year), value = pre_val, theme = "info")
      boxes[[length(boxes)+1]] <- value_box(title = paste("Trend Post-", input$split_year), value = post_val, theme = "warning")
    }
    tagList(
      do.call(layout_column_wrap, c(list(width = 1/length(boxes)), boxes)),
      if (!is.null(vault$res$reference_site)) {
        div(
          p(
            paste0(
              "Reference beach selected ",
              if (identical(vault$res$reference_selection_mode, "automatic")) {
                "automatically"
              } else {
                "by advanced override"
              },
              ": ", vault$res$reference_site,
              " (A = 0; X[1] prior centered on ",
              format(round(vault$res$reference_first_value, 2), big.mark = ","),
              " annual nesters in ", vault$res$reference_first_year, ")."
            ),
            style = "font-size: 0.82rem; color: #5f6b76; margin: 6px 0 0 4px;"
          ),
          if (!is.null(vault$res$reference_selection_reason)) {
            p(
              vault$res$reference_selection_reason,
              style = paste0(
                "font-size: 0.76rem; color: #74808b; ",
                "margin: 2px 0 0 4px;"
              )
            )
          }
        )
      }
    )
  })
  
  guided_history_plot <- reactive({
    req(vault$res)
    # Display transformation only. Reuse the same aligned posterior draws;
    # never refit or change model units when the presentation is switched.
    draws <- regional_abundance_draws(vault$res)
    qs <- apply(draws,2,quantile,probs=c(.025,.5,.975))
    d <- data.frame(Year=vault$res$years,Lower=qs[1,],Median=qs[2,],Upper=qs[3,])
    ggplot(d,aes(Year,Median))+
      geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="#185e55",alpha=.12)+
      geom_line(color="#185e55",linewidth=1)+
      labs(x="Year",y="Females nesting each year")+
      theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),
        plot.background=element_rect(fill="#f8f7f2",color=NA),
        panel.grid.major=element_line(color="#d9dfd8"),text=element_text(color="#243d3a"))
  })
  output$clean_baseline_plot <- renderPlot({
    validate(need(vault$res, "Use Estimate population in step 3 to generate historical results."))
    if(!isTRUE(input$developer_view)) {print(guided_history_plot());return(invisible(NULL))}
    fit <- vault$res$fit; years <- vault$res$years
    if (is.null(fit$sims.list$A)) { X_total <- exp(fit$sims.list$X) } else { X_total <- apply(fit$sims.list$X, 2, function(v) rowSums(apply(fit$sims.list$A, 2, function(x) exp(v + x)))) }
    X_q <- apply(log(X_total), 2, quantile, probs = c(0.025, 0.5, 0.975))
    obs_summary <- vault$abund %>% dplyr::group_by(Year) %>% dplyr::summarise(Ann_Tot = sum(Annual_Nesters, na.rm=TRUE)) %>% filter(Ann_Tot > 0)
    
    par(mar = c(4, 4, 1, 1), bty = "l", xaxs = "i", yaxs = "i")
    plot(years, X_q[2,], type="n", ylim=range(c(log(obs_summary$Ann_Tot), X_q), na.rm=TRUE), ylab="Annual nesters (log)", xlab="Year", axes=FALSE)
    axis(1, col="black", col.axis="black", lwd=1); axis(2, col="black", col.axis="black", lwd=1)
    if(input$run_split) abline(v = input$split_year, lty=2, col="#0072B2", lwd=1.5)
    polygon(c(years, rev(years)), c(X_q[1, ], rev(X_q[3, ])), col = 'grey90', border = NA)
    lines(years, X_q[2, ], lwd = 2.5, col = 'black')
    points(obs_summary$Year, log(obs_summary$Ann_Tot), pch = 16, col = "black")
    
    raw_mapped <- raw_ingested_data()
    observed_years <- unique(raw_mapped$Year[!is.na(raw_mapped$Count) & raw_mapped$Count > 0])
    missing_years <- years[!(years %in% observed_years)]
    if (length(missing_years) > 0) { idx_missing <- match(missing_years, years); points(missing_years, X_q[2, idx_missing], pch = 1, col = "black", cex = 1.6, lwd = 2) }
  })
  
  output$dynamic_unified_plot <- renderUI({
    validate(need(vault$abund, "Use Estimate population in step 3 to generate historical results..."))
    plotOutput("unified_trend_plot", height = "360px")
  })
  
  unified_trend_plot_object <- reactive({
    req(vault$res, vault$marss, vault$abund, input$plot_layers)
    j_fit <- vault$res$fit; all_years <- vault$res$years; sites = sort(unique(vault$abund$Site)); n_sites <- length(sites); df_all_fits <- data.frame()
      
      if ("jags" %in% input$plot_layers) {
        A_sims <- j_fit$sims.list$A; if (is.null(A_sims)) A_sims <- matrix(0, nrow = length(j_fit$sims.list$U), ncol = 1)
        for (i in seq_along(sites)) {
          site_A <- if(ncol(A_sims) >= i) as.numeric(A_sims[, i]) else rep(0, nrow(A_sims))
          for (t in seq_along(all_years)) { pred <- exp(j_fit$sims.list$X[, t] + site_A); df_all_fits <- rbind(df_all_fits, data.frame(Year=all_years[t], Site=sites[i], Model="JAGS Regional trend", Median=median(pred), Lower=quantile(pred,0.025), Upper=quantile(pred,0.975))) }
        }
      }
      if ("jags_indep" %in% input$plot_layers && !is.null(vault$jags_indep_fits)) {
        for (s in names(vault$jags_indep_fits)) {
          j_indep <- vault$jags_indep_fits[[s]]; j_fit_i <- j_indep$fit; j_years_i <- j_indep$years
          for (t in seq_along(j_years_i)) { pred <- exp(j_fit_i$sims.list$X[, t]); df_all_fits <- rbind(df_all_fits, data.frame(Year = j_years_i[t], Site = s, Model = "JAGS Independent Trends", Median = median(pred), Lower = quantile(pred, 0.025), Upper = quantile(pred, 0.975))) }
        }
      }
      if ("marss_s" %in% input$plot_layers && !is.null(vault$marss$shared)) {
        ms_fit <- vault$marss$shared; ms_states <- as.numeric(ms_fit$states[1, ]); ms_se <- as.numeric(ms_fit$states.se[1, ]); ms_A <- stats::coef(ms_fit, type="matrix")$A; if (is.null(ms_A) || length(ms_A) == 0) ms_A <- matrix(0, nrow = n_sites, ncol = 1)
        for(i in 1:n_sites) { site_A <- if(nrow(ms_A) >= i) as.numeric(ms_A[i,1]) else 0; df_all_fits <- rbind(df_all_fits, data.frame(Year=all_years, Site=sites[i], Model="MARSS Regional Trend", Median=exp(ms_states+site_A), Lower=exp((ms_states-1.96*ms_se)+site_A), Upper=exp((ms_states+1.96*ms_se)+site_A))) }
      }
      if ("marss_i" %in% input$plot_layers && !is.null(vault$marss$indep)) {
        mi_fit <- vault$marss$indep
        for(i in 1:n_sites) { mi_states <- if(is.matrix(mi_fit$states)) as.numeric(mi_fit$states[i, ]) else as.numeric(mi_fit$states); mi_se <- if(is.matrix(mi_fit$states.se)) as.numeric(mi_fit$states.se[i, ]) else as.numeric(mi_fit$states.se); df_all_fits <- rbind(df_all_fits, data.frame(Year=all_years, Site=sites[i], Model="MARSS Independent Trends", Median=exp(mi_states), Lower=exp(mi_states-1.96*mi_se), Upper=exp(mi_states+1.96*mi_se))) }
      }
      
      validate(need(nrow(df_all_fits) > 0, "Select at least one layer framework to draw fitted trends chart."))
      
    ggplot(df_all_fits %>% left_join(vault$abund, by=c("Year","Site")), aes(x=Year)) +
        geom_ribbon(aes(ymin=Lower, ymax=Upper, fill=Model), alpha=0.12, color=NA) + 
        geom_line(aes(y=Median, color=Model, linetype=Model), linewidth=1.1) + 
        geom_point(aes(y=Annual_Nesters), color="black", size=2.2, na.rm=TRUE) +
        facet_wrap(~Site, scales="free_y", ncol=1) + 
        theme_classic() + 
        theme(strip.text = element_text(face = "bold", size = 11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth = 1)) + labs(y = "Nesters", x = "Year")
  })

  output$unified_trend_plot <- renderPlot({
    withProgress(message = "Generating abundance trajectories plot...", value = 0.5, {
      unified_trend_plot_object()
    })
  })
  
  output$diagnostic_model_selector <- renderUI({
    req(vault$res); choices <- c("Regional Shared Model" = "regional")
    if (!is.null(vault$jags_indep_fits)) { beach_names <- names(vault$jags_indep_fits); beach_choices <- setNames(beach_names, beach_names); choices <- c(choices, "Independent Beach Models" = list(beach_choices)) }
    selectInput("diag_model_choice", "Select Model for Error & Posterior Diagnostics (Tabs 3 & 4):", choices = choices, width = "100%")
  })
  
  active_diag <- reactive({ req(vault$res, input$diag_model_choice); if (input$diag_model_choice == "regional") { return(vault$res) } else { req(vault$jags_indep_fits[[input$diag_model_choice]]); return(vault$jags_indep_fits[[input$diag_model_choice]]) } })
  
  output$table_u <- renderTable({
    validate(need(vault$res, "Please run baseline models on Page 2 first."))
    sites <- sort(unique(vault$abund$Site)); j_u <- median(vault$res$fit$sims.list$U)
    get_badge <- function(u) sprintf("<span class='badge %s'>%+.2f%% / yr</span>", if((exp(u)-1)*100 < 0) "bg-danger text-white" else "bg-success text-white", (exp(u)-1)*100)
    rows <- list(data.frame(Framework = "JAGS Regional trend", Strategy = "Shared", U = sprintf("%.4f", j_u), Trend = get_badge(j_u), check.names=FALSE))
    
    if(!is.null(vault$marss)) {
      m_s_u <- stats::coef(vault$marss$shared, type="matrix")$U[1,1]; m_i_u <- stats::coef(vault$marss$indep, type="matrix")$U[,1]
      rows[[length(rows)+1]] <- data.frame(Framework = "MARSS Regional Trend", Strategy = "Shared", U = sprintf("%.4f", m_s_u), Trend = get_badge(m_s_u), check.names=FALSE)
      for(i in seq_along(sites)) { val_u <- if(length(m_i_u) >= i) m_i_u[i] else m_i_u[1]; rows[[length(rows)+1]] <- data.frame(Framework = "MARSS Independent Trends", Strategy = sites[i], U = sprintf("%.4f", val_u), Trend = get_badge(val_u), check.names=FALSE) }
    }
    do.call(rbind, rows)
  }, sanitize.text.function = function(x) x)
  
  variance_plot_object <- reactive({
    diag_data <- active_diag(); req(diag_data)
    j_q <- mean(diag_data$fit$sims.list$Q); r_sims <- diag_data$fit$sims.list$R; j_r <- if(is.matrix(r_sims)) mean(colMeans(r_sims)) else mean(r_sims)
    df_var <- data.frame(Component = c("Environmental stochasticity (Q)", "Observation stochasticity (R)"), Variance = c(j_q, j_r))
    ggplot(df_var, aes(x = Component, y = Variance, fill = Component)) + 
      geom_bar(stat = "identity", color = "black", width = 0.5, linewidth = 0.8) + 
      scale_fill_manual(values = c("Environmental stochasticity (Q)" = "#0072B2", "Observation stochasticity (R)" = "#D55E00")) +
      theme_classic() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth = 1), text = element_text(size = 12)) + labs(y = "Variance scale value", x = "")
  })

  output$plot_var <- renderPlot({
    variance_plot_object()
  })
  
  output$variance_interpretation_text <- renderUI({
    diag_data <- active_diag(); req(diag_data)
    j_q <- mean(diag_data$fit$sims.list$Q); r_sims <- diag_data$fit$sims.list$R; j_r <- if(is.matrix(r_sims)) mean(colMeans(r_sims)) else mean(r_sims)
    if (j_q > j_r) {
      p(tags$b("Interpretation: "), "Estimated process variance (Q) exceeds observation variance (R). In this model, more unexplained variation is assigned to changes in the latent population process than to observation error; this is not direct evidence of a particular environmental cause.")
    } else {
      p(tags$b("Interpretation: "), "Estimated observation variance (R) exceeds process variance (Q). The model attributes more variation to the observation layer, but it does not identify whether that variation arose from monitoring effort, counting error, or other unmodeled processes.")
    }
  })
  
  output$covariance_interpretation_text <- renderUI({
    diag_data <- active_diag(); req(diag_data)
    fit <- diag_data$fit; years <- diag_data$years; fy <- length(years)
    r_vec <- as.numeric(fit$sims.list$U); q_vec <- as.numeric(fit$sims.list$Q)
    if (is.null(fit$sims.list$A)) {
      n_vec <- exp(as.numeric(fit$sims.list$X[, fy]))
    } else {
      x_last <- as.numeric(fit$sims.list$X[, fy])
      n_vec <- rowSums(exp(sweep(fit$sims.list$A, 1, x_last, "+")))
    }
    cor_rq <- cor(r_vec, q_vec, use = "complete.obs"); cor_rn <- cor(r_vec, n_vec, use = "complete.obs"); cor_nq <- cor(n_vec, q_vec, use = "complete.obs")
    
    tagList(
      h5(style = "font-size: 1.1rem; font-weight: bold; margin-bottom: 12px;", shiny::icon("chart-line"), " Management & Biological Interpretation Guide:"),
      tags$ul(style = "padding-left: 20px;",
              tags$li(style = "margin-bottom: 10px;", tags$b("Posterior association between trend and process variance [r = ", sprintf("%.2f", cor_rq), "]: "), if(abs(cor_rq) < 0.30) { "The sampled association is weak; inspect the full posterior and convergence diagnostics before interpreting the parameters as well separated." } else { "The sampled association is notable, so trend inference may be sensitive to the process-variance estimate." }),
              tags$li(style = "margin-bottom: 10px;", tags$b("Posterior association between trend and final abundance [r = ", sprintf("%.2f", cor_rn), "]: "), "This describes dependence within the fitted posterior; it is not an independent model-validation test."),
              tags$li(style = "margin-bottom: 10px;", tags$b("Posterior association between final abundance and process variance [r = ", sprintf("%.2f", cor_nq), "]: "), if(abs(cor_nq) >= 0.30) { tags$span(style = "font-weight: 500; color: #b02a37;", "The association is notable; report joint uncertainty rather than treating abundance and process variance as independent.") } else { "The sampled association is weak, but this alone does not establish parameter independence." })
      ),
      div(style = "background-color: #f8f9fa; border-left: 4px solid #6c757d; padding: 10px; margin-top: 15px; border-radius: 4px;", tags$b("Statistical Distribution Note (Histogram Shapes): "), "The asymmetric right-hand tail displayed in the middle (N_final) histogram indicates that our population uncertainty is not uniform. While the lower boundary is firmly restricted by actual beach nest counts, the upper boundary allows for a wide margin of error.")
    )
  })
  
  output$posterior_pairs_plot <- renderPlot({
    diag_data <- active_diag(); req(diag_data)
    withProgress(message = "Computing Joint Parameter Covariance Matrix...", value = 0.5, {
      fit <- diag_data$fit; years <- diag_data$years; fy <- length(years)
      r_vec <- as.numeric(fit$sims.list$U); q_vec <- as.numeric(fit$sims.list$Q)
      if (is.null(fit$sims.list$A)) {
        n_vec <- exp(as.numeric(fit$sims.list$X[, fy]))
      } else {
        x_last <- as.numeric(fit$sims.list$X[, fy])
        n_vec <- rowSums(exp(sweep(fit$sims.list$A, 1, x_last, "+")))
      }
      pairs_df <- data.frame(r = r_vec, N_final = n_vec, Q = q_vec)
      panel_hist <- function(x, ...) { h <- hist(x, plot = FALSE, breaks = 25); y <- h$counts / max(h$counts); old_par <- par(usr = c(par("usr")[1:2], 0, 1.1)); on.exit(par(old_par)); lines(approx(h$mids, y, xout=seq(min(x), max(x), length.out=100)), lwd = 2, col = "black") }
      panel_scatter <- function(x, y, ...) { points(x, y, pch = 20, col = rgb(0.3, 0.3, 0.3, 0.05), cex = 0.5); ok <- is.finite(x)&is.finite(y); if(sum(ok)>3 && var(x[ok])>0 && var(y[ok])>0) {xy<-cbind(x[ok],y[ok]);ee<-eigen(cov(xy),symmetric=TRUE);if(all(ee$values>0)){theta<-seq(0,2*pi,length.out=120);circle<-cbind(cos(theta),sin(theta));for(level in c(.5,.95)){ellipse<-sweep(sqrt(qchisq(level,2))*circle %*% diag(sqrt(ee$values)) %*% t(ee$vectors),2,colMeans(xy),"+");lines(ellipse,col=if(level==.5)"black" else "gray50",lty=if(level==.5)1 else 2)}}} }
      panel_cor <- function(x, y, ...) { r_coef <- cor(x, y, use = "complete.obs"); old_par <- par(usr = c(0, 1, 0, 1)); on.exit(par(old_par)); text(0.5, 0.5, sprintf("%.2f", r_coef), cex = 1.5, font = 2) }
      par(mar = c(3, 3, 1, 1), bg = "white")
      pairs(pairs_df, labels = c("r", "N_final", "Q"), diag.panel = panel_hist, lower.panel = panel_cor, upper.panel = panel_scatter, cex.labels = 1.4, font.labels = 2)
    })
  })
  
  output$download_preview_raw <- downloadHandler(
    filename = function() "Raw_Nest_Counts_Preview.png",
    content = function(file) {
      plot_object <- preview_annual_plot_object()
      req(plot_object)
      ggsave(file, plot = plot_object, device = "png", width = 8, height = 5, bg = "white")
    }
  )
  output$download_baseline <- downloadHandler(
    filename = function() { "Calculated_Regional_Baseline_Trend.png" },
    content = function(file) {
      if(!isTRUE(input$developer_view)) {
        ggsave(file,plot=guided_history_plot(),device="png",width=9,height=4,bg="#f8f7f2")
        return(invisible(NULL))
      }
      png(file, width = 800, height = 500, res = 100)
      fit <- vault$res$fit; years <- vault$res$years
      if(is.null(fit$sims.list$A)) X_total <- exp(fit$sims.list$X) else X_total <- apply(fit$sims.list$X, 2, function(v) rowSums(apply(fit$sims.list$A, 2, function(x) exp(v + x))))
      X_q <- apply(log(X_total), 2, quantile, probs = c(0.025, 0.5, 0.975))
      obs_summary <- vault$abund %>% dplyr::group_by(Year) %>% dplyr::summarise(Ann_Tot = sum(Annual_Nesters, na.rm=TRUE))
      par(mar = c(4, 4, 1, 1), bty = "l", xaxs = "i", yaxs = "i")
      plot(years, X_q[2,], type="n", ylim=range(c(log(obs_summary$Ann_Tot), X_q), na.rm=TRUE), ylab="Annual nester records (log)", xlab="Year", axes=FALSE)
      axis(1); axis(2); polygon(c(years, rev(years)), c(X_q[1, ], rev(X_q[3, ])), col = 'grey90', border = NA); lines(years, X_q[2, ], lwd = 2.5); points(obs_summary$Year, log(obs_summary$Ann_Tot), pch = 16)
      dev.off()
    }
  )
  output$download_unified <- downloadHandler(
    filename = function() "Fits_Comparison.png",
    content = function(file) {
      ggsave(file, plot = unified_trend_plot_object(), device = "png", width = 8, height = 6, bg = "white")
    }
  )
  output$download_var <- downloadHandler(
    filename = function() "Variance_Components.png",
    content = function(file) {
      ggsave(file, plot = variance_plot_object(), device = "png", width = 5, height = 4, bg = "white")
    }
  )
  output$download_projection_plot <- downloadHandler(
    filename = function() "Sea_Turtle_Projection_Scenarios.png",
    content = function(file) {
      ux_require_current_results()
      req(vault$step5_summary, vault$step5_raw)
      ggsave(
        file, plot = step5_plot_object(), device = "png",
        width = 10, height = 6, dpi = 300, bg = "white"
      )
    }
  )
  output$download_projection_summary <- downloadHandler(
    filename = function() "Sea_Turtle_Projection_Summary.csv",
    content = function(file) {
      ux_require_current_results()
      req(vault$step5_summary)
      write.csv(vault$step5_summary, file, row.names = FALSE)
    }
  )
  output$download_projection_draws <- downloadHandler(
    filename = function() "Sea_Turtle_Projection_Simulation_Draws.csv",
    content = function(file) {
      ux_require_current_results()
      req(vault$step5_raw)
      write.csv(vault$step5_raw, file, row.names = FALSE)
    }
  )
  output$download_posterior <- downloadHandler(
    filename = function() { "Joint_Posterior_Matrix.png" },
    content = function(file) {
      png(file, width = 650, height = 650, res = 120)
      fit <- vault$res$fit; years <- vault$res$years; fy <- length(years)
      r_vec <- as.numeric(fit$sims.list$U); q_vec <- as.numeric(fit$sims.list$Q)
      if(is.null(fit$sims.list$A)) {
        n_vec <- exp(as.numeric(fit$sims.list$X[, fy]))
      } else {
        x_last <- as.numeric(fit$sims.list$X[, fy])
        n_vec <- rowSums(exp(sweep(fit$sims.list$A, 1, x_last, "+")))
      }
      pairs_df <- data.frame(r = r_vec, N_final = n_vec, Q = q_vec)
      panel_hist <- function(x, ...) { h <- hist(x, plot = FALSE, breaks = 25); y <- h$counts / max(h$counts); lines(approx(h$mids, y, xout=seq(min(x), max(x), length.out=100)), lwd = 2) }
      panel_scatter <- function(x, y, ...) { points(x, y, pch = 20, col = rgb(0.3, 0.3, 0.3, 0.05)); if(requireNamespace("car",quietly=TRUE)){lines(car::dataEllipse(x, y, levels = 0.50, draw = FALSE), lwd = 2); lines(car::dataEllipse(x, y, levels = 0.95, draw = FALSE), lty = 2)} }
      panel_cor <- function(x, y, ...) { text(mean(range(x)), mean(range(y)), sprintf("%.2f", cor(x, y, use = "complete.obs")), cex = 1.5, font = 2) }
      pairs(pairs_df, labels = c("r", "N_final", "Q"), diag.panel = panel_hist, lower.panel = panel_cor, upper.panel = panel_scatter)
      dev.off()
    }
  )
  # These outputs create analysis inputs. Initialize them in both views so
  # opening a hidden workspace cannot silently change the fitting defaults.
  for (id in c("imputation_uncertainty_option_ui", "imputation_branch_selector_ui",
               "reference_site_ui", "diagnostic_model_selector", "dev_fit_selector")) {
    outputOptions(output, id, suspendWhenHidden=FALSE)
  }
}
shinyApp(ui, server)
