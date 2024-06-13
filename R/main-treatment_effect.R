#' Draw one subjects from each group
#'
#' @export
#' @param data original data frame
#' @param id_var name of the subject identifier variable
#' @param group_var grouping variable
#' @return new data frame
sample_subjects <- function(data, id_var = "id", group_var = "arm") {
  N <- 1 # had to be fixed to 1 due to id issues
  id_fac <- data[[id_var]]
  checkmate::assert_factor(id_fac)
  sample_rows_by_group(data, N = N, group_var = group_var)
}

#' Analyse treatment effect
#'
#' @export
#' @param fit_post Posterior fit.
#' @param fit_prior Prior fit. Needs to have same number of draws as
#' \code{fit_post}.
#' @param t_pred Vector of time points (in days). Chosen based on data range
#' if \code{t_pred=NULL}. Chosen to be an evenly spaced vector from
#' 0 to at most 2 years every 8 weeks if \code{t_pred="auto"}. If
#' \code{t_pred=NULL}, it will be changed to \code{"auto"} if \code{br=TRUE}.
#' @param time_var name of the time variable
#' @param group_var grouping variable
#' @param method Used method for analyzing the treatment effect. Can be either
#' \itemize{
#'  \item \code{new_sub}: Simulating new imaginary subjects
#'  \item \code{group_est}: Using group-level parameter estimates
#' }
#' @param br Use best response definition of depth? See the \code{dor}
#' method of the \code{FunctionDraws} class.
treatment_effect <- function(fit_post, fit_prior, time_var = "t",
                             t_pred = NULL, group_var = "arm",
                             method = "new_sub", br = FALSE) {
  # Validate input
  checkmate::assert_choice(method, choices = c("new_sub", "group_est"))
  if (is.null(t_pred) && br) {
    t_pred <- "auto"
  }

  # "New" subjects
  if (method == "new_sub") {
    zero_log_z <- FALSE
  } else if (method == "group_est") {
    zero_log_z <- TRUE
  } else {
    stop("invalid method")
  }
  p_new <- predict_new_subjects(
    fit_post, fit_prior, time_var, t_pred, group_var, zero_log_z
  )
  traj <- p_new$function_draws(data_scale = F)$zero_offset_by_id(group_var)

  # Return
  traj_log <- traj
  traj <- traj$exp()
  list(
    traj_log = traj_log,
    traj = traj,
    depth = traj$dor(group_var = group_var, br = br),
    duration = traj$dur(group_var = group_var, time_var = time_var),
    duration_responding = traj$dur(
      group_var = group_var,
      time_var = time_var,
      only_responding = TRUE
    ),
    p_new = p_new
  )
}



#' Generate predictions for new imaginary subjects from each group
#'
#' @export
#' @param fit_post posterior fit
#' @param fit_prior prior fit
#' @param t_pred Vector of time points. Chosen automatically if \code{NULL}.
#' @param time_var time variable
#' @param group_var grouping variable
#' @param zero_log_z Use zero for all \code{log_z_*} parameters. This
#' effectively takes group means of all \code{HierOffsetTerm}s.
#' @param prior_param_names Names of parameters for which prior draws should
#' be used. If \code{NULL} (default), these are attempted be to detected
#' automatically. These are usually all individual-specific parameters whose
#' values we do not know for "new" subjects so we draw them from the prior.
predict_new_subjects <- function(fit_post, fit_prior,
                                 time_var, t_pred = NULL,
                                 group_var = "arm",
                                 zero_log_z = FALSE,
                                 prior_param_names = NULL) {
  checkmate::assert_class(fit_post, "TSModelFit")
  checkmate::assert_class(fit_prior, "TSModelFit")
  t_data <- fit_post$get_data("LON")[[time_var]]
  t_pred <- t_pred_auto(t_pred, t_data)

  # Pick one subject from each group
  newdat <- group_pred_input(fit_post, group_var, t_pred, time_var)

  # Take draws
  if (is.null(prior_param_names)) {
    prior_param_names <- id_specific_param_names(fit_post)
  }
  d <- fit_post$draws()
  d_prior <- fit_prior$draws()
  for (par in prior_param_names) {
    message("using prior draws of ", hl_string(par))
    d[[par]] <- d_prior[[par]]
  }

  # Edit draws
  if (zero_log_z) {
    d <- new_draws(d, fit_post$log_z_pars(), 0)
  }

  # Predict using partly posterior and partly prior draws
  fit_post$predict(newdat, fitted_params = d)
}

# Helper
id_specific_param_names <- function(fit) {
  id_var <- fit$get_model("LON")$id_var
  pn <- names(fit$draws())
  pat <- c(
    paste0("^log_z_offset_", id_var), # HierOffsetTerm
    paste0("^offset_", id_var), # GroupedOffsetTerm
    paste0("^xi_.*", id_var, ".*") # GPTerm
  )
  pat <- paste(pat, collapse = "|")
  idx <- grepl(pn, pattern = pat)
  pn[which(idx)]
}


# Helper
group_pred_input <- function(fit, group_var, t_pred, time_var) {
  model <- fit$get_model("LON")
  id_var <- model$id_var
  newsubs <- sample_subjects(fit$get_data("LON"), id_var, group_var)
  newdat <- extend_df(newsubs, t_pred, time_var)
  newdat[[model$y_var]] <- 1 # value has no effect
  newdat
}




#' Visualize depth of response
#'
#' @export
#' @param te A list returned by \code{\link{treatment_effect}}.
#' @inheritParams plot_rvar_df
#' @inheritParams plot_dur
#' @param ... arguments passed to ggdist
plot_dor <- function(te, halfeye = FALSE, control_group_name = NULL,
                     ...) {
  depth <- te$depth
  if (!is.null(control_group_name)) {
    depth <- dor_diff(depth, control_group_name)
    lab <- "value - control"
  } else {
    lab <- "value"
  }
  plot_rvar_arm(depth, "depth",
    horizontal = FALSE,
    halfeye = halfeye, ...
  ) + ylab(lab)
}

#' Visualize duration of response
#'
#' @export
#' @param te A list returned by \code{\link{treatment_effect}}.
#' @param only_responding Compute only for the trajectories that have response.
#' Currently not implemented!
#' @inheritParams plot_rvar_df
#' @param metric Can be one of \code{"duration"}, \code{"start"}, \code{"end"},
#' or \code{"TTB120"}.
#' @param control_group_name If given, the values are differences to control,
#' or proportion of samples larger than control (out of non-equal samples) if
#' \code{gt = TRUE}.
#' @param ... arguments passed to \code{ggdist}.
plot_dur <- function(te, only_responding = FALSE,
                     metric = "duration", halfeye = FALSE,
                     control_group_name = NULL, ...) {
  gtr <- function(x) unlist_rvar(lapply(x, function(x) x[metric]))
  ok <- c("duration", "start", "end", "TTB120")
  checkmate::assert_choice(metric, choices = ok)
  if (only_responding) {
    stop("not implemented")
    # does not work
    # plotter <- function(x) plot_rvar_arm(gtr(x$duration_responding), "duration")
  }
  dur <- te$duration
  if (!is.null(control_group_name)) {
    dur <- dur_diff(dur, control_group_name)
    lab <- "value - control"
  } else {
    lab <- "value"
  }
  plot_rvar_arm(gtr(dur), metric,
    horizontal = TRUE, halfeye = halfeye, ...
  ) + xlab(lab)
}

#' Get duration differences to control
#'
#' @export
#' @param dur List of duration \code{rvar} vectors.
#' @param control_group_name Control arm name.
#' @param gt See the \code{control_group_name} argument.
#' @return \itemize{
#' \item A list with \code{rvar}s that are differences to control, or
#' \item if \code{gt = TRUE}, matrix with number of rows equal to number of
#' arms, where elements are for the proportion of samples larger than control
#' (out of non-equal samples)
#' }
dur_diff <- function(dur, control_group_name = NULL, gt = FALSE) {
  if (is.null(control_group_name)) {
    control_group_name <- names(dur)[1]
  }
  dur_diff <- list()
  if (gt) {
    dur_diff <- c()
  }
  for (d in 1:length(dur)) {
    a <- dur[[d]]
    b <- dur[[control_group_name]]
    if (gt) {
      r <- sum(a > b) / sum(a != b)
    } else {
      r <- a - b
    }
    if (gt) {
      dur_diff <- rbind(dur_diff, r)
    } else {
      dur_diff[[d]] <- r
    }
  }
  if (gt) {
    rownames(dur_diff) <- names(dur)
  } else {
    names(dur_diff) <- names(dur)
  }
  dur_diff
}

#' Get duration differences to control
#'
#' @export
#' @param dor Vector of depth \code{rvar}s.
#' @inheritParams dur_diff
dor_diff <- function(dor, control_group_name = NULL, gt = FALSE) {
  if (is.null(control_group_name)) {
    control_group_name <- names(dor)[1]
  }
  a <- dor
  b <- dor[[control_group_name]]
  if (gt) {
    diff <- sum(a > b) / sum(a != b)
  } else {
    diff <- a - b
  }
  diff
}

#' Summarize duration of response
#'
#' @export
#' @inheritParams plot_dur
summarize_dur <- function(te, only_responding = FALSE) {
  a <- te$duration
  if (only_responding) {
    a <- te$duration_responding
  }
  df <- data.frame(sapply(a, function(x) as.character(x)))
  colnames(df) <- names(a)
  df
}
