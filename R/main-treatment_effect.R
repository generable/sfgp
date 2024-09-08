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
#' 0 to at most 2 years every 8 weeks if \code{t_pred="auto"}.
#' @param time_var name of the time variable
#' @param group_var grouping variable
#' @param method Used method for analyzing the treatment effect. Can be either
#' \itemize{
#'  \item \code{new_sub}: Simulating new imaginary subjects
#'  \item \code{group_est}: Using group-level parameter estimates
#' }
treatment_effect <- function(fit_post, fit_prior, time_var = "t",
                             t_pred = NULL, group_var = "arm",
                             method = "new_sub") {
  # Validate input
  checkmate::assert_choice(method, choices = c("new_sub", "group_est"))

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
    metrics = trajectory_metrics(traj, group_var, time_var),
    p_new = p_new
  )
}


#' Compute trajectory metrics.
#'
#' @export
#' @param trajectories An object of class \code{\link{FunctionDraws}}.
#' @param id_var Name of id variable.
#' @param time_var Name of the time variable.
#' @return The following metrics are computed for each draw for each id.
#' \itemize{
#'  \item \code{depth} = depth of response (RECIST)
#'  \item \code{depth_br} = depth of response (best response definition, can
#'  be negative)
#'  \item \code{duration} = duration of response (RECIST)
#'  \item \code{start} = start time of response (RECIST)
#'  \item \code{end} = end time of response (RECIST)
#'  \item \code{exists} = whether response exists (RECIST)
#'  \item \code{endures} = whether response endures whole time interval (RECIST)
#'  \item \code{ttb120} = Time to 1.2 x baseline (different duration metric).
#' }
#' See the \code{treatment-effect} vignette for definition of the metrics.
trajectory_metrics <- function(trajectories, id_var, time_var) {
  checkmate::assert_class(trajectories, "FunctionDraws")
  checkmate::assert_character(id_var, len = 1)
  # use first provided time_var
  time_var <- time_var[[1]]
  checkmate::assert_character(time_var, len = 1)
  df <- trajectories$as_data_frame_long()
  df %>%
    dplyr::group_by(!!sym(id_var), .draw_idx) %>%
    dplyr::summarize(
      depth = depth_of_response(value),
      depth_br = best_response(value),
      duration = response_dur(!!sym(time_var), value),
      start = response_start(!!sym(time_var), value),
      end = response_end(!!sym(time_var), value),
      exists = response_exists(!!sym(time_var), value),
      endures = response_endures(!!sym(time_var), value),
      ttb120 = ttb120(!!sym(time_var), value),
      .groups = "drop"
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
  t_data <- fit_post$get_data("LON")[[time_var[[1]]]]
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

#' Example implementation of plotting trajectory metrics.
#'
#' @export
#' @param df A full long data frame returned by
#' \code{\link{trajectory_metrics}}.
#' @param group_var Name of grouping variable.
#' @param metric Name of the metric. See possible metric names in
#' documentation of \code{\link{trajectory_metrics}}.
plot_metric <- function(df, metric = "depth", group_var = "arm") {
  ggplot(df, aes(x = !!sym(group_var), y = !!sym(metric))) +
    ggdist::stat_slabinterval()
}
