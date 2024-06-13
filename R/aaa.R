# MAIN DOCUMENTATION PAGE -------------------------------------------------

#' The 'sfgp' package.
#'
#' @description SF+GP modeling using 'Stan'.
#' @author Juho Timonen (first.last at iki.fi)
#' @keywords tumor GP Stan Bayesian
#'
#' @section Getting started:
#' See the following \code{R6} classes.
#' \itemize{
#'  \item \code{\link{TSModel}}: Main model class.
#'  \item \code{\link{TSModelFit}}: Fit class.
#'  \item \code{\link{TermList}}: Class describing model terms.
#'  \item \code{\link{FunctionDraws}}: Class to hold fitted function
#'  distributions.
#' }
#'
#' @section Data:
#' The data that you wish to analyze with 'sfgp' should be in an \R
#' \code{data.frame} where columns correspond to measured variables and rows
#' correspond to observations.
#'
#' @section Vignettes:
#' \itemize{
#'  \item Simple GP model: How to create a fit a simple GP model.
#' }
#'
#' @name sfgp-package
#' @aliases sfgp
#' @import ggplot2 stats
#' @importFrom R6 R6Class
#' @importFrom posterior as_draws_rvars
#'
#'
"_PACKAGE"

#' A small artificial test data
#'
#' @format A data frame with 176 rows and 6 variables:
#' \describe{
#'   \item{id}{factor, individual id}
#'   \item{time}{numeric, a continuous variable}
#'   \item{arm}{factor, three levels}
#'   \item{sex}{factor, two levels}
#'   \item{weight}{numeric, a continuous variable}
#'   \item{y}{numeric, a continuous variable}
#'   \item{f_true}{numeric, the true signal used to generate the data}
#' }
#' @family built-in datasets
"testdata"


#' Run an example
#'
#' @description Fits a model to simple simulated data.
#' @export
#' @param num_bf Number of basis functions.
#' @param scale_bf Basis function domain scale.
#' @param formula The model formula.
#' @param ... Other arguments to the \code{$fit()} method of
#' \code{\link{TSModel}}.
#' @return An \code{\link{TSModelFit}} object.
example <- function(num_bf = 32, scale_bf = 1.5, formula = "y ~ gp(x)",
                    ...) {
  form <- stats::as.formula(formula)
  m <- TSModel$new(form)
  xx <- seq(1, 10, by = 0.15)
  ff <- 20 + 5 * sin(xx) + 2 * sin(5 * xx) + xx
  yy <- ff + stats::rnorm(n = length(xx), mean = 0, sd = 1)
  a <- data.frame(x = xx, y = yy)
  tc1 <- list(num_bf = num_bf, scale_bf = scale_bf)
  tc <- list(f_gp_x = tc1)
  m$fit(data = a, term_confs = tc, ...)
}

#' Run an example
#'
#' @description Fits a model to \code{testdata}.
#' @export
#' @param num_bf Number of basis functions.
#' @param scale_bf Basis function domain scale.
#' @param formula The model formula.
#' @param ... Other arguments to the \code{$fit()} method of
#' \code{\link{TSModel}}.
#' @return An \code{\link{TSModelFit}} object.
example2 <- function(num_bf = 24, scale_bf = 1.5,
                     formula = "y ~ gp(time) + gp(time,arm) + gp(time,id)",
                     ...) {
  form <- stats::as.formula(formula)
  dat <- testdata
  m <- TSModel$new(form)
  dat <- add_sff_input(dat, m)
  tc1 <- list(num_bf = num_bf, scale_bf = scale_bf)
  tc <- list(gp_x = tc1)
  m$fit(data = dat, term_confs = tc, ...)
}

#' Extend data frame
#'
#' @export
#' @param df original data frame
#' @param t vector of new time values
#' @param time_var name of time variable
#' @return new data frame
extend_df <- function(df, t, time_var) {
  u <- df_unique_factor_rows(df)
  df <- df_replicate_rows(u, length(t))
  df[[time_var]] <- rep(t, nrow(u))
  df
}

#' Visualize 'rvar' in data frame
#'
#' @export
#' @param df the data frame
#' @param group_var name of the grouping variable
#' @param rvar_name name of the 'rvar' column
#' @param horizontal orient horizontally?
#' @param halfeye use halfeye plot instead of pointinterval
#' @param ... arguments passed to ggdist
plot_rvar_df <- function(df, rvar_name, group_var, horizontal = FALSE,
                         halfeye = FALSE, ...) {
  # create plot
  if (horizontal) {
    plt <- ggplot2::ggplot(df, aes(
      y = !!sym(group_var),
      xdist = !!sym(rvar_name)
    ))
  } else {
    plt <- ggplot2::ggplot(df, aes(
      x = !!sym(group_var),
      ydist = !!sym(rvar_name)
    ))
  }

  # update
  if (halfeye) {
    plt <- plt + ggdist::stat_halfeye(...)
  } else {
    plt <- plt + ggdist::stat_pointinterval(...)
  }

  plt <- plt + ggtitle(rvar_name)
  if (horizontal) {
    plt <- plt + xlab("value")
  } else {
    plt <- plt + ylab("value")
  }

  # return
  plt
}

#' Add input for 'FormulaSFTerm'
#'
#' @description
#' Duplicates original data columns, adding new columns with the \code{_kg} and
#' \code{_ks} suffixes if these are missing.
#'
#' @export
#' @param df the data frame
#' @param model an object of class \code{\link{TSModel}}
add_sff_input <- function(df, model) {
  t <- get_sff_term(model)
  vars_ks <- as.vector(t$term_list_ks$input_vars())
  vars_kg <- as.vector(t$term_list_kg$input_vars())
  vars <- c(vars_kg, vars_ks)
  for (v in vars) {
    if (v %in% colnames(df)) {
      message("column ", hl_string(v), " already exists in the data frame")
    } else {
      parts <- strsplit(v, split = "_", fixed = TRUE)[[1]]
      v_base <- paste(parts[1:(length(parts) - 1)], sep = "_")
      new_col <- df[[v_base]]
      if (is.null(new_col)) {
        stop(hl_string(v_base), " not found in the data frame")
      }
      message(
        "copying column ", hl_string(v_base), " to new column ", hl_string(v)
      )
      df[[v]] <- new_col
    }
  }
  df
}

#' Create a dense grid until event time for each subject (for JointModel)
#'
#' @export
#' @param df_tte time-to-event data frame
#' @param df_lon the longitudinal data frame
#' @param id_var id variable
#' @param time_var time variable
#' @param num_points number of grid points
#' @param even_spacing space grid points evenly?
create_jm_grid <- function(df_tte, df_lon, num_points = 30, id_var = "id",
                           time_var = "t", even_spacing = FALSE) {
  t_min <- 1e-6 # because hazard might not be defined at t = 0
  id <- df_tte[[id_var]]
  time <- df_tte[[time_var]]
  checkmate::assert_factor(id)
  checkmate::assert_numeric(time)
  ids <- as.numeric(levels(id))
  df_grid <- NULL
  t_max <- max(time)
  if (even_spacing) {
    t_seq <- seq(0, t_max, length.out = num_points)
  } else {
    t_seq <- exp(seq(0, 3, length.out = num_points)) - 1
    t_seq <- t_max / max(t_seq) * t_seq
    t_seq[t_seq <= t_min] <- t_min
  }
  for (lev in ids) {
    idx_tte <- which(id == lev)
    idx_tte <- idx_tte[length(idx_tte)]
    idx_lon <- which(df_lon[[id_var]] == lev)[1]
    row_tte <- df_tte[idx_tte, , drop = FALSE] # take last row of the subject
    row_lon <- df_lon[idx_lon, , drop = FALSE] # take first row of the subject
    row_lon[[time_var]] <- NULL # remove old time variable
    row <- cbind(row_tte, row_lon)
    df_j <- row[rep(1, num_points), ] # repeat row max_num_points times
    df_j[[time_var]] <- t_seq
    df_grid <- rbind(df_grid, df_j)
  }
  df_grid
}


#' Name of a term in the formula syntax to its Stan code name
#'
#' @export
#' @param term A string.
term_to_code <- function(term) {
  checkmate::assert_character(term)
  form <- as.formula(paste0("y~", term))
  a <- create_termlist(form, NULL)
  a$terms[[1]]$stanname_base()
}

#' Expose functions in a given Stan file to R
#'
#' @description
#' Does not seem to work.
#'
#' @param stan_file Path to the Stan file relative to the directory
#' \code{inst/functions}, without the \code{.stan} suffix.
#' @param ... Optional arguments to \code{cmdstanr::cmdstan_model()}.
expose_stan_functions <- function(stan_file, ...) {
  checkmate::assert_character(stan_file)
  sc <- read_stan_function(stan_file)
  sc <- paste("functions {", sc, "}", sep = "\n")
  sf <- cmdstanr::write_stan_file(sc)
  a <- cmdstanr::cmdstan_model(sf, ...)
  a$expose_functions()
  a
}
