#' The time-to-event model fit class
#'
#' @export
JointModelFit <- R6::R6Class("JointModelFit",
  inherit = TSModelFit,
  public = list(

    #' @description Print object description.
    print = function() {
      cat("An R6 object of class JointModelFit.\n")
    },

    #' @description
    #' Get the entire joint model or a particular submodel.
    #'
    #' @param name Name of the submodel. If \code{NULL} (default), the entire
    #' joint model object is returned.
    get_model = function(name = NULL) {
      if (is.null(name)) {
        return(private$model)
      }
      private$model[[name]]
    },

    #' @description
    #' Extract one TTE-related function as a \code{\link{FunctionDraws}} object.
    #'
    #' @param component Base name of the function in the 'Stan' code.
    #' @param dataname Data where the function was evaluated.
    function_draws_tte = function(component = "inst_haz", dataname = NULL) {
      f_name <- component
      model <- self$get_model("tte")
      if (!is.null(dataname)) {
        if (dataname == "LON") {
          model <- self$get_model("lon")
        }
        f_name <- paste0(f_name, "_", dataname)
      } else {
        dataname <- "TTE"
      }
      id_name <- model$id_var
      dat <- self$get_data(dataname)
      x <- dat[, c(id_name, "t"), drop = FALSE]
      f <- self$draws(name = f_name)
      FunctionDraws$new(x, f, f_name)
    },

    #' @description
    #' Extract baseline hazard as a \code{\link{FunctionDraws}} object.
    #' @param dataname Data where the function was evaluated.
    h0 = function(dataname = "GRD") {
      self$function_draws_tte("h0", dataname)
    },

    #' @description
    #' Extract instant hazard as a \code{\link{FunctionDraws}} object.
    #' @param dataname Data where the function was evaluated.
    instant_hazard = function(dataname = NULL) {
      self$function_draws_tte("inst_haz", dataname)
    },

    #' @description
    #' Extract cumulative hazard as a \code{\link{FunctionDraws}} object.
    #' @param dataname Data where the function was evaluated.
    cumulative_hazard = function(dataname = NULL) {
      self$function_draws_tte("cum_haz", dataname)
    },

    #' @description
    #' Extract survival function as a \code{\link{FunctionDraws}} object.
    #' @param dataname Data where the function was evaluated.
    survival_function = function(dataname = NULL) {
      s <- self$cumulative_hazard(dataname)$scale(a = -1)$exp()
      s$rename("survival_function")
      s
    },

    #' @description
    #' Visualize model fit (instant hazard).
    #'
    #' @param ... Arguments passed to the plot method of
    #' \code{\link{FunctionDraws}}.
    #' @param plot_events Should the event occurrences be plotted?
    plot_hazard = function(plot_events = TRUE, ...) {
      fd <- self$instant_hazard(dataname = "GRD")
      plt <- fd$plot(...) + ggtitle("Instant hazard")
      dat <- self$get_data("GRD")
      dat_tte <- self$get_data("TTE")
      id_var <- self$get_model()$lon$id_var
      if (plot_events) {
        map <- aes(x = t, y = 0, group = !!sym(id_var), pch = event)
        plt <- plt + geom_point(dat_tte, mapping = map, inherit.aes = F)
      }
      plt
    },

    #' Compute predictions at test points
    #'
    #' @param data_lon Data frame. If \code{NULL}, the data used in fitting
    #' is used.
    #' @param data_tte Data frame. If \code{NULL}, the data used in fitting
    #' is used.
    #' @param data_grid Data frame. If \code{NULL}, the data used in fitting
    #' is used.
    #' @param fitted_params Parameters or model fit.
    #' @param ... Other arguments to \code{generate_quantities()}.
    #' @return A new \code{\link{JointModelFit}} object.
    predict = function(data_lon = NULL, data_tte = NULL, data_grid = NULL,
                       fitted_params = NULL, ...) {
      if (is.null(data_lon)) {
        data_lon <- self$get_data("LON")
      }
      if (is.null(data_tte)) {
        data_tte <- self$get_data("TTE")
      }
      if (is.null(data_grid)) {
        data_grid <- self$get_data("GRD")
      }
      model <- self$get_model()

      # Prepare 'Stan' data
      dat <- model$create_standata(
        data_lon = data_lon,
        data_tte = data_tte,
        data_grid = data_grid,
        term_confs = self$term_confs,
        num_bf = NULL,
        scale_bf = NULL,
        skip_transform = NULL,
        prior_only = FALSE,
        set_transforms = FALSE
      )

      # Call 'Stan'
      gq <- self$gq(stan_data = dat$stan_data, fitted_params = fitted_params, ...)

      # Original data
      data_orig <- list(
        LON = dat$data_lon_orig,
        TTE = data_tte,
        GRD = data_grid
      )

      # Return
      JointModelFit$new(
        model, gq, data_orig, dat$stan_data, dat$full_term_confs
      )
    },

    #' Compute predictions at given time points
    #'
    #' @param t_test The test time points
    #' @param t_var Name of the time variable
    #' @param t_var_copy Names of other continuous covariates to which
    #' \code{t_test} should be copied.
    #' @param y_test Test y
    #' @return A new fit object
    predict_time = function(t_test, t_var = "t", y_test = NULL,
                            t_var_copy = NULL) {
      dat_lon <- self$get_data("LON")
      y_name <- self$get_model("lon")$y_var
      df_test <- create_df_predict_time(
        dat_lon, t_test, t_var, t_var_copy,
        y_test, y_name
      )
      self$predict(data_lon = df_test)
    }
  )
)
