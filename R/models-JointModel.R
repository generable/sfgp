#' The joint model class (R6 class)
#'
#' @export
#' @field tte The time-to-event model, has class \code{\link{TTEModel}}.
#' @field lon The longitudinal model, has class \code{\link{TSModel}}.
JointModel <- R6::R6Class("JointModel",
  inherit = StanModel,
  private = list(
    loglik_suffix = "joint",
    stanfiles_functions_impl = function() {
      c(
        self$tte$stanfiles_functions(),
        self$lon$stanfiles_functions()
      )
    },
    stancode_data_impl = function(datanames) {
      dn <- c("LON", "TTE", "GRD")
      c1 <- self$lon$stancode_data(dn)
      c2 <- self$tte$stancode_data(dn)
      paste(c1, c2, sep = "\n")
    },
    stancode_tdata_impl = function(datanames) {
      dn <- c("LON", "TTE", "GRD")
      c1 <- self$tte$stancode_tdata(dn)
      c2 <- self$lon$stancode_tdata(dn)
      paste(c1, c2, sep = "\n")
    },
    stancode_pars_impl = function() {
      c1 <- self$tte$stancode_pars()
      c2 <- self$lon$stancode_pars()
      paste(c1, c2, sep = "\n")
    },
    stancode_tpars_impl = function(datanames) {
      dn <- c("TTE", "GRD")
      c1 <- self$lon$stancode_tpars(c("LON", dn))
      c2 <- self$tte$stancode_tpars(dn)
      paste(c1, c2, sep = "\n")
    },
    stancode_model_impl = function() {
      c1 <- self$lon$stancode_model()
      c2 <- self$tte$stancode_model()
      paste(c1, c2, sep = "\n")
    },
    stancode_gq_impl = function() {
      c1 <- self$lon$stancode_gq()
      c2 <- self$tte$stancode_gq()
      paste(c1, c2, sep = "\n")
    }
  ),
  # PUBLIC
  public = list(
    tte = NULL,
    lon = NULL,

    #' @description
    #' Create model
    #'
    #' @param lon An object of class \code{\link{StanModel}}.
    #' @param h0 An object of class \code{\link{BaselineHazard}}.
    #' @param compile Should the 'Stan' model code be created and compiled.
    initialize = function(lon, h0 = NULL, compile = TRUE) {
      if (is.null(h0)) {
        h0 <- WeibullHazard$new()
      }
      checkmate::assert_class(h0, "BaselineHazard")
      checkmate::assert_class(lon, "StanModel")
      tte <- TTEModel$new(h0 = h0, link_name = "f_sum")
      tte$id_var <- lon$id_var
      self$tte <- tte
      self$lon <- lon
      super$initialize(compile)
    },

    #' @description
    #' The model description as a string
    string = function() {
      str <- paste0(
        class_name_hl(self), ":\n  ",
        "* Longitudinal submodel = ", self$lon$string(), "\n  ",
        "* Time-to-event submodel = ", self$tte$string(), "\n"
      )
      str
    },

    #' @description
    #' Create the 'Stan' data list from data frames.
    #'
    #' @param data_lon Longitudinal data, a data frame.
    #' @param data_tte The events data, a data frame.
    #' @param data_grid The integration grid data, a data frame.
    #' @param term_confs A list that specifies configuration of model terms.
    #' If name of any term is not found from the list, \code{$default_conf()}
    #' of that \code{FormulaTerm} is used.
    #' @param num_bf If not \code{NULL}, configurations of all
    #' \code{GPTerm}s are updated with this value.
    #' @param scale_bf If not \code{NULL}, configurations of all
    #' \code{GPTerm}s are updated with this value.
    #' @param skip_transform Term names whose input transform should be
    #' skipped.
    #' @param set_transforms If longitudinal data transforms should be set
    #' based on the given \code{data_lon}. This should be \code{TRUE} when
    #' fitting a model, and \code{FALSE} when computing predictions using GQ.
    #' @param prior_only Sample from prior only?
    #' @return A list.
    create_standata = function(data_lon,
                               data_tte,
                               data_grid = NULL,
                               term_confs = NULL,
                               num_bf = NULL,
                               scale_bf = NULL,
                               skip_transform = NULL,
                               prior_only = FALSE,
                               set_transforms = TRUE) {
      # Stan input of longitudinal model at longitudinal observations
      dl_lon <- self$lon$create_standata(
        data_lon, term_confs, num_bf, scale_bf,
        skip_transform, prior_only,
        set_transforms
      )

      # Stan input of longitudinal model at event observations
      full_term_confs <- dl_lon$full_term_confs
      dl_tte <- self$lon$term_list$create_standata(
        data_tte, "TTE", full_term_confs
      )

      # Stan input of longitudinal model at integration grid
      dl_grid <- self$lon$term_list$create_standata(
        data_grid, "GRD", full_term_confs
      )

      # Stan input of hazard function and TTE model
      dl_tte_haz <- list(
        t_hazard_TTE = get_x_from_data(data_tte, "t"),
        t_hazard_GRD = get_x_from_data(data_grid, "t")
      )
      dl_tte_events <- self$tte$create_standata(data_tte$event)


      # Combine and take unique fields
      sd <- c(
        dl_lon$stan_data, dl_tte,
        dl_grid, dl_tte_haz, dl_tte_events
      )
      sd <- sd[unique(names(sd))]

      # Return
      list(
        stan_data = sd,
        full_term_confs = full_term_confs,
        data_lon_orig = dl_lon$orig_data
      )
    },


    #' @description
    #' Fit the model.
    #'
    #' @param data_lon Longitudinal data, a data frame.
    #' @param data_tte The events data, a data frame.
    #' @param data_grid The integration grid data, a data frame. Can be
    #' created using \code{\link{create_jm_grid}}.
    #' @param term_confs A list that specifies configuration of model terms.
    #' If name of any term is not found from the list, \code{$default_conf()}
    #' of that \code{FormulaTerm} is used.
    #' @param num_bf If not \code{NULL}, configurations of all
    #' \code{GPTerm}s are updated with this value.
    #' @param scale_bf If not \code{NULL}, configurations of all
    #' \code{GPTerm}s are updated with this value.
    #' @param skip_transform Term names whose input transform should be
    #' skipped.
    #' @param prior_only Sample from prior only.
    #' @param ... Arguments passed to \code{sample} method of the
    #' 'CmdStanR' model.
    #' @return A \code{\link{JointModelFit}} object.
    fit = function(data_lon,
                   data_tte,
                   data_grid,
                   term_confs = NULL,
                   num_bf = NULL,
                   scale_bf = NULL,
                   skip_transform = NULL,
                   prior_only = FALSE,
                   ...) {
      # Get Stan model object
      stan_model <- self$get_stanmodel()
      if (is.null(stan_model)) {
        stop("Stan model does not exist, you need to call compile()!")
      }

      # Create Stan input list
      dat <- self$create_standata(
        data_lon, data_tte, data_grid, term_confs, num_bf, scale_bf,
        skip_transform, prior_only, TRUE
      )

      # Call 'Stan'
      stan_fit <- stan_model$sample(data = dat$stan_data, ...)

      # Create the fit object
      data_orig <- list(
        LON = dat$data_lon_orig,
        TTE = data_tte,
        GRD = data_grid
      )
      JointModelFit$new(
        self, stan_fit, data_orig, dat$stan_data, dat$full_term_confs
      )
    }
  )
)
