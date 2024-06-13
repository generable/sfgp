#' The time-to-event model class (R6 class)
#'
#' @export
#' @field h0 An object of class \code{\link{BaselineHazard}}.
#' @field id_var name of id variable
TTEModel <- R6::R6Class("TTEModel",
  inherit = StanModel,
  private = list(
    link_name = NULL,
    loglik_suffix = "tte",
    stanfiles_functions_impl = function() {
      c(
        self$h0$stanfiles_functions(),
        "hazard/integrate",
        "hazard/inst_hazard"
      )
    },
    stancode_data_impl = function(datanames) {
      datanames <- c("TTE", "GRD")
      c1 <- self$h0$stancode_data(datanames)
      c2a <- "  // Events"
      c2 <- paste0("  int<lower=0> num_events;")
      c3 <- paste0("  array[num_events] int<lower=1> event_ids;")
      paste(c1, c2a, c2, c3, sep = "\n")
    },
    stancode_tdata_impl = function(datanames) {
      dn <- c("GRD", "TTE")
      c1 <- self$h0$stancode_tdata(dn)
      c2 <- "  // Indices (TTEModel)"
      id_var <- self$id_var
      c3 <- paste0(
        "  array[G_", id_var, ", 2] int inds_", dn, " = inds_array(",
        self$stanname_id(dn), ", G_", id_var, "); ",
        collapse = "\n"
      )
      paste(c1, c2, c3, sep = "\n")
    },
    stancode_pars_impl = function() {
      paste(
        self$h0$stancode_pars(),
        "  real assoc; // association parameter",
        sep = "\n"
      )
    },
    stancode_tpars_impl = function(datanames) {
      dn <- c("TTE", "GRD")
      c0 <- self$h0$stancode_tpars(dn)

      c1a <- "  // Instant and cumulative hazard"
      t_out <- self$h0$stanname_t(dn[1])
      t_grid <- self$h0$stanname_t(dn[2])
      y_grid <- paste0("inst_haz_", dn[2])
      inds_out <- paste0("inds_", dn[1])
      inds_grid <- paste0("inds_", dn[2])
      c1 <- stancall_inst_haz(dn, private$link_name)
      c2 <- paste0(
        "  vector[n_", dn[1], "] cum_haz_", dn[1],
        " = integrate_1d_many(", t_out, ",", t_grid, ", ", y_grid, ", ",
        inds_out, ", ", inds_grid, ");"
      )
      c3 <- "\n  // Log likelihood"
      c4 <- paste0(
        "  real log_lik_", private$loglik_suffix,
        " = sum(-cum_haz_TTE) + ",
        "sum(log(inst_haz_TTE[event_ids]));"
      )
      paste(c0, c1a, c1, c2, c3, c4, sep = "\n")
    },
    stancode_model_impl = function() {
      c1 <- self$h0$stancode_model()
      c2 <- "  assoc ~ normal(0, 1);"
      c3 <- stancode_loglik(private$loglik_suffix)
      paste(c1, c2, c3, sep = "\n")
    },
    stancode_gq_impl = function() {
      dn <- "TTE"
      c1 <- self$h0$stancode_gq()
      c2 <- paste0("  vector[n_", dn, "] inst_haz = inst_haz_", dn, ";")
      c3 <- paste0("  vector[n_", dn, "] cum_haz = cum_haz_", dn, ";")
      paste(c1, c2, c3, sep = "\n")
    }
  ),
  # PUBLIC
  public = list(
    h0 = NULL,
    id_var = "id",

    #' @description Get name of id variable in 'Stan' code.
    #' @param datanames names of data sets
    stanname_id = function(datanames) {
      paste0("dat_", self$id_var, "_", datanames)
    },

    #' @description
    #' Create model
    #'
    #' @param h0 An object of class \code{\link{BaselineHazard}}.
    #' @param link_name Base name of the link variable in 'Stan' code.
    initialize = function(h0 = NULL, link_name = "f_sum") {
      private$link_name <- link_name
      if (is.null(h0)) {
        h0 <- WeibullHazard$new()
      }
      checkmate::assert_class(h0, "BaselineHazard")
      self$h0 <- h0
      super$initialize(compile = FALSE)
    },

    #' @description
    #' Create the 'Stan' data list from a data vector.
    #'
    #' @param event a logical vector
    #' @return A list.
    create_standata = function(event) {
      checkmate::assert_logical(event)

      # Return
      list(
        event_ids = which(event),
        num_events = sum(event)
      )
    }
  )
)


# Stan code line that computes instant hazard
stancall_inst_haz <- function(dn, link_name) {
  paste0("  vector[n_", dn, "] inst_haz_", dn,
    " = compute_inst_hazard(h0_", dn, ", assoc, ",
    link_name, "_", dn, ", y_loc, y_scale, delta);",
    collapse = "\n"
  )
}
