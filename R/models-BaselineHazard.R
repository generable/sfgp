#' Hazard function (R6 class)
BaselineHazard <- R6::R6Class("BaselineHazard",
  inherit = StanCodeCreator,
  private = list(
    t_name = "t_hazard",
    stancode_data_impl = function(datanames) {
      paste0(" vector[n_", datanames, "]", self$stanname_t(datanames), ";",
        collapse = "\n"
      )
    }
  ),
  public = list(


    #' @description
    #' Get name of the hazard time variable in 'Stan' code
    #'
    #' @param datanames names of data sets
    stanname_t = function(datanames) {
      paste0(private$t_name, "_", datanames)
    }
  )
)


#' Weibull hazard function
#'
#' @export
WeibullHazard <- R6::R6Class("WeibullHazard",
  inherit = BaselineHazard,
  private = list(
    prior_lambda = NULL,
    prior_gamma = NULL,
    stanfiles_functions_impl = function() {
      c("hazard/weibull")
    },
    stancode_pars_impl = function() {
      paste(
        "  real<lower=0> h0_lambda; // hazard fun param (scale)",
        "  real<lower=0> h0_gamma; // hazard fun param (shape)",
        sep = "\n"
      )
    },
    stancode_tpars_impl = function(datanames) {
      c1 <- paste0(
        "  vector[n_", datanames, "] h0_", datanames, " = weibull_haz(",
        self$stanname_t(datanames), ", h0_lambda, h0_gamma);"
      )
      paste(c1, collapse = "\n")
    },
    stancode_model_impl = function() {
      l1 <- paste0("  h0_lambda ~ ", private$prior_lambda, "; // weibull param")
      l2 <- paste0("  h0_gamma ~ ", private$prior_gamma, "; // weibull param")
      paste(l1, l2, sep = "\n")
    }
  ),

  # PUBLIC
  public = list(


    #' @description
    #' Create hazard
    #'
    #' @param prior_lambda prior of Weibull scale param
    #' @param prior_gamma prior of Weibull shape param
    initialize = function(prior_lambda = "gamma(5, 5)",
                          prior_gamma = "inv_gamma(3, 6)") {
      private$prior_lambda <- prior_lambda
      private$prior_gamma <- prior_gamma
    }
  )
)
