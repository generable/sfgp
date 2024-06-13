# Define the EmaxTerm class
EmaxTerm <- R6::R6Class("EmaxTerm",
  inherit = FormulaTerm,
  lock_class = TRUE,
  private = list(
    latex_type = "EM",
    latex_param_names = function() {
      c("\text{ED}_{50}", "E_{\text{MAX}}", "gamma")
    },
    stanfiles_functions_impl = function() {
      c("emax")
    }
  ),
  public = list(
    prior_log_ed50 = "student_t(4, 5, 2)",
    prior_gamma = "inv_gamma(5, 2)",
    prior_emax = "student_t(4, 0, 1)",
    initialize = function(x_name, z_name) {
      checkmate::assert_character(x_name, any.missing = FALSE)
      checkmate::assert_character(z_name, any.missing = FALSE)
      private$typename <- "emax"
      private$suffix <- paste0(x_name, "X", z_name)
      self$x_name <- x_name
      self$z_name <- z_name
    },
    stancode_data = function(used_names, datanames) {
      x <- self$stanlines_data_x(datanames)
      z <- self$stanlines_data_z(datanames)
      lines <- c(x$lines, z$lines)
      names <- c(x$stannames, z$stannames)
      build_stancode_data(lines, names, used_names)
    },
    stancode_pars = function() {
      gn <- paste0("G_", self$z_name)
      line1 <- paste0("  vector[", gn, "] log_ED50;\n")
      line2 <- paste0("  real<lower=0> Emax", "; // max effect")
      line3 <- paste0("  real<lower=1> gamma", "; // shape")
      paste0(line1, "\n", line2, "\n", line3, "\n")
    },
    stancode_tpars = function(datanames) {
      sn <- self$stanname(datanames)
      sx <- self$stanname_x(datanames)
      nn <- paste0("n_", datanames)
      dnz <- self$stanname_z(datanames)
      fc <- paste0("emax(", sx, ", log_ED50[", dnz, "], Emax, gamma)")
      c1 <- paste0("  vector[", nn, "] ", sn, " = ", fc, ";")
      paste(c1, collapse = "\n")
    },
    stancode_model = function() {
      line1 <- paste0("  log_ED50 ~ ", self$prior_log_ed50, ";")
      line2 <- paste0("  Emax ~ ", self$prior_emax, ";")
      line3 <- paste0("  gamma ~ ", self$prior_gamma, ";")
      paste(line1, line2, line3, sep = "\n")
    },
    standata = function(datasets, conf) {
      conf <- self$ensure_conf(conf)
      c(
        self$create_standata_x(datasets),
        self$create_standata_z(datasets)
      )
    }
  )
)
