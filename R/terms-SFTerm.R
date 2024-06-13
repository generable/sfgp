# Used by FormulaParser
create_sfterm <- function(covariates, hierarchy) {
  x_name <- covariates[1]
  z_name <- covariates[2]
  h_name <- hierarchy
  if (is.na(z_name)) {
    out <- SFTerm$new(x_name = x_name)
  } else {
    formula_def <- c(
      paste0("kg ~ offset(", z_name, "_kg)"),
      paste0("ks ~ offset(", z_name, "_ks)")
    )
    if (!is.null(h_name)) {
      formula_def <- c(
        paste0("kg ~ offset(", z_name, "_kg | ", h_name, "_kg)"),
        paste0("ks ~ offset(", z_name, "_ks | ", h_name, "_ks)")
      )
    }
    out <- create_sffterm(covariates, formula_def)
  }
  out
}

# Used by FormulaParser
create_sffterm <- function(covariates, formula_def) {
  if (length(formula_def) != 2) {
    stop("error parsing formula definition for sf parameters")
  }
  x_name <- covariates[1]
  for (f in formula_def) {
    ff <- as.formula(f)
    y_var <- FormulaParser$new(ff)$get_y_name()
    if (y_var == "ks") {
      form_ks <- ff
    } else if (y_var == "kg") {
      form_kg <- ff
    } else {
      stop("sff formula must have either kg or ks on the RHS")
    }
  }
  FormulaSFTerm$new(
    x_name = x_name, formula_kg = form_kg, formula_ks = form_ks
  )
}

# Define the SFTerm class
SFTerm <- R6::R6Class("SFTerm",
  inherit = FormulaTerm,
  lock_class = TRUE,
  private = list(
    latex_type = "log-SF",
    latex_param_names = function() {
      c("k_{g}", "k_{s}")
    },
    stanfiles_functions_impl = function() {
      c("sf")
    }
  ),
  public = list(
    prior_kg = "student_t(4, 0, 1)",
    prior_ks = "student_t(4, 0, 1)",
    x_transform = MaxScaleTransform$new(),
    initialize = function(x_name) {
      checkmate::assert_character(x_name, any.missing = FALSE)
      if (inherits(self, "FormulaSFTerm")) {
        private$typename <- "log_sff"
      } else {
        private$typename <- "log_sf"
      }
      private$suffix <- x_name
      self$x_name <- x_name
    },
    stancode_data = function(used_names, datanames) {
      x <- self$stanlines_data_x(datanames)
      build_stancode_data(x$lines, x$stannames, used_names)
    },
    stancode_pars = function() {
      scb <- stancode_bounds(lower = 0, upper = NULL)
      line1 <- paste0("  real", scb, " kg", "; // growth rate")
      line2 <- paste0("  real", scb, " ks", "; // shrinkage rate")
      paste0(line1, "\n", line2, "\n")
    },
    stancode_tpars = function(datanames) {
      sn <- self$stanname(datanames)
      sx <- self$stanname_x(datanames)
      nn <- paste0("n_", datanames)
      fc <- paste0("log_sf_norm(", sx, ", kg, ks)")
      c1 <- paste0("  vector[", nn, "] ", sn, " = ", fc, ";")
      paste(c1, collapse = "\n")
    },
    stancode_model = function() {
      line1 <- paste0("  kg ~ ", self$prior_kg, ";")
      line2 <- paste0("  ks ~ ", self$prior_ks, ";")
      paste(line1, line2, sep = "\n")
    },
    standata = function(datasets, conf) {
      conf <- self$ensure_conf(conf)
      self$create_standata_x(datasets)
    }
  )
)

# Define the FormulaSFTerm class
FormulaSFTerm <- R6::R6Class("FormulaSFTerm",
  inherit = SFTerm,
  lock_class = TRUE,
  private = list(
    latex_param_names = function() {
      c("\\mathbf{k}_{g}", "\\mathbf{k}_{s}")
    },
    stanfiles_functions_impl = function() {
      f1 <- self$term_list_kg$stanfiles_functions()
      f2 <- self$term_list_ks$stanfiles_functions()
      c("sf", f1, f2)
    }
  ),
  public = list(
    term_list_kg = NULL,
    term_list_ks = NULL,
    prior_log_C_kg = "normal(-4, 4)",
    prior_log_C_ks = "normal(-3, 4)",
    prior_kg = NULL,
    prior_ks = NULL,
    initialize = function(x_name, formula_kg, formula_ks) {
      self$term_list_kg <- create_termlist(formula_kg, self$prior_kg)
      self$term_list_ks <- create_termlist(formula_ks, self$prior_ks)
      self$term_list_kg$fsum_name <- "log_kg"
      self$term_list_ks$fsum_name <- "log_ks"
      super$initialize(x_name)
    },
    input_vars = function() {
      v1 <- super$input_vars()
      vars_ks <- as.vector(self$term_list_ks$input_vars())
      vars_kg <- as.vector(self$term_list_kg$input_vars())
      vars_no_ks <- replicate_without_suffix(vars_ks, "ks")
      vars_no_kg <- replicate_without_suffix(vars_kg, "kg")
      unique(c(v1, vars_ks, vars_kg, vars_no_kg, vars_no_ks))
    },
    stancode_data = function(used_names, datanames) {
      x <- self$stanlines_data_x(datanames)
      lines <- c(x$lines)
      dn <- c(x$stannames)
      out <- build_stancode_data(lines, dn, used_names)
      sc1 <- self$term_list_kg$stancode_data(datanames)
      sc2 <- self$term_list_ks$stancode_data(datanames)
      list(
        code = paste(out$code, sc1, sc2, sep = "\n"),
        stannames = out$stannames
      )
    },
    stancode_tdata = function(datanames) {
      paste(
        self$term_list_kg$stancode_tdata(datanames),
        self$term_list_ks$stancode_tdata(datanames),
        sep = "\n"
      )
    },
    stancode_pars = function() {
      code_C <- paste0("  real log_C_kg;\n  real log_C_ks;")
      paste(
        self$term_list_kg$stancode_pars(),
        self$term_list_ks$stancode_pars(),
        code_C,
        sep = "\n"
      )
    },
    stancode_tpars = function(datanames) {
      code <- paste(
        self$term_list_kg$stancode_tpars(datanames),
        self$term_list_ks$stancode_tpars(datanames),
        sep = "\n"
      )
      head <- "  // FormulaSFTerm"
      l1 <- par_to_nat_scale(self$term_list_kg$fsum_name, "kg", datanames)
      l2 <- par_to_nat_scale(self$term_list_ks$fsum_name, "ks", datanames)
      sn <- self$stanname(datanames)
      xn <- self$stanname_x(datanames)
      l3 <- paste0("  vector[n_", datanames, "] ", sn, " = ",
        "log_sf_norm_vec(", xn, ", kg_", datanames, ", ks_",
        datanames, ");",
        collapse = "\n"
      )
      paste(code, head, l1, l2, l3, sep = "\n")
    },
    stancode_model = function() {
      code1 <- self$term_list_kg$stancode_model()
      code2 <- self$term_list_ks$stancode_model()
      code3 <- paste0("  log_C_kg ~ ", self$prior_log_C_kg, ";")
      code4 <- paste0("  log_C_ks ~ ", self$prior_log_C_ks, ";")
      paste0(code1, code2, code3, code4)
    },
    stancode_gq = function() {
      code1 <- self$term_list_kg$stancode_gq()
      code2 <- self$term_list_ks$stancode_gq()
      code3 <- super$stancode_gq()
      code4 <- paste0("\n  // SF parameters")
      code5 <- paste0("  vector[n_LON] kg = kg_LON;")
      code6 <- paste0("  vector[n_LON] ks = ks_LON;")
      paste(code1, code2, code3, code4, code5, code6, sep = "\n")
    },
    standata = function(datasets, conf) {
      if (length(datasets) != 1) {
        stop("standata() for multiple datasets not implemntd for FormulaSFTerm")
      }
      data <- datasets[[1]]
      dataname <- names(datasets)[1]
      conf <- self$ensure_conf(conf)
      sd_x <- self$create_standata_x(datasets)
      sd_1 <- self$term_list_kg$create_standata(data, dataname, conf$kg)
      sd_2 <- self$term_list_ks$create_standata(data, dataname, conf$ks)
      c(sd_x, sd_1, sd_2)
    }
  )
)

# Stan code for raw parameters
raw_params_stancode <- function(param_name, G_hier, G_id) {
  line1 <- paste0("  vector[", G_hier, "] log_mu_", param_name, ";")
  line2 <- paste0("  vector[", G_hier, "] log_sigma_", param_name, ";")
  line3 <- paste0("  vector[", G_id, "] log_z_", param_name, ";")
  paste(line1, line2, line3, "\n", sep = "\n")
}

# Hierarchical parameters to vectors with length equal to data
param_vecs_log_hier <- function(param_name, h_name, z_name) {
  pn_z <- param_name
  list(
    z = paste0("log_z_", pn_z, "[", z_name, "]"),
    mu = paste0("log_mu_", param_name, "[", h_name, "]"),
    sigma = paste0("log_sigma_", param_name, "[", h_name, "]")
  )
}

# Hierarchical prior on log scale
hier_prior_stancode <- function(param_name, prior) {
  z <- paste0("log_z_", param_name)
  m <- paste0("log_mu_", param_name)
  s <- paste0("log_sigma_", param_name)
  line1 <- paste0("  ", z, " ~ ", prior$log_z, ";")
  line2 <- paste0("  ", m, " ~ ", prior$log_mu, ";")
  line3 <- paste0("  ", s, " ~ ", prior$log_sigma, ";")
  paste(line1, line2, line3, "\n", sep = "\n")
}

# kg/ks parameter on natural scale
par_to_nat_scale <- function(f_sum_name, pn, datanames) {
  paste0("  vector[n_", datanames, "] ", pn, "_",
    datanames, " = exp(log_C_", pn, " + ",
    f_sum_name, "_", datanames, ");",
    collapse = "\n"
  )
}
