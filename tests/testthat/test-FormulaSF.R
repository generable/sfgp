test_that("creating an SFF term with formula works", {
  dat <- sfsim()
  form <- y ~ sff(t | ks ~ offset(id_ks) + gp(t_ks), kg ~ offset(id_kg) + gp(t_kg))
  m <- TSModel$new(form)
  dat <- add_sff_input(dat, m)
  term <- m$term_list$terms$f_log_sff_t
  num_bf <- sample(10:20, 1)
  fit <- m$fit(dat, chains = 1, num_bf = num_bf)
  expect_true(inherits(term, "FormulaSFTerm"))

  # Check that term configuration works for terms inside FormulaSFTerm
  expect_equal(fit$get_stan_data()$B_t_ks, num_bf)

  # TODO: Not working!
  # conf_sff <- list(kg = list(f_gp_t_kg = list(scale_bf = 2.3)))
  # conf2 <- list(f_log_sff_t = conf_sff)
  # fit2 <- m$fit(dat, chains = 1, term_confs = conf2, scale_bf = NULL)

  plt <- fit$plot()
  expect_s3_class(plt, "ggplot")

  # Plot inferred ks, kg as a function of time
  f_t_ks <- fit$create_functiondraws(c("id_ks", "t_ks"), "ks")
  f_t_kg <- fit$create_functiondraws(c("id_kg", "t_kg"), "kg")
  plt_ks <- f_t_ks$log()$plot() # should depend on time
  plt_kg <- f_t_kg$log()$plot() # should depend on time
  expect_s3_class(plt_ks, "ggplot")
  expect_s3_class(plt_kg, "ggplot")
})
