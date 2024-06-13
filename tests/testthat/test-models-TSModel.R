test_that("generating only the Stan code works", {
  code1 <- stancode_ts(y ~ gp(x), print = FALSE)
  expect_match(code1, "Term f_gp_x")
  expect_gt(nchar(code1), 20)
})


test_that("creating a TSModel works", {
  m <- TSModel$new(hello ~ gp(foo) + gp(foo, bar))
  expect_output(m$print(), "TSModel")
  expect_output(cat(m$term_list$latex()), "alpha")
  expect_output(m$term_list$terms[[1]]$print())
  expect_equal(m$y_var, "hello")
  expect_equal(m$term_list$length(), 3)
  tn <- m$term_list$terms[[1]]$name_long()
  expect_true(is.character(tn))
  expect_gt(nchar(tn), 3)
})

test_that("creating the Stan data works", {
  m <- TSModel$new(hello ~ gp(foo) + gp(foo, bar), compile = F, id_var = "bar")
  a <- data.frame(
    foo = c(-100, 2, 3, 4),
    hello = c(0, 3, 2, 1),
    bar = as.factor(c(1, 1, 2, 2))
  )
  expect_error(m$fit(), "you need to call compile")
  sbf <- 2.7
  confs <- m$term_list$fill_term_confs(num_bf = 2, scale_bf = sbf)
  m$term_list$set_transforms(a)
  sd <- m$create_standata(a, term_confs = confs)
  sd <- sd$stan_data
  expect_equal(sd$n_LON, 4)
  expect_equal(sd$B_foo, 2)
  expect_equal(max(sd$dat_foo_unit_LON), 1)
  expect_equal(min(sd$dat_foo_unit_LON), -1)
  expect_equal(sd$L_foo, sbf)
  expect_equal(m$id_var, "bar")
})


test_that("fitting a model and plotting function draws work", {
  m <- TSModel$new(hello ~ gp(foo) + gp(foo, bar), delta = 0.05)
  a <- data.frame(
    foo = c(1, 2, 3, 4),
    hello = c(0, 3, 2, 1),
    bar = as.factor(c(1, 1, 2, 2))
  )
  num_bf <- 8
  fit <- m$fit(data = a, refresh = 0, num_bf = num_bf)
  B_foo <- fit$term_confs[["f_gp_foo"]]$num_bf
  expect_equal(B_foo, num_bf)
  expect_equal(m$get_delta(), 0.05)
  fd <- fit$function_draws()
  f1 <- fit$function_draws("f_gp_foo")
  f2 <- fit$function_draws("f_gp_fooXbar")
  f3 <- f1 + f2
  f4 <- f3 - f1
  expect_output(fd$print(), "FunctionDraws")
  expect_output(f1$print(), "FunctionDraws")
  expect_output(f2$print(), "FunctionDraws")
  expect_output(f3$print(), "FunctionDraws")
  expect_output(f4$print(), "FunctionDraws")
  expect_s3_class(fd$plot(), "ggplot")
  expect_s3_class(f1$plot(), "ggplot")
  expect_s3_class(f2$plot(), "ggplot")
  expect_s3_class(f3$plot(), "ggplot")
  expect_s3_class(f4$plot(), "ggplot")
  expect_s3_class(fit$plot(), "ggplot")
  plt <- fit$plot(f_reference = rep(1, 4))
  ll <- fit$loglik()
  ep <- fit$measurement_error()
  expect_equal(length(ll), nrow(a))
  expect_equal(length(ep), nrow(a))
  expect_s3_class(plt, "ggplot")

  df <- fit$fit_quality_summary()
  expect_equal(nrow(df), 1)
  expect_equal(ncol(df), 2) # colnames c("id", "error_perc")
})

test_that("the simplest sf() example works", {
  # Simplest example
  fit <- example(formula = "y~sf(x)")
  expect_s3_class(fit$plot(), "ggplot")
})


test_that("a more complex example works", {
  r <- example(
    formula = "y ~ sf(x) + gp(x)",
    iter_warmup = 500, iter_sampling = 500, chains = 1
  )
  p1 <- r$plot()
  p2 <- (r$function_draws(data_scale = FALSE) - r$function_draws("f_gp_x"))$plot()
  p3 <- (r$function_draws("f_gp_x") + r$function_draws("f_baseline_id"))$plot()
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  p4 <- r$predict()$plot(predictive = FALSE, capped = FALSE)
  expect_s3_class(p4, "ggplot")
})

test_that("simplest model with empty formula (only grouped offset) works", {
  a <- TSModel$new(y ~ .)
  r <- a$fit(testdata)
  expect_equal(ncol(r$function_draws()$get_input()), 1)
})

test_that("example2() works with grouped sf and prior makes sense", {
  a <- example2(formula = "y ~ sf(time, id)", prior_only = TRUE, chains = 2)
  expect_s3_class(a$plot(), "ggplot")
  f <- a$function_draws(data_scale = F)
  fmax <- max(max(f$get_output()))
  expect_lt(fmax, expected = Inf) # can still be huge
  plt <- a$plot(filter_by = "arm", kept_vals = 1, group_by = "id", facet_by = "id")
  expect_s3_class(plt, "ggplot")
})

test_that("grouped sf models work and predicting works", {
  dat <- sfsim()
  m1 <- TSModel$new(y ~ sf(t, id) + gp(t))
  dat1 <- add_sff_input(dat, m1)

  expect_output(cat(m1$term_list$latex()), "alpha")
  f1_alt <- m1$fit(dat1, chains = 1, skip_transform = "f_log_sff_t")
  f1 <- m1$fit(dat1, chains = 1)

  m2 <- TSModel$new(y ~ sf(t, id | arm))
  dat2 <- add_sff_input(dat, m2)
  f2 <- m2$fit(dat2, chains = 1)

  expect_s3_class(f1$plot(), "ggplot")
  expect_s3_class(f1_alt$plot(), "ggplot")
  expect_s3_class(f2$plot(), "ggplot")
  expect_error(m1$term_list$get_term("foo"), "is not a Stan variable name")
  expect_error(m2$term_list$get_term("foo"), "is not a Stan variable name")

  pred1 <- f1$predict_time(seq(0, 58, by = 1))
  pred2 <- f1$predict_time(seq(50, 58, by = 1))
  p1 <- pred1$plot(plot_y = FALSE)
  p2 <- pred2$plot(plot_y = FALSE)

  p1_alt <- p1 + ggplot2::xlim(50, 58)
  expect_s3_class(p1_alt, "ggplot") # these two should
  expect_s3_class(p2, "ggplot") # look the same but with randomness

  # Check Stan data
  t_range <- range(dat$t)
  sd <- f1$get_stan_data()
  expect_equal(max(sd$dat_t_unit_LON), 1)
  expect_equal(min(sd$dat_t_unit_LON), -1)
  sd_pred1 <- pred1$get_stan_data()
  sd_pred2 <- pred2$get_stan_data()
  expect_gt(max(sd_pred1$dat_t_unit_LON), 1)
  expect_gt(max(sd_pred2$dat_t_unit_LON), 1)
})

test_that("treatment effect estimation works", {
  dat <- sfsim()
  model <- TSModel$new(y ~ sf(t, id | arm) + gp(t) + gp(t, arm))
  dat <- add_sff_input(dat, model)

  fit_post <- model$fit(dat, chains = 1)
  fit_prior <- model$fit(dat, chains = 1, prior_only = TRUE)

  # Treatment effect methods
  te <- treatment_effect(fit_post, fit_prior, time_var = "t", group_var = "arm")
  te2 <- treatment_effect(fit_post, fit_prior,
    time_var = "t", group_var = "arm", method = "group_est"
  )
  te3 <- treatment_effect(fit_post, fit_prior,
    time_var = "t", group_var = "arm", br = TRUE
  )
  expect_s3_class(te$traj$plot(), "ggplot")
  expect_s3_class(te2$traj$plot(), "ggplot")
  expect_s3_class(te3$traj$plot(), "ggplot")
  sumr1 <- summarize_dur(te)
  sumr2 <- summarize_dur(te, only_responding = TRUE)
  expect_equal(ncol(sumr1), 3)
  expect_equal(ncol(sumr2), 3)

  # SF Trajectories
  sf1 <- te$p_new$function_draws("f_log_sff_t")
  sf2 <- te2$p_new$function_draws("f_log_sff_t")
  expect_s3_class(sf1$plot(), "ggplot") # these should be similar but
  expect_s3_class(sf2$plot(), "ggplot") # this should have less variation


  # Depth of response
  plt_dor <- plot_dor(te, halfeye = T, .width = c(0.99, 0.95))
  expect_s3_class(plt_dor, "ggplot")

  # Duration of response
  plt_dur <- plot_dur(te)
  expect_s3_class(plt_dur, "ggplot")

  # Probability of difference
  d1 <- dor_diff(te$depth, gt = T)
  d2 <- dur_diff(te$duration, gt = T)
  expect_equal(length(d1), 3)
  expect_equal(nrow(d2), 3)
})


test_that("prior_sigma_informed() works", {
  pr <- prior_sigma_informed()
  expect_equal(pr, "normal(0.06509, 0.01506)")
})
