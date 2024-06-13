test_that("creating and fitting a JointModel works", {
  # Create model
  tsm <- TSModel$new(y ~ gp(t), compile = F, id_var = "subject")
  m <- JointModel$new(lon = tsm)
  expect_output(print(m), "Longitudinal submodel")

  # Create data
  t <- c(55, 55, 60, 50, 60, 60)
  event <- as.logical(c(1, 1, 0, 1, 0, 0))
  subject <- as.factor(c(1, 2, 3, 4, 5, 6))
  dat_tte <- data.frame(subject, t, event)
  dat_ts <- sfsim(n1 = 2, n2 = 2, n3 = 2)
  colnames(dat_ts)[1] <- "subject"
  dat_grid <- create_jm_grid(dat_tte, dat_ts,
    num_points = 20,
    id_var = "subject"
  )

  # Fit
  fit <- m$fit(
    data_lon = dat_ts, data_tte = dat_tte, data_grid = dat_grid, chains = 1
  )
  expect_output(print(fit), "JointModelFit")

  # Tests
  plt1 <- fit$plot()
  expect_s3_class(plt1, "ggplot")
  expect_s3_class(fit$h0()$plot(), "ggplot")

  cum_haz <- fit$cumulative_hazard()
  expect_equal(cum_haz$get_name(), "cum_haz")

  expect_s3_class(fit$plot_hazard(), "ggplot")

  # Predict
  p <- fit$predict()
  expect_output(print(p), "JointModelFit")
  plt2 <- p$plot()
  expect_s3_class(plt2, "ggplot")

  # Predict dense in time
  t_test <- seq(0, 50, by = 2)
  pp <- fit$predict_time(t_test)

  # Smooth plot
  plt3 <- pp$function_draws()$plot()
  expect_s3_class(plt3, "ggplot")
})
