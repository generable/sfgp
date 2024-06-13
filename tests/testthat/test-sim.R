test_that("data simulation works", {
  dat <- sfsim()
  expect_equal(ncol(dat), 9)
  plt <- ggplot(dat, aes(x = t, y = f, group = id, color = arm)) +
    geom_line() +
    geom_point(aes(x = t, y = y, group = id), inherit.aes = F) +
    facet_wrap(. ~ id)
  expect_s3_class(plt, "ggplot")
})

test_that("sf-only makes sense for sim data of different scales", {
  model <- TSModel$new(y ~ sf(t, id | arm))
  scales <- c(0.001, 0.1, 10, 10000)
  plots_prior <- list()
  plots_post <- list()
  j <- 0
  for (s in scales) {
    j <- j + 1
    set.seed(123)
    tt <- paste0("ts0_scale = ", s)
    dat <- sfsim(ts0_scale = s)
    dat <- add_sff_input(dat, model)
    fit_prior <- model$fit(dat, prior_only = TRUE, refresh = 0, chains = 1)
    p1 <- fit_prior$plot(color_by = "arm_kg") + scale_y_log10() + ggtitle(tt)
    fit_post <- model$fit(dat, refresh = 0, chains = 1)
    p2 <- fit_post$plot(color_by = "arm_kg") + ggtitle(tt)

    plots_prior[[j]] <- p1
    plots_post[[j]] <- p2
  }
  J <- length(scales)
  expect_equal(length(plots_post), J)
  expect_equal(length(plots_prior), J)
})
