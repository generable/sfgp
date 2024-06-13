test_that("Emax term works", {
  sim_emax <- function(x, e0, emax, ed50, gamma, sigma) {
    xg <- x^gamma
    f <- e0 + (emax * xg) / (ed50^gamma + xg)
    y <- exp(f + sigma * rnorm(length(f)))
    data.frame(x = x, f = f, y = y)
  }

  gamma <- 10
  e0 <- 3
  x <- seq(0, 400, by = 4)
  df1 <- sim_emax(x, e0, 3, 160, gamma, 0.25)
  df2 <- sim_emax(x, e0, 3, 260, gamma, 0.25)
  df3 <- sim_emax(x, e0, 3, 160, gamma, 0.25)
  df4 <- sim_emax(x, e0, 3, 260, gamma, 0.25)

  df1$dose <- df1$x
  df2$dose <- df2$x
  df3$dose <- df3$x
  df4$dose <- df4$x

  df <- rbind(df1, df2, df3, df4)
  df$id <- as.factor(rep(1:4, each = nrow(df1)))
  df$mutation <- as.factor(rep(rep(0:1, each = nrow(df1)), 2))

  # Fitting a model
  m <- TSModel$new(y ~ emax(x, mutation))
  fit <- m$fit(data = df, chains = 1)

  # Studying the fit
  plt1 <- fit$function_draws()$log()$plot(facet_by = "id")
  plt2 <- fit$plot(capped = FALSE, facet_by = "id")
  expect_s3_class(plt1, "ggplot")
  expect_s3_class(plt2, "ggplot")
})
