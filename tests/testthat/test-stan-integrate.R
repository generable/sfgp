# test_that("expose_stan_functions() works", {
#   # Expose Stan function to R
#   b <- expose_stan_functions("hazard/integrate")
#
#   # Test
#   h <- 0.01
#   xx <- seq(0, 10, by = h)
#   yy <- sin(xx) + 0.4 * sin(2.5 * xx)
#   x <- seq(0, 10, by = 0.2)
#   I_analytic <- -cos(x) - 0.16 * cos(2.5 * x) + 1.16
#   I_numeric <- b$functions$integrate_1d(x, xx, yy)
#   mae <- max(abs(I_analytic - I_numeric))
#   expect_lt(mae, 0.001)
#
#   # Integration for 3 subjects
#   t_grid <- c(xx, xx, xx)
#   t_out <- c(x, x, x)
#   y_grid <- c(yy, yy + 2, yy + 1)
#   G <- length(xx)
#   N <- length(x)
#   io <- list(c(1, N), c(N + 1, 2 * N), c(2 * N + 1, 3 * N))
#   ig <- list(c(1, G), c(G + 1, 2 * G), c(2 * G + 1, 3 * G))
#   I_numeric <- b$functions$integrate_1d_many(t_out, t_grid, y_grid, io, ig)
#   i2 <- io[[2]]
#   I2 <- I_numeric[i2[1]:i2[2]]
#
#   # Test that integral is correct for subject 2
#   I2_analytic <- 2 * x - cos(x) - 0.16 * cos(2.5 * x) + 1.16
#   mae2 <- max(abs(I2_analytic - I2))
#   expect_lt(mae, 0.001)
# })
