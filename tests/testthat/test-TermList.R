test_that("a TermList can be created", {
  a <- create_termlist(y ~ gp(x) + sf(x), NULL)
  expect_equal(a$length(), 2)
  expect_output(a$print())
  expect_error(a$get_term("foo"), "is not a Stan variable name")
  sn <- a$stan_names()
  t <- a$get_term(sn[2])
  expect_match(t$stanname_base(), "f_log_sf_x")
  expect_equal(length(a$stanfiles_functions()), 5)
})

test_that("creating Stan code from TermList works", {
  a <- create_termlist(y ~ gp(x) + sf(x), NULL)
  expect_match(a$stancode_data("OBS"), "f_gp_x")
  code_tp <- a$stancode_tpars(c("OBS", "PRED"))
  expect_output(cat(code_tp), "f_log_sf_x_PRED")
  expect_output(cat(a$stancode_gq()), "f_sum")
})

test_that("creating Stan data for TermList works", {
  a <- create_termlist(y ~ gp(x) + sf(x), NULL)
  d <- data.frame(x = c(1, 2, 3, 4))
  expect_message(a$check_transforms(d), "[1, 4]")
  a$set_transforms(d)
  expect_message(a$check_transforms(d), "[-1, 1]")
  expect_message(a$check_transforms(d), "[0.25, 1]")

  # Should work
  r <- a$create_standata(data = d, dataname = "ASD")
  expect_true("L_x" %in% names(r))

  # cannot predict outside [-L, L]
  d2 <- data.frame(x = c(1, 2, 3, 4, 10))
  expect_error(
    a$create_standata(data = d2, dataname = "TEST"),
    "GP approximation not valid"
  )
})

test_that("A custom prior can be given", {
  prior_sf <- list(prior_ks = "normal(3, 4)")
  prior_invalid <- list(prior_asd = "normal(3, 4)")
  a <- create_termlist(y ~ sf(x), list(f_log_sf_x = prior_sf))
  expect_output(cat(a$stancode_model()), "3, 4")
  expect_error(
    create_termlist(y ~ sf(x), list(f_log_sf_x = prior_invalid)),
    "cannot add bindings to a locked environment"
  )
  expect_error(
    create_termlist(y ~ sf(x), list(FOO = prior_sf)),
    "is not a Stan variable name"
  )
})
