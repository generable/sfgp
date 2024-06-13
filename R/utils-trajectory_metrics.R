# Best response after baseline
best_response <- function(y) {
  checkmate::assert_numeric(y, lower = 0)
  checkmate::assert_vector(y, min.len = 2)
  baseline <- y[1]
  y_after_baseline <- y[2:length(y)]
  A <- baseline
  C <- baseline - min(y_after_baseline)
  C / A
}

# Time to 120% of baseline
ttb120 <- function(t, y) {
  checkmate::assert_numeric(y, lower = 0)
  inds <- which(y >= 1.2 * y[1])
  if (length(inds) == 0) {
    return(max(t))
  }
  t_idx <- min(inds)
  t[t_idx]
}

# Depth of response from a single trajectory
depth_of_response <- function(y) {
  checkmate::assert_numeric(y, lower = 0)
  baseline <- y[1]
  A <- baseline
  C <- baseline - min(y)
  C / A
}

# Depth of response from an rvars trajectory
depth_of_response_rvar <- function(y, br = FALSE) {
  ymat <- posterior::as_draws_matrix(posterior::merge_chains(y))
  D <- posterior::ndraws(ymat)
  dor <- rep(0, D)
  for (d in seq_len(D)) {
    x <- as.vector(ymat[d, ])
    if (br) {
      dor[d] <- best_response(x) # can be negative, max 1
    } else {
      dor[d] <- depth_of_response(x) # between [0,1]
    }
  }
  posterior::rvar(dor)
}

# Depth of response from a FunctionDraws
depth_of_response_functiondraws <- function(fd, id_var, br = FALSE) {
  a <- split_to_subjects(fd, id_var)
  getr <- function(x) {
    depth_of_response_rvar(x$get_output(), br = br)
  }
  rvar_list <- lapply(a, getr)
  unlist_rvar(rvar_list)
}

# Duration of response from a single trajectory
duration_of_response <- function(t, y) {
  L <- length(t)
  checkmate::assert_numeric(t, len = L)
  checkmate::assert_numeric(y, len = L)
  baseline <- y[1]
  during_response <- 0
  response_start <- 0
  response_end <- 0
  response_exists <- FALSE
  min_to_date <- baseline
  for (it in seq_len(L)) {
    # update min to date with each step
    if (y[it] < min_to_date) {
      min_to_date <- y[it]
    }
    # response starts at SLD <= 30% below baseline value
    if (y[it] <= 0.7 * baseline) {
      if (during_response == 0 && response_start == 0) {
        response_start <- t[it]
        response_exists <- TRUE
        during_response <- 1
      }
    }
    if (during_response == 1) {
      # response ends at SLD >= 20% above baseline value
      if (y[it] >= 1.2 * min_to_date) {
        response_end <- t[it]
        during_response <- 0
      }
    }
    if (it == L) {
      # trajectory ends with subject still in response
      if (during_response == 1) {
        response_end <- t[it]
      }
    }
  }

  # return
  list(
    start = response_start,
    end = response_end,
    response_endures = as.logical(during_response),
    response_exists = response_exists
  )
}

# Duration of response from an rvars trajectory
duration_of_response_rvars <- function(t, y, only_responding = FALSE) {
  ymat <- posterior::as_draws_matrix(posterior::merge_chains(y))
  D <- posterior::ndraws(ymat)
  exists <- rep(0, D)
  endures <- rep(0, D)
  start <- rep(0, D)
  end <- rep(0, D)
  TTB120 <- rep(0, D)
  for (d in seq_len(D)) {
    x <- as.vector(ymat[d, ])
    ret <- duration_of_response(t, x)
    exists[d] <- as.numeric(ret$response_exists)
    endures[d] <- as.numeric(ret$response_endures)
    start[d] <- ret$start
    end[d] <- ret$end
    TTB120[d] <- ttb120(t, x)
  }
  duration <- end - start

  # filter and return
  mat <- cbind(exists, endures, duration, start, end, TTB120)
  if (only_responding) {
    idx_rows <- which(mat[, 1] == 1)
    mat <- mat[idx_rows, , drop = F]
  }
  posterior::rvar(mat)
}

# Duration of response from a functiondraws
duration_of_response_functiondraws <- function(fd, id_var, time_var,
                                               only_responding) {
  a <- split_to_subjects(fd, id_var)
  getr <- function(x) {
    t <- x$get_input()[[time_var]]
    y <- x$get_output()
    duration_of_response_rvars(t, y, only_responding)
  }
  rvar_list <- lapply(a, getr)
  rvar_list
}


# Plot trajectory and indicate response start and end
trajectory_plot <- function(t, y) {
  dur <- duration_of_response(t, y)
  TTB <- ttb120(t, y)
  dor <- depth_of_response(y)
  duration <- dur$end - dur$start
  plt <- ggplot(data.frame(t = t, y = y), aes(t, y)) +
    geom_line() +
    geom_point() +
    ggtitle(paste0(
      "Depth = ", round(dor, 3),
      ", duration = ", round(duration, 3),
      ", resp. exists = ", dur$response_exists,
      ", resp. endures = ", dur$response_endures,
      ", TTB120 = ", round(TTB, 3)
    ))
  plt +
    geom_hline(yintercept = y[1], color = "gray20", lty = 2) +
    geom_vline(xintercept = dur$start, color = "steelblue3") +
    geom_vline(xintercept = dur$end, color = "firebrick") +
    geom_vline(xintercept = TTB, color = "orange")
}

# need to do this because 'unlist()' doesn't work for a list of rvars
unlist_rvar <- function(rvar_list) {
  nams <- names(rvar_list)
  J <- length(rvar_list)
  out <- rvar_list[[1]]
  if (J > 1) {
    for (j in 2:J) {
      out <- c(out, rvar_list[[j]])
    }
  }
  names(out) <- nams
  out
}
