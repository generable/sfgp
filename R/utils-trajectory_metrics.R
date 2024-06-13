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

# Duration
response_dur <- function(t, y) {
  d <- duration_of_response(t, y)
  d$end - d$start
}

# Duration
response_end <- function(t, y) {
  duration_of_response(t, y)$end
}

# Duration
response_start <- function(t, y) {
  duration_of_response(t, y)$start
}

# Duration
response_exists <- function(t, y) {
  duration_of_response(t, y)$response_exists
}

# Duration
response_endures <- function(t, y) {
  duration_of_response(t, y)$response_endures
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
