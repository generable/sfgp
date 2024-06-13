# Union of data frames
df_union <- function(x, y) {
  checkmate::assert_class(x, "data.frame")
  checkmate::assert_class(y, "data.frame")
  checkmate::assert_true(nrow(x) == nrow(y))
  common_names <- intersect(names(x), names(y))
  for (name in common_names) {
    if (!all(x[[name]] == y[[name]], na.rm = TRUE)) {
      stop(paste0("Columns '", name, "' are not identical!"))
    }
  }
  rem <- !(names(y) %in% common_names)
  rem_names <- colnames(y)[rem]
  all_names <- c(colnames(x), rem_names)
  combined_df <- cbind(x, y[, rem])
  colnames(combined_df) <- all_names
  return(combined_df)
}

# Remove certain suffixes from strings
replicate_without_suffix <- function(strings, sfx) {
  regex <- paste0("_", sfx, "$")
  # Identify strings that end with suffix
  has_suffix <- grepl(regex, strings)

  # Remove suffix and replicate the string if it has the suffix, else
  # return original string
  ifelse(has_suffix, gsub(regex, "", strings), strings)
}

# Colorize string
colorize_string <- function(x, col) {
  if (interactive()) {
    x <- paste0(col, x, "\u001b[0m")
  }
  x
}

# Number string
number_string <- function(x) {
  col <- "\u001b[34;1m" # bold blue
  colorize_string(x, col)
}

# Highlight string
hl_string <- function(x) {
  col <- "\u001b[33m" # orange
  colorize_string(x, col)
}

# Highlight class
hl_class <- function(x) {
  col <- "\u001b[31m" # red
  colorize_string(x, col)
}

# Indices of data frame columns that are factors
df_factor_cols <- function(df) {
  which(sapply(df, inherits, "factor"))
}

# Only unique rows of the factors of a data frame
df_unique_factor_rows <- function(df) {
  col_inds <- df_factor_cols(df)
  if (is.null(col_inds)) {
    stop("no factor columns found in data frame")
  }
  df <- unique(df[, col_inds, drop = FALSE])
  rownames(df) <- NULL
  df
}

# Replicate df rows
df_replicate_rows <- function(df, n) {
  return(df[rep(seq_len(nrow(df)), each = n), , drop = FALSE])
}

# Base R data frame row filter
df_filter_rows <- function(df, filter_by = NULL, kept_vals = NULL) {
  if (!is.null(filter_by)) {
    rows_keep <- which(df[[filter_by]] %in% kept_vals)
    df <- df[rows_keep, , drop = FALSE]
  }
  df
}

# Sample N data rows from each group
sample_rows_by_group <- function(data, N = 1, group_var = "arm") {
  checkmate::assert_integerish(N, lower = 1)
  group_fac <- data[[group_var]]
  checkmate::assert_factor(group_fac)
  levs <- unique(levels(group_fac))
  new_dat <- NULL
  for (lev in levs) {
    inds <- which(group_fac == lev)
    new_rows <- pick_random_rows(data[inds, , drop = F], size = N)
    new_dat <- rbind(new_dat, new_rows)
  }
  new_dat
}


# Pick random row of a data frame
pick_random_rows <- function(df, size = 1) {
  idx <- sample(nrow(df), size = size, replace = FALSE)
  df[idx, , drop = FALSE]
}

# Visualize an 'rvar' vector by arm
plot_rvar_arm <- function(rvar, quantity_name = "value",
                          horizontal = FALSE,
                          halfeye = FALSE, ...) {
  df <- create_df_by_arm(rvar, quantity_name)
  plot_rvar_df(df, quantity_name, "arm", horizontal,
    halfeye = halfeye, ...
  )
}

# Helper
create_df_by_arm <- function(rvar, name) {
  arm <- names(rvar)
  if (length(arm) == 0) {
    arm <- paste0("Arm ", seq_len(length(rvar)))
  }
  df <- data.frame(arm, rvar)
  colnames(df) <- c("arm", name)
  df$arm <- as.factor(arm)
  df
}


# automatic selection of prediction time points
t_pred_auto <- function(t_pred = NULL, t_data = NULL) {
  if (!is.null(t_pred)) {
    if (length(t_pred) == 1) {
      if (t_pred == "auto") {
        t_max <- min(2 * 365, max(t_data))
        message(
          "t_pred was 'auto', setting it to 10 evenly spaced points with",
          " t_max = ", t_max
        )
        t_pred <- seq(0, t_max, length.out = 10)
      } else {
        stop("invalid t_pred argument")
      }
    } else {
      checkmate::assert_numeric(t_pred, min.len = 2)
    }
  } else {
    message("t_pred was NULL, setting it based on data")
    t_range <- range(t_data)
    t_min <- min(0, t_range[1])
    t_max <- 1.1 * t_range[2]
    t_pred <- seq(t_min, t_max, length.out = 40)
  }
  t_pred
}

# Finalizer
finalize_stan_data <- function(stan_data, prior_only) {
  stan_data <- stan_data[unique(names(stan_data))] # remove duplicates
  checkmate::assert_logical(prior_only)
  stan_data$prior_only <- as.numeric(prior_only)
  stan_data
}

# Main class name
class_name <- function(x) {
  class(x)[1]
}

# Main class name
class_name_hl <- function(x) {
  hl_class(class_name(x))
}
