#' Simulate from SF model (3 arms)
#'
#' @export
#' @param n1 number of subjects in arm 1
#' @param n2 number of subjects in arm 2
#' @param n3 number of subjects in arm 3
#' @param ts0_scale initial tumor sizes will be sampled from the uniform
#' distribution \code{ts0_scale * U(1,3)}
#' @param mu_ks Mean \code{k_s} parameters for the three arms.
#' @param mu_kg Mean \code{k_g} parameters for the three arms.
#' @param times Time points where measurements are made. Scatter is added to
#' these.
#' @param kg_sd Jitter added to kg.
#' @param ks_sd Jitter added to ks.
#' @param sigma Noise param value.
#' @param x_effect Magnitude of effect of x on ks.
#' @param log_C_kg offset for log kg
#' @param log_C_ks offset for log ks
#' @param t_sd Jitter added to time points.
#' @return a data frame
sfsim <- function(n1 = 4, n2 = 4, n3 = 4, ts0_scale = 20,
                  mu_ks = c(0.01, 2, 4),
                  mu_kg = c(3, 1.5, 0.05),
                  times = c(0.5, 3, 6, 12, 24, 48),
                  kg_sd = 0.2,
                  ks_sd = 0.2,
                  sigma = 0.2,
                  x_effect = 0,
                  log_C_kg = -7,
                  log_C_ks = -5,
                  t_sd = 0.2) {
  df <- sim_id_and_arm(n1, n2, n3)
  df_final <- NULL
  MU_KG <- mu_kg
  MU_KS <- mu_ks
  for (j in 1:nrow(df)) {
    t_j <- times
    t_j <- sort(t_j + t_sd * rnorm(length(t_j)))
    n_j <- length(t_j)
    ID <- df$id[j]
    ARM <- df$arm[j]
    X <- stats::runif(1)
    id_j <- rep(ID, n_j)
    arm_j <- rep(ARM, n_j)
    x_j <- rep(X, n_j)
    ts0 <- ts0_scale * runif(1, 1, 3)
    # f_x_ks <- x_effect * (sin(6 * X) + 0.5)
    f_x_ks <- sim_emax(X, x_effect, 0.35, 9)
    kg <- exp(log_C_kg + rnorm(1, mean = MU_KG[ARM], sd = kg_sd))
    ks <- exp(log_C_ks + rnorm(1, mean = MU_KS[ARM], sd = ks_sd) + f_x_ks)
    kg_vec <- rep(kg, n_j)
    ks_vec <- rep(ks, n_j)
    f_j <- sfsim_one(ts0, kg_vec, ks_vec, t_j)
    y_j <- rlnorm(n = length(f_j), meanlog = log(f_j), sdlog = sigma)
    df_j <- data.frame(
      id = id_j, arm = arm_j, t = t_j, x = x_j, f_x_ks = f_x_ks, f = f_j,
      y = y_j, true_ks = ks, true_kg = kg
    )
    df_final <- rbind(df_final, df_j)
  }
  df_final
}

# Emax
sim_emax <- function(x, emax, ed50, gamma) {
  xg <- x^gamma
  (emax * xg) / (ed50^gamma + xg)
}


# Simulate one trajectory from SF
sfsim_one <- function(ts0, kg, ks, t) {
  ts0 * (exp(kg * t) + exp(-ks * t) - 1.0)
}


# Simulate id and arm factors
sim_id_and_arm <- function(n1, n2, n3) {
  arm1 <- rep(1, n1)
  arm2 <- rep(2, n2)
  arm3 <- rep(3, n3)
  arm <- c(arm1, arm2, arm3)
  id <- c(1:n1, n1 + (1:n2), n1 + n2 + (1:n2))
  df <- data.frame(id, arm)
  df$id <- as.factor(df$id)
  df$arm <- as.factor(df$arm)
  df
}
