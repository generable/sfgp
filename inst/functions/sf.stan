
  // Log Stein Fojo (without ts0)
  vector log_sf_norm(vector t, real kg, real ks) {
    return(log_diff_exp(log_sum_exp(kg .* t, -ks .* t), 0.0));
  }

  // Vectorized log Stein Fojo (without ts0)
  // - all arguments should have the same length
  vector log_sf_norm_vec(vector t, vector kg, vector ks) {
    return(log_diff_exp(log_sum_exp(kg .* t, -ks .* t), 0.0));
  }
