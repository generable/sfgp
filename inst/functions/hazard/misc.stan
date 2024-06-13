
  // Cumulative mean of vector
  vector cumulative_mean(vector x) {
    int n = num_elements(x);
    return(cumulative_sum(x) ./ seq_len(n));
  }

  // Linear 1d interpolation
  vector interp_1d(vector y, vector x, vector x_out) {
    int left = 1;
    int right = 1;
    real w = 1.0;
    int N_in = num_elements(x);
    int N_out = num_elements(x_out);
    vector[N_out] y_out;
    for (j in 1 : N_out) {
      while (x[right] < x_out[j]) {
        right = right + 1;
      }
      while (x[left + 1] < x_out[j]) {
        left = left + 1;
      }
      w = (x[right] - x_out[j]) / (x[right] - x[left]);
      y_out[j] = w * y[left] + (1 - w) * y[right];
    }
    return(y_out);
  }
