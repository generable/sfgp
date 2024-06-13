
  // Linear interpolation
  // - y = vector to interpolate, length N_in
  // - x = incresing vector, length N_in
  // - x_out = increasing vector, length N_out, values must be in
  //   (min(x),max(x)]
  // Returns a vector, length N_out, corresponding to interpolated
  // values y_out
  vector interp_1d_linear(vector y, vector x, vector x_out) {
    int left = 1;
    int right = 1;
    real w = 1.0;
    real dx = 0.0;
    int N_in = num_elements(x);
    int N_out = num_elements(x_out);
    vector[N_out] y_out = rep_vector(y[1], N_out);
    for (j in 1 : N_out) {
      while (x[right] < x_out[j]) {
        right = right + 1;
      }
      while (x[left + 1] < x_out[j]) {
        left = left + 1;
      }
      dx = x[right] - x[left];
      if(dx <= 0){
        y_out[j] = y[left];
      } else {
        w = (x[right] - x_out[j]) / dx;
        y_out[j] = w * y[left] + (1 - w) * y[right];
      }
    }
    return y_out;
  }

  // Trapezoidal rule
  vector trapezoid(vector x, vector y) {
    int N = num_elements(x);
    vector[N-1] dx = x[2:N] - x[1:(N-1)];
    vector[N-1] ym = 0.5*(y[2:N] + y[1:(N-1)]);
    return(cumulative_sum(dx .* ym));
  }

  // Integrate from 0 to all t_out
  // Using the trapezoidal rule
  // Assumes that v has equidistant evaluations of the hazard at t = 0, ... , T
  vector integrate_1d(vector t_out, vector t_grid, vector y_grid) {
    int N_grid = num_elements(t_grid);
    vector[N_grid] Y_grid = rep_vector(0.0, N_grid);
    Y_grid[2:N_grid] = trapezoid(t_grid, y_grid);
    return(interp_1d_linear(Y_grid, t_grid, t_out));
  }

  // Compute numeric integral for each individual
  vector integrate_1d_many(vector t_out, vector t_grid, vector y_grid,
         array[,] int inds_out, array[,] int inds_grid) {
      int G = size(inds_out);
      int N_out = num_elements(t_out);
      vector[N_out] Y_int;
      for(g in 1:G){
        int iso = inds_out[g, 1];
        int ieo = inds_out[g, 2];
        int isg = inds_grid[g, 1];
        int ieg = inds_grid[g, 2];
        Y_int[iso:ieo] = integrate_1d(
          t_out[iso:ieo], t_grid[isg:ieg], y_grid[isg:ieg]
        );
      }
      return(Y_int);
  }
