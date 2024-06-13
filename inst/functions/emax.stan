
  // EMAX (without E0)
  vector emax(vector x, vector log_ED50, real Emax, real gamma) {
    int N = num_elements(x);
    vector[N] xg = x^gamma;
    return((Emax * xg) ./ (exp(log_ED50)^gamma + xg));
  }
