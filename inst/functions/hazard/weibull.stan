
  // revised Weibull hazard function
  vector weibull_haz(vector x, real lambda, real gamma) {
    return(gamma * lambda * x ^ (gamma - 1.0));
  }
