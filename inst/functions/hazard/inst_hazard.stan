
  // new instant hazard function
  vector compute_inst_hazard(vector h0, real assoc, vector link,
      real loc, real scale, real delta) {
    int n = num_elements(h0);
    vector[n] link_norm = (exp(link) + delta - loc) / scale;
    return(h0 .* exp(assoc * link_norm));
  }

  // old instant hazard function
  //vector compute_inst_hazard(vector h0, real assoc, vector link) {
  //  return(h0 .* exp(assoc * link));
  //}
