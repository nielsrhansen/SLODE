// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


// Problem with this is that it does not constrain parameters
arma::vec ridge(arma::mat &x, arma::vec &y, double &Uy_n2, arma::vec x_n, const double thr = 1e-7) {
  /*
   * Finds ridge-limit lambda->0 of beta using SVD
   * x_n is column norms of x (should be calculated once upfront)
   * Uy_n2 is squared norm of U'y in SVD (recall rss = sq_norm(y) - sq_norm(U'y), with near-zero singular values are removed as columns of U)
   */
  arma::mat U;
  arma::vec s;
  arma::mat V;

  arma::svd(U, s, V, x);

  // Truncate
  arma::uvec keepers = arma::find(arma::abs(s) % x_n / x.n_rows > thr);
  U = U.cols(keepers);
  V = V.cols(keepers);
  s = s.elem(keepers);

  // Get U'y
  arma::vec Uy = U.t() * y;

  // Its norm
  Uy_n2 = arma::sum(arma::square(Uy));

  // Manipulat to get \hat beta (recal \hat\beta = V(D^-1)U'y, where if D has zeros, the entries in D^-1 are also zero)
  Uy /= s;
  return V * Uy;
}

arma::sp_mat modelfits(arma::mat x, arma::vec y, arma::umat models) {
  /*
   * x is full design matrix  (nxp)
   * y is observations        (n)
   * models is 0/1 matrix     (mxp) rows indicating what to include in each model
   *
   * return is sparse mx(p+1)-matrix first col is rss, remainding are parameters estimates
   */

  unsigned int n = y.n_elem, p = models.n_cols, m = models.n_rows;
  if (n != x.n_rows) Rcpp::stop("Number of elements in y does not match number of rows in x.");
  if (p != x.n_cols) Rcpp::stop("Number of columns in models does not match number of columns in x.");

  // For filling sp_mat return
  arma::uvec fmodels = arma::find(models);
  unsigned int n_betas = fmodels.n_elem;
  arma::umat locations(2, n_betas + m);
  arma::vec values(n_betas + m);

  // For getting rss
  double Uy_n2 = 0.0, y_n2 = arma::sum(arma::square(y));

  // x colnorms
  arma::vec x_n = arma::sqrt(arma::sum(arma::square(x), 0).t());

  unsigned int ind = 0;
  for (unsigned int i_model = 0; i_model < m; i_model++) {
    // Location (for each model (=row) the first column holds rss, the next coef estimates)
    arma::uvec i_loc  = arma::find(models.row(i_model));
    arma::mat x_sub   = x.cols(i_loc);
    arma::vec x_n_sub = x_n.elem(i_loc);
    arma::vec beta    = ridge(x_sub, y, Uy_n2, x_n_sub);// old one: arma::solve(x_sub, y);

    // Fill model location (i.e., fix row)
    locations(0, arma::span(ind, ind + i_loc.n_elem)).fill(i_model);

    // Fill in location and value for rss
    locations(1, ind) = 0;
    values(ind) = y_n2 - Uy_n2;

    // Fill in location and value for beta
    if (!i_loc.is_empty()) {
      locations(1, arma::span(ind + 1, ind + i_loc.n_elem)) = i_loc.t() + 1;
      values(arma::span(ind + 1, ind + i_loc.n_elem)) = beta;
    }

    // Advance index
    ind += i_loc.n_elem + 1;
  }

  arma::sp_mat ret(locations, values, m, p + 1);

  return ret;
}
