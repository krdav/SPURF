#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


mat logistic(const mat& t, double M, double B) {
  return 1 / (1 + exp(-B * (t - M)));
}

// [[Rcpp::export]]
mat compute_kmeans_clust_profs(const mat& dff_kmeans_centers,
                               const mat& dff2_subsamp) {
  // 1) compute distances from subsampled profiles to all cluster centroids
  // 2) assign cluster profile estimates to subsampled profiles
  mat clust_profs(size(dff2_subsamp), fill::none);

  for (unsigned int i = 0; i < dff2_subsamp.n_rows; i++) {
    vec dist_vec(dff_kmeans_centers.n_rows, fill::none);

    for (unsigned int j = 0; j < dff_kmeans_centers.n_rows; j++) {
      dist_vec(j) = norm(dff2_subsamp.row(i) - dff_kmeans_centers.row(j), 1);
    }

    clust_profs.row(i) =
        dff_kmeans_centers.rows(find(dist_vec == min(dist_vec)));
  }

  return clust_profs;
}


///////////////////////// PREDICT PROFILES


// [[Rcpp::export]]
mat predict_profs(const rowvec& alpha, const urowvec& alpha_cols_assign,
                  const mat& dff2_subsamp, const List& dff2_profs,
                  unsigned int alpha_groups) {
  mat alpha_mat(alpha.begin(), alpha_groups, dff2_profs.size());
  mat pred_profs(size(dff2_subsamp), fill::zeros);

  for (uword i = 0; i < dff2_profs.size(); i++) {
    mat prof_i = dff2_profs[i];
    pred_profs +=
        prof_i.each_row() % alpha_mat(alpha_cols_assign - 1, uvec({i})).t();
  }
  pred_profs += dff2_subsamp.each_row() %
                (1 - sum(alpha_mat.rows(alpha_cols_assign - 1), 1).t());

  return pred_profs;
}


///////////////// F


// [[Rcpp::export]]
double f(const rowvec& alpha, const urowvec& alpha_cols_assign,
         const mat& dff2_subsamp, const mat& dff2_full, const List& dff2_profs,
         unsigned int alpha_groups) {
  mat pred_profs = predict_profs(alpha, alpha_cols_assign, dff2_subsamp,
                                 dff2_profs, alpha_groups);
  return accu(square(dff2_full - pred_profs)) / (149 * dff2_full.n_rows);
}


///////////////// PENALIZED F


// [[Rcpp::export]]
double f_penal(const rowvec& alpha, const urowvec& alpha_cols_assign,
               const mat& dff2_subsamp, const mat& dff2_full,
               const List& dff2_profs, unsigned int alpha_groups, int d,
               double lambda1, double lambda2) {
  double err = f(alpha, alpha_cols_assign, dff2_subsamp, dff2_full, dff2_profs,
                 alpha_groups);
  double coef_penal = norm(alpha, 1);
  double fusion_penal = 0;
  for (uword i = 0; i < dff2_profs.size(); i++) {
    fusion_penal += norm(
        diff(alpha(span(i * alpha_groups, (i + 1) * alpha_groups - 1)), d), 1);
  }

  return err + lambda1 * coef_penal + lambda2 * fusion_penal;
}


///////////// JACCARD F


// [[Rcpp::export]]
mat calculate_jacc(const mat& pred_profs, const mat& dff2_full, double cutoff,
                   double B) {
  mat pred_sets = logistic(pred_profs, cutoff, B);
  mat full_sets = logistic(dff2_full, cutoff, B);
  mat jacc_mat(dff2_full.n_rows, 149, fill::none);

  for (uword i = 0; i < dff2_full.n_rows; i++) {
    for (uword j = 0; j < 149; j++) {
      if (all(dff2_full(uvec({i}), regspace<uvec>(21 * j, 21 * (j + 1) - 1)) ==
              rowvec({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 1}))) {
        jacc_mat(i, j) = 0;
        continue;
      }

      const rowvec& pred_set =
          pred_sets(uvec({i}), regspace<uvec>(21 * j, 21 * (j + 1) - 1));
      const rowvec& full_set =
          full_sets(uvec({i}), regspace<uvec>(21 * j, 21 * (j + 1) - 1));

      double num = dot(pred_set, full_set);
      rowvec denom_vec = pred_set + full_set;
      denom_vec(find(denom_vec > 1)).ones();
      double denom = sum(denom_vec);

      if (num == 0 && denom == 0) {
        jacc_mat(i, j) = 1;
      } else {
        jacc_mat(i, j) = num / denom;
      }
    }
  }

  return -jacc_mat;
}


// [[Rcpp::export]]
double f_jacc(const rowvec& alpha, const urowvec& alpha_cols_assign,
              const mat& dff2_subsamp, const mat& dff2_full,
              const List& dff2_profs, unsigned int alpha_groups, double cutoff,
              double B) {
  mat pred_profs = predict_profs(alpha, alpha_cols_assign, dff2_subsamp,
                                 dff2_profs, alpha_groups);
  return accu(calculate_jacc(pred_profs, dff2_full, cutoff, B)) /
         (149 * dff2_full.n_rows);
}


/////////////// PENALIZED JACCARD F


// [[Rcpp::export]]
double f_penal_jacc(const rowvec& alpha, const urowvec& alpha_cols_assign,
                    const mat& dff2_subsamp, const mat& dff2_full,
                    const List& dff2_profs, unsigned int alpha_groups,
                    double cutoff, double B, int d, double lambda1,
                    double lambda2) {
  double err = f_jacc(alpha, alpha_cols_assign, dff2_subsamp, dff2_full,
                      dff2_profs, alpha_groups, cutoff, B);
  double coef_penal = norm(alpha, 1);
  double fusion_penal = 0;
  for (uword i = 0; i < dff2_profs.size(); i++) {
    fusion_penal += norm(
        diff(alpha(span(i * alpha_groups, (i + 1) * alpha_groups - 1)), d), 1);
  }

  return err + lambda1 * coef_penal + lambda2 * fusion_penal;
}


///////////// PREDICT NOGAP PROFILES


// [[Rcpp::export]]
mat predict_nogap_profs(const rowvec& alpha, const urowvec& alpha_cols_assign,
                        const mat& dff2_subsamp, const List& dff2_profs,
                        const umat& dff2_nogap, unsigned int alpha_groups) {
  mat alpha_mat(alpha.begin(), alpha_groups, dff2_profs.size());
  mat pred_profs(size(dff2_subsamp), fill::zeros);
  pred_profs.cols(regspace<uvec>(20, 21, pred_profs.n_cols - 1)).fill(1);

  for (unsigned int i = 0; i < dff2_subsamp.n_rows; i++) {
    uvec col_inds = find(dff2_nogap.row(i));
    pred_profs(uvec({i}), col_inds).zeros();

    for (uword j = 0; j < dff2_profs.size(); j++) {
      mat prof_j = dff2_profs[j];
      pred_profs(uvec({i}), col_inds) +=
          prof_j(uvec({i}), col_inds) %
          alpha_mat(alpha_cols_assign - 1, uvec({j})).eval()(col_inds).t();
    }
    pred_profs(uvec({i}), col_inds) +=
        dff2_subsamp(uvec({i}), col_inds) %
        (1 - sum(alpha_mat.rows(alpha_cols_assign - 1).eval().rows(col_inds), 1)
                 .t());
  }

  return pred_profs;
}


//////////////// NOGAP F


// [[Rcpp::export]]
double f_nogap(const rowvec& alpha, const urowvec& alpha_cols_assign,
               const mat& dff2_subsamp, const mat& dff2_full,
               const List& dff2_profs, const umat& dff2_nogap,
               unsigned int alpha_groups) {
  mat pred_profs = predict_nogap_profs(alpha, alpha_cols_assign, dff2_subsamp,
                                       dff2_profs, dff2_nogap, alpha_groups);
  return accu(square(dff2_full - pred_profs)) / (accu(dff2_nogap) / 21);
}


////////////////// NOGAP JACCARD F


// [[Rcpp::export]]
double f_nogap_jacc(const rowvec& alpha, const urowvec& alpha_cols_assign,
                    const mat& dff2_subsamp, const mat& dff2_full,
                    const List& dff2_profs, const umat& dff2_nogap,
                    unsigned int alpha_groups, double cutoff, double B) {
  mat pred_profs = predict_nogap_profs(alpha, alpha_cols_assign, dff2_subsamp,
                                       dff2_profs, dff2_nogap, alpha_groups);
  return accu(calculate_jacc(pred_profs, dff2_full, cutoff, B)) /
         (accu(dff2_nogap) / 21);
}
