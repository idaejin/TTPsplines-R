// TT P-spline hot paths (general d) — RcppArmadillo
// Compile via: Rcpp::sourceCpp("src/tt_pspline_nd.cpp")
//
// Core layout matches R: array dim c(r_left, p, r_right), vec order (a, j, b)
// with a fastest.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

namespace {

inline arma::mat core_slice_j(const arma::cube& core, int j) {
  // core: rl x p x rr  -> Cj: rl x rr at basis index j
  const int rl = core.n_rows;
  const int rr = core.n_slices;
  arma::mat Cj(rl, rr);
  for (int b = 0; b < rr; ++b) {
    for (int a = 0; a < rl; ++a) {
      Cj(a, b) = core(a, j, b);
    }
  }
  return Cj;
}

inline double ridge_scale_arma(const arma::mat& xtx, double multiplier = 1e-7) {
  double scale = arma::mean(xtx.diag());
  if (!std::isfinite(scale) || scale <= 0.0) scale = 1.0;
  return multiplier * scale;
}

}  // namespace

// [[Rcpp::export]]
arma::mat tt_design_core_d_cpp(const arma::mat& Left,
                               const arma::mat& Right,
                               const arma::mat& Bk) {
  const int n = Bk.n_rows;
  const int rl = Left.n_cols;
  const int rr = Right.n_cols;
  const int p = Bk.n_cols;
  arma::mat X(n, rl * p * rr);
  int col = 0;
  for (int b = 0; b < rr; ++b) {
    for (int j = 0; j < p; ++j) {
      for (int a = 0; a < rl; ++a) {
        X.col(col++) = Left.col(a) % Bk.col(j) % Right.col(b);
      }
    }
  }
  return X;
}

// [[Rcpp::export]]
arma::mat contract_left_step_cpp(const arma::mat& left,
                                 const arma::cube& core,
                                 const arma::mat& Bk) {
  const int n = Bk.n_rows;
  const int p = core.n_cols;
  const int rr = core.n_slices;
  arma::mat out(n, rr, arma::fill::zeros);
  for (int j = 0; j < p; ++j) {
    arma::mat tmp = left * core_slice_j(core, j);
    tmp.each_col() %= Bk.col(j);
    out += tmp;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat contract_right_step_cpp(const arma::mat& right,
                                  const arma::cube& core,
                                  const arma::mat& Bk) {
  const int n = Bk.n_rows;
  const int p = core.n_cols;
  const int rl = core.n_rows;
  arma::mat out(n, rl, arma::fill::zeros);
  for (int j = 0; j < p; ++j) {
    arma::mat Cj = core_slice_j(core, j);          // rl x rr
    arma::mat tmp = right * Cj.t();                // n x rl
    tmp.each_col() %= Bk.col(j);
    out += tmp;
  }
  return out;
}

// [[Rcpp::export]]
arma::vec gaussian_core_update_cpp(const arma::vec& yc,
                                   const arma::mat& Left,
                                   const arma::mat& Right,
                                   const arma::mat& Bk,
                                   const arma::mat& penalty,
                                   double lambda) {
  arma::mat X = tt_design_core_d_cpp(Left, Right, Bk);
  arma::mat xtx = X.t() * X;
  const double ridge = ridge_scale_arma(xtx);
  arma::mat system = xtx + lambda * penalty;
  system.diag() += ridge;
  arma::vec rhs = X.t() * yc;
  return arma::solve(system, rhs, arma::solve_opts::likely_sympd);
}

// Build left interfaces: list of length d; L[[0]] = ones(n,1)
static std::vector<arma::mat> build_left_interfaces(
    const std::vector<arma::cube>& cores,
    const std::vector<arma::mat>& basis) {
  const int d = cores.size();
  const int n = basis[0].n_rows;
  std::vector<arma::mat> L(d);
  L[0] = arma::ones(n, 1);
  arma::mat cur = L[0];
  for (int k = 0; k < d - 1; ++k) {
    cur = contract_left_step_cpp(cur, cores[k], basis[k]);
    L[k + 1] = cur;
  }
  return L;
}

static std::vector<arma::mat> build_right_interfaces(
    const std::vector<arma::cube>& cores,
    const std::vector<arma::mat>& basis) {
  const int d = cores.size();
  const int n = basis[0].n_rows;
  std::vector<arma::mat> R(d);
  R[d - 1] = arma::ones(n, 1);
  arma::mat cur = R[d - 1];
  for (int k = d - 1; k >= 1; --k) {
    cur = contract_right_step_cpp(cur, cores[k], basis[k]);
    R[k - 1] = cur;
  }
  return R;
}

static std::vector<arma::cube> list_to_cubes(const List& cores_list) {
  const int d = cores_list.size();
  std::vector<arma::cube> cores(d);
  for (int k = 0; k < d; ++k) {
    cores[k] = as<arma::cube>(cores_list[k]);
  }
  return cores;
}

static std::vector<arma::mat> list_to_mats(const List& mat_list) {
  const int d = mat_list.size();
  std::vector<arma::mat> out(d);
  for (int k = 0; k < d; ++k) {
    out[k] = as<arma::mat>(mat_list[k]);
  }
  return out;
}

static List cubes_to_list(const std::vector<arma::cube>& cores) {
  List out(cores.size());
  for (std::size_t k = 0; k < cores.size(); ++k) {
    out[k] = cores[k];
  }
  return out;
}

// [[Rcpp::export]]
arma::vec tt_contraction_d_cpp(const List& cores_list, const List& basis_list) {
  std::vector<arma::cube> cores = list_to_cubes(cores_list);
  std::vector<arma::mat> basis = list_to_mats(basis_list);
  const int n = basis[0].n_rows;
  arma::mat cur = arma::ones(n, 1);
  for (std::size_t k = 0; k < cores.size(); ++k) {
    cur = contract_left_step_cpp(cur, cores[k], basis[k]);
  }
  return cur.col(0);
}

//' Full TT-ALS fit (Gaussian, fixed anisotropic λ).
//' Returns cores, mu, optional jacobian, block sizes.
// [[Rcpp::export]]
List tt_fit_d_cpp(const arma::vec& y,
                  const List& basis_list,
                  const List& init_cores,
                  const arma::vec& lambda,
                  const List& penalties_list,
                  int sweeps,
                  bool return_jacobian = true) {
  std::vector<arma::mat> basis = list_to_mats(basis_list);
  std::vector<arma::cube> cores = list_to_cubes(init_cores);
  std::vector<arma::mat> penalties = list_to_mats(penalties_list);
  const int d = cores.size();
  const int n = y.n_elem;
  if (static_cast<int>(basis.size()) != d) {
    stop("basis and cores length mismatch");
  }
  if (static_cast<int>(lambda.n_elem) != d) {
    stop("lambda length must equal d");
  }

  const double intercept = arma::mean(y);
  arma::vec yc = y - intercept;

  for (int sweep = 0; sweep < sweeps; ++sweep) {
    // Incremental: rebuild L once per sweep start; update after each core
    std::vector<arma::mat> L = build_left_interfaces(cores, basis);
    std::vector<arma::mat> R = build_right_interfaces(cores, basis);
    for (int k = 0; k < d; ++k) {
      arma::vec coef = gaussian_core_update_cpp(
          yc, L[k], R[k], basis[k], penalties[k], lambda[k]);
      const int rl = cores[k].n_rows;
      const int p = cores[k].n_cols;
      const int rr = cores[k].n_slices;
      cores[k] = arma::cube(coef.memptr(), rl, p, rr);
      // refresh left interface for next cores
      if (k + 1 < d) {
        L[k + 1] = contract_left_step_cpp(L[k], cores[k], basis[k]);
      }
      // refresh right interfaces for earlier cores (if needed later in sweep)
      // When going left-to-right, R[k] for current/future does not depend on core k.
      // So no R update needed within forward sweep.
    }
  }

  arma::vec mu = intercept + tt_contraction_d_cpp(cubes_to_list(cores), basis_list);

  List out = List::create(
      _["cores"] = cubes_to_list(cores),
      _["intercept"] = intercept,
      _["mu"] = mu,
      _["lambda"] = lambda,
      _["d"] = d,
      _["sweeps"] = sweeps);

  if (return_jacobian) {
    std::vector<arma::mat> L = build_left_interfaces(cores, basis);
    std::vector<arma::mat> R = build_right_interfaces(cores, basis);
    int m = 0;
    IntegerVector block_sizes(d);
    for (int k = 0; k < d; ++k) {
      block_sizes[k] = cores[k].n_rows * cores[k].n_cols * cores[k].n_slices;
      m += block_sizes[k];
    }
    arma::mat J(n, m);
    int offset = 0;
    for (int k = 0; k < d; ++k) {
      arma::mat Xk = tt_design_core_d_cpp(L[k], R[k], basis[k]);
      J.cols(offset, offset + Xk.n_cols - 1) = Xk;
      offset += Xk.n_cols;
    }
    out["jacobian"] = J;
    out["block_sizes"] = block_sizes;
  }

  return out;
}

//' Effective DF: tr[(J'J + P)^{-1} J'J]
// [[Rcpp::export]]
double effective_df_cpp(const arma::mat& jacobian, const arma::mat& penalty) {
  arma::mat xtx = jacobian.t() * jacobian;
  const double ridge = ridge_scale_arma(xtx, 1e-9);
  arma::mat system = xtx + penalty;
  system.diag() += ridge;
  arma::mat infl = arma::solve(system, xtx, arma::solve_opts::likely_sympd);
  return arma::trace(infl);
}

// =====================================================================
// Weighted ALS + conditional GCV (Gaussian / Bernoulli / Poisson)
// =====================================================================

namespace {

arma::vec solve_ridge_escalating(const arma::mat& A, const arma::vec& b,
                                 double base_ridge) {
  const int m = A.n_rows;
  const double scale = (base_ridge > 0.0) ? base_ridge
                                          : ridge_scale_arma(A, 1e-4);
  const double facs[8] = {1.0, 10.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e8};
  for (int i = 0; i < 8; ++i) {
    arma::mat M = A;
    M.diag() += facs[i] * scale;
    arma::mat R;
    if (arma::chol(R, M)) {
      arma::vec y = arma::solve(arma::trimatl(R.t()), b);
      arma::vec g = arma::solve(arma::trimatu(R), y);
      if (g.is_finite()) return g;
    }
  }
  arma::mat M = A;
  M.diag() += std::max(1e-2, 1e8 * scale);
  arma::vec g = arma::solve(M, b, arma::solve_opts::fast);
  return g;
}

arma::mat solve_ridge_mat(const arma::mat& A, const arma::mat& B,
                          double base_ridge) {
  const double scale = (base_ridge > 0.0) ? base_ridge
                                          : ridge_scale_arma(A, 1e-4);
  const double facs[7] = {1.0, 10.0, 1e2, 1e3, 1e4, 1e5, 1e6};
  for (int i = 0; i < 7; ++i) {
    arma::mat M = A;
    M.diag() += facs[i] * scale;
    arma::mat R;
    if (arma::chol(R, M)) {
      arma::mat Y = arma::solve(arma::trimatl(R.t()), B);
      arma::mat out = arma::solve(arma::trimatu(R), Y);
      if (out.is_finite()) return out;
    }
  }
  arma::mat M = A;
  M.diag() += std::max(1e-2, 1e6 * scale);
  return arma::solve(M, B, arma::solve_opts::fast);
}

double edf_S(const arma::mat& S, const arma::mat& P, double lambda) {
  const double ridge = ridge_scale_arma(S, 1e-6);
  arma::mat A = S + lambda * P;
  arma::mat MinvS = solve_ridge_mat(A, S, ridge);
  return arma::trace(MinvS);
}

arma::vec solve_gb(const arma::mat& S, const arma::mat& P,
                   const arma::vec& b, double lambda) {
  const double ridge = ridge_scale_arma(S, 1e-6);
  return solve_ridge_escalating(S + lambda * P, b, ridge);
}

struct GcvEval {
  double value;
  arma::vec g;
  double ed;
  double rss;
};

GcvEval conditional_gcv_eval(const arma::vec& yw,
                             const arma::mat& Xw,
                             const arma::mat& S,
                             const arma::mat& P,
                             const arma::vec& b,
                             double lambda) {
  GcvEval out;
  const int n = yw.n_elem;
  out.g = solve_gb(S, P, b, lambda);
  arma::vec resid = yw - Xw * out.g;
  out.rss = arma::dot(resid, resid);
  out.ed = edf_S(S, P, lambda);
  const double denom = (n - out.ed) * (n - out.ed);
  out.value = (!std::isfinite(denom) || denom < 1e-12)
                  ? R_PosInf
                  : (n * out.rss / denom);
  return out;
}

// Golden-section minimize of f on [a,b]
template <typename Fun>
double golden_minimize(Fun f, double a, double b, double tol) {
  const double resphi = 2.0 - (1.0 + std::sqrt(5.0)) / 2.0;
  double x1 = a + resphi * (b - a);
  double x2 = b - resphi * (b - a);
  double f1 = f(x1);
  double f2 = f(x2);
  int guard = 0;
  while (std::abs(b - a) > tol && guard++ < 80) {
    if (f1 < f2) {
      b = x2;
      x2 = x1;
      f2 = f1;
      x1 = a + resphi * (b - a);
      f1 = f(x1);
    } else {
      a = x1;
      x1 = x2;
      f1 = f2;
      x2 = b - resphi * (b - a);
      f2 = f(x2);
    }
  }
  return 0.5 * (a + b);
}

struct OptLam {
  double lambda;
  arma::vec g;
  double value;
  double ed;
  int n_eval;
};

OptLam optimize_lambda_gcv(const arma::vec& yw,
                           const arma::mat& Xw,
                           const arma::mat& S,
                           const arma::mat& P,
                           const arma::vec& b,
                           double lambda_min,
                           double lambda_max,
                           double tol) {
  OptLam out;
  out.n_eval = 0;
  const double lo = std::log(lambda_min);
  const double hi = std::log(lambda_max);
  auto obj = [&](double ll) {
    ++out.n_eval;
    return conditional_gcv_eval(yw, Xw, S, P, b, std::exp(ll)).value;
  };
  const double ll_hat = golden_minimize(obj, lo, hi, tol);
  const double lam = std::exp(ll_hat);
  GcvEval fit = conditional_gcv_eval(yw, Xw, S, P, b, lam);
  ++out.n_eval;
  out.lambda = lam;
  out.g = fit.g;
  out.value = fit.value;
  out.ed = fit.ed;
  return out;
}

void fill_cube_from_vec(arma::cube& C, const arma::vec& g) {
  const int rl = C.n_rows;
  const int p = C.n_cols;
  const int rr = C.n_slices;
  // R/array column-major: a fastest, then j, then b — matches arma cube
  for (int b = 0; b < rr; ++b) {
    for (int j = 0; j < p; ++j) {
      for (int a = 0; a < rl; ++a) {
        const int idx = a + j * rl + b * rl * p;
        C(a, j, b) = g(idx);
      }
    }
  }
}

void glm_working(const std::string& family,
                 const arma::vec& y,
                 const arma::vec& eta,
                 arma::vec& w,
                 arma::vec& z,
                 arma::vec& mu) {
  const int n = y.n_elem;
  w.set_size(n);
  z.set_size(n);
  mu.set_size(n);
  if (family == "bernoulli" || family == "binomial") {
    for (int i = 0; i < n; ++i) {
      double m = 1.0 / (1.0 + std::exp(-eta(i)));
      if (m < 1e-5) m = 1e-5;
      if (m > 1.0 - 1e-5) m = 1.0 - 1e-5;
      mu(i) = m;
      const double var = m * (1.0 - m);
      w(i) = std::max(var, 1e-4);
      z(i) = eta(i) + (y(i) - m) / var;
    }
  } else if (family == "poisson") {
    for (int i = 0; i < n; ++i) {
      double et = eta(i);
      if (et < -20.0) et = -20.0;
      if (et > 20.0) et = 20.0;
      double m = std::exp(et);
      if (m < 1e-8) m = 1e-8;
      mu(i) = m;
      w(i) = m;
      z(i) = et + (y(i) - m) / m;
    }
  } else {  // gaussian
    mu = eta;
    w.fill(1.0);
    z = y;  // caller should pass centered target via intercept handling
  }
}

double glm_deviance(const std::string& family,
                    const arma::vec& y,
                    const arma::vec& mu) {
  const int n = y.n_elem;
  double dev = 0.0;
  if (family == "bernoulli" || family == "binomial") {
    for (int i = 0; i < n; ++i) {
      double m = mu(i);
      if (m < 1e-12) m = 1e-12;
      if (m > 1.0 - 1e-12) m = 1.0 - 1e-12;
      dev += -2.0 * (y(i) * std::log(m) + (1.0 - y(i)) * std::log(1.0 - m));
    }
  } else if (family == "poisson") {
    for (int i = 0; i < n; ++i) {
      double m = std::max(mu(i), 1e-12);
      const double yi = y(i);
      const double term = (yi > 0.0) ? (yi * std::log(yi / m)) : 0.0;
      dev += 2.0 * (term - (yi - m));
    }
  } else {
    arma::vec r = y - mu;
    dev = arma::dot(r, r);
  }
  return dev;
}

arma::vec glm_invlink(const std::string& family, const arma::vec& eta) {
  const int n = eta.n_elem;
  arma::vec mu(n);
  if (family == "bernoulli" || family == "binomial") {
    for (int i = 0; i < n; ++i) mu(i) = 1.0 / (1.0 + std::exp(-eta(i)));
  } else if (family == "poisson") {
    for (int i = 0; i < n; ++i) {
      double et = eta(i);
      if (et < -20.0) et = -20.0;
      if (et > 20.0) et = 20.0;
      mu(i) = std::exp(et);
    }
  } else {
    mu = eta;
  }
  return mu;
}

double glm_init_intercept(const std::string& family, const arma::vec& y) {
  const double m = arma::mean(y);
  if (family == "bernoulli" || family == "binomial") {
    double p = m;
    if (p < 0.05) p = 0.05;
    if (p > 0.95) p = 0.95;
    return std::log(p / (1.0 - p));
  }
  if (family == "poisson") {
    return std::log(std::max(m, 0.1));
  }
  return m;
}

}  // namespace

//' Build weighted Gram S=X'WX and RHS b=X'W zc (materializes X).
// [[Rcpp::export]]
List weighted_core_system_cpp(const arma::vec& zc,
                              const arma::mat& Left,
                              const arma::mat& Right,
                              const arma::mat& Bk,
                              const arma::vec& weight) {
  arma::mat X = tt_design_core_d_cpp(Left, Right, Bk);
  arma::vec sw = arma::sqrt(weight);
  for (arma::uword i = 0; i < sw.n_elem; ++i) {
    if (!std::isfinite(sw(i)) || sw(i) < 0.0) sw(i) = 0.0;
  }
  arma::mat Xw = X;
  Xw.each_col() %= sw;
  arma::vec zw = sw % zc;
  arma::mat S = Xw.t() * Xw;
  arma::vec b = Xw.t() * zw;
  return List::create(
      _["S"] = S,
      _["b"] = b,
      _["X"] = X,
      _["Xw"] = Xw,
      _["yw"] = zw);
}

//' Conditional GCV value / coef for fixed λ (weighted or unweighted via Xw,yw).
// [[Rcpp::export]]
List conditional_gcv_cpp(const arma::vec& yw,
                         const arma::mat& Xw,
                         const arma::mat& S,
                         const arma::mat& P,
                         const arma::vec& b,
                         double lambda) {
  GcvEval e = conditional_gcv_eval(yw, Xw, S, P, b, lambda);
  return List::create(
      _["value"] = e.value,
      _["g"] = e.g,
      _["ed"] = e.ed,
      _["rss"] = e.rss);
}

//' Brent/golden optimize λ on conditional GCV (log scale).
// [[Rcpp::export]]
List optimize_lambda_gcv_cpp(const arma::vec& yw,
                             const arma::mat& Xw,
                             const arma::mat& S,
                             const arma::mat& P,
                             const arma::vec& b,
                             double lambda_min = 1e-2,
                             double lambda_max = 1e2,
                             double tol = 1e-3) {
  OptLam o = optimize_lambda_gcv(yw, Xw, S, P, b, lambda_min, lambda_max, tol);
  return List::create(
      _["lambda"] = o.lambda,
      _["g"] = o.g,
      _["value"] = o.value,
      _["ed"] = o.ed,
      _["n_eval"] = o.n_eval);
}

//' Weighted core update with optional cGCV λ (criterion: "fixed" | "gcv").
// [[Rcpp::export]]
List weighted_core_update_cgcv_cpp(const arma::vec& zc,
                                   const arma::mat& Left,
                                   const arma::mat& Right,
                                   const arma::mat& Bk,
                                   const arma::vec& weight,
                                   const arma::mat& penalty,
                                   double lambda,
                                   std::string criterion = "gcv",
                                   double lambda_min = 1e-2,
                                   double lambda_max = 1e2,
                                   double tol = 1e-3) {
  List sys = weighted_core_system_cpp(zc, Left, Right, Bk, weight);
  arma::mat S = as<arma::mat>(sys["S"]);
  arma::vec b = as<arma::vec>(sys["b"]);
  arma::mat Xw = as<arma::mat>(sys["Xw"]);
  arma::vec yw = as<arma::vec>(sys["yw"]);
  int n_eval = 0;
  double lam = lambda;
  arma::vec g;
  double value = NA_REAL;
  double ed = NA_REAL;
  if (criterion == "gcv") {
    OptLam o = optimize_lambda_gcv(yw, Xw, S, penalty, b,
                                   lambda_min, lambda_max, tol);
    lam = o.lambda;
    g = o.g;
    value = o.value;
    ed = o.ed;
    n_eval = o.n_eval;
  } else {
    g = solve_gb(S, penalty, b, lam);
    GcvEval e = conditional_gcv_eval(yw, Xw, S, penalty, b, lam);
    value = e.value;
    ed = e.ed;
    n_eval = 1;
  }
  return List::create(
      _["g"] = g,
      _["lambda"] = lam,
      _["value"] = value,
      _["ed"] = ed,
      _["n_eval"] = n_eval);
}

//' Gaussian TT-ALS with per-core conditional GCV (Rcpp).
// [[Rcpp::export]]
List tt_cgcv_fit_cpp(const arma::vec& y,
                     const List& basis_list,
                     const List& init_cores,
                     const List& penalties_list,
                     arma::vec lambda_init,
                     int sweeps = 12,
                     double lambda_min = 1e-2,
                     double lambda_max = 1e2,
                     double tol = 1e-3,
                     double tol_lambda = 1e-3,
                     bool return_jacobian = false) {
  std::vector<arma::mat> basis = list_to_mats(basis_list);
  std::vector<arma::cube> cores = list_to_cubes(init_cores);
  std::vector<arma::mat> penalties = list_to_mats(penalties_list);
  const int d = cores.size();
  const int n = y.n_elem;
  if (static_cast<int>(lambda_init.n_elem) != d) {
    const double lam0 = lambda_init.n_elem > 0 ? lambda_init(0) : 1.0;
    lambda_init = arma::vec(d);
    lambda_init.fill(lam0);
  }
  arma::vec lambda = lambda_init;
  const double intercept = arma::mean(y);
  arma::vec yc = y - intercept;

  int n_eval = 0;
  int n_core_visits = 0;
  int n_sweeps_done = 0;
  arma::vec prev_lam = lambda;

  for (int sweep = 0; sweep < sweeps; ++sweep) {
    std::vector<arma::mat> L = build_left_interfaces(cores, basis);
    std::vector<arma::mat> R = build_right_interfaces(cores, basis);
    for (int k = 0; k < d; ++k) {
      arma::mat X = tt_design_core_d_cpp(L[k], R[k], basis[k]);
      arma::mat S = X.t() * X;
      arma::vec b = X.t() * yc;
      OptLam o = optimize_lambda_gcv(yc, X, S, penalties[k], b,
                                     lambda_min, lambda_max, tol);
      fill_cube_from_vec(cores[k], o.g);
      lambda(k) = o.lambda;
      n_eval += o.n_eval;
      ++n_core_visits;
      if (k + 1 < d) {
        L[k + 1] = contract_left_step_cpp(L[k], cores[k], basis[k]);
      }
    }
    ++n_sweeps_done;
    double dlog = 0.0;
    for (int k = 0; k < d; ++k) {
      const double a = std::log(std::max(lambda(k), 1e-12));
      const double b = std::log(std::max(prev_lam(k), 1e-12));
      dlog = std::max(dlog, std::abs(a - b));
    }
    if (dlog < tol_lambda) break;
    prev_lam = lambda;
  }

  arma::vec mu = intercept + tt_contraction_d_cpp(cubes_to_list(cores), basis_list);
  List out = List::create(
      _["cores"] = cubes_to_list(cores),
      _["intercept"] = intercept,
      _["mu"] = mu,
      _["lambda"] = lambda,
      _["d"] = d,
      _["n_sweeps"] = n_sweeps_done,
      _["n_core_visits"] = n_core_visits,
      _["n_criterion_evals"] = n_eval,
      _["family"] = "gaussian",
      _["method"] = "TT-cGCV",
      _["backend"] = "Rcpp");

  if (return_jacobian) {
    std::vector<arma::mat> L = build_left_interfaces(cores, basis);
    std::vector<arma::mat> R = build_right_interfaces(cores, basis);
    int m = 0;
    IntegerVector block_sizes(d);
    for (int k = 0; k < d; ++k) {
      block_sizes[k] = cores[k].n_rows * cores[k].n_cols * cores[k].n_slices;
      m += block_sizes[k];
    }
    arma::mat J(n, m);
    int offset = 0;
    for (int k = 0; k < d; ++k) {
      arma::mat Xk = tt_design_core_d_cpp(L[k], R[k], basis[k]);
      J.cols(offset, offset + Xk.n_cols - 1) = Xk;
      offset += Xk.n_cols;
    }
    out["jacobian"] = J;
    out["block_sizes"] = block_sizes;
  }
  return out;
}

//' GLM PIRLS + weighted TT-ALS with per-core conditional GCV (Rcpp).
//' family: "bernoulli" | "poisson" (gaussian → use tt_cgcv_fit_cpp).
// [[Rcpp::export]]
List tt_glm_pirls_cgcv_cpp(const arma::vec& y,
                           const List& basis_list,
                           const List& init_cores,
                           const List& penalties_list,
                           std::string family,
                           arma::vec lambda_init,
                           int pirls_iter = 12,
                           int als_sweeps = 4,
                           double lambda_min = 1e-2,
                           double lambda_max = 1e2,
                           double tol = 1e-3,
                           double tol_dev = 1e-6,
                           bool select_lambda = true) {
  if (family == "binomial") family = "bernoulli";
  if (family != "bernoulli" && family != "poisson") {
    stop("tt_glm_pirls_cgcv_cpp: family must be bernoulli or poisson");
  }
  std::vector<arma::mat> basis = list_to_mats(basis_list);
  std::vector<arma::cube> cores = list_to_cubes(init_cores);
  std::vector<arma::mat> penalties = list_to_mats(penalties_list);
  const int d = cores.size();
  const int n = y.n_elem;
  if (static_cast<int>(lambda_init.n_elem) != d) {
    const double lam0 = lambda_init.n_elem > 0 ? lambda_init(0) : 1.0;
    lambda_init = arma::vec(d);
    lambda_init.fill(lam0);
  }
  arma::vec lambda = lambda_init;

  // Shrink random init
  for (int k = 0; k < d; ++k) cores[k] *= 0.05;

  double intercept = glm_init_intercept(family, y);
  arma::vec f = tt_contraction_d_cpp(cubes_to_list(cores), basis_list);
  arma::vec eta = intercept + f;
  arma::vec mu = glm_invlink(family, eta);
  double dev = glm_deviance(family, y, mu);

  int n_eval = 0;
  int n_core_visits = 0;
  int n_pirls_done = 0;
  std::vector<double> hist_pirls;
  std::vector<double> hist_dev;
  std::vector<double> hist_maxeta;

  for (int it = 0; it < pirls_iter; ++it) {
    arma::vec w, z, mu_w;
    glm_working(family, y, eta, w, z, mu_w);

    for (int sw = 0; sw < als_sweeps; ++sw) {
      arma::vec zc = z - intercept;
      std::vector<arma::mat> L = build_left_interfaces(cores, basis);
      std::vector<arma::mat> R = build_right_interfaces(cores, basis);
      for (int k = 0; k < d; ++k) {
        arma::mat X = tt_design_core_d_cpp(L[k], R[k], basis[k]);
        arma::vec swt = arma::sqrt(w);
        for (arma::uword i = 0; i < swt.n_elem; ++i) {
          if (!std::isfinite(swt(i)) || swt(i) < 0.0) swt(i) = 0.0;
        }
        arma::mat Xw = X;
        Xw.each_col() %= swt;
        arma::vec yw = swt % zc;
        arma::mat S = Xw.t() * Xw;
        arma::vec b = Xw.t() * yw;
        double lam = lambda(k);
        arma::vec g;
        if (select_lambda) {
          OptLam o = optimize_lambda_gcv(yw, Xw, S, penalties[k], b,
                                         lambda_min, lambda_max, tol);
          lam = o.lambda;
          g = o.g;
          n_eval += o.n_eval;
        } else {
          g = solve_gb(S, penalties[k], b, lam);
          ++n_eval;
        }
        fill_cube_from_vec(cores[k], g);
        lambda(k) = lam;
        ++n_core_visits;
        if (k + 1 < d) {
          L[k + 1] = contract_left_step_cpp(L[k], cores[k], basis[k]);
        }
      }
      f = tt_contraction_d_cpp(cubes_to_list(cores), basis_list);
      intercept = arma::dot(w, z - f) / std::max(arma::accu(w), 1e-12);
    }

    f = tt_contraction_d_cpp(cubes_to_list(cores), basis_list);
    eta = intercept + f;
    mu = glm_invlink(family, eta);
    const double dev_new = glm_deviance(family, y, mu);
    ++n_pirls_done;
    hist_pirls.push_back(static_cast<double>(it + 1));
    hist_dev.push_back(dev_new);
    hist_maxeta.push_back(arma::max(arma::abs(eta)));
    if (!std::isfinite(dev_new)) break;
    const double rel = std::abs(dev - dev_new) / std::max(1.0, std::abs(dev));
    dev = dev_new;
    if (it > 2 && rel < tol_dev) break;
  }

  List history = List::create(
      _["pirls"] = hist_pirls,
      _["deviance"] = hist_dev,
      _["max_abs_eta"] = hist_maxeta);

  return List::create(
      _["cores"] = cubes_to_list(cores),
      _["intercept"] = intercept,
      _["mu"] = mu,
      _["eta"] = eta,
      _["lambda"] = lambda,
      _["deviance"] = dev,
      _["d"] = d,
      _["n_pirls"] = n_pirls_done,
      _["n_core_visits"] = n_core_visits,
      _["n_criterion_evals"] = n_eval,
      _["history"] = history,
      _["family"] = family,
      _["method"] = select_lambda ? "TT-cGCV-PIRLS" : "TT-PIRLS-fixed",
      _["backend"] = "Rcpp");
}
