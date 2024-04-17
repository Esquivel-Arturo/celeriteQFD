//#ifndef _CELERITE2_STANINTERFACE_HPP_DEFINED_
//#define _CELERITE2_STANINTERFACE_HPP_DEFINED_

#include <Eigen/Core>
#include <iostream>
using namespace std;


/* #include <stan/model/model_header.hpp>

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math; */



#define UNUSED(x) (void)(x)

#define CAST_BASE(TYPE, VAR) Eigen::MatrixBase<TYPE> &VAR = const_cast<Eigen::MatrixBase<TYPE> &>(VAR##_out)

#define CAST_VEC(TYPE, VAR, ROWS)                                                                                                                    \
  CAST_BASE(TYPE, VAR);                                                                                                                              \
  VAR.derived().resize(ROWS)

#define CAST_MAT(TYPE, VAR, ROWS, COLS)                                                                                                              \
  CAST_BASE(TYPE, VAR);                                                                                                                              \
  VAR.derived().resize(ROWS, COLS)

const int THE_WORKSPACE_VARIABLE_MUST_BE_ROW_MAJOR = 0;
#define ASSERT_ROW_MAJOR(TYPE) EIGEN_STATIC_ASSERT((TYPE::ColsAtCompileTime == 1) || TYPE::IsRowMajor, THE_WORKSPACE_VARIABLE_MUST_BE_ROW_MAJOR)

namespace internal {

template <bool do_update = true>
struct update_workspace {
  template <typename A, typename B>
  static void apply(Eigen::Index n, const Eigen::MatrixBase<A> &a, Eigen::MatrixBase<B> const &b_out) {
    CAST_BASE(B, b);
    b.row(n) = a;
  }
};

template <>
struct update_workspace<false> {
  template <typename A, typename B>
  static void apply(Eigen::Index n, const Eigen::MatrixBase<A> &a, Eigen::MatrixBase<B> const &b_out) {
    UNUSED(n);
    UNUSED(a);
    UNUSED(b_out);
  }
};

template <bool is_solve = false>
struct update_f {
  template <typename A, typename B, typename C, typename D>
  static void apply(const Eigen::MatrixBase<A> &a, const Eigen::MatrixBase<B> &b, const Eigen::MatrixBase<C> &c, Eigen::MatrixBase<D> const &d_out) {
    CAST_BASE(D, d);
    d.noalias() += a * b;
    UNUSED(c);
  }

  template <typename A, typename B, typename C, typename D, typename E, typename F, typename G>
  static void reverse(const Eigen::MatrixBase<A> &a, const Eigen::MatrixBase<B> &b, const Eigen::MatrixBase<C> &c, const Eigen::MatrixBase<D> &d,
                      Eigen::MatrixBase<E> const &e_out, Eigen::MatrixBase<F> const &f_out, Eigen::MatrixBase<G> const &g_out) {
    CAST_BASE(E, e);
    CAST_BASE(F, f);
    e.noalias() += b * d.transpose();
    f.noalias() += a * d;
    UNUSED(c);
    UNUSED(g_out);
  }
};

template <>
struct update_f<true> {
  template <typename A, typename B, typename C, typename D>
  static void apply(const Eigen::MatrixBase<A> &a, const Eigen::MatrixBase<B> &b, const Eigen::MatrixBase<C> &c, Eigen::MatrixBase<D> const &d_out) {
    CAST_BASE(D, d);
    d.noalias() += a * c;
    UNUSED(b);
  }

  template <typename A, typename B, typename C, typename D, typename E, typename F, typename G>
  static void reverse(const Eigen::MatrixBase<A> &a, const Eigen::MatrixBase<B> &b, const Eigen::MatrixBase<C> &c, const Eigen::MatrixBase<D> &d,
                      Eigen::MatrixBase<E> const &e_out, Eigen::MatrixBase<F> const &f_out, Eigen::MatrixBase<G> const &g_out) {
    CAST_BASE(E, e);
    CAST_BASE(G, g);
    e.noalias() += c * d.transpose();
    g.noalias() += a * d;
    UNUSED(b);
    UNUSED(f_out);
  }
};

template <bool is_solve = false>
struct update_z {
  template <typename A, typename B>
  static void apply(const Eigen::MatrixBase<A> &a, Eigen::MatrixBase<B> const &b_out) {
    CAST_BASE(B, b);
    b.noalias() += a;
  }
};

template <>
struct update_z<true> {
  template <typename A, typename B>
  static void apply(const Eigen::MatrixBase<A> &a, Eigen::MatrixBase<B> const &b_out) {
    CAST_BASE(B, b);
    b.noalias() -= a;
  }
};

template <bool is_solve, bool do_update = true, typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut,
          typename Work>
void forward(const Eigen::MatrixBase<Input> &t,                // (N,)
             const Eigen::MatrixBase<Coeffs> &c,               // (J,)
             const Eigen::MatrixBase<LowRank> &U,              // (N, J)
             const Eigen::MatrixBase<LowRank> &V,              // (N, J)
             const Eigen::MatrixBase<RightHandSide> &Y,        // (N, Nrhs)
             Eigen::MatrixBase<RightHandSideOut> const &Z_out, // (N, Nrhs)
             Eigen::MatrixBase<Work> const &F_out              // (N, J * Nrhs)
) {
  ASSERT_ROW_MAJOR(Work);

  typedef typename LowRank::Scalar Scalar;
  typedef typename Eigen::internal::plain_row_type<RightHandSide>::type RowVector;
  typedef typename Eigen::internal::plain_col_type<Coeffs>::type CoeffVector;
  typedef typename Eigen::Matrix<Scalar, LowRank::ColsAtCompileTime, RightHandSide::ColsAtCompileTime> Inner;

  Eigen::Index N = U.rows(), J = U.cols(), nrhs = Y.cols();
  CAST_BASE(RightHandSideOut, Z); // Must already be the right shape
  CAST_BASE(Work, F);
  if (do_update) {
    F.derived().resize(N, J * nrhs);
    F.row(0).setZero();
  }

  CoeffVector p(J);
  Inner Fn(J, nrhs);
  Eigen::Map<typename Eigen::internal::plain_row_type<Work>::type> ptr(Fn.data(), 1, J * nrhs);

  // This will track the previous row allowing for inplace operations
  RowVector tmp = Y.row(0);

  Fn.setZero();
  for (Eigen::Index n = 1; n < N; ++n) {
    p = exp(c.array() * (t(n - 1) - t(n)));
    update_f<is_solve>::apply(V.row(n - 1).transpose(), tmp, Z.row(n - 1), Fn);
    tmp = Y.row(n);
    update_workspace<do_update>::apply(n, ptr, F);
    Fn = p.asDiagonal() * Fn;
    update_z<is_solve>::apply(U.row(n) * Fn, Z.row(n));
  }
}

template <bool is_solve, bool do_update = true, typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut,
          typename Work>
void backward(const Eigen::MatrixBase<Input> &t,                // (N,)
              const Eigen::MatrixBase<Coeffs> &c,               // (J,)
              const Eigen::MatrixBase<LowRank> &U,              // (N, J)
              const Eigen::MatrixBase<LowRank> &V,              // (N, J)
              const Eigen::MatrixBase<RightHandSide> &Y,        // (N, Nrhs)
              Eigen::MatrixBase<RightHandSideOut> const &Z_out, // (N, Nrhs)
              Eigen::MatrixBase<Work> const &F_out              // (N, J * Nrhs)
) {
  ASSERT_ROW_MAJOR(Work);

  typedef typename LowRank::Scalar Scalar;
  typedef typename Eigen::internal::plain_row_type<RightHandSide>::type RowVector;
  typedef typename Eigen::internal::plain_col_type<Coeffs>::type CoeffVector;
  typedef typename Eigen::Matrix<Scalar, LowRank::ColsAtCompileTime, RightHandSide::ColsAtCompileTime> Inner;

  Eigen::Index N = U.rows(), J = U.cols(), nrhs = Y.cols();
  CAST_BASE(RightHandSideOut, Z); // Must already be the right shape
  CAST_BASE(Work, F);
  if (do_update) {
    F.derived().resize(N, J * nrhs);
    F.row(N - 1).setZero();
  }

  CoeffVector p(J);
  Inner Fn(J, nrhs);
  Eigen::Map<typename Eigen::internal::plain_row_type<Work>::type> ptr(Fn.data(), 1, J * nrhs);

  // This will track the previous row allowing for inplace operations
  RowVector tmp = Y.row(N - 1);

  Fn.setZero();
  for (Eigen::Index n = N - 2; n >= 0; --n) {
    p = exp(c.array() * (t(n) - t(n + 1)));
    update_f<is_solve>::apply(U.row(n + 1).transpose(), tmp, Z.row(n + 1), Fn);
    tmp = Y.row(n);
    update_workspace<do_update>::apply(n, ptr, F);
    Fn = p.asDiagonal() * Fn;
    update_z<is_solve>::apply(V.row(n) * Fn, Z.row(n));
  }
}

template <bool is_solve, typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOptional, typename Work,
          typename RightHandSideInternal, typename InputOut, typename CoeffsOut, typename LowRankOut, typename RightHandSideOut>
void forward_rev(const Eigen::MatrixBase<Input> &t,                      // (N,)
                 const Eigen::MatrixBase<Coeffs> &c,                     // (J,)
                 const Eigen::MatrixBase<LowRank> &U,                    // (N, J)
                 const Eigen::MatrixBase<LowRank> &V,                    // (N, J)
                 const Eigen::MatrixBase<RightHandSideOptional> &Y,      // (N, Nrhs)
                 const Eigen::MatrixBase<RightHandSide> &Z,              // (N, Nrhs)
                 const Eigen::MatrixBase<Work> &F,                       // (N, J * Nrhs)
                 Eigen::MatrixBase<RightHandSideInternal> const &bZ_out, // (N, Nrhs)
                 Eigen::MatrixBase<InputOut> const &bt_out,              // (N,)
                 Eigen::MatrixBase<CoeffsOut> const &bc_out,             // (J,)
                 Eigen::MatrixBase<LowRankOut> const &bU_out,            // (N, J)
                 Eigen::MatrixBase<LowRankOut> const &bV_out,            // (N, J)
                 Eigen::MatrixBase<RightHandSideOut> const &bY_out       // (N, Nrhs)  -  Must be the right shape already (and zeroed)
) {
  ASSERT_ROW_MAJOR(Work);

  typedef typename LowRank::Scalar Scalar;
  typedef typename Eigen::internal::plain_col_type<Coeffs>::type CoeffVector;
  typedef typename Eigen::Matrix<Scalar, LowRank::ColsAtCompileTime, RightHandSide::ColsAtCompileTime> Inner;

  Eigen::Index N = U.rows(), J = U.cols(), nrhs = Y.cols();
  CAST_VEC(InputOut, bt, N);
  CAST_VEC(CoeffsOut, bc, J);
  CAST_MAT(LowRankOut, bU, N, J);
  CAST_MAT(LowRankOut, bV, N, J);
  CAST_BASE(RightHandSideOut, bY);
  CAST_BASE(RightHandSideInternal, bZ);

  Scalar dt, factor;
  CoeffVector p(J), bp(J);
  Inner Fn(J, nrhs), bF(J, nrhs);
  Eigen::Map<typename Eigen::internal::plain_row_type<Work>::type> ptr(Fn.data(), 1, J * nrhs);
  bF.setZero();
  for (Eigen::Index n = N - 1; n >= 1; --n) {
    dt  = t(n - 1) - t(n);
    p   = exp(c.array() * dt);
    ptr = F.row(n);

    // Reverse: update_z<is_solve>::apply(U.row(n) * Fn, Z.row(n));
    update_z<is_solve>::apply(bZ.row(n) * (p.asDiagonal() * Fn).transpose(), bU.row(n));
    update_z<is_solve>::apply(U.row(n).transpose() * bZ.row(n), bF);

    // Reverse: Fn = P.row(n - 1).asDiagonal() * Fn;
    bp.array() = (Fn * bF.transpose()).diagonal().array() * p.array();
    bc.noalias() += dt * bp;
    factor = (c.array() * bp.array()).sum();
    bt(n) -= factor;
    bt(n - 1) += factor;
    bF = p.asDiagonal() * bF;

    // Reverse: update_f<is_solve>::apply(V.row(n - 1).transpose(), Y.row(n - 1), Z.row(n - 1), Fn);
    update_f<is_solve>::reverse(V.row(n - 1), Y.row(n - 1), Z.row(n - 1), bF, bV.row(n - 1), bY.row(n - 1), bZ.row(n - 1));
  }
}

template <bool is_solve, typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename Work, typename RightHandSideInternal,
          typename InputOut, typename CoeffsOut, typename LowRankOut, typename RightHandSideOut>
void backward_rev(const Eigen::MatrixBase<Input> &t,                      // (N,)
                  const Eigen::MatrixBase<Coeffs> &c,                     // (J,)
                  const Eigen::MatrixBase<LowRank> &U,                    // (N, J)
                  const Eigen::MatrixBase<LowRank> &V,                    // (N, J)
                  const Eigen::MatrixBase<RightHandSide> &Y,              // (N, Nrhs)
                  const Eigen::MatrixBase<RightHandSide> &Z,              // (N, Nrhs)
                  const Eigen::MatrixBase<Work> &F,                       // (N, J * Nrhs)
                  Eigen::MatrixBase<RightHandSideInternal> const &bZ_out, // (N, Nrhs)
                  Eigen::MatrixBase<InputOut> const &bt_out,              // (N,)
                  Eigen::MatrixBase<CoeffsOut> const &bc_out,             // (J,)
                  Eigen::MatrixBase<LowRankOut> const &bU_out,            // (N, J)
                  Eigen::MatrixBase<LowRankOut> const &bV_out,            // (N, J)
                  Eigen::MatrixBase<RightHandSideOut> const &bY_out       // (N, Nrhs)  -  Must be the right shape already (and zeroed)
) {
  ASSERT_ROW_MAJOR(Work);

  typedef typename LowRank::Scalar Scalar;
  typedef typename Eigen::internal::plain_col_type<Coeffs>::type CoeffVector;
  typedef typename Eigen::Matrix<Scalar, LowRank::ColsAtCompileTime, RightHandSide::ColsAtCompileTime> Inner;

  Eigen::Index N = U.rows(), J = U.cols(), nrhs = Y.cols();
  CAST_VEC(InputOut, bt, N);
  CAST_VEC(CoeffsOut, bc, J);
  CAST_MAT(LowRankOut, bU, N, J);
  CAST_MAT(LowRankOut, bV, N, J);
  CAST_BASE(RightHandSideOut, bY);
  CAST_BASE(RightHandSideInternal, bZ);

  Scalar dt, factor;
  CoeffVector p(J), bp(J);
  Inner Fn(J, nrhs), bF(J, nrhs);
  Eigen::Map<typename Eigen::internal::plain_row_type<Work>::type> ptr(Fn.data(), 1, J * nrhs);
  bF.setZero();
  for (Eigen::Index n = 0; n <= N - 2; ++n) {
    dt  = t(n) - t(n + 1);
    p   = exp(c.array() * dt);
    ptr = F.row(n);

    // Reverse: update_z<is_solve>::apply(V.row(n) * Fn, Z.row(n));
    update_z<is_solve>::apply(bZ.row(n) * (p.asDiagonal() * Fn).transpose(), bV.row(n));
    update_z<is_solve>::apply(V.row(n).transpose() * bZ.row(n), bF);

    // Reverse: Fn = P.row(n).asDiagonal() * Fn;
    bp.array() = (Fn * bF.transpose()).diagonal().array() * p.array();
    bc.noalias() += dt * bp;
    factor = (c.array() * bp.array()).sum();
    bt(n + 1) -= factor;
    bt(n) += factor;
    bF = p.asDiagonal() * bF;

    // Reverse: update_f<is_solve>::apply(U.row(n + 1).transpose(), Y.row(n + 1), Z.row(n + 1), Fn);
    update_f<is_solve>::reverse(U.row(n + 1), Y.row(n + 1), Z.row(n + 1), bF, bU.row(n + 1), bY.row(n + 1), bZ.row(n + 1));
  }
}

}

/// original name space forward
namespace forward{
/**
 * \brief Get the dense representation of a celerite matrix
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param a     (N,): The diagonal component
 * @param U     (N, J): The first low rank matrix
 * @param V     (N, J): The second low rank matrix
 * @param K_out (N, N): The dense matrix
 */
template <typename Input, typename Coeffs, typename Diag, typename LowRank, typename Dense>
void to_dense(const Eigen::MatrixBase<Input> &t,    // (N,)
              const Eigen::MatrixBase<Coeffs> &c,   // (J,)
              const Eigen::MatrixBase<Diag> &a,     // (N,)
              const Eigen::MatrixBase<LowRank> &U,  // (N, J)
              const Eigen::MatrixBase<LowRank> &V,  // (N, J)
              Eigen::MatrixBase<Dense> const &K_out // (N,)
) {
  typedef typename Eigen::internal::plain_row_type<LowRank>::type RowVector;

  Eigen::Index N = U.rows(), J = U.cols();
  CAST_MAT(Dense, K, N, N);

  RowVector p(1, J);
  for (Eigen::Index m = 0; m < N; ++m) {
    p.setConstant(1.0);
    K(m, m) = a(m);
    for (Eigen::Index n = m + 1; n < N; ++n) {
      p.array() *= exp(-c.array() * (t(n) - t(n - 1)));
      K(n, m) = (U.row(n).array() * V.row(m).array() * p.array()).sum();
      K(m, n) = K(n, m);
    }
  }
}

/**
 * \brief Compute the Cholesky factorization of the system
 *
 * This computes `d` and `W` such that:
 *
 * `K = L*diag(d)*L^T`
 *
 * where `K` is the celerite matrix and
 *
 * `L = 1 + tril(U*W^T)`
 *
 * This can be safely applied in place: `d_out` can point to `a` and `W_out` can
 * point to `V`, and the memory will be reused. In this particular case, the
 * `celerite2::core::factor_rev` function doesn't use `a` and `V`, but this
 * won't be true for all `_rev` functions.
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param a     (N,): The diagonal component
 * @param U     (N, J): The first low rank matrix
 * @param V     (N, J): The second low rank matrix
 * @param d_out (N,): The diagonal component of the Cholesky factor
 * @param W_out (N, J): The second low rank component of the Cholesky factor
 * @param S_out (N, J*J): The cached value of the S matrix at each step
 */
template <bool update_workspace = true, typename Input, typename Coeffs, typename Diag, typename LowRank, typename DiagOut, typename LowRankOut,
          typename Work>
Eigen::Index factor(const Eigen::MatrixBase<Input> &t,          // (N,)
                    const Eigen::MatrixBase<Coeffs> &c,         // (J,)
                    const Eigen::MatrixBase<Diag> &a,           // (N,)
                    const Eigen::MatrixBase<LowRank> &U,        // (N, J)
                    const Eigen::MatrixBase<LowRank> &V,        // (N, J)
                    Eigen::MatrixBase<DiagOut> const &d_out,    // (N,)
                    Eigen::MatrixBase<LowRankOut> const &W_out, // (N, J)
                    Eigen::MatrixBase<Work> const &S_out        // (N, J*J)
) {
  ASSERT_ROW_MAJOR(Work);

  typedef typename Diag::Scalar Scalar;
  typedef typename Eigen::internal::plain_row_type<LowRank>::type RowVector;
  typedef typename Eigen::internal::plain_col_type<Coeffs>::type CoeffVector;

  Eigen::Index N = U.rows(), J = U.cols();
  CAST_VEC(DiagOut, d, N);
  CAST_MAT(LowRankOut, W, N, J);
  CAST_BASE(Work, S);
  if (update_workspace) {
    S.derived().resize(N, J * J);
    S.row(0).setZero();
  }

  // This is a temporary vector used to minimize computations internally
  RowVector tmp;
  CoeffVector p;

  // This holds the accumulated value of the S matrix at each step
  Eigen::Matrix<Scalar, LowRank::ColsAtCompileTime, LowRank::ColsAtCompileTime, Eigen::ColMajor> Sn(J, J);

  // This is a flattened pointer to Sn that is used for copying the data
  Eigen::Map<typename Eigen::internal::plain_row_type<Work>::type> ptr(Sn.data(), 1, J * J);

  // First row
  Sn.setZero();
  d(0)               = a(0);
  W.row(0).noalias() = V.row(0) / d(0);

  // The rest of the rows
  for (Eigen::Index n = 1; n < N; ++n) {
    p = exp(c.array() * (t(n - 1) - t(n)));

    // Update S_n = diag(P) * (S_n-1 + d*W*W.T) * diag(P)
    Sn.noalias() += d(n - 1) * W.row(n - 1).transpose() * W.row(n - 1);
    Sn = p.asDiagonal() * Sn;

    // Save the current value of Sn to the workspace
    // Note: This is actually `diag(P) * (S + d*W*W.T)` without the final `* diag(P)`
    internal::update_workspace<update_workspace>::apply(n, ptr, S);

    // Incorporate the second diag(P) that we didn't include above for bookkeeping
    Sn *= p.asDiagonal();

    // Update d = a - U * S * U.T
    tmp  = U.row(n) * Sn;
    d(n) = a(n) - tmp * U.row(n).transpose();
    if (d(n) <= 0.0) return n;

    // Update W = (V - U * S) / d
    W.row(n).noalias() = (V.row(n) - tmp) / d(n);
  }

  return 0;
}

/**
 * \brief Apply a strictly lower matrix multiply
 *
 * This computes:
 *
 * `Z += tril(U * V^T) * Y`
 *
 * where `tril` is the strictly lower triangular function.
 *
 * Note that this will *update* the value of `Z`.
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param U     (N, J): The first low rank matrix
 * @param W     (N, J): The second low rank matrix
 * @param Y     (N, Nrhs): The matrix to be multiplied
 * @param Z_out (N, Nrhs): The matrix to be updated
 * @param F_out (N, J*Nrhs): The workspace
 */
template <bool update_workspace = true, typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut,
          typename Work>
void solve_lower(const Eigen::MatrixBase<Input> &t,                // (N,)
                 const Eigen::MatrixBase<Coeffs> &c,               // (J,)
                 const Eigen::MatrixBase<LowRank> &U,              // (N, J)
                 const Eigen::MatrixBase<LowRank> &W,              // (N, J)
                 const Eigen::MatrixBase<RightHandSide> &Y,        // (N, nrhs)
                 Eigen::MatrixBase<RightHandSideOut> const &Z_out, // (N, nrhs)
                 Eigen::MatrixBase<Work> const &F_out              // (N, J*nrhs)
) {
  ASSERT_ROW_MAJOR(Work);
  CAST_BASE(RightHandSideOut, Z);
  Z = Y;
  internal::forward<true, update_workspace>(t, c, U, W, Y, Z, F_out);
}

/**
 * \brief Compute the solution of a upper triangular linear equation
 *
 * This computes `Z` such that:
 *
 * `Y = L^T * Y`
 *
 * where
 *
 * `L = 1 + tril(U*W^T)`
 *
 * This can be safely applied in place.
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param U     (N, J): The first low rank matrix
 * @param W     (N, J): The second low rank matrix
 * @param Y     (N, Nrhs): The right hand side
 * @param Z_out (N, Nrhs): The solution of this equation
 * @param F_out (N, J*Nrhs): The workspace
 */
template <bool update_workspace = true, typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut,
          typename Work>
void solve_upper(const Eigen::MatrixBase<Input> &t,                // (N,)
                 const Eigen::MatrixBase<Coeffs> &c,               // (J,)
                 const Eigen::MatrixBase<LowRank> &U,              // (N, J)
                 const Eigen::MatrixBase<LowRank> &W,              // (N, J)
                 const Eigen::MatrixBase<RightHandSide> &Y,        // (N, nrhs)
                 Eigen::MatrixBase<RightHandSideOut> const &Z_out, // (N, nrhs)
                 Eigen::MatrixBase<Work> const &F_out              // (N, J*nrhs)
) {
  ASSERT_ROW_MAJOR(Work);
  CAST_BASE(RightHandSideOut, Z);
  Z = Y;
  internal::backward<true, update_workspace>(t, c, U, W, Y, Z, F_out);
}

/**
 * \brief Apply a strictly lower matrix multiply
 *
 * This computes:
 *
 * `Z += tril(U * V^T) * Y`
 *
 * where `tril` is the strictly lower triangular function.
 *
 * Note that this will *update* the value of `Z`.
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param U     (N, J): The first low rank matrix
 * @param V     (N, J): The second low rank matrix
 * @param Y     (N, Nrhs): The matrix to be multiplied
 * @param Z_out (N, Nrhs): The matrix to be updated
 * @param F_out (N, J*Nrhs): The workspace
 */
template <bool update_workspace = true, typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut,
          typename Work>
void matmul_lower(const Eigen::MatrixBase<Input> &t,                // (N,)
                  const Eigen::MatrixBase<Coeffs> &c,               // (J,)
                  const Eigen::MatrixBase<LowRank> &U,              // (N, J)
                  const Eigen::MatrixBase<LowRank> &V,              // (N, J)
                  const Eigen::MatrixBase<RightHandSide> &Y,        // (N, nrhs)
                  Eigen::MatrixBase<RightHandSideOut> const &Z_out, // (N, nrhs)
                  Eigen::MatrixBase<Work> const &F_out              // (N, J*nrhs)
) {
  internal::forward<false, update_workspace>(t, c, U, V, Y, Z_out, F_out);
}

/**
 * \brief Apply a strictly upper matrix multiply
 *
 * This computes:
 *
 * `Z += triu(V * U^T) * Y`
 *
 * where `triu` is the strictly lower triangular function.
 *
 * Note that this will *update* the value of `Z`.
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param U     (N, J): The first low rank matrix
 * @param V     (N, J): The second low rank matrix
 * @param Y     (N, Nrhs): The matrix to be multiplied
 * @param Z_out (N, Nrhs): The matrix to be updated
 * @param F_out (N, J*Nrhs): The workspace
 */
template <bool update_workspace = true, typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut,
          typename Work>
void matmul_upper(const Eigen::MatrixBase<Input> &t,                // (N,)
                  const Eigen::MatrixBase<Coeffs> &c,               // (J,)
                  const Eigen::MatrixBase<LowRank> &U,              // (N, J)
                  const Eigen::MatrixBase<LowRank> &V,              // (N, J)
                  const Eigen::MatrixBase<RightHandSide> &Y,        // (N, nrhs)
                  Eigen::MatrixBase<RightHandSideOut> const &Z_out, // (N, nrhs)
                  Eigen::MatrixBase<Work> const &F_out              // (N, J*nrhs)
) {
  internal::backward<false, update_workspace>(t, c, U, V, Y, Z_out, F_out);
}

/**
 * \brief The general lower-triangular dot product of a rectangular celerite system
 *
 * @param t1     (N,): The left input coordinates (must be sorted)
 * @param t2     (M,): The right input coordinates (must be sorted)
 * @param c      (J,): The transport coefficients
 * @param U      (N, J): The first low rank matrix
 * @param V      (M, J): The second low rank matrix
 * @param Y      (M, Nrhs): The matrix that will be left multiplied by the celerite model
 * @param Z_out  (N, Nrhs): The result of the operation
 * @param F_out  (M, J*Nrhs): The workspace
 */
template <bool do_update = true, typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut,
          typename Work>
void general_matmul_lower(const Eigen::MatrixBase<Input> &t1,               // (N,)
                          const Eigen::MatrixBase<Input> &t2,               // (M,)
                          const Eigen::MatrixBase<Coeffs> &c,               // (J,)
                          const Eigen::MatrixBase<LowRank> &U,              // (N, J)
                          const Eigen::MatrixBase<LowRank> &V,              // (M, J)
                          const Eigen::MatrixBase<RightHandSide> &Y,        // (M, nrhs)
                          Eigen::MatrixBase<RightHandSideOut> const &Z_out, // (N, nrhs)
                          Eigen::MatrixBase<Work> const &F_out              // (M, J*nrhs)
) {
  ASSERT_ROW_MAJOR(Work);

  typedef typename LowRank::Scalar Scalar;
  typedef typename Eigen::internal::plain_col_type<Coeffs>::type CoeffVector;
  typedef typename Eigen::Matrix<Scalar, LowRank::ColsAtCompileTime, RightHandSide::ColsAtCompileTime> Inner;

  Eigen::Index N = t1.rows(), M = t2.rows(), J = c.rows(), nrhs = Y.cols();

  CAST_BASE(RightHandSideOut, Z);
  CAST_BASE(Work, F);
  if (do_update) {
    F.derived().resize(M, J * nrhs);
    F.row(0).setZero();
  }

  CoeffVector p(J);
  Inner Fm = V.row(0).transpose() * Y.row(0);
  Eigen::Map<typename Eigen::internal::plain_row_type<Work>::type> ptr(Fm.data(), 1, J * nrhs);
  internal::update_workspace<do_update>::apply(0, ptr, F);

  Scalar tn = t2(0);
  Eigen::Index n, m = 1;
  for (n = 0; n < N; ++n)
    if (t1(n) >= tn) break;
  for (; n < N; ++n) {
    tn = t1(n);
    while (m < M && t2(m) <= tn) {
      p  = exp(c.array() * (t2(m - 1) - t2(m)));
      Fm = p.asDiagonal() * Fm;
      Fm.noalias() += V.row(m).transpose() * Y.row(m);
      internal::update_workspace<do_update>::apply(m, ptr, F);
      m++;
    }
    p = exp(c.array() * (t2(m - 1) - tn));
    Z.row(n).noalias() += U.row(n) * p.asDiagonal() * Fm;
  }
}

/**
 * \brief The general upper-triangular dot product of a rectangular celerite system
 *
 * @param t1     (N,): The left input coordinates (must be sorted)
 * @param t2     (M,): The right input coordinates (must be sorted)
 * @param c      (J,): The transport coefficients
 * @param U      (N, J): The first low rank matrix
 * @param V      (M, J): The second low rank matrix
 * @param Y      (M, Nrhs): The matrix that will be left multiplied by the celerite model
 * @param Z_out  (N, Nrhs): The result of the operation
 * @param F_out  (M, J*Nrhs): The workspace
 */
template <bool do_update = true, typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut,
          typename Work>
void general_matmul_upper(const Eigen::MatrixBase<Input> &t1,               // (N,)
                          const Eigen::MatrixBase<Input> &t2,               // (M,)
                          const Eigen::MatrixBase<Coeffs> &c,               // (J,)
                          const Eigen::MatrixBase<LowRank> &U,              // (N, J)
                          const Eigen::MatrixBase<LowRank> &V,              // (M, J)
                          const Eigen::MatrixBase<RightHandSide> &Y,        // (M, nrhs)
                          Eigen::MatrixBase<RightHandSideOut> const &Z_out, // (N, nrhs)
                          Eigen::MatrixBase<Work> const &F_out              // (M, J*nrhs)
) {
  ASSERT_ROW_MAJOR(Work);

  typedef typename LowRank::Scalar Scalar;
  typedef typename Eigen::internal::plain_col_type<Coeffs>::type CoeffVector;
  typedef typename Eigen::Matrix<Scalar, LowRank::ColsAtCompileTime, RightHandSide::ColsAtCompileTime> Inner;

  Eigen::Index N = t1.rows(), M = t2.rows(), J = c.rows(), nrhs = Y.cols();

  CAST_BASE(RightHandSideOut, Z);
  CAST_BASE(Work, F);
  if (do_update) {
    F.derived().resize(M, J * nrhs);
    F.row(0).setZero();
  }

  CoeffVector p(J);
  Inner Fm = V.row(M - 1).transpose() * Y.row(M - 1);
  Eigen::Map<typename Eigen::internal::plain_row_type<Work>::type> ptr(Fm.data(), 1, J * nrhs);

  Scalar tn = t2(M - 1);
  Eigen::Index n, m = M - 2;
  for (n = N - 1; n >= 0; --n)
    if (t1(n) < tn) break;
  for (; n >= 0; --n) {
    tn = t1(n);
    while (m >= 0 && t2(m) > tn) {
      p  = exp(c.array() * (t2(m) - t2(m + 1)));
      Fm = p.asDiagonal() * Fm;
      Fm.noalias() += V.row(m).transpose() * Y.row(m);
      internal::update_workspace<do_update>::apply(m, ptr, F);
      m--;
    }
    p = exp(c.array() * (tn - t2(m + 1)));
    Z.row(n).noalias() += U.row(n) * p.asDiagonal() * Fm;
  }
}

} 
// interface



#define MakeEmptyWork(BaseType)                                                                                                                      \
  typedef typename BaseType::Scalar Scalar;                                                                                                          \
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Empty;                                                              \
  Empty

namespace interfaces {



/**
 * \brief Compute the Cholesky factorization of the system
 *
 * This computes `d` and `W` such that:
 *
 * `diag(a) + tril(U*V^T) + triu(V*U^T) = L*diag(d)*L^T`
 *
 * where
 *
 * `L = 1 + tril(U*W^T)`
 *
 * This can be safely applied in place: `d_out` can point to `a` and `W_out` can
 * point to `V`, and the memory will be reused.
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param a     (N,): The diagonal component
 * @param U     (N, J): The first low rank matrix
 * @param V     (N, J): The second low rank matrix
 * @param d_out (N,): The diagonal component of the Cholesky factor
 * @param W_out (N, J): The second low rank component of the Cholesky factor
 */
template <typename Input, typename Coeffs, typename Diag, typename LowRank, typename DiagOut, typename LowRankOut>
Eigen::Index factor(const Eigen::MatrixBase<Input> &t,         // (N,)
                    const Eigen::MatrixBase<Coeffs> &c,        // (J,)
                    const Eigen::MatrixBase<Diag> &a,          // (N,)
                    const Eigen::MatrixBase<LowRank> &U,       // (N, J)
                    const Eigen::MatrixBase<LowRank> &V,       // (N, J)
                    Eigen::MatrixBase<DiagOut> const &d_out,   // (N,)
                    Eigen::MatrixBase<LowRankOut> const &W_out // (N, J)
) {
  MakeEmptyWork(Diag) S;
  return forward::factor<false>(t, c, a, U, V, d_out, W_out, S);
}

/**
 * \brief Compute the solution of a lower triangular linear equation
 *
 * This computes `Z` such that:
 *
 * `Y = L * Y`
 *
 * where
 *
 * `L = 1 + tril(U*W^T)`
 *
 * This can be safely applied in place.
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param U     (N, J): The first low rank matrix
 * @param W     (N, J): The second low rank matrix
 * @param Y     (N, Nrhs): The right hand side
 * @param Z_out (N, Nrhs): The solution of this equation
 */
template <typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut>
void solve_lower(const Eigen::MatrixBase<Input> &t,               // (N,)
                 const Eigen::MatrixBase<Coeffs> &c,              // (J,)
                 const Eigen::MatrixBase<LowRank> &U,             // (N, J)
                 const Eigen::MatrixBase<LowRank> &W,             // (N, J)
                 const Eigen::MatrixBase<RightHandSide> &Y,       // (N, nrhs)
                 Eigen::MatrixBase<RightHandSideOut> const &Z_out // (N, nrhs)
) {
  MakeEmptyWork(Input) F;
  forward::solve_lower<false>(t, c, U, W, Y, Z_out, F);
}

/**
 * \brief Compute the solution of a upper triangular linear equation
 *
 * This computes `Z` such that:
 *
 * `Y = L^T * Y`
 *
 * where
 *
 * `L = 1 + tril(U*W^T)`
 *
 * This can be safely applied in place.
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param U     (N, J): The first low rank matrix
 * @param W     (N, J): The second low rank matrix
 * @param Y     (N, Nrhs): The right hand side
 * @param Z_out (N, Nrhs): The solution of this equation
 */
template <typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut>
void solve_upper(const Eigen::MatrixBase<Input> &t,               // (N,)
                 const Eigen::MatrixBase<Coeffs> &c,              // (J,)
                 const Eigen::MatrixBase<LowRank> &U,             // (N, J)
                 const Eigen::MatrixBase<LowRank> &W,             // (N, J)
                 const Eigen::MatrixBase<RightHandSide> &Y,       // (N, nrhs)
                 Eigen::MatrixBase<RightHandSideOut> const &Z_out // (N, nrhs)
) {
  MakeEmptyWork(Input) F;
  forward::solve_upper<false>(t, c, U, W, Y, Z_out, F);
}

/**
 * \brief Apply a strictly lower matrix multiply
 *
 * This computes:
 *
 * `Z += tril(U * V^T) * Y`
 *
 * where `tril` is the strictly lower triangular function.
 *
 * Note that this will *update* the value of `Z`.
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param U     (N, J): The first low rank matrix
 * @param V     (N, J): The second low rank matrix
 * @param Y     (N, Nrhs): The matrix to be multiplied
 * @param Z_out (N, Nrhs): The matrix to be updated
 */
template <typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut>
void matmul_lower(const Eigen::MatrixBase<Input> &t,               // (N,)
                  const Eigen::MatrixBase<Coeffs> &c,              // (J,)
                  const Eigen::MatrixBase<LowRank> &U,             // (N, J)
                  const Eigen::MatrixBase<LowRank> &V,             // (N, J)
                  const Eigen::MatrixBase<RightHandSide> &Y,       // (N, nrhs)
                  Eigen::MatrixBase<RightHandSideOut> const &Z_out // (N, nrhs)
) {
  MakeEmptyWork(Input) F;
  forward::matmul_lower<false>(t, c, U, V, Y, Z_out, F);
}

/**
 * \brief Apply a strictly upper matrix multiply
 *
 * This computes:
 *
 * `Z += triu(V * U^T) * Y`
 *
 * where `triu` is the strictly lower triangular function.
 *
 * Note that this will *update* the value of `Z`.
 *
 * @param t     (N,): The input coordinates (must be sorted)
 * @param c     (J,): The transport coefficients
 * @param U     (N, J): The first low rank matrix
 * @param V     (N, J): The second low rank matrix
 * @param Y     (N, Nrhs): The matrix to be multiplied
 * @param Z_out (N, Nrhs): The matrix to be updated
 */
template <typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut>
void matmul_upper(const Eigen::MatrixBase<Input> &t,               // (N,)
                  const Eigen::MatrixBase<Coeffs> &c,              // (J,)
                  const Eigen::MatrixBase<LowRank> &U,             // (N, J)
                  const Eigen::MatrixBase<LowRank> &V,             // (N, J)
                  const Eigen::MatrixBase<RightHandSide> &Y,       // (N, nrhs)
                  Eigen::MatrixBase<RightHandSideOut> const &Z_out // (N, nrhs)
) {
  MakeEmptyWork(Input) F;
  forward::matmul_upper<false>(t, c, U, V, Y, Z_out, F);
}

/**
 * \brief The general lower-triangular dot product of a rectangular celerite system
 *
 * @param t1     (N,): The left input coordinates (must be sorted)
 * @param t2     (M,): The right input coordinates (must be sorted)
 * @param c      (J,): The transport coefficients
 * @param U      (N, J): The first low rank matrix
 * @param V      (M, J): The second low rank matrix
 * @param Y      (M, Nrhs): The matrix that will be multiplied
 * @param Z_out  (N, Nrhs): The result of the operation
 */
template <typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut>
void general_matmul_lower(const Eigen::MatrixBase<Input> &t1,              // (N,)
                          const Eigen::MatrixBase<Input> &t2,              // (M,)
                          const Eigen::MatrixBase<Coeffs> &c,              // (J,)
                          const Eigen::MatrixBase<LowRank> &U,             // (N, J)
                          const Eigen::MatrixBase<LowRank> &V,             // (M, J)
                          const Eigen::MatrixBase<RightHandSide> &Y,       // (M, nrhs)
                          Eigen::MatrixBase<RightHandSideOut> const &Z_out // (N, nrhs)
) {
  MakeEmptyWork(Input) F;
  forward::general_matmul_lower<false>(t1, t2, c, U, V, Y, Z_out, F);
}

/**
 * \brief The general upper-triangular dot product of a rectangular celerite system
 *
 * @param t1     (N,): The left input coordinates (must be sorted)
 * @param t2     (M,): The right input coordinates (must be sorted)
 * @param c      (J,): The transport coefficients
 * @param U      (N, J): The first low rank matrix
 * @param V      (M, J): The second low rank matrix
 * @param Y      (M, Nrhs): The matrix that will be multiplied
 * @param Z_out  (N, Nrhs): The result of the operation
 */
template <typename Input, typename Coeffs, typename LowRank, typename RightHandSide, typename RightHandSideOut>
void general_matmul_upper(const Eigen::MatrixBase<Input> &t1,              // (N,)
                          const Eigen::MatrixBase<Input> &t2,              // (M,)
                          const Eigen::MatrixBase<Coeffs> &c,              // (J,)
                          const Eigen::MatrixBase<LowRank> &U,             // (N, J)
                          const Eigen::MatrixBase<LowRank> &V,             // (M, J)
                          const Eigen::MatrixBase<RightHandSide> &Y,       // (M, nrhs)
                          Eigen::MatrixBase<RightHandSideOut> const &Z_out // (N, nrhs)
) {
  MakeEmptyWork(Input) F;
  forward::general_matmul_upper<false>(t1, t2, c, U, V, Y, Z_out, F);
}

}

namespace terms{

#ifndef CELERITE_MAX_WIDTH
#define CELERITE_MAX_WIDTH 32
#endif

struct dimension_mismatch : public std::exception {
  const char *what() const throw() { return "dimension mismatch"; }
};

template <int J1, int J2>
struct sum_width {
  constexpr static int value = (J1 == Eigen::Dynamic || J2 == Eigen::Dynamic) ? Eigen::Dynamic : (J1 + J2);
};

/**
 * The abstract base class from which terms should inherit
 */
template <typename T, int J_ = Eigen::Dynamic>
class Term {
  protected:
  constexpr static int Width = ((0 < J_) && (J_ <= CELERITE_MAX_WIDTH)) ? J_ : Eigen::Dynamic;
  static constexpr int Order = (Width != 1) ? Eigen::RowMajor : Eigen::ColMajor;

  public:
  /**
   * \typedef Scalar
   * The underlying scalar type of this `Term` (should probably always be `double`)
   */
  typedef T Scalar;

  /**
   * \typedef Vector
   * An `Eigen` vector with data type `Scalar`
   */
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

  /**
   * \typedef LowRank
   * The `Eigen` type for the low-rank matrices used internally
   */
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Width, Order> LowRank;

  /**
   * \typedef CoeffVector
   * The `Eigen` type for a fixed width vector of coefficients
   */
  typedef Eigen::Matrix<Scalar, Width, 1> CoeffVector;

  /**
   * \typedef Coeffs
   * A tuple of vectors giving the coefficients for the celerite model
   */
  typedef std::tuple<Vector, Vector, Vector, Vector, Vector, Vector> Coeffs;

  /**
   * \typedef Matrices
   * A tuple of matrices representing this celerite process
   */
  typedef std::tuple<CoeffVector, Vector, LowRank, LowRank> Matrices;

  Term(){};

  int get_width() const { return Width; }

  /**
   * Set the coefficients of the term
   *
   * @param ar     (J_real,): The real amplitudes.
   * @param cr     (J_real,): The real exponential.
   * @param ac     (J_comp,): The complex even amplitude.
   * @param bc     (J_comp,): The complex odd amplitude.
   * @param cc     (J_comp,): The complex exponential.
   * @param dc     (J_comp,): The complex frequency.
   */
  void set_coefficients(const Vector &ar, const Vector &cr, const Vector &ac, const Vector &bc, const Vector &cc, const Vector &dc) {
    Eigen::Index nr = ar.rows(), nc = ac.rows();

    ar_.resize(nr);
    cr_.resize(nr);
    ac_.resize(nc);
    bc_.resize(nc);
    cc_.resize(nc);
    dc_.resize(nc);

    ar_ << ar;
    cr_ << cr;
    ac_ << ac;
    bc_ << bc;
    cc_ << cc;
    dc_ << dc;
  }

  /**
   * Get the coefficients of the term as a tuple
   */
  Coeffs get_coefficients() const { return std::make_tuple(ar_, cr_, ac_, bc_, cc_, dc_); }

  /**
   * Get the matrices required to represent the celerite process
   *
   * @param x    (N,): The independent coordinates of the data.
   * @param diag (N,): The diagonal variance of the process.
   */
  Matrices get_celerite_matrices(const Vector &x, const Vector &diag) const {
    Eigen::Index N = x.rows();
    if (diag.rows() != N) throw dimension_mismatch();

    Eigen::Index nr = ar_.rows();
    Eigen::Index nc = ac_.rows();
    Eigen::Index J  = nr + 2 * nc;
    if (Width != Eigen::Dynamic && Width != J) throw dimension_mismatch();

    CoeffVector c(J);
    Vector a = diag.array() + (ar_.sum() + ac_.sum());
    LowRank U(N, J), V(N, J);

    c << cr_, cc_, cc_;

    U.block(0, 0, N, nr).rowwise() = ar_.transpose();
    V.block(0, 0, N, nr).setConstant(Scalar(1));

    auto arg                   = (x * dc_.transpose()).array().eval();
    auto ca                    = cos(arg).eval();
    auto sa                    = sin(arg).eval();
    U.block(0, nr, N, nc)      = ca.array().rowwise() * ac_.transpose().array() + sa.array().rowwise() * bc_.transpose().array();
    U.block(0, nr + nc, N, nc) = sa.array().rowwise() * ac_.transpose().array() - ca.array().rowwise() * bc_.transpose().array();
    V.block(0, nr, N, nc)      = ca;
    V.block(0, nr + nc, N, nc) = sa;

    return std::make_tuple(c, a, U, V);
  }

  /**
   * Adding two terms builds a new term where the coefficients have been concatenated
   *
   * @param other (Term): The term to add to this one.
   */
  template <typename Other>
  Term<typename std::common_type<Scalar, typename Other::Scalar>::type, sum_width<Width, Other::Width>::value> operator+(const Other &other) const {
    typedef typename std::common_type<Scalar, typename Other::Scalar>::type NewScalar;

    auto coeffs = other.get_coefficients();

    Eigen::Index nr = ar_.rows() + std::get<0>(coeffs).rows();
    Eigen::Index nc = ac_.rows() + std::get<2>(coeffs).rows();

    Eigen::Matrix<NewScalar, Eigen::Dynamic, 1> ar(nr), cr(nr), ac(nc), bc(nc), cc(nc), dc(nc);

    ar << ar_, std::get<0>(coeffs);
    cr << cr_, std::get<1>(coeffs);
    ac << ac_, std::get<2>(coeffs);
    bc << ac_, std::get<3>(coeffs);
    cc << ac_, std::get<4>(coeffs);
    dc << ac_, std::get<5>(coeffs);

    Term<NewScalar, sum_width<Width, Other::Width>::value> new_term;
    new_term.set_coefficients(ar, cr, ac, bc, cc, dc);

    return new_term;
  }

  private:
  Vector ar_, cr_, ac_, bc_, cc_, dc_;
};

/**
 * \class RealTerm
 * The simplest celerite model
 *
 * @param a: The amplitude of the term.
 * @param c: The exponent of the term.
 */
template <typename T>
class RealTerm : public Term<T, 1> {
  public:
  /**
   * \typedef Scalar
   * The underlying scalar type of this `Term` (should probably always be `double`)
   */
  typedef T Scalar;
  constexpr static int Width = 1;
  using typename Term<Scalar, 1>::Vector;
  using typename Term<Scalar, 1>::LowRank;
  RealTerm(const Scalar &a, const Scalar &c) {
    Vector ar(1), cr(1), ac, bc, cc, dc;
    ar << a;
    cr << c;
    this->set_coefficients(ar, cr, ac, bc, cc, dc);
  };
};

/**
 * \class ComplexTerm
 * A general celerite model
 *
 * @param a: The real part of the amplitude.
 * @param b: The complex part of the amplitude.
 * @param c: The real part of the exponent.
 * @param d: The complex part of the exponent.
 */
template <typename T>
class ComplexTerm : public Term<T, 2> {
  public:
  /**
   * \typedef Scalar
   * The underlying scalar type of this `Term` (should probably always be `double`)
   */
  typedef T Scalar;
  constexpr static int Width = 2;
  using typename Term<Scalar, 2>::Vector;
  using typename Term<Scalar, 2>::LowRank;
  ComplexTerm(const Scalar &a, const Scalar &b, const Scalar &c, const Scalar &d) {
    Vector ar, cr, ac(1), bc(1), cc(1), dc(1);
    ac << a;
    bc << b;
    cc << c;
    dc << d;
    this->set_coefficients(ar, cr, ac, bc, cc, dc);
  };
};

/**
 * \class SHOTerm
 * A term representing a stochastically-driven, damped harmonic oscillator
 *
 * @param S0:  The power at `omega = 0`.
 * @param w0:  The undamped angular frequency.
 * @param Q:   The quality factor.
 * @param eps: A regularization parameter used for numerical stability.
 */
template <typename T>
class SHOTerm : public Term<T, 2> {
  public:
  /**
   * \typedef Scalar
   * The underlying scalar type of this `Term` (should probably always be `double`)
   */
  typedef T Scalar;
  constexpr static int Width = 2;
  using typename Term<Scalar, 2>::Vector;
  using typename Term<Scalar, 2>::LowRank;
  SHOTerm(const Scalar &S0, const Scalar &w0, const Scalar &Q, const Scalar &eps = 1e-5) {
    Vector ar, cr, ac, bc, cc, dc;
    if (Q < 0.5) {
      ar.resize(2);
      cr.resize(2);
      //auto f = stan::math::sqrt(stan::math::max(1.0 - 4.0 * Q * Q, eps));
      // need to make stan happy by using stan's math rather than std::sqrt
      auto f = stan::math::sqrt(1.0 - 4.0 * Q * Q >= eps ? 1.0 - 4.0 * Q * Q : eps);
      auto a = 0.5 * S0 * w0 * Q;
      auto c = 0.5 * w0 / Q;
      ar(0)  = a * (1 + 1 / f);
      ar(1)  = a * (1 - 1 / f);
      cr(0)  = c * (1 - f);
      cr(1)  = c * (1 + f);
    } else {
      ac.resize(1);
      bc.resize(1);
      cc.resize(1);
      dc.resize(1);
      auto f = stan::math::sqrt((4.0 * Q * Q - 1 >= eps ? 4.0 * Q * Q - 1 : eps));
      //auto f = stan::math::sqrt(4.0 * Q * Q - 1);
      auto a = S0 * w0 * Q;
      auto c = 0.5 * w0 / Q;
      ac(0)  = a;
      bc(0)  = a / f;
      cc(0)  = c;
      dc(0)  = c * f;
    }
    this->set_coefficients(ar, cr, ac, bc, cc, dc);
  };
};

}

// Calculate the log Likelihood function of SHO kernel, see SHOTerm for more details on parameters
template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__, typename T6__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__>::type>::type
logLikSHO(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
              const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& y,
              const T2__& S0,
              const T3__& w0,
              const T4__& Q,
              const T5__& eps,
              const Eigen::Matrix<T6__, Eigen::Dynamic, 1>& diag, std::ostream* pstream__){
  typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__>::type>::type local_scalar_t__;
  typedef local_scalar_t__ fun_return_scalar_t__;

  // local copy of data, for whatever reason it is necessary for the solver
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> tloc;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> yloc;
  tloc = t;
  yloc = y;

  local_scalar_t__ epsloc;
  epsloc = eps;
    // declear a SHO term
  terms::SHOTerm<local_scalar_t__> curr_SHOTerm(S0, w0, Q, epsloc);
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> c;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> a;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> U;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> V;

  std::tie (c, a, U, V) = curr_SHOTerm.get_celerite_matrices(tloc,diag);// get those magic matrixes
  Eigen::Index flag;
  flag = interfaces::factor(t, c, a, U, V, a, V); // to reuse memory

  // // now a is d (the diagonal of the Chol decomposition), and V is the W
  // // s.t. `K = L*diag(d)*L^T`, `L = 1 + tril(U*W^T)`
  // // lower solver, see python lib celerite2/celerite2.py#L272
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> Z;
  interfaces::solve_lower(tloc, c, U, V, yloc, Z);
  Eigen::Array<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> temp;
  temp = a.array();
  temp = temp.inverse() * Z.array();
  temp = temp * Z.array();
  return(-0.5 * sum(log(a)) - 0.5 * sum(temp));
}



// Calculate the log Likelihood function of Rotation kernel, see RoatationTerm for more details on parameters
template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__, typename T6__, typename T7__, typename T8__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__, T7__, typename boost::math::tools::promote_args<T8__>::type>::type>::type
logLikRotation(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
                   const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& y,
                   const T2__& sigma,
                   const T3__& period,
                   const T4__& Q0,
                   const T5__& dQ,
                   const T6__& f,
                   const T7__& eps,
                   const Eigen::Matrix<T8__, Eigen::Dynamic, 1>& diag, 
                   std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__, T7__, typename boost::math::tools::promote_args<T8__>::type>::type>::type  local_scalar_t__;
    //typedef double local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;

  // local copy of data, for whatever reason it is necessary for the solver
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> tloc;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> yloc;
  tloc = t;
  yloc = y;

  // two component of SHO:
  local_scalar_t__ S1, w1, Q1, S2, w2, Q2, amp, epsloc;
  amp =  sigma * sigma/ (1 + f);
  Q1 = 0.5 + Q0 + dQ;
  w1 = 4 * 3.1415926 * Q1 / (period * stan::math::sqrt(4 * Q1 * Q1 - 1));
  S1 = amp / (w1 * Q1);
  Q2 = 0.5 + Q0;
  w2 = 8 * 3.1415926 * Q2 / (period * stan::math::sqrt(4 * Q2 * Q2 - 1));
  S2 = f * amp / (w2 * Q2);
  epsloc = eps;
   
  //  // get the two terms and add them up
  terms::SHOTerm<local_scalar_t__> curr_SHOTerm1(S1, w1, Q1, epsloc);
  //terms::SHOTerm<local_scalar_t__> curr_SHOTerm2(S2, w2, Q2, epsloc);
  terms::SHOTerm<local_scalar_t__> Rotation(S2, w2, Q2, epsloc);
  Rotation + curr_SHOTerm1;


  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> c;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> a;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> U;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> V;

  std::tie (c, a, U, V) = Rotation.get_celerite_matrices(t,diag);// get those magic matrixes

  Eigen::Index flag;

  flag = interfaces::factor(t, c, a, U, V, a, V); // to reuse memory
  
  // // now a is d (the diagonal of the Chol decomposition), and V is the W
  // // s.t. `K = L*diag(d)*L^T`, `L = 1 + tril(U*W^T)`
  // // lower solver, see python lib celerite2/celerite2.py#L272
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> Z;
  interfaces::solve_lower(tloc, c, U, V, yloc, Z);
  Eigen::Array<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> temp;
  temp = a.array();
  temp = temp.inverse() * Z.array();
  
  temp = temp * Z.array();
  return(-0.5 * sum(log(a)) - 0.5 * sum(temp));

  
}

// Apply cholesky decomposition of the covariance matrix of Rotation kernel to a vector y
template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__, typename T6__, typename T7__, typename T8__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__, T7__, typename boost::math::tools::promote_args<T8__>::type>::type>::type, Eigen::Dynamic, 1>
dotCholRotation(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
                    const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& y,
                    const T2__& sigma,
                    const T3__& period,
                    const T4__& Q0,
                    const T5__& dQ,
                    const T6__& f,
                    const T7__& eps,
                    const Eigen::Matrix<T8__, Eigen::Dynamic, 1>& diag, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__, T7__, typename boost::math::tools::promote_args<T8__>::type>::type>::type  local_scalar_t__;
    //typedef double local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;

  // local copy of data, for whatever reason it is necessary for the solver
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> tloc;
  //Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> yloc;
  tloc = t;
  //yloc = y;

  // two component of SHO:
  local_scalar_t__ S1, w1, Q1, S2, w2, Q2, amp, epsloc;
  amp =  sigma * sigma/ (1 + f);
  Q1 = 0.5 + Q0 + dQ;
  w1 = 4 * 3.1415926 * Q1 / (period * stan::math::sqrt(4 * Q1 * Q1 - 1));
  S1 = amp / (w1 * Q1);
  Q2 = 0.5 + Q0;
  w2 = 8 * 3.1415926 * Q2 / (period * stan::math::sqrt(4 * Q2 * Q2 - 1));
  S2 = f * amp / (w2 * Q2);
  epsloc = eps;
   
  //  // get the two terms and add them up
  terms::SHOTerm<local_scalar_t__> curr_SHOTerm1(S1, w1, Q1, epsloc);
  //terms::SHOTerm<local_scalar_t__> curr_SHOTerm2(S2, w2, Q2, epsloc);
  terms::SHOTerm<local_scalar_t__> Rotation(S2, w2, Q2, epsloc);
  Rotation + curr_SHOTerm1;


  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> c;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> a;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> U;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> V;

  std::tie (c, a, U, V) = Rotation.get_celerite_matrices(t,diag);// get those magic matrixes

  Eigen::Index flag;

  flag = interfaces::factor(t, c, a, U, V, a, V); // to reuse memory
  
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> Z;
  Z = y;
  Eigen::Array<local_scalar_t__, Eigen::Dynamic, 1> temp;
  temp = a.array();
  temp = sqrt(temp);
  temp = temp * y.array();
  Z = temp.matrix();
  interfaces::matmul_lower(tloc, c, U, V, Z, Z);
  return(Z);
  //return(0.5);
  
}


// Apply cholesky decomposition of the covariance matrix of SHO kernel to a vector y

template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__, typename T6__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__>::type>::type,Eigen::Dynamic, 1>
dotCholSHO(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
              const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& y,
              const T2__& S0,
              const T3__& w0,
              const T4__& Q,
              const T5__& eps,
              const Eigen::Matrix<T6__, Eigen::Dynamic, 1>& diag, std::ostream* pstream__){
  typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__>::type>::type local_scalar_t__;
  typedef local_scalar_t__ fun_return_scalar_t__;

  // local copy of data, for whatever reason it is necessary for the solver
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> tloc;
  //Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> yloc;
  tloc = t;
  //yloc = y;

  local_scalar_t__ epsloc;
  epsloc = eps;
    // declear a SHO term
  terms::SHOTerm<local_scalar_t__> curr_SHOTerm(S0, w0, Q, epsloc);
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> c;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> a;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> U;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> V;

  std::tie (c, a, U, V) = curr_SHOTerm.get_celerite_matrices(tloc,diag);// get those magic matrixes
  Eigen::Index flag;
  flag = interfaces::factor(tloc, c, a, U, V, a, V); // to reuse memory
  
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> Z;
  Z = y;
  Eigen::Array<local_scalar_t__, Eigen::Dynamic, 1> temp;
  temp = a.array();
  temp = sqrt(temp);
  temp = temp * y.array();
  Z = temp.matrix();
  interfaces::matmul_lower(tloc, c, U, V, Z, Z);
  return(Z);
}


template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__, typename T5__, typename T6__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__>::type>::type,Eigen::Dynamic, 1>
dotCholQuasiPeriod(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& t,
              const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& y,
              const T2__& B,
              const T3__& L,
              const T4__& P,
              const T5__& C,
              const Eigen::Matrix<T6__, Eigen::Dynamic, 1>& diag, std::ostream* pstream__){
  typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__, T5__, T6__>::type>::type local_scalar_t__;
  typedef local_scalar_t__ fun_return_scalar_t__;

  // local copy of data, for whatever reason it is necessary for the solver
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> tloc;
  //Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> yloc;
  tloc = t;
  //yloc = y;

  local_scalar_t__ a1,b1,c1,d1,a2,b2,c2,d2;

  a1 = B /(2 + C);
  b1 = 0;
  c1 = 1/L;
  d1 = 2.*3.1415926535/P;

  a2 = a1*(1+C);
  b2 = 0;
  c2 = c1;
  d2 = 0;

  terms::ComplexTerm<local_scalar_t__> quasiperiodicterm(a1,b1,c1,d1);
  terms::ComplexTerm<local_scalar_t__> term2(a2,b2,c2,d2);
  quasiperiodicterm + term2;

  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> c;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> a;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> U;
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> V;

  std::tie (c, a, U, V) = quasiperiodicterm.get_celerite_matrices(tloc,diag);// get those magic matrixes
  Eigen::Index flag;
  flag = interfaces::factor(tloc, c, a, U, V, a, V); // to reuse memory
  
  Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> Z;
  Z = y;
  Eigen::Array<local_scalar_t__, Eigen::Dynamic, 1> temp;
  temp = a.array();
  temp = sqrt(temp);
  temp = temp * y.array();
  Z = temp.matrix();
  interfaces::matmul_lower(tloc, c, U, V, Z, Z);
  return(Z);
}





//#endif