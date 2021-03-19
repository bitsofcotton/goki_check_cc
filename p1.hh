/*
BSD 3-Clause License

Copyright (c) 2020-2021, kazunobu watatsu
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#if !defined(_P1_)

// Get invariant structure that
// \[0, &alpha;\[ register computer with deterministic calculation
// that whole in/out is shown in the varlen.
// Please refer bitsofcotton/randtools .
// [[A I], [I A^-1]] [x y], if rank A = dim x = dim y,
// rank [[A I], [I A^-1]] == rank A.
// so with this, U [[B], [O]] [x y] == 0.
// this returns one of u_k that concerns 0.
// skip attemps maximum skip status number if original data has a noised ones.
template <typename T> SimpleVector<T> invariantP1(const SimpleVector<T>& in, const int& varlen, const int& skip = 0, const bool& computer = false) {
  assert(varlen < in.size());
#if defined(_FLOAT_BITS_)
  static const auto epsilon(T(1) >> int64_t(mybits - 2));
#else
  static const auto epsilon(std::numeric_limits<T>::epsilon());
#endif
  static const auto threshold_inner(sqrt(epsilon));
  SimpleMatrix<T> A((in.size() - varlen) * 2 + 1, varlen + 1);
  SimpleVector<T> fvec(A.cols());
  SimpleVector<T> one(A.rows());
  const auto nin(sqrt(in.dot(in)));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < fvec.size(); i ++)
    A(A.rows() - 1, i) = fvec[i] = T(0);
  for(int i = 0; i < in.size() - varlen; i ++) {
    for(int j = 0; j < varlen; j ++)
      A(i, j) = in[i + j] / nin;
    A(i, varlen) = T(1) / sqrt(T(A.rows() * A.cols()));
    if(computer) {
      auto mul(A(i, 0));
      for(int j = 1; j < A.cols(); j ++)
        mul *= A(i, j);
      if(mul == T(0)) return fvec;
      assert(mul != T(0));
      A.row(i) /= pow(abs(mul), T(1) / T(A.cols()));
      for(int j = 0; j < A.cols(); j ++)
        A(i, j) = T(1) / A(i, j);
    }
    A.row(i + A.rows() / 2) = - A.row(i);
  }
  auto denom(A.row(0).dot(A.row(0)));
  for(int i = 1; i < A.rows() / 2; i ++)
    denom = max(denom, A.row(i).dot(A.row(i)));
  A /= sqrt(denom);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < A.rows(); i ++)
    one[i] = T(1);
  one /= sqrt(one.dot(one));
  A(A.rows() - 1, varlen - 1) = T(1);
  SimpleMatrix<T> Pt(A.cols(), A.rows());
  for(int i = 0; i < Pt.rows(); i ++)
    for(int j = 0; j < Pt.cols(); j ++)
      Pt(i, j) = T(0);
  for(int i = 0; i < A.cols(); i ++) {
    const auto Atrowi(A.col(i));
    const auto work(Atrowi - Pt.projectionPt(Atrowi));
    Pt.row(i) = work / sqrt(work.dot(work));
    const auto npt(Pt.row(i).dot(Pt.row(i)));
    if(! (isfinite(npt) && ! isnan(npt)))
      return fvec;
  }
  const auto R(Pt * A);
  const auto mone(Pt * (- one));
  SimpleVector<T> on;
  SimpleVector<T> orth;
  for(int n_fixed = 0 ; n_fixed < Pt.rows() - 1; n_fixed ++) {
    on = Pt.projectionPt(- one);
    std::vector<pair<T, int> > onM;
    onM.reserve(on.size());
    for(int j = 0; j < on.size(); j ++)
      onM.emplace_back(std::make_pair(on[j], j));
    std::sort(onM.begin(), onM.end());
    auto fidx(- 1);
    for(int i = onM.size() - 1 - skip; i < onM.size(); i ++)
      if(threshold_inner < onM[i].first) {
        // XXX: select one of these.
        // fidx = onM[onM.size() - 1].second;
        // fidx = onM[i].second;
        fidx = onM[(onM.size() - 1 + i) / 2].second;
        break;
      }
    if(fidx < 0)
      break;
    orth = Pt.col(fidx);
    const auto norm2orth(orth.dot(orth));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / norm2orth);
  }
  fvec = R.solve(mone - Pt * on);
  const auto nfv(fvec.dot(fvec));
  return (isfinite(nfv) ? (nfv == T(0) ? fvec : fvec / sqrt(nfv))
    : fvec * T(0));
}


template <typename T> class P1I {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P1I();
  inline P1I(const int& stat, const int& var);
  inline ~P1I();
  inline T next(const T& in, const int& skip = 0, const T& computer = T(0));
  std::vector<Vec> invariant;
private:
  inline const T& sgn(const T& x) const;
  Vec buf;
  int varlen;
  int t;
};

template <typename T> inline P1I<T>::P1I() {
  varlen = t = 0;
}

template <typename T> inline P1I<T>::P1I(const int& stat, const int& var) {
  buf.resize(stat + (varlen = var) * 2 - 1);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
  t = 0;
}

template <typename T> inline P1I<T>::~P1I() {
  ;
}

template <typename T> inline const T& P1I<T>::sgn(const T& x) const {
  static const T zero(0);
  static const T one(1);
  static const T mone(- one);
  return x != zero ? (x < zero ? mone : one) : zero;
}

template <typename T> T P1I<T>::next(const T& in, const int& skip, const T& computer) {
  for(int i = 0; i < buf.size() - 1; i ++)
    buf[i]  = buf[i + 1];
  buf[buf.size() - 1] = in + computer;
  if(t ++ < buf.size()) return computer;
  // N.B. to compete with cheating, we calculate long term same invariant each.
  const auto invariant0(invariantP1<T>(buf, varlen, abs(skip), computer != T(0)));
  if(invariant0.dot(invariant0) == T(0)) return computer;
  if(int(invariant.size()) < skip) {
    invariant.emplace_back(invariant0);
    return computer;
  } else {
    if(! invariant.size()) invariant.emplace_back(invariant0);
    for(int i = 0; i < invariant.size() - 1; i ++)
      std::swap(invariant[i], invariant[i + 1]);
    invariant[invariant.size() - 1] = invariant0;
  }
  auto avg(invariant[0]);
  for(int i = 1; i < invariant.size(); i ++)
    avg += invariant[i];
  auto res(avg[varlen] * sqrt(T((buf.size() - varlen) * 2 + 1) * T(varlen + 1)));
  for(int i = 0; i < varlen - 1; i ++)
    res += avg[i] / buf[i - varlen + buf.size() + 1];
  return res = avg[varlen - 1] / res;
}

#define _P1_
#endif

