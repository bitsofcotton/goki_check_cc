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
template <typename T> const SimpleVector<T> invariantP1(const SimpleVector<T>& in, const int& varlen, const int& skip = 0) {
#if defined(_FLOAT_BITS_)
  const auto epsilon(T(1) >> int64_t(mybits - 2));
#else
  const auto epsilon(std::numeric_limits<T>::epsilon());
#endif
  const auto threshold_inner(sqrt(epsilon));
  SimpleMatrix<T> A((in.size() - varlen) * 2 + 1, varlen);
  SimpleVector<T> fvec(A.cols());
  SimpleVector<T> one(A.rows());
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < A.rows() / 2; i ++) {
    for(int j = 0; j < A.cols(); j ++)
      A(i, j) = in[i + j];
    const auto norm(sqrt(A.row(i).dot(A.row(i))));
    if(norm != T(0))
      A.row(i) /= norm;
    A.row(i + A.rows() / 2) = - A.row(i);
  }
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < A.rows(); i ++)
    one[i] = T(1);
  one /= sqrt(one.dot(one));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < fvec.size() - 1; i ++)
    A(A.rows() - 1, i) = fvec[i] = T(0);
  A(A.rows() - 1, A.cols() - 1) = T(1);
  fvec[fvec.size() - 1] = T(0);
  SimpleMatrix<T> Pt(A.cols(), A.rows());
  for(int i = 0; i < Pt.rows(); i ++)
    for(int j = 0; j < Pt.cols(); j ++)
      Pt(i, j) = T(0);
  for(int i = 0; i < A.cols(); i ++) {
    const auto Atrowi(A.col(i));
    const auto work(Atrowi - Pt.projectionPt(Atrowi));
    Pt.row(i) = work / sqrt(work.dot(work));
  }
  const auto R(Pt * A);
  SimpleVector<T> on;
  SimpleVector<T> orth;
  for(int n_fixed = 0 ; n_fixed < Pt.rows() - 1; n_fixed ++) {
    on = Pt.projectionPt(- one);
    std::vector<pair<T, int> > onM;
    onM.reserve(on.size());
    for(int j = 0; j < on.size(); j ++)
      onM.emplace_back(std::make_pair(on[j], j));
    std::sort(onM.begin(), onM.end());
    const auto& fidx(onM[onM.size() - 1 - skip].second);
    assert(isfinite(on[fidx]));
    if(on[fidx] <= threshold_inner)
      break;
    orth = Pt.col(fidx);
    const auto norm2orth(orth.dot(orth));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / norm2orth);
  }
  fvec = R.solve(Pt * on);
  const auto nfv(fvec.dot(fvec));
  return (isfinite(nfv) ? (nfv == T(0) ? fvec : fvec / sqrt(nfv))
    : fvec * T(0));
}


template <typename T> class P1 {
public:
  typedef SimpleVector<T> Vec;
  inline P1();
  inline P1(const int& stat, const int& var);
  inline ~P1();
  inline T next(const T& in, const int& skip = 0);
private:
  Vec buf;
  int varlen;
  int t;
};

template <typename T> inline P1<T>::P1() {
  varlen = t = 0;
}

template <typename T> inline P1<T>::P1(const int& stat, const int& var) {
  buf.resize(stat + (varlen = var) - 1);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
  t = 0;
}

template <typename T> inline P1<T>::~P1() {
  ;
}

template <typename T> inline T P1<T>::next(const T& in, const int& skip) {
  for(int i = 0; i < buf.size() - 1; i ++)
    buf[i] = buf[i + 1];
  buf[buf.size() - 1] = in;
  if(t ++ < buf.size()) return T(0);
  const auto fvec(invariantP1<T>(buf, varlen, skip));
  SimpleVector<T> p(fvec.size() - 1);
  for(int i = 0; i < p.size(); i ++)
    p[i] = fvec[i] / fvec[fvec.size() - 1];
  p /= sqrt(p.dot(p));
  T res(0);
  for(int i = 0; i < p.size(); i ++)
    res += buf[i - p.size() + buf.size()] * p[i];
  return isfinite(res) ? res : T(0);
}

#define _P1_
#endif

