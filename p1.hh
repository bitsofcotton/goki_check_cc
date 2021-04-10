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
// \[- &alpha, &alpha;\[ register computer with deterministic calculation
// that whole in/out is shown in the varlen.
// Please refer bitsofcotton/randtools .
// [[A I], [I A^-1]] [x y], if rank A = dim x = dim y,
// rank [[A I], [I A^-1]] == rank A.
// so with this, U [[B], [O]] [x y] == 0.
// this returns one of u_k that concerns 0.
// suppose input stream is in [-1,1]^k.
// with noised ones, we can use catgp instead of this.
template <typename T> SimpleVector<T> invariantP1(const SimpleVector<T>& in, const int& varlen, const int& ratio = 0) {
  SimpleMatrix<T> A((in.size() - varlen + 1) * 2 - 1, varlen + 2);
  assert(A.cols() <= (A.rows() + 1) / 2);
  SimpleVector<T> fvec(A.cols());
  SimpleVector<T> one(A.rows());
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < fvec.size(); i ++)
    fvec[i] = T(0);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < A.rows(); i ++)
    one[i] = T(1);
  one /= sqrt(one.dot(one));
  for(int i = 0; i < (A.rows() + 1) / 2; i ++) {
    for(int j = 0; j < varlen; j ++)
      A(i, j) = in[i + j];
    A(i, varlen) = T(1) / sqrt(T(A.rows() * A.cols()));
    A(i, varlen + 1) = T(i + 1) / T((A.rows() + 1) / 2 + 1) / sqrt(T(A.rows() * A.cols()));
    T pd(0);
    for(int j = 0; j < A.cols(); j ++)
      pd += log(A(i, j) + T(2));
    A.row(i) *= exp((ratio ? T(ratio) * pd : - pd) / T(A.cols()));
    assert(isfinite(A.row(i).dot(A.row(i))));
    if(i + (A.rows() + 1) / 2 < A.rows())
      A.row(i + (A.rows() + 1) / 2) = - A.row(i);
  }
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
  for(int n_fixed = 0; n_fixed < Pt.rows() - 1; n_fixed ++) {
    const auto on(Pt.projectionPt(- one));
    auto fidx(0);
    for( ; fidx < on.size(); fidx ++)
      if(T(0) < on[fidx]) break;
    for(int i = 1; i < on.size(); i ++)
      if(T(0) < on[i] && on[i] <= on[fidx]) fidx = i;
    if(on[fidx] <= T(0)) break;
    const auto orth(Pt.col(fidx));
    const auto norm2orth(orth.dot(orth));
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
    for(int j = 0; j < Pt.cols(); j ++)
      Pt.setCol(j, Pt.col(j) - orth * Pt.col(j).dot(orth) / norm2orth);
  }
  fvec = R.solve(- Pt * one);
  return isfinite(fvec.dot(fvec)) ? fvec : fvec * T(0);
}


template <typename T> class P1I {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline P1I();
  inline P1I(const int& stat, const int& var);
  inline ~P1I();
  inline T next(const T& in, const int& ratio = 0);
private:
  Vec buf;
  int varlen;
  int t;
};

template <typename T> inline P1I<T>::P1I() {
  varlen = t = 0;
}

template <typename T> inline P1I<T>::P1I(const int& stat, const int& var) {
  assert(0 <= stat && 1 <= var);
  buf.resize(stat + (varlen = var) * 2 + 1);
  for(int i = 0; i < buf.size(); i ++)
    buf[i] = T(0);
  t = 0;
}

template <typename T> inline P1I<T>::~P1I() {
  ;
}

template <typename T> T P1I<T>::next(const T& in, const int& ratio) {
  for(int i = 0; i < buf.size() - 1; i ++)
    buf[i]  = buf[i + 1];
  buf[buf.size() - 1] = in;
  if(t ++ < buf.size()) return T(0);
  // N.B. to compete with noise, we calculate each.
  //      we can use catgp on worse noised ones.
  const auto nin(sqrt(buf.dot(buf)));
  const auto invariant(invariantP1<T>(buf / nin, varlen, ratio));
        auto work(invariant);
  for(int i = 1; i < varlen; i ++)
    work[i - 1] = buf[i - varlen + buf.size()] / nin;
  work[varlen - 1] = work[varlen - 2];
  work[varlen + 1] = work[varlen] =
    T(1) / sqrt(T((buf.size() - varlen + 1) * 2) * T(varlen + 2));
  return (invariant.dot(work) - invariant[varlen - 1] * work[varlen - 1]) / invariant[varlen - 1] * nin;
}

#define _P1_
#endif

