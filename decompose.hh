/*
BSD 3-Clause License

Copyright (c) 2020, kazunobu watatsu
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
#if !defined(_DECOMPOSE_)

template <typename T> class Decompose {
public:
  typedef SimpleVector<T> Vec;
  typedef SimpleMatrix<T> Mat;
  inline Decompose();
  inline Decompose(const int& size);
  inline ~Decompose();
         Vec mimic(const Vec& dst, const Vec& src, const T& intensity = T(1));
  inline Vec mimic0(const Vec& dst, const Vec& src);
         Vec next(const Vec& in0);
  Mat complementMat(const Vec& in);
  T   lasterr;
  std::vector<Mat> bA;
  Mat A;
};

template <typename T> inline Decompose<T>::Decompose() {
  ;
}

template <typename T> inline Decompose<T>::Decompose(const int& size) {
  P0<T> p0;
  A.resize(size, size);
  for(int i = 0; i < A.rows(); i ++)
    for(int j = 0; j < A.cols(); j ++)
      A(i, j) = T(i == j ? 1 : 0);
  bA.emplace_back(A);
  for(int i = 1; i < size; i ++) {
    SimpleMatrix<T> AA(size, size);
    for(int j = 0; j < A.rows(); j ++) {
      const auto jj(T(j) * T(i + 1) / T(size - 1));
      AA.row(j) = p0.taylor(A.cols(), (jj - floor(jj)) * T(size - 1));
    }
    bA.emplace_back(AA);
    A += AA;
  }
}

template <typename T> inline Decompose<T>::~Decompose() {
  ;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::mimic(const Vec& dst, const Vec& src, const T& intensity) {
  const int  size(bA.size());
  const auto dcnt(int(dst.size() + size - 1) / size - 1);
  const auto scnt(int(src.size() + size - 1) / size - 1);
  assert(0 < dcnt && 0 < scnt);
  Vec dd(size);
  Vec ss(size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++) {
    dd[i] = ss[i] = T(0);
    const auto dbegin(i * dcnt);
    const auto dend(i < size - 1 ? min((i + 1) * dcnt, int(dst.size())) : int(dst.size()));
    const auto sbegin(i * scnt);
    const auto send(i < size - 1 ? min((i + 1) * scnt, int(src.size())) : int(src.size()));
    for(int j = dbegin; j < dend; j ++)
      dd[i] += dst[j];
    for(int j = sbegin; j < send; j ++)
      ss[i] += src[j];
    dd[i] /= T(dend - dbegin);
    ss[i] /= T(send - sbegin);
  }
  const auto ddd(mimic0(dd, ss) * intensity + dd * (T(1) - intensity));
  auto res(dst);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < ddd.size(); i ++) {
    const auto ratio(ddd[i] / dd[i]);
    for(int j = i * dcnt;
            j < (i < ddd.size() - 1 ? min((i + 1) * dcnt, int(dst.size())) : int(dst.size()));
            j ++)
      res[j] *= ratio;
  }
  return res;
}

template <typename T> inline typename Decompose<T>::Vec Decompose<T>::mimic0(const Vec& dst, const Vec& src) {
  return complementMat(next(dst)) * complementMat(next(src)).solve(dst);
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::next(const Vec& in0) {
  assert(A.rows() == in0.size());
  auto in(in0);
  auto avg(in[0]);
  for(int i = 1; i < in.size(); i ++)
    avg += in[i];
  avg /= T(in.size());
  for(int i = 0; i < in.size(); i ++)
    in[i] -= avg;
  auto f(A.solve(in));
  return f /= f.dot(f);
}

template <typename T> typename Decompose<T>::Mat Decompose<T>::complementMat(const Vec& in) {
  Mat B(A.rows(), A.cols());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < bA.size(); i ++)
    B.setCol(i, bA[i] * in);
  return B;
}

#define _DECOMPOSE_
#endif

