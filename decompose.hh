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
         Vec mimic(const Vec& dst, const Vec& src, const T& intensity = T(1)) const;
         Vec emphasis(const Vec& dst, const T& intensity = T(1)) const;
  inline Vec mimic0(const Vec& dst, const Vec& src) const;
         Vec next(const Vec& in0) const;
  Mat complementMat(const Vec& in) const;
  T   lasterr;
  std::vector<Mat> bA;
  Mat A;
         Vec prepare(const Vec& in, const int& idx = 0) const;
         Vec apply(const Vec& v, const Vec& dst, const Vec& src, const int& idx = 0) const;
};

template <typename T> inline Decompose<T>::Decompose() {
  ;
}

template <typename T> inline Decompose<T>::Decompose(const int& size) {
  static P0<T> p0;
  A.resize(size, size);
  for(int i = 0; i < A.rows(); i ++)
    for(int j = 0; j < A.cols(); j ++)
      A(i, j) = T(0);
  for(int i = 0; i < size; i ++) {
    SimpleMatrix<T> AA(size, size);
    for(int j = 0; j < A.rows(); j ++) {
      const auto jj(T(j) * T(i + 1) / T(size));
      AA.row(j) = p0.taylor(A.cols(), (jj - floor(jj)) * T(A.cols()));
    }
    bA.emplace_back(AA);
    A += AA;
  }
}

template <typename T> inline Decompose<T>::~Decompose() {
  ;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::mimic(const Vec& dst, const Vec& src, const T& intensity) const {
  const int  size(bA.size());
  const auto size2(dst.size() / size);
  const auto size3(src.size() / size);
        auto res(dst);
  for(int i = 0; i < size2; i ++) {
    const auto dd(prepare(dst, i));
    const auto m0(apply(dst, mimic0(dd,
                    prepare(src, i * size3 / size2) * intensity +
                    dd * (T(1) - intensity)), dd, i));
    for(int j = 0; j < size; j ++)
      res[(j * size2 + i + size2 / 2) % res.size()] =
        m0[(j * size2 + i + size2 / 2) % m0.size()];
  }
  return res;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::emphasis(const Vec& dst, const T& intensity) const {
  const int  size(bA.size());
  const auto size2(dst.size() / size);
        auto res(dst);
  std::vector<SimpleVector<T> > dd;
  std::vector<SimpleVector<T> > ndd;
  std::vector<SimpleVector<T> > freq;
  dd.reserve(size2);
  ndd.reserve(size2);
  freq.reserve(size2);
  T normfreq(0);
  for(int i = 0; i < size2; i ++) {
    dd.emplace_back(prepare(dst, i));
    ndd.emplace_back(next(dd[i]));
    freq.emplace_back(complementMat(ndd[i]).solve(dd[i]));
    normfreq += freq[i].dot(freq[i]);
  }
  normfreq /= size2;
  normfreq  = sqrt(normfreq);
  for(int i = 0; i < size2; i ++) {
    for(int j = 0; j < freq[i].size() - 1; j ++)
      freq[i][j] += intensity * T(j + 1) / T(freq[i].size() - 1) * normfreq;
    const auto m0(apply(dst, complementMat(ndd[i]) * freq[i], dd[i], i));
    for(int j = 0; j < size; j ++)
      res[(j * size2 + i + size2 / 2) % res.size()] =
        m0[(j * size2 + i + size2 / 2) % m0.size()];
  }
  return res;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::prepare(const Vec& in, const int& idx) const {
  const int  size(bA.size());
  const auto cnt(int(in.size() + size - 1) / size - 1);
  assert(0 < cnt);
  Vec res(size);
#if defined(_OPENMP)
#pragma omp parallel
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++) {
    res[i] = T(0);
    const auto ibegin(i * cnt);
    const auto iend(i < size - 1 ? min((i + 1) * cnt, int(in.size())) : int(in.size()));
    for(int j = ibegin; j < iend; j ++)
      res[i] += in[(j + idx) % in.size()];
    res[i] /= T(iend - ibegin);
  }
  return res;
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::apply(const Vec& v, const Vec& dst, const Vec& src, const int& idx) const {
  const int  size(bA.size());
  assert(dst.size() == size && src.size() == size);
  const auto cnt(int(v.size() + size - 1) / size - 1);
  assert(0 < cnt);
  auto res(v);
#if defined(_OPENMP)
#pragma omp for schedule(static, 1)
#endif
  for(int i = 0; i < size; i ++) {
    const auto ratio(dst[i] / src[i]);
    for(int j = i * cnt;
            j < (i < res.size() - 1 ? min((i + 1) * cnt, int(res.size())) : int(res.size()));
            j ++)
      res[(j + idx) % res.size()] *= ratio;
  }
  return res;
}

template <typename T> inline typename Decompose<T>::Vec Decompose<T>::mimic0(const Vec& dst, const Vec& src) const {
  return complementMat(next(dst)) * complementMat(next(src)).solve(dst);
}

template <typename T> typename Decompose<T>::Vec Decompose<T>::next(const Vec& in0) const {
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

template <typename T> typename Decompose<T>::Mat Decompose<T>::complementMat(const Vec& in) const {
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

